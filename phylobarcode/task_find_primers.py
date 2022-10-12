#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import re, numpy as np, pandas as pd

# legacy code, no need to create a separate logger
#log_format = logging.Formatter(fmt='phylobarcode_primer %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

def find_primers (fastafile = None, primer_opt_size = 20, border = 400, num_return = 100, taxon = None, output = None):
    if primer_opt_size is None: primer_opt_size = 20
    if border is None:          border = 400
    if num_return is None:      num_return = 100
    if output is None:
        output = "primers." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, writing to file {output}")

    taxon_df = get_taxonomy (taxon) # may be None is taxon is None
    fas = read_fasta_as_list (fastafile)
    ldic = {}
    rdic = {}
    for i, seqfasta in enumerate(fas):
        if i and not i%100:
            logger.info (f"Processing sequence {i}")
        left, right = get_primers (seqfasta.seq, seqfasta.id, border = border, num_return = num_return, taxon_df = taxon_df)
        if left is not None:
            for x in left:
                if x[0] not in ldic: ldic[x[0]] = [x[1:]]
                else:                ldic[x[0]].append(x[1:])
        if right is not None:
            for x in right:
                if x[0] not in rdic: rdic[x[0]] = [x[1:]]
                else:                rdic[x[0]].append(x[1:])
    save_primers_to_file (ldic, rdic, output)

def get_primers (sequence, seqname, primer_opt_size = 20, border = 400, num_return=100, taxon_df = None):
    seqlen = len(sequence)
    primer_task = "pick_primer_list" # "pick_sequencing_primers" "generic" "pick_primer_list"
    primer_min_size = 14
    primer_max_size = 28
    if primer_opt_size < primer_min_size:
        primer_opt_size = primer_min_size
    if primer_opt_size > primer_max_size:
        primer_opt_size = primer_max_size
    max_n_accepted = 1
    if border < primer_opt_size:
        border = primer_opt_size
    if num_return < 1:
        num_return = 1
    if seqlen - 2 * border - 2 * primer_opt_size < seqlen * 0.1:
        border = int(seqlen * 0.45 - primer_opt_size)
        if border < primer_opt_size:
            logger.info (f"Border is too small for sequence {seqname}")
            return None, None;
    min_product_size = seqlen - 2 * border - 2 * primer_opt_size;
    if (min_product_size > 0.8 * seqlen):
        min_product_size = int(0.8 * seqlen)

    arguments = (     # sequence must be in one line, and there must be no extra newlines or spaces
        "PRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_INTERNAL_OLIGO=0\nPRIMER_PICK_RIGHT_PRIMER=1\n"
        f"SEQUENCE_ID={seqname}\nPRIMER_TASK={primer_task}\n"
        f"PRIMER_MIN_SIZE={primer_min_size}\nPRIMER_MAX_SIZE={primer_max_size}\n"
        f"PRIMER_MAX_NS_ACCEPTED={max_n_accepted}\n"
        f"PRIMER_PRODUCT_SIZE_RANGE={min_product_size}-{seqlen}\n" # 5000-9500
        f"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,{border},{seqlen-border},{border}\n"  # left, left_len, right, right_len : 0,500,9000,554
        f"PRIMER_NUM_RETURN={num_return}\nPRIMER_OPT_SIZE={primer_opt_size}\n"
        f"SEQUENCE_TEMPLATE={sequence}\n=\n"
        )
    child = subprocess.Popen(" primer3_core " , stdin=subprocess.PIPE, stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
    output = child.communicate(input=arguments)[0].split("\n")

    # seqname is used in warnings and to map to taxon; seqlen to calculate distance from border
    return extract_primer_from_output (output, seqname, seqlen, taxon_df) 

def find_primers_parallel (fastafile = None, primer_opt_size = 20, border = 400, num_return = 100, taxon = None, 
        output = None, nthreads = 2):
    if primer_opt_size is None: primer_opt_size = 20
    if border is None:          border = 400
    if num_return is None:      num_return = 100
    if output is None:
        output = "findprimer." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, writing to files {output}_<suffix>")
    from multiprocessing import Pool
    from functools import partial

    taxon_df = get_taxonomy (taxon) # may be None is taxon is None

    fas = read_fasta_as_list (fastafile)
    if nthreads > len(fas): nthreads = len(fas)
    chunk_size = len(fas)//nthreads + 1 # no garantee that all threads will be used, specially if chumk_size is small like 2 or 3
    chunk_ids = [i for i in range(0, len(fas), chunk_size)] + [len(fas)] # from [i] to [i+1] excusive
    nthreads = len(chunk_ids) - 1
    logger.info (f"Using {nthreads} threads from those available")

    with Pool(nthreads) as p:
        results = p.map(partial(
            pool_get_primers_parallel, border=border, num_return=num_return, fasta=fas, ids=chunk_ids, taxon_df = taxon_df), 
            [i for i in range(nthreads)])
    
    ldic = {}
    rdic = {}
    for r_threaded in results:
        for [left,right] in r_threaded:
            if left is not None:
                for x in left:
                    if x[0] not in ldic: ldic[x[0]] = [x[1:]]
                    else:                ldic[x[0]].append(x[1:])
            if right is not None:
                for x in right:
                    if x[0] not in rdic: rdic[x[0]] = [x[1:]]
                    else:                rdic[x[0]].append(x[1:])
    save_primers_to_file (ldic, rdic, output)

def pool_get_primers_parallel (thread_number, border = 400, num_return=100, fasta=None, ids=None, taxon_df = None):
    res = []
    for i in range(ids[thread_number], ids[thread_number+1]):
        rec = fasta[i]
        left, right = get_primers (sequence = rec.seq, seqname = rec.id, border = border, num_return = num_return, taxon_df = taxon_df)
        res.append([left, right])
    return res

def save_primers_to_file (ldic, rdic, output): # below, each v = dict[k] = [[penalty1, distance1,taxon], [penalty2,distance2,taxon], ...] 
    def get_uniq (x,j): # x[2] is taxon name , x[3] is genus name
        return len(set([i[j] for i in x]))
    def get_fun (func, x, j):
        return func([i[j] for i in x])

    llist = [[k, get_uniq(v,3), get_uniq(v,2), len(v), 
        round(get_fun(max,v,0),4), int(get_fun(min,v,1)), int(get_fun(max,v,1))] for k, v in ldic.items()]
    rlist = [[k, get_uniq(v,3), get_uniq(v,2), len(v), 
        round(get_fun(max,v,0),4), int(get_fun(min,v,1)), int(get_fun(max,v,1))] for k, v in rdic.items()]
    ## this numpy below doesnt work with strings 
    #llist = [[k, len(v), np.amax(v,axis=0)[0], int(np.amin(v,axis=0)[1]), int(np.amax(v,axis=0)[1])] for k, v in ldic.items()]
    #rlist = [[k, len(v), np.amax(v,axis=0)[0], int(np.amin(v,axis=0)[1]), int(np.amax(v,axis=0)[1])] for k, v in rdic.items()]
    llist.sort(key=lambda x: (x[1],x[2],x[3]), reverse=True)
    rlist.sort(key=lambda x: (x[1],x[2],x[3]), reverse=True)

    with open_anyformat (f"{output}_l.tsv.xz", "w") as f:
        f.write (str("primer\tgenus_diversity\ttaxon_diversity\tfrequency\tpenalty\tmin_distance\tmax_distance\n").encode())
        for x in llist:
            row = "\t".join([str(i) for i in x]) + "\n"
            f.write (str(row).encode()) # encode() is needed in case output is compressed
    with open_anyformat (f"{output}_r.tsv.xz", "w") as f:
        f.write (str("primer\tgenus_diversity\ttaxon_diversity\tfrequency\tpenalty\tmin_distance\tmax_distance\n").encode())
        for x in rlist:
            row = "\t".join([str(i) for i in x]) + "\n"
            f.write (str(row).encode()) # encode() is needed in case output is compressed
    logger.info (f"Saved primers to {output}_l.tsv.xz and {output}_r.tsv.xz")

def extract_primer_from_output (output, seqname, seqlen, taxon_df = None):
    '''
    Extracts the primers from the output of primer3_core. Returns lists of left and right primers, with
    [sequence,penalty, position] per primer.
    '''
    taxonname = "unknown"
    if taxon_df is not None and seqname in taxon_df["seqid"].values:
        taxonname = str(taxon_df.loc[taxon_df['seqid'] == seqname, 'taxonid'].iloc[0])
    genusname = "unknown"
    if taxon_df is not None and seqname in taxon_df["seqid"].values:
        genusname = str(taxon_df.loc[taxon_df['seqid'] == seqname, 'genus'].iloc[0])

    def xtract_left_or_right (out, side, seqlen):
        y=[re.match(f"PRIMER_{side}_(\d+)_SEQUENCE=(.*)",x) for x in out] # most will be "None" since won't match
        pseq = [[int(i.group(1)),i.group(2)] for i in y if i is not None]
        y=[re.match(f"PRIMER_{side}_(\d+)_PENALTY=(.*)",x) for x in out] # most will be "None" since won't match
        ppen = [[int(i.group(1)),float(i.group(2))] for i in y if i is not None]
        y=[re.match(f"PRIMER_{side}_(\d+)=(\d+),(\d+)",x) for x in out] ## =start,length (but we dont use length)
        pstart = [[int(i.group(1)),int(i.group(2))] for i in y if i is not None]
        # pseq=[idx, sequence]; ppen=[idx, penalty]; pstart=[idx, start]; we want [sequence,penalty, start]
        if not len(pseq):
            logger.warning (f"No {side} primers found for sequence {seqname}")
            return None
        idx = [i[0] for i in pseq] + [i[0] for i in ppen] + [i[0] for i in pstart] # flattened version of two index lists (e.g. [0,1,2...,0,1,2...])
        min_idx = min(idx) ## make sure it starts from 0 (which is already the case in current version of primer3...)
        spl = [[None,None,None, taxonname, genusname] for i in range (max(idx) + 1 - min_idx)] ## length = max between two lists flattened (base_zero -> max +1 is length)
        for i in pseq:
            spl[i[0]-min_idx][0] = str(i[1]).upper()
        for i in ppen:
            spl[i[0]-min_idx][1] = i[1]
        if seqlen: # "start" of right primer is actually distance from end of sequence
            for i in pstart:
                spl[i[0]-min_idx][2] = seqlen - i[1]
        else: # start is distance from beginning of sequence
            for i in pstart:
                spl[i[0]-min_idx][2] = i[1]
        return [i for i in spl if i[0] is not None] # each elem is [sequence,penalty,dist_from_border,taxonname,genusname]

    left = xtract_left_or_right (output, "LEFT", None)
    right = xtract_left_or_right (output, "RIGHT", seqlen)
    return left, right

def get_taxonomy (taxon):
    if taxon is None:
        logger.info(f"No taxonomic information provided, will not calculate taxonomic representativity")
        return None
    if not os.path.isfile(taxon):
        logger.error(f"Taxonomic information file {taxon} not found")
        return None
    taxon_df = pd.read_csv (taxon, compression="infer", sep="\t", dtype='unicode')
    # use taxonid from gff file and NOT from GTDB taxonomy (since we may have fewer); however gtdb_taxonomy has e.g.
    # d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_H;g__Priestia;s__Priestia megaterium
    taxon_df = taxon_df[["seqid","gff_taxonid", "gtdb_taxonomy"]]
    taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, drop_gtdb_column = True) # True = drop original column
    taxon_df = taxon_df.rename(columns={"gff_taxonid":"taxonid"})
    taxon_df.fillna("unknown", inplace=True)
    logger.info(f"Read {len(taxon_df)} entries with taxonomic information from file {taxon}")
    return taxon_df

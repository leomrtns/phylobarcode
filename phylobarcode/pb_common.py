import os, logging, xxhash
from Bio import Seq, SeqIO, AlignIO
import random, datetime, sys, re, glob, collections, subprocess, itertools, pathlib, base64, string
import lzma, gzip, bz2

# legacy code, now every module shares the same parent logger
#log_format = logging.Formatter(fmt='phylobarcode_common %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

base62 = string.digits + string.ascii_letters + '_.-+~@'  # 66 elements actually (62 is alphanum only)
len_base62 = len (base62)

def remove_prefix_suffix (strlist):
    def all_same(x):  # https://stackoverflow.com/a/6719272/204903
        return all(x[0] == y for y in x)
    if isinstance(strlist, str): return ""
    if len(strlist) == 1: return ""
    char_tuples = zip(*strlist)
    fix_tuples  = itertools.takewhile(all_same, char_tuples)
    prefix = ''.join(x[0] for x in fix_tuples)
    inverse = [x[::-1] for x in strlist]
    char_tuples = zip(*inverse)
    fix_tuples  = itertools.takewhile(all_same, char_tuples)
    suffix = ''.join(x[0] for x in fix_tuples)
    suffix = suffix[::-1]

    l_pre = len(prefix) ## we could skip "prefix" and store lenght of fix_tuples but this is more readable
    l_suf = len(suffix)
    return [x[l_pre:len(x)-l_suf] for x in strlist] # does not work well for 'lefT' and 'righT' 

def split_gtdb_taxonomy_from_dataframe (taxon_df, gtdb_column = "gtdb_taxonomy", drop_gtdb_column = True, replace = None):
    '''
    Splits the GTDB taxonomy string into a list of taxonomic ranks: 
    d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_H;g__Priestia;s__Priestia megaterium
    '''
    linneus = {"phylum}":1, "class":2, "order":3, "family":4, "genus":5, "species":6}
    for k,v in linneus.items():
        taxon_df[k] = taxon_df[gtdb_column].str.split(";").str[v].str.split("__").str[1]
    if replace is not None:
        for k in linneus.keys():
            taxon_df[k] = taxon_df[k].fillna(replace)
    if (drop_gtdb_column): taxon_df.drop(columns = [gtdb_column], inplace = True)
    return taxon_df

def seq_to_base62 (seq): # https://stackoverflow.com/questions/1119722/base-62-conversion
    '''
    Convert a sequence to a base62 string of its integer representation.
    '''
    integer = xxhash.xxh128_intdigest(seq) # alias of xxh3_128_intdigest()
    if integer == 0:
        return base62[0]
    ret = ''
    while integer != 0:
        ret = base62[integer % len_base62] + ret
        integer //= len_base62
    return ret

def read_fasta_as_list (filename, clean_sequence=True, substring=None):
    unaligned = []
    with open_anyformat (filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if clean_sequence:
                record.seq  = Seq.Seq(str(record.seq.upper()).replace(".","N"))
            if substring is None or any([x in record.description for x in substring]):
                unaligned.append(record)
    logger.debug("Read %s sequences from file %s", str(len(unaligned)), filename)
    return unaligned

def read_fasta_headers_as_list (filename):
    seqnames = []
    with open_anyformat (filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqnames.append(record.description)
    logger.debug("Read %s sequence headers from file %s", str(len(seqnames)), filename)
    return seqnames

def mafft_align_seqs (sequences=None, infile = None, outfile = None, prefix = None, nthreads = 1): # list not dict
    if (sequences is None) and (infile is None):
        logger.error("You must give me a fasta object or a file")
        return None
    if prefix is None: prefix = "./"
    hash_name = '%012x' % random.randrange(16**12)
    if infile is None: ifl = f"{prefix}/mafft_{hash_name}.fasta"
    else: ifl = infile # if both infile and sequences are present, it will save (overwrite) infile
    if outfile is None: ofl = f"{prefix}/mafft_{hash_name}.aln"
    else: ofl = outfile # in this case it will not exclude_reference
    if sequences: SeqIO.write(sequences, ifl, "fasta") ## else it should be present in infile
    if nthreads < 1: nthreads = -1 # mafft default to use all available threads

    runstr = f"mafft --auto --ep 0.23 --leavegappyregion --thread {nthreads} {ifl} > {ofl}"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    aligned = AlignIO.read(ofl, "fasta")

    if infile is None:  os.remove(ifl)
    if outfile is None: os.remove(ofl)
    return aligned

def cdhit_cluster_seqs (sequences=None, infile = None, outfile = None, prefix = None, nthreads = 1, 
        id = 0.9, fast = True): # list not dict
    def read_clstr_file (clstr_file):
        clusters = {}
        with open(clstr_file, "r") as clstr:
            for line in clstr:
                if line.startswith(">"):
                    cluster_id = line.strip().split()[1] # e.g. ">Cluster 0"
                    clusters[cluster_id] = []
                else:
                    x = line.strip().split()[2].strip(">")
                    clusters[cluster_id].append(x[:x.index("...")]) # e.g. "0	100nt, >seq1... at 100%"
        clusters = [v for v in clusters.values()] # dict of lists to list of lists (we don't need cluster ids)
        return clusters

    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me a fasta object or a file")
    if prefix is None: prefix = "./"
    hash_name = '%012x' % random.randrange(16**12)
    if infile is None: ifl = f"{prefix}/cdhit_{hash_name}.fasta"
    else: ifl = infile # if both infile and sequences are present, it will save (overwrite) infile
    if outfile is None: ofl = f"{prefix}/cdhit_{hash_name}.reps.fasta"
    else: ofl = outfile # in this case it will not exclude_reference
    if sequences: SeqIO.write(sequences, ifl, "fasta") ## else it should be present in infile
    if nthreads < 1: nthreads = 0 # cdhit default to use all available threads
    if fast is True: algo = "0" # sequence is clustered to the first cluster that meet the threshold
    else: algo = "1" # sequence is clustered to the most similar cluster that meet the threshold

    runstr = f"cd-hit -i {ifl} -o {ofl} -c {id} -M 0 -T {nthreads} -d 0 -aS 0.5 -aL 0.5 -g {algo} -s 0.5 -p 0"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    representatives = SeqIO.parse(ofl, "fasta")
    clusters = read_clstr_file(f"{ofl}.clstr")

    if infile is None:  os.remove(ifl)
    if outfile is None: os.remove(ofl)
    os.remove(f"{ofl}.clstr") ## always remove the clstr file
    return representatives, clusters

def calc_freq_N_from_string (genome):
    l = len(genome)
    if (l):
        number_Ns = sum([genome.upper().count(nuc) for nuc in ["N", "-"]])
        return number_Ns / l
    else: return 1.

def calc_freq_ACGT_from_string (genome):
    l = len(genome)
    if (l):
        number_ACGTs = sum([genome.upper().count(nuc) for nuc in ["A", "C", "G", "T"]])
        return number_ACGTs / l
    else: return 0.

def remove_duplicated_sequences_list (sequences): # input is list, returns a dict
    uniq_seqs = {}
    uniq_qual = {}
    duplicates = []
    for x in sequences:
        seq = str(x.seq)
        quality = len(seq) - sum([seq.upper().count(nuc) for nuc in ["N", "-"]])
        if x.id in uniq_seqs.keys(): # sequence name has been seen before 
            if uniq_qual[x.id] < quality: # replaces if better quality
                uniq_qual[x.id] = quality
                uniq_seqs[x.id] = x
            duplicates.append(x.id)
        else: # first time we see this sequence
            uniq_qual[x.id] = quality
            uniq_seqs[x.id] = x
    if len(duplicates)>0:
        logger.warning ("%s duplicate (i.e. with same name) sequences were resolved by choosing the one with highest quality", len(duplicates))
        duplicates = list(set(duplicates))
        logger.debug ("And the sequence names are:\n%s\n", "\n".join(duplicates))
    else:
        logger.info ("Checked for duplicates but all sequences have distinct names")
    return uniq_seqs, uniq_qual

def save_sequence_dict_to_file (seqs, fname=None, use_seq_id = False):
    if fname is None: fname = "tmp." + '%012x' % random.randrange(16**12) + ".aln.xz"
    logger.info(f"Saving sequences to file {fname}")
    with open_anyformat (fname, "w") as fw: 
        for name, rec in seqs.items():
            if use_seq_id is True: # default is to use dict key
                name = rec.id
            if rec:  ## missing/query sequences
                seq = str(rec.seq)
                fw.write(str(f">{name}\n{seq}\n").encode())
                rec.id = name ## make sure alignment will have same names
    logger.info(f"Finished saving sequences")
    return os.path.basename(fname)

def open_anyformat (fname, mode = "r"):
    if (mode == "r"): openmode = "rt"
    else:             openmode = "wb"
    if   fname.endswith(".bz2"): this_open = bz2.open #if "bz2" in filename[-5:]: this_open = bz2.open
    elif fname.endswith(".gz"):  this_open = gzip.open
    elif fname.endswith(".xz"):  this_open = lzma.open
    else:  
        this_open = open
#      if (mode == "w"): openmode = "w"  ## only raw file for writting doesn't need "wb"
    return this_open (fname, openmode) 


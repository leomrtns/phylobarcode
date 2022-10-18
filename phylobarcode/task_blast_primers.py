#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import pandas as pd, numpy as np
import itertools, io, multiprocessing
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

# legacy code (no need to create a distinct stream)
#log_format = logging.Formatter(fmt='phylobarcode__blast %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

def blast_primers_from_tsv (tsv=None, output=None, database = None, evalue=1, task="blastn", max_target_seqs = 1000,
        taxon=None, nthreads=1):
    if tsv is None: 
        logger.error("No tsv file provided")
        return
    if output is None:
        output = "blastprimers." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, using {output}")
    if database is None:
        logger.error("I need a blast database (full path with prefix of DB in blast format;" +
            "usually filename without '.nal' or '.nsq')")
        return
    if evalue > 10000: evalue = 10000
    if evalue < 1e-3: evalue = 1e-3
    if max_target_seqs > 1e9: max_target_seqs = 1e9
    if max_target_seqs < 1: max_target_seqs = 1
    if task not in ["blastn-short", "blastn"]:
        logger.error("I need a task = 'blastn-short' or 'blastn'; I'll use 'blastn'")

    df_l = pd.read_csv (tsv, compression="infer", sep="\t", dtype='unicode')
    primers_l = df_l["primer"].tolist()
    logger.info(f"Read {len(primers_l)} primers from file {tsv}")
    if taxon:
        taxon_df = pd.read_csv (taxon, compression="infer", sep="\t", dtype='unicode')
        # use taxonid from gff file and NOT from GTDB taxonomy (since we may have fewer); however gtdb_taxonomy has e.g.
        # d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_H;g__Priestia;s__Priestia megaterium
        taxon_df = taxon_df[["seqid","gff_taxonid", "gtdb_taxonomy"]]
        taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, drop_gtdb_column = True) # True = drop original column
        taxon_df = taxon_df.rename(columns={"seqid":"sseqid", "gff_taxonid":"taxonid"})
        logger.info(f"Read {len(taxon_df)} entries with taxonomic information from file {taxon}")
    else:
        taxon_df = None
        logger.info(f"No taxonomic information provided, will not calculate taxonomic representativity")

    if nthreads < 2: # single BLAST job, still this one job will use all available CPUs
        ncpus = multiprocessing.cpu_count()
        if ncpus < 1: ncpus = 1
        logger.info(f"Using single thread with {ncpus} CPUs per blast job")
        logger.info(f"Running blast on left primers from file {tsv}")
        blast_l = query_primers_blastn_on_database (primers_l, database=database, 
            evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus)
        blast_l, df_l = calculate_stats_and_write_files (blast_l, output, df_l, taxon_df)
    else: # user wants several BLAST jobs
        try:
            from multiprocessing import Pool
            from functools import partial
        except ImportError:
            ncpus = multiprocessing.cpu_count()
            logger.error(f"Multiprocessing not available for python, will run one blast process with {ncpus} CPUs")
            logger.info(f"Running blast on left primers from file {tsv}")
            blast_l = query_primers_blastn_on_database (primers_l, database=database, 
                evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus)
            blast_l, df_l = calculate_stats_and_write_files (blast_l, output, df_l, taxon_df)

        else:
            logger.info(f"Running blast on left primers from file {tsv}")
            blast_l = query_primers_blastn_on_database_parallel (primers_l, database=database, 
                evalue=evalue, task=task, max_target_seqs=max_target_seqs, nthreads=nthreads)
            blast_l, df_l = calculate_stats_and_write_files (blast_l, output, df_l, taxon_df)
    
def blast_primers_from_tsvORIGINAL (tsv=None, output=None, database = None, evalue=1, task="blastn", max_target_seqs = 1000,
        taxon=None, nthreads=1):
    if tsv is None: 
        logger.error("Not a single tsv file provided, and I need two (one for for 5'/left and one for 3'/right primers)")
        return
    if len(tsv) !=2:
        logger.error("I need two tsv files, one with left primers and one for right primers")
        return
    if output is None or len(output) != 2: ## this should not happen if function called from main script
        prefix = "blastprimers." + '%012x' % random.randrange(16**12) 
        output = [prefix + "_left", prefix + "_right"]
        logger.warning (f"No output file specified, using {output[0]} and {output[1]}")
    if database is None:
        logger.error("I need a blast database (full path with prefix of DB in blast format;" +
            "usually filename without '.nal' or '.nsq')")
        return
    if evalue > 10000: evalue = 10000
    if evalue < 1e-3: evalue = 1e-3
    if max_target_seqs > 1e9: max_target_seqs = 1e9
    if max_target_seqs < 1: max_target_seqs = 1
    if task not in ["blastn-short", "blastn"]:
        logger.error("I need a task = 'blastn-short' or 'blastn'; I'll use 'blastn'")

    df_l = pd.read_csv (tsv[0], compression="infer", sep="\t", dtype='unicode')
    primers_l = df_l["primer"].tolist()
    logger.info(f"Read {len(primers_l)} primers from file {tsv[0]}")
    df_r = pd.read_csv (tsv[1], compression="infer", sep="\t", dtype='unicode')
    primers_r = df_r["primer"].tolist()
    logger.info(f"Read {len(primers_r)} primers from file {tsv[1]}; Running blast now")
    if taxon:
        taxon_df = pd.read_csv (taxon, compression="infer", sep="\t", dtype='unicode')
        # use taxonid from gff file and NOT from GTDB taxonomy (since we may have fewer); however gtdb_taxonomy has e.g.
        # d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_H;g__Priestia;s__Priestia megaterium
        taxon_df = taxon_df[["seqid","gff_taxonid", "gtdb_taxonomy"]]
        taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, drop_gtdb_column = True) # True = drop original column
        taxon_df = taxon_df.rename(columns={"seqid":"sseqid", "gff_taxonid":"taxonid"})
        logger.info(f"Read {len(taxon_df)} entries with taxonomic information from file {taxon}")
    else:
        taxon_df = None
        logger.info(f"No taxonomic information provided, will not calculate taxonomic representativity")

    if nthreads < 2: # single BLAST job, still this one job will use all available CPUs
        ncpus = multiprocessing.cpu_count()
        if ncpus < 1: ncpus = 1
        logger.info(f"Using single thread with {ncpus} CPUs per blast job")
        logger.info(f"Running blast on left primers from file {tsv[0]}")
        blast_l = query_primers_blastn_on_database (primers_l, database=database, 
            evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus)
        blast_l, df_l = calculate_stats_and_write_files (blast_l, output[0], df_l, taxon_df)
        logger.info(f"Running blast on right primers from file {tsv[1]}")
        blast_r = query_primers_blastn_on_database (primers_r, database=database, 
            evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus)
        blast_r, df_r = calculate_stats_and_write_files (blast_r, output[1], df_r, taxon_df)
    else: # user wants several BLAST jobs
        try:
            from multiprocessing import Pool
            from functools import partial
        except ImportError:
            ncpus = multiprocessing.cpu_count()
            logger.error(f"Multiprocessing not available for python, will run one blast process with {ncpus} CPUs")
            logger.info(f"Running blast on left primers from file {tsv[0]}")
            blast_l = query_primers_blastn_on_database (primers_l, database=database, 
                evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus)
            blast_l, df_l = calculate_stats_and_write_files (blast_l, output[0], df_l, taxon_df)

            logger.info(f"Running blast on right primers from file {tsv[1]}")
            blast_r = query_primers_blastn_on_database (primers_r, database=database, 
                evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus)
            blast_r, df_r = calculate_stats_and_write_files (blast_r, output[1], df_r, taxon_df)

        else:
            logger.info(f"Running blast on left primers from file {tsv[0]}")
            blast_l = query_primers_blastn_on_database_parallel (primers_l, database=database, 
                evalue=evalue, task=task, max_target_seqs=max_target_seqs, nthreads=nthreads)
            blast_l, df_l = calculate_stats_and_write_files (blast_l, output[0], df_l, taxon_df)

            logger.info(f"Running blast on right primers from file {tsv[1]}")
            blast_r = query_primers_blastn_on_database_parallel (primers_r, database=database, 
                evalue=evalue, task=task, max_target_seqs=max_target_seqs, nthreads=nthreads)
            blast_r, df_r = calculate_stats_and_write_files (blast_r, output[1], df_r, taxon_df)
    

def query_primers_blastn_on_database_parallel (primer_list, database=None, evalue = 1, task="blastn-short", max_target_seqs=1000, nthreads=1):
    n_primers = len(primer_list)
    if nthreads > n_primers: nthreads = len(primer_list)
    chunk_size = n_primers//nthreads + 1 # no garantee that all threads will be used, specially if chumk_size is small like 2 or 3
    ## see extract_riboproterins, ids below might not be needed since python handles properly indices exceeding list length
    chunk_ids = [i for i in range(0, n_primers, chunk_size)] + [n_primers] # from [i] to [i+1] excusive
    nthreads = len(chunk_ids) - 1 
    primers_per_thread = [primer_list[chunk_ids[i]:chunk_ids[i+1]] for i in range(nthreads)]
    ncpus = multiprocessing.cpu_count() // nthreads
    if ncpus < 1: ncpus = 1
    logger.info (f"Using {nthreads} threads from those available for blast, and {ncpus} CPUs per blast job")

    from multiprocessing import Pool
    from functools import partial
    with Pool(nthreads) as p:
        results = p.map(
                partial(query_primers_blastn_on_database, database=database, evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus), 
                primers_per_thread)
    results = [i for i in results if i is not None]
    if len(results): return pd.concat(results)
    else:            return None

def query_primers_blastn_on_database (primer_list, database=None, evalue = 1, task="blastn-short", max_target_seqs=1000, ncpus=1):
    """ 
    task = "blastn-short" (more precise) or task = "blastn" 
    nthreads = number of threads for blast query
    """
    if database is None:
        logger.error(f"I need a database (full path with prefix of DB in blast format; usually filename without '.nal' \
                or '.nsq');");
        return 
    ## output is format 6 plus aligned query plus aligned subject
    cline = f"blastn -query - -db {database} -evalue {evalue} -task \"{task}\" -word_size 7 -out - \
              -outfmt \"6 std qseq sseq\"  -max_target_seqs {max_target_seqs} -num_threads {ncpus}"
    child = subprocess.Popen(cline, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                            universal_newlines=True, shell=(sys.platform!="win32"))
    records = [SeqRecord(Seq.Seq(seq), id=seq) for seq in primer_list] ## head is primer itself (like "> ACGTGT") 
    SeqIO.write(records, child.stdin, "fasta")
    child.stdin.close()
    # print(child.stdout.read()) ## can be used only once since child is closed
    try:
        df = pd.read_table (io.StringIO(child.stdout.read()), header=None)
    except pd.errors.EmptyDataError:
        logger.error("No data returned from blast; If you are sure the database location and name are correct (I guess"
                "they are not!), then try increasing evalue and not using `accurate` task")
        return None
    default_outfmt6_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq ssqeq'.strip().split(' ')
    df.columns = default_outfmt6_cols
    return df

def calculate_stats_and_write_files (blast, output, df, taxon_df=None):
    ## add taxonomic information to raw blast results ##
    if taxon_df is not None: blast = blast.merge(taxon_df, on="sseqid", how="left")
    # difference between full primer and matched region
    blast["len_match_diff"] = blast["qseqid"].str.len() - blast["length"]

    ofname = output + ".raw.tsv.xz"
    blast.to_csv (ofname, sep="\t", index=False)
    logger.info(f"Wrote {len(blast)} raw BLAST hits to file {ofname}")
    df0 = stats_merge_blast_and_primers (blast, df)
    ofname = output + ".stats.tsv.xz"
    df0.to_csv (ofname, sep="\t", index=False)
    logger.info(f"Wrote stats about {len(df0)} primers to file {ofname}")
    return blast, df0

def stats_merge_blast_and_primers (blast_df, primer_df):
    """
    add columns to primer tables with information from blast results (and taxonomic info if available)
    """
    def summarise_blast_df (blast_df, suffix=""):
        df = blast_df.groupby("qseqid").agg({
            "sseqid":["count","nunique"], 
            "len_match_diff":[("mean", lambda x: round(np.mean(x),3))],
            "evalue": [("mean", lambda x: round(np.mean(x),3))], # tuple = (name, function); must be in a list
            "bitscore":"max", 
            "qseq":longest_mode, 
            "mismatch":[
                ("mean", lambda x: round(np.mean(x),3)), 
                ("n_mismatches",lambda x: x.gt(0).sum()),
                ("perfect_matches", lambda x: x.eq(0).sum())
                ] # list of tuples 
            })
        df.columns = ["_".join(i) for i in df.columns] # multiindex to single index
        df = df.reset_index() # qseqid is index, not column
        df["sseqid_count"] = round(df["sseqid_count"]/df["sseqid_nunique"], 2) # average number of hits per sequence
        if "taxonid" in blast_df.columns: ## taxonomic information is available
            df2 = blast_df.groupby("qseqid").agg({
                "taxonid":"nunique", 
                "family":"nunique",
                "genus":"nunique",
                "species":"nunique"
                }).reset_index() # qseqid is index, not column
            df = df.merge(df2, on="qseqid", how="left")
            
        df = df.astype({ # nullable integer type https://pandas.pydata.org/pandas-docs/stable/user_guide/integer_na.html (notice capital Int)
            "sseqid_count":'Int64', # unline int, can hold NaN values without converting to float 
            "sseqid_nunique":'Int64',
            "mismatch_n_mismatches":'Int64',
            "mismatch_perfect_matches":'Int64',
            "taxonid":"Int64",
            "family":"Int64",
            "genus":"Int64",
            "species":"Int64"
            }, errors="ignore")
        df.rename(columns={
            "qseqid":"primer", 
            "sseqid_count":"avge_hits" + suffix,
            "sseqid_nunique":"unique_hits" + suffix,
            "len_match_diff_mean":"avge_match_diff" + suffix,
            "evalue_mean":"avge_evalue" + suffix,
            "bitscore_max":"max_bitscore" + suffix,
            "mismatch_mean":"avge_mismatches" + suffix,
            "mismatch_n_mismatches": "n_mismatches" + suffix,
            "mismatch_perfect_matches": "perfect_matches" + suffix,
            "taxonid": "unique_taxonid" + suffix,
            "family": "unique_family" + suffix,
            "genus": "unique_genus" + suffix,
            "species": "unique_species" + suffix,
            "qseq_longest_mode":"aligned_primer" + suffix
            }, inplace=True)
        return df

    def longest_mode (x):
        x = pd.Series.mode(x)
        if isinstance(x, str): return x
        return sorted(x, key=len, reverse=True)[0]
    
    df = summarise_blast_df (blast_df, "_all")

    ## now on good matches only 
    blast_df = blast_df[blast_df["len_match_diff"] < 2] # almost perfect matches
    df2 = summarise_blast_df (blast_df, "")
    df = df.merge(df2, on="primer", how="left")

    primer_df = primer_df.merge(df, on="primer", how="left")

    return sort_primers_by_performance (primer_df)

sort_columns_dict = {
        "clust_idx": {"sort":True, "type":"string"},
        "unique_hits": {"sort":False, "type":"Int64"},
        "avge_hits": {"sort":True, "type":"float"},
        "unique_family": {"sort":False, "type":"Int64"},
        "unique_genus": {"sort":False, "type":"Int64"},
        "unique_species": {"sort":False, "type":"Int64"},
        "avge_hits_all": {"sort":True, "type":"float"},
        "perfect_matches": {"sort":False, "type":"Int64"},
        "avge_mismatches": {"sort":True, "type":"float"},
        "n_mismatches": {"sort":True, "type":"Int64"},
        "taxon_diversity": {"sort":False, "type":"Int64"}, # how many taxonid this primer was generated from
        "frequency": {"sort":False, "type":"Int64"}, # how many sequences this primer was generated from
        }

def sort_primers_by_performance (df):
    sort_columns = [k for k in sort_columns_dict.keys() if k in df.columns] # taxon info may be missing
    for k in sort_columns: # must use "Int64" type for nullable integer columns
        df[k] = df[k].astype(sort_columns_dict[k]["type"])
    return df.sort_values(sort_columns, ascending=[sort_columns_dict[k]["sort"] for k in sort_columns], ignore_index=True)

def select_primers_from_tsv (tsv = None, subsample = None, n_elements = None, output = None):
    '''
    First, select subsample percent of primers from each variable, and then select n_elements primers from the resulting
    set.  Also, delete unused columns from "blast_primers"
    '''
    if tsv is None: raise ValueError("tsv file is required")

    df = pd.read_csv(tsv, sep="\t", compression="infer", dtype="unicode")
    logger.info(f"Read {len(df)} primers from file {tsv}")

    if subsample is None: subsample = 90
    if n_elements is None or n_elements > len(df): n_elements = len(df)
    subsample = int(len(df) * (subsample / 100))
    if subsample > len(df): subsample = len(df)
    if subsample < 1: subsample = 1
    if n_elements < 1: n_elements = int (n_elements * len(df)) ## can be a fraction instead of an integer
    if n_elements < 1: n_elements = 1
    n_elements = int(n_elements)

    if output is None:
        output = "selected." + '%012x' % random.randrange(16**12) 
        logger.warning(f"No output file name was provided, using {output}.tsv")

    # remove unused columns (also makes it clear to user which columns are used)
    to_prune = [i for i in df.columns if i.endswith("_all")]
    df.drop(columns=to_prune, inplace=True)
    logger.info(f"Removed {len(to_prune)} columns (unrelated to best hits); will now select {n_elements} primers " + 
            f"after univariate selection of {subsample} primers from each variable")
    # select best primers according the each variable
    sel_primers = []
    cols_present = [x for x in sort_columns_dict.keys() if x in df.columns]
    for col in cols_present:
        df1 = df.sort_values(col, ascending=sort_columns_dict[col]["sort"], ignore_index=True)
        sel_primers.append(set(df1.iloc[:subsample]["primer"]))
    sel_primers = set.union(*sel_primers)
    df = df[df["primer"].isin(sel_primers)]
    logger.info(f"Selected {len(sel_primers)} primers based on univariate subsampling")
    # sort and avoid consecutive primers from same cluster
    df = sort_primers_by_performance(df)
    if "cluster" in df.columns:
        df["clust_idx"] = df.groupby("cluster").cumcount() # idx = 1 best, idx = 2 second best, etc.
        df = sort_primers_by_performance(df) # sort again, separating primers from same cluster
        df.drop(columns=["clust_idx"], inplace=True) # remove temporary column
        logger.info(f"Reordered to avoid consecutive primers from same cluster; {len(df['cluster'].unique())} clusters found")

    # select best n_elements primers according to sorting
    df = df.head(int(n_elements))
    df.to_csv(f"{output}.tsv", sep="\t", index=False)
    logger.info(f"Selected {len(df)} primers and saved to {output}.tsv")
    return 

def cross_hits_between_primer_sets (df_l, df_r): ## TODO: unfinished
    """
    original "refseq_riboprotein" compared all primer pairs for a given genome, keeping 
    only primers with a single hit in the genome. 
    """
    def remove_close_numbers (x0, closeness=2):
        """
        from set of numbers closer than `closeness`, keep average between extremes
        """
        x = sorted(x0)
        l = len(x)
        b1 = [x[0]] + [x[i]     for i in range(1,l) if x[  i] - x[  i-1] > closeness]
        b2 =          [x[l-i-1] for i in range(1,l) if x[l-i] - x[l-i-1] > closeness][::-1] + [x[l-1]]
        return [(g+h)/2 for g,h in zip(b1,b2)]

    return df_l, df_r


##    legacy code    ##

def query_string_blastn_on_database (seqstring, database=None, evalue = 1, task="blastn", max_target_seqs=1000, nthreads=1):
    """ 
    task = "blastn-short" (more precise) or task = "blastn" 
    this is the original function using biopython NCBIXML (you need to see refseq_riboprotein.py for how to parse output)
    """
    if database is None:
        database = database_fallback
        logger.error(f"I need a database (full path with prefix of DB in blast format; usually filename without '.nal' \
                or '.nsq'); I will try to use {database} but there are no guarantees it will work")
    cline = f"blastn -query - -db {database} -evalue {evalue} -task \"{task}\" -word_size 7 -out - \
              -outfmt 5 -max_target_seqs {max_target_seqs} -num_threads {nthreads}"
    child = subprocess.Popen(cline, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                            universal_newlines=True, shell=(sys.platform!="win32"))
    SeqIO.write(SeqRecord(Seq.Seq(seqstring)), child.stdin, "fasta")
    child.stdin.close()
    return NCBIXML.read(child.stdout)



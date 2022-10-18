#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import pandas as pd, numpy as np
import itertools, io, multiprocessing
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("phylobarcode_global_logger")

def pair_primers_from_raw_tsv (tsv=None, output=None, taxon=None, nthreads=1):
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
    
    if taxon is not None:
        taxon_df = pd.read_csv (taxon, compression="infer", sep="\t", dtype='unicode')
        taxon_df = taxon_df[["seqid","gff_taxonid", "gtdb_taxonomy"]]
        taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, drop_gtdb_column = True) # True = drop original column
        # use taxonid from gff file and NOT from GTDB taxonomy (since we may have fewer); use "genus" etc from GTDB
        taxon_df = taxon_df.rename(columns={"seqid":"sseqid", "gff_taxonid":"taxonid"})
        logger.info(f"Read {len(taxon_df)} entries with taxonomic information from file {taxon}")
    else: 
        taxon_df = None
        logger.info(f"No taxonomic information provided, I expect the raw blast files to have it already then")

    df_l = read_raw_blast_and_add_taxonomy (tsv[0], taxon_df)
    df_r = read_raw_blast_and_add_taxonomy (tsv[1], taxon_df)
    genomes_l = extract_location_from_genomes (df_l)
    genomes_r = extract_location_from_genomes (df_r)
    ## STOPHERE

def read_raw_blast_and_add_taxonomy (tsv, taxon_df=None):
    df = pd.read_csv (tsv, compression="infer", sep="\t", dtype='unicode')
    if taxon_df is None and "genus" not in df.columns:
        logger.error(f"Taxonomic table (from `merge_fasta_gff`) not provided, and no columns `taxonid` or `genus` in raw blast file {tsv}")
        sys.exit(1)
    if taxon_df is None: 
        logger.info(f"Columns `taxonid` and `genus` from raw blast file {tsv} will be used since no taxonomic table (from `merge_fasta_gff`) was provided")
        return df
    # we have taxon_df and perhaps taxon information from raw blast file
    replace_cols = ["taxonid", "phylum", "class", "order", "family", "genus", "species"]
    replace_cols = [c for c in replace_cols if c in taxon_df.columns] # perhaps the columns changed? anyway this should not happen
    replace_cols = [i for i in replace_cols if i in df.columns]
    if len(replace_cols) > 0:
        logger.info(f"Taxonomic info from raw blast file will be replaced with information from taxonomic table")
        df.drop (replace_cols, axis=1, inplace=True)
    df = df.merge(taxon_df, on="sseqid", how="left")
    logger.info(f"Read {len(df)} entries from raw blast file {tsv} and added taxonomic information")
    return df

def extract_location_from_genomes (df):
    df_start = df.groupby("sseqid")["sstart"].apply(list).reset_index(name='location')
    df_end = df.groupby("sseqid")["send"].apply(list).reset_index(name='location') # each row will have "sseqid" and list of locations
    locs = {}
    for row in df_start.itertuples(index=False):
        locs[row.sseqid] = remove_close_numbers (row.location + df_end[df_end.sseqid == row.sseqid].iloc[0].location)
    ## STOPHERE


def remove_close_numbers (x0, closeness=2):
    """
    from set of numbers closer than `closeness`, keep average between extremes
    """
    x = sorted(x0)
    l = len(x)
    b1 = [x[0]] + [x[i]     for i in range(1,l) if x[  i] - x[  i-1] > closeness]
    b2 =          [x[l-i-1] for i in range(1,l) if x[l-i] - x[l-i-1] > closeness][::-1] + [x[l-1]]
    return [(g+h)/2 for g,h in zip(b1,b2)]

def cross_hits_between_primer_sets (df_l, df_r): ## TODO: unfinished
    """
    original "refseq_riboprotein" compared all primer pairs for a given genome, keeping 
    only primers with a single hit in the genome. 
    """

    return df_l, df_r

## note: 
## dist = [abs(x.sbjct_start - y) for x in hitaln.hsps for y in loc_r_genome[genome_sample_id]]
## dist = [x if x < hitaln.length/2 else hitaln.length - x for x in dist] # distance to closest end of genome

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

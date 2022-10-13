#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
from phylobarcode.pb_kmer import * 
import pandas as pd, numpy as np
import itertools, pathlib, shutil, gzip
from Bio import pairwise2
from sklearn import cluster

# legacy code, no need to create a distinct stream
# log_format = logging.Formatter(fmt='phylobarcode_clustr %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

# FIXME: use parasail.sg_stats_striped_16().matches and len(seq) instead of scores 
# (see https://manpages.debian.org/stretch/vsearch/vsearch.1 `--iddef` for possibilities)

def cluster_flanks_from_fasta (fastafile = None, border = 400, output = None, identity = None, min_samples = 2, scratch = None, nthreads=1):
    if border is None: border = 400
    if output is None:
        output = "flanking." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, writing to file {output}")
    if identity and not isinstance(identity, float):
        logger.warning ("Identity threshold must be a float, reverting to default=0.95")
        identity = 0.95
    if identity and identity < 0.1:
        logger.warning (f"Identity threshold must be between 0.1 and 1, value {identity} is too low, setting to 0.1")
        identity = 0.1
    if identity and identity > 1.0:
        logger.warning (f"Identity threshold must be between 0.1 and 1, value {identity} is too high, setting to 1.0")
        identity = 1.0

    fas = read_fasta_as_list (fastafile)
    flank_l = []
    flank_r = []
    ofile = [f"{output}_l.fasta.gz",f"{output}_r.fasta.gz"]
    with open_anyformat (ofile[0], "w") as f_l, open_anyformat (ofile[1], "w") as f_r:
        for seqfasta in fas:
            seqfasta.seq = str(seqfasta.seq) # it's a SeqRecord
            seqlen = len(seqfasta.seq)
            this_b = border
            if seqlen - 2 * border < seqlen * 0.1:
                this_b = int (seqlen * 0.45)
            f_l.write (str(f">{seqfasta.id}\n{seqfasta.seq[:this_b]}\n").encode())
            f_r.write (str(f">{seqfasta.id}\n{seqfasta.seq[seqlen-this_b:]}\n").encode())
            flank_l.append (seqfasta.seq[:this_b])
            flank_r.append (seqfasta.seq[seqlen-this_b:])
    logger.info (f"Flanking regions from all samples saved to {ofile[0]} and {ofile[1]};")

    if (identity is not None):
        ofile_centroid = [f"{output}_centroid_l.fasta.gz",f"{output}_centroid_r.fasta.gz"]
        if scratch: 
            pathlib.Path(scratch).mkdir(parents=True, exist_ok=True)
        logger.info (f"Clustering and excluding redundant left sequences with vsearch")
        find_centroids_from_file_vsearch (fastafile=ofile[0], output=ofile_centroid[0], identity=identity, nthreads=nthreads, scratch=scratch)
        logger.info (f"Clustering and excluding redundant right sequences with vsearch")
        find_centroids_from_file_vsearch (fastafile=ofile[1], output=ofile_centroid[1], identity=identity, nthreads=nthreads, scratch=scratch)
        if scratch: # delete scratch subdirectory
            shutil.rmtree(pathlib.Path(scratch))
    else:
        ofile_centroid = [f"{output}_clustered_l.fasta.gz",f"{output}_clustered_r.fasta.gz"]
        seqnames = [i.id for i in fas]
        logger.info (f"Clustering and chosing representative left sequences with OPTICS")
        find_representatives_from_sequences_optics (flank_l, names=seqnames, output=ofile_centroid[0], min_samples=min_samples, nthreads=nthreads)
        logger.info (f"Clustering and chosing representative right sequences with OPTICS")
        find_representatives_from_sequences_optics (flank_r, names=seqnames, output=ofile_centroid[1], min_samples=min_samples, nthreads=nthreads)
    
    logger.info (f"Finished. Reduced sequence sets saved to files {ofile_centroid[0]} and {ofile_centroid[1]};")
    return
    
def cluster_primers_from_tsv (tsv = None, output = None, min_samples = 10, subsample=100, kmer_length = 5, 
        threshold = 0.7, n_best = -1, nthreads = 1):
    if tsv is None: 
        logger.error("No tsv file provided")
        return
    if output is None:
        output = "clusters." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, writing to file {output}.tsv.xz")
    if min_samples < 3: min_samples = 3
    if subsample < 1e-3: subsample = 1e-3
    if subsample > 100: subsample = 100
    if kmer_length < 3: kmer_length = 3
    if kmer_length > 18: kmer_length = 18
    if threshold < 0.01: threshold = 0.01
    if threshold > 1: threshold = 1

    df = pd.read_csv (tsv, compression="infer", sep="\t", dtype='unicode')
    logger.info(f"Read {len(df)} primers from file {tsv}")
    if "cluster" in df.columns:
        logger.warning (f"Column 'cluster' already exists in file {tsv}, overwriting")
        df.drop(columns=["cluster"], inplace=True)
    
    df = subsample_primers (df, subsample=subsample)

    df = cluster_primers_with_vsearch (df)
    logger.info (f"Found {len(df['v_cluster'].unique())} clusters with vsearch, writing to file {output}.tsv.xz")
    df.to_csv (f"{output}.tsv.xz", sep="\t", index=False)
    return

    df = kmer_clustering_dataframe (df, kmer_length = kmer_length, threshold = threshold, nthreads=nthreads)
    logger.info(f"Rough kmer-based pre-clustering done, writing to file {output}.tsv.xz")
    df.to_csv (f"{output}.tsv.xz", sep="\t", index=False)

    df = batch_cluster_primers_from_kmer_cluster (df, min_samples=min_samples, nthreads=nthreads)
    if n_best > 0: df = df.groupby("cluster").head(n_best)
        
    logger.info(f"Total of {df['cluster'].nunique()} clusters found, writing to file {output}.tsv.xz")
    df.to_csv (f"{output}.tsv.xz", sep="\t", index=False)
    return

def batch_cluster_primers_from_kmer_cluster (df, min_samples = 5, nthreads = 1):
    df1 = df[df.groupby("kmer_cluster")["kmer_cluster"].transform("count") >= min_samples]
    df1 = df1[["primer", "kmer_cluster"]] # only primer sequence and current cluster
    n_big_clusters = len(df1.kmer_cluster.unique())
    logger.info(f"Clustering {len(df1)} primers from {n_big_clusters} large kmer clusters")
    df_list = []
    for i, cluster_id in enumerate(df1.kmer_cluster.unique()):
        if i and i % (n_big_clusters//10) == 0: logger.info(f"Clustering {i}th kmer cluster")
        primers = df1[df1["kmer_cluster"] == cluster_id]["primer"].tolist()
        primers = [str(i) for i in primers]
        score_mat = create_NW_score_matrix_parallel (primers, nthreads=nthreads)
        distmat = score_to_distance_matrix_fraction (score_mat, mafft=True)
        with np.errstate(divide='ignore'): # silence OPTICS warning (https://stackoverflow.com/a/59405142/204903)
            cl = cluster.OPTICS(min_samples=2, min_cluster_size=2, metric="precomputed", n_jobs=nthreads).fit(distmat)
        inc = itertools.count(max(cl.labels_)) ## counter
        cl.labels_ = [f"{cluster_id}_{next(inc)}" if i == -1 else f"{cluster_id}_{i}" for i in cl.labels_] ## transform noise samples into their own clusters
        df_local = pd.DataFrame({"primer":primers, "cluster":cl.labels_}) # new column "cluster"
        df_list.append(df_local)
    df1 = pd.concat(df_list) ## looks like df1, but has column "cluster" only
    df = df.merge(df1, on="primer", how="left") # merge with original dataframe (which now has "kmer_cluster" and "cluster" columns)
    logger.debug(df)
    df["cluster"].fillna(df["kmer_cluster"], inplace=True) # now "cluster" is a refined version of "kmer_cluster"
    df.drop(columns=["kmer_cluster"], inplace=True)
    # spread out clusters (so that member of same cluster do not appear consecutively)
    df["clust_idx"] = df.groupby("cluster").cumcount() # idx = 1 best, idx = 2 second best, etc.
    df = df.sort_values(
            by=["clust_idx","genus_diversity","taxon_diversity","frequency","max_distance","penalty"], 
            ascending=[True,False,False,False,True,True])
    df.drop(columns=["clust_idx"], inplace=True) # remove temporary column
    return df

def kmer_clustering_dataframe (df, kmer_length = None, threshold = None, nthreads=1):
    if kmer_length is None: kmer_length = 5
    if threshold is None: threshold = 0.7
    primers = df["primer"].tolist()
    df["genus_diversity"] = df["genus_diversity"].astype(int)
    df["taxon_diversity"] = df["taxon_diversity"].astype(int)
    df["frequency"] = df["frequency"].astype(int)
    df["max_distance"] = df["max_distance"].astype(int)
    df["penalty"] = df["penalty"].astype(float)
    df = df.sort_values(
        by=["genus_diversity","taxon_diversity","frequency","max_distance","penalty"],
        ascending=[False,False,False,True,True])

    primers = [str(i) for i in primers]
    kmers = [single_kmer(i,k = kmer_length) for i in primers]
    logger.debug (f"Number of kmers: {len(kmers)}")
    cluster_1,_,_ = cluster_single_kmers_parallel (kmers, threshold=threshold, jaccard=True, nthreads=nthreads)
    n_clusters = len(set(cluster_1))
    logger.debug (f"Number of clusters: {n_clusters}")
    df["kmer_cluster"] = cluster_1
    df["kmer_cluster"] = df["kmer_cluster"].astype(str)
    
    df["clust_idx"] = df.groupby("kmer_cluster").cumcount() # idx = 1 best, idx = 2 second best, etc.
    df = df.sort_values(
            by=["clust_idx","genus_diversity","taxon_diversity","frequency","max_distance","penalty"], 
            ascending=[True,False,False,False,True,True])
    df.drop(columns=["clust_idx"], inplace=True) # remove temporary column
    return df

def subsample_primers (df, subsample=100):
    if subsample >= 100: return df
    logger.info (f"Subsampling {len(df)} primers to {subsample:.2f}% of the original set over frequencies, distance, and penalty (from primer3)")
    subsample /= 100 # pandas percentile uses 0-1 scale
    df["genus_diversity"] = df["genus_diversity"].astype(int)
    df["taxon_diversity"] = df["taxon_diversity"].astype(int)
    df["frequency"] = df["frequency"].astype(int)
    df["max_distance"] = df["max_distance"].astype(int)
    df["penalty"] = df["penalty"].astype(float)
    threshold = {
        "genus_diversity": df["genus_diversity"].quantile(1. - subsample), # always smaller than taxon_diversity
        "taxon_diversity": df["taxon_diversity"].quantile(1. - subsample),
        "frequency": df["frequency"].quantile(1. - subsample), # quantile goes from min to max, thus we want 100-x% smallest value
        "max_distance": df["max_distance"].quantile(subsample),
        "penalty": df["penalty"].quantile(subsample)
        }
    max_div = df["taxon_diversity"].max() # if GTDB file was not given then all have {taxon/genus}_diversity = 1
    if threshold["taxon_diversity"] > 1: # all primers have freq at least one thus all would be chosen  
        df = df[(df["taxon_diversity"] >= threshold["taxon_diversity"]) | 
                (df["genus_diversity"] >= threshold["genus_diversity"]) | # genus_diversity is too strict, used only here
                (df["frequency"] >= threshold["frequency"]) | 
                (df["max_distance"] < threshold["max_distance"]) |
                (df["penalty"] < threshold["penalty"])]
    elif max_div <= 1: # taxonomic info probably missing when finding primers; must neglect this column
        df = df[(df["frequency"] >= threshold["frequency"]) |
                (df["max_distance"] < threshold["max_distance"]) |
                (df["penalty"] < threshold["penalty"])]
    else: # max diversity is > 1 however threshold is 1, thus all primers would have been chosen
        df = df[(df["taxon_diversity"] > 1) | # OR large diversity OR taxon_div=1 but small distance or penalty
                ((df["frequency"] >= threshold["frequency"]) & 
                 (df["max_distance"] < threshold["max_distance"]) &
                 (df["penalty"] < threshold["penalty"]))]
    logger.info (f"Subsampling done, {len(df)} primers kept using thresholds {threshold}")
    df = df.sort_values(
            by=["genus_diversity","taxon_diversity","frequency","max_distance","penalty"], 
            ascending=[False,False,False,True,True])
    return df

def create_NW_score_matrix (seqlist, use_parasail = True, band_size = 0): ## seqs don't need to be aligned, must be strings
    if use_parasail:
        try:  
            import parasail
        except ImportError: 
            logger.warning("Parasail module not installed, reverting to Bio.pairwise2 from Biopython")
            use_parasail = False
    logger.info(f"Calculating pairwise alignment scores for {size} sequences using a single thread")
    size = len(seqlist)
    scoremat = np.zeros((size, size))
    if use_parasail is True and band_size == 0:
        mymat = parasail.pam10
        for i in range(size): 
            scoremat[i,i] = parasail.sg_striped_16(str(seqlist[i]), str(seqlist[i]), 10, 2, mymat).score # nw (sg doenst penalise begin and end )
        for i,j in itertools.combinations(range(size),2): #parasail._stats_ also gives x.length, x.score
            #scoremat[i,j] = scoremat[j,i] = parasail.nw_stats_striped_16(str(seqlist[i]), str(seqlist[j]), 11,1, parasail.blosum30).matches
            scoremat[i,j] = scoremat[j,i] = parasail.sg_striped_16(str(seqlist[i]), str(seqlist[j]), 10, 2, mymat).score
    elif use_parasail is True and isinstance(band_size, int): # banded: not semi-global but "full" NW, with simple penalty matrix
        mymat = parasail.matrix_create("ACGT", 3, -1)
        for i in range(size):
            scoremat[i,i] = parasail.nw_banded(str(seqlist[i]), str(seqlist[i]), 10, 2, band_size, mymat).score # global Needleman-Wunsch 
        for i,j in itertools.combinations(range(size),2): #parasail._stats_ also gives x.length, x.score
            scoremat[i,j] = scoremat[j,i] = parasail.nw_banded(str(seqlist[i]), str(seqlist[j]), 10, 2, band_size, mymat).score
    else:
        for i in range(size): 
            scoremat[i,i] = float(len(seqlist[i]))  # diagonals have sequence lengths (=best possible score!)
        for i,j in itertools.combinations(range(size),2): 
            scoremat[i,j] = scoremat[j,i] = pairwise2.align.globalxx(seqlist[i], seqlist[j], score_only=True)
    return scoremat

def create_NW_score_matrix_parallel (seqlist, use_parasail = True, band_size = 0, nthreads = 1): ## seqs don't need to be aligned, must be strings
    """
    Create a score matrix for a list of pairs of sequences, using the Needleman-Wunsch algorithm.
    """
    if nthreads == 1: return create_NW_score_matrix (seqlist, use_parasail, band_size)
    from multiprocessing import Pool
    from functools import partial
    if use_parasail:
        try:  
            import parasail
        except ImportError: 
            logger.warning("Parasail module not installed, reverting to Bio.pairwise2 from Biopython")
            use_parasail = False
    size = len(seqlist)
    logger.debug(f"Calculating pairwise alignment scores for {size} sequences using {nthreads} threads")
    pairs = list(itertools.combinations(range(size),2)) + [(i,i) for i in range(size)]
    pairchunks = [pairs[i::nthreads] for i in range(nthreads)]

    with Pool(len(pairchunks)) as p:
        scorelist = p.map(
            partial(nw_score_matrix_for_pairlist, seqlist=seqlist, use_parasail=use_parasail, band_size=band_size), 
            pairchunks)
    scorelist = [item for onethread in scorelist for item in onethread] ## flatten list of lists 
    pairs = [item for onethread in pairchunks for item in onethread] ## flatten list of lists
    scoremat = np.zeros((size, size))
    for (i,j), d_ij in zip(pairs, scorelist):
        scoremat[i,j] = scoremat[j,i] = d_ij
    return scoremat

def nw_score_matrix_for_pairlist (pairlist, seqlist, use_parasail, band_size):
    scorelist = []
    if use_parasail is True and band_size == 0:
        import parasail
        mymat = parasail.pam10
        for i, j in pairlist:
            d_ij = parasail.sg_striped_16(seqlist[i], seqlist[j], 10, 2, mymat).score
            scorelist.append(d_ij)
    elif use_parasail is True and isinstance(band_size, int): # banded: not semi-global but "full" NW, with simple penalty matrix
        import parasail
        mymat = parasail.matrix_create("ACGT", 3, -1)
        for i, j in pairlist:
            d_ij = parasail.nw_banded(seqlist[i], seqlist[j], 10, 2, band_size, mymat).score
            scorelist.append(d_ij)
    else:
        for i, j in pairlist:
            if i == j: 
                d_ij = float(len(seqlist[i])) # diagonals have lengths (=best possible score!)
            else:
                d_ij = pairwise2.align.globalxx(seqlist[i], list[j], score_only=True)
            scorelist.append(d_ij)
    return scorelist

def score_to_distance_matrix_fraction (scoremat, mafft = False):
    """
    receives a score matrix (from create_NW_score_matrix) and returns a distance matrix as fraction of indels or Satoh's method (mafft = True)
    """
    distmat = np.zeros(scoremat.shape)
    offset = scoremat.min() - 1.
    scoremat -= offset
    if mafft: # distance as Satoh in MAFFT 
        for i,j in itertools.combinations(range(distmat.shape[0]),2):
            distmat[i,j] = distmat[j,i] = 1. - scoremat[i,j]/min(scoremat[i,i],scoremat[j,j])
    else: # distance = fraction of indels
        for i,j in itertools.combinations(range(distmat.shape[0]),2):
            distmat[i,j] = distmat[j,i] = (scoremat[i,i] + scoremat[j,j]) / scoremat[i,j] - 2. # avoids division by zero
    return distmat
    
def find_centroids_from_file_vsearch (fastafile=None, output=None, identity=0.95, nthreads=1, scratch=None):
    if fastafile is None:
        logger.error("No fasta file provided")
        return
    if output is None:
        output = "centroids." + '%012x' % random.randrange(16**12) + ".fasta.gz"
        logger.warning (f"No output file specified, writing to file {output}")
    if scratch is None:
        scratch = "."
        logger.warning (f"No scratch directory provided, writing to current directory")
    if nthreads < 1: nthreads = 1
    tmpfile = os.path.join(scratch, "tmp." + '%012x' % random.randrange(16**12) + ".fasta")
    # run vsearch and store centroids into unzipped tmpfile
    runstr = f"vsearch --cluster_fast {fastafile} --id {identity} --centroids {tmpfile} --threads {nthreads}"
    proc_run = subprocess.check_output(runstr, shell=(sys.platform!="win32"), universal_newlines=True)
    # gzip tmpfile into output file
    with open(tmpfile, 'rb') as f_in, gzip.open(output, 'wb') as f_out: f_out.writelines(f_in)
    # delete tmp file
    pathlib.Path(tmpfile).unlink()

# affinity propagation returns representatives, using similiarity matrix as input; birch needs features
def find_representatives_from_sequences_optics (sequences=None, names=None, output=None, min_samples=2, nthreads=1):
    if sequences is None:
        logger.error("No sequences provided to OPTICS")
        return
    if names is None:
        names = [f"seq{i}" for i in range(len(sequences))]
    if output is None:
        output = "representatives." + '%012x' % random.randrange(16**12) + ".fasta.gz"
        logger.warning (f"No output file specified, writing to file {output}")
    if (min_samples > len(sequences)//3): min_samples = len(sequences)//3
    if (min_samples < 2): min_samples = 2

    score_matrix = create_NW_score_matrix_parallel (sequences, nthreads=nthreads)
    distmat = score_to_distance_matrix_fraction (score_matrix, maftt=True)
    logger.info(f"Finished calculating pairwise scores; will now calculate OPTICS clusterings")
    with np.errstate(divide='ignore'): # silence OPTICS warning (https://stackoverflow.com/a/59405142/204903)
        cl = cluster.OPTICS(min_samples=min_samples, min_cluster_size=2, metric="precomputed", n_jobs=nthreads).fit(distmat)

    idx = [i for i,j in enumerate(cl.labels_) if j < 0] # all noisy points (negative labels) are representatives
    cl = [[i,j,k,l] for i,(j,k,l) in enumerate(zip(cl.labels_, cl.reachability_, names))] # we have [index, label, reachability, name]
    cl = [x for x in cl if x[1] >= 0] # excluding noisy seqs; we need one per cluster, with min reachability distance
    cl = sorted(cl, key=lambda x: (x[1], x[2])) # sort by cluster label, breaking ties with reachability
#    print ("\n".join(["\t".join(map(str,i)) for i in cl])) # DEBUG
    cl = [list(v)[0] for k,v in itertools.groupby(cl, key=lambda x: x[1])] # groupby x[1] i.e. cluster label and return first element of each group
    idx += [x[0] for x in cl] # add all cluster representatives to indx_1

    # write representatives to file
    with open_anyformat(output, "w") as f:
        for i in idx:
            f.write(str(f">{names[i]}\n{sequences[i]}\n").encode())
    logger.info(f"Wrote {len(idx)} representatives to file {output}")


def cluster_primers_with_vsearch (df, maxdiffs=8, maxgaps=10, identity=0.8, scratch = None, nthreads=1):
    '''
    still tends to oversplit (from 241k to 65k clusters with maxgaps=maxdiffs=30 and id=0.2)
    '''
    def df_from_consensus_fasta (fasfile):
        falist = read_fasta_as_list (fasfile)
        logger.info(f"Read {len(falist)} sequences from consensus file {fasfile}")
        #seqnames are  ">centroid=CGGCCTTCTCGACCGC;seqs=1"
        falist = [[re.search('centroid=(.*);seqs', s.id).group(1), s.seq] for s in falist]
        df = pd.DataFrame(falist, columns=['primer', 'consensus'])
        return df
    def df_from_profile (profile):
        lines = [l.strip() for l in open(profile, 'r').readlines() if l.strip() != ""]
        logger.info(f"Read {len(lines)} lines from profile file {profile}")
        seqname = None
        prof = {}
        sizes = {}
        for l in lines:
            if l.startswith(">"): # >centroid=AAAGTATCATAATTGACAACTTGTCCAT;seqs=35
                seqname = re.search('centroid=(.*);seqs', l).group(1)
                sizes[seqname] = int(re.search(';seqs=(\d+)', l).group(1))
                prof[seqname] = []
            else: prof[seqname].append(l)
        del (lines)
        ACGT = "ACGT"
        prolist = []
        for seqname in [k for k in prof.keys() if sizes[k] > 1]: 
            tbl = [l.split("\t") for l in prof[seqname]]
            tbl = sorted(tbl, key=lambda x: int(x[0]))
            seq = [ACGT[ np.argmax(x[2:6]) ] for x in tbl]
            ambiguous = [np.where(np.array(x[2:6])==max(x[2:6]))[0].shape[0] > 1 for x in tbl]
            seq = [x if not y else x.lower() for x,y in zip(seq, ambiguous)]
            prolist.append([seqname, "".join(seq)])
        df = pd.DataFrame(prolist, columns=['primer', 'profile'])
        logger.info(f"Finished creating profile dataframe with {len(df)} rows")
        return df

    hash_name = '%012x' % random.randrange(16**12)
    if scratch is None:
        scratch = f"scratch_vsearch.{hash_name}"
    if not os.path.exists(scratch):
        pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # create scratch subdirectory
        scratch_exists = False
    else: scratch_exists = True

    fastafile = os.path.join(scratch, f"primers.fasta")
    ucfile = os.path.join(scratch, f"primers.uc")
    consfile = os.path.join(scratch, f"consensus.fasta")
    profile = os.path.join(scratch, f"profile.txt")
    if nthreads < 1: nthreads = 1
    
    with open_anyformat(fastafile, "w") as f:
        for seq in df["primer"].tolist():
            f.write(str(f">{seq}\n{seq}\n").encode())

    # run vsearch and store centroids into unzipped tmpfile
    runstr = f"vsearch --cluster_fast {fastafile} --minseqlength 8 --maxdiffs {maxdiffs} --maxgaps {maxgaps} " + \
             f"--id {identity} --threads {nthreads} --uc {ucfile} --consout {consfile} --profile {profile}"
    proc_run = subprocess.check_output(runstr, shell=(sys.platform!="win32"), universal_newlines=True)
    # https://manpages.debian.org/stretch/vsearch/vsearch.1 : no headers, and first column is a hits (H), centroid (S), or cluster (C) info
    vclus = pd.read_csv(ucfile, sep="\t", names=["rectype","v_cluster","primer"], usecols=[0,1,8]) 
    vclus = vclus[vclus["rectype"].isin(["S","H"])] # keep only centroids and non-centroids (Hits), excluding summary "C"
    vclus.drop("rectype", axis=1, inplace=True)
    ## add (shorter) consensus info
    vcons = df_from_consensus_fasta (consfile)
    vclus = pd.merge(vclus, vcons, on="primer", how="left")
    ## add (longer consensus) profile info
    vprof = df_from_profile (profile)
    vclus = pd.merge(vclus, vprof, on="primer", how="left")

    # export consensus to all from same cluster
    vclus["consensus"] = vclus.groupby("v_cluster")["consensus"].transform(lambda x: x.mode()[0] if x.notna().any() else np.nan) 
    # export profile to all from same cluster (some clusters do not have a profile)
    vclus["profile"] = vclus.groupby("v_cluster")["profile"].transform(lambda x: x.mode()[0] if x.notna().any() else np.nan)
    vclus["profile"] = vclus["profile"].fillna(vclus["consensus"])

    # delete tmp files
    for file in [fastafile, ucfile, consfile, profile]:
        pathlib.Path(file).unlink()
    if not scratch_exists:
        shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory
    orig_cols = df.columns.tolist()
    vclus = vclus.merge(df, on="primer", how="left")
    vclus = vclus[orig_cols + ["v_cluster", "consensus", "profile"]]
    return vclus


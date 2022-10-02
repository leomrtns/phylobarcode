#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
from phylobarcode.pb_kmer import short_kmer
import pandas as pd, numpy as np
import itertools, pathlib, shutil, gzip
from Bio import pairwise2
from sklearn import cluster

# legacy code, no need to create a distinct stream
# log_format = logging.Formatter(fmt='phylobarcode_clustr %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

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
    
def cluster_primers_from_tsv (tsv = None, output = None, min_samples = 2, subsample=100, nthreads = 1):
    if tsv is None: 
        logger.error("No tsv file provided")
        return
    if output is None:
        output = "clusters." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, writing to file {output}.tsv")
    if subsample < 1e-7: subsample = 1e-7
    if subsample > 100: subsample = 100

    df = pd.read_csv (tsv, compression="infer", sep="\t", dtype='unicode')
    #df.set_index("primer", drop=False, inplace=True) # keep column with primer sequences
    logger.info(f"Read {len(df)} primers from file {tsv}")
    df = subsample_primers (df, subsample=subsample)
    primers = df["primer"].tolist()
    primers = [str(i) for i in primers]
    score_mat = create_NW_score_matrix_parallel (primers, nthreads=nthreads)
    logger.info(f"Pairwise distances calculated; will now cluster primers")
    distmat = score_to_distance_matrix_fraction (score_mat, mafft=True)
    with np.errstate(divide='ignore'): # silence OPTICS warning (https://stackoverflow.com/a/59405142/204903)
        cl = cluster.OPTICS(min_samples=min_samples, min_cluster_size=2, metric="precomputed", n_jobs=nthreads).fit(distmat)
    df["cluster"] = cl.labels_
    logger.info(f"Clustering done, writing to file {output}.tsv")
    df.to_csv (f"{output}.tsv", sep="\t", index=False)

def subsample_primers (df, subsample=100):
    if subsample > 100: return df
    logger.info (f"Subsampling {len(df)} primers to {subsample:.2f}% of the original set over frequency, distance, and penalty (from primer3)")
    subsample /= 100 # pandas percentile uses 0-1 scale
    df["frequency"] = df["frequency"].astype(int)
    df["max_distance"] = df["max_distance"].astype(int)
    df["penalty"] = df["penalty"].astype(float)
    threshold = {
            "frequency": df["frequency"].quantile(1. - subsample), # quantile goes from min to max, thus we want 100-x% smallest value
            "max_distance": df["max_distance"].quantile(subsample),
            "penalty": df["penalty"].quantile(subsample)
            }
    if threshold["frequency"] > 1: # all primers have freq at least one thus all would be chosen  
        df = df[(df["frequency"] >= threshold["frequency"]) | 
                (df["max_distance"] < threshold["max_distance"]) |
                (df["penalty"] < threshold["penalty"])]
    else:
        df = df[(df["frequency"] > threshold["frequency"]) | # OR large frequency OR freq=1 but small distance or penalty
                ((df["max_distance"] < threshold["max_distance"]) &
                 (df["penalty"] < threshold["penalty"]))]
    logger.info (f"Subsampling done, {len(df)} primers kept using thresholds {threshold}")
    df = df.sort_values(by=["frequency","max_distance","penalty"], ascending=[False,True,True])
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
        for i in range(size): 
            scoremat[i,i] = parasail.sg_striped_16(str(seqlist[i]), str(seqlist[i]), 9,1, parasail.blosum30).score # nw (sg doenst penalise begin and end )
        for i,j in itertools.combinations(range(size),2): #parasail._stats_ also gives x.length, x.score
            #scoremat[i,j] = scoremat[j,i] = parasail.nw_stats_striped_16(str(seqlist[i]), str(seqlist[j]), 11,1, parasail.blosum30).matches
            scoremat[i,j] = scoremat[j,i] = parasail.sg_striped_16(str(seqlist[i]), str(seqlist[j]), 9,1, parasail.blosum30).score
    elif use_parasail is True and isinstance(band_size, int): # banded: not semi-global but "full" NW, with simple penalty matrix
        mymat = parasail.matrix_create("ACGT", 2, -1)
        for i in range(size):
            scoremat[i,i] = parasail.nw_banded(str(seqlist[i]), str(seqlist[i]), 8, 1, band_size, mymat).score # global Needleman-Wunsch 
        for i,j in itertools.combinations(range(size),2): #parasail._stats_ also gives x.length, x.score
            scoremat[i,j] = scoremat[j,i] = parasail.nw_banded(str(seqlist[i]), str(seqlist[j]), 8, 1, band_size, mymat).score
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
    logger.info(f"Calculating pairwise alignment scores for {size} sequences using {nthreads} threads")
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
        for i, j in pairlist:
            d_ij = parasail.sg_striped_16(seqlist[i], seqlist[j], 9,1, parasail.blosum30).score
            scorelist.append(d_ij)
    elif use_parasail is True and isinstance(band_size, int): # banded: not semi-global but "full" NW, with simple penalty matrix
        import parasail
        mymat = parasail.matrix_create("ACGT", 2, -1)
        for i, j in pairlist:
            d_ij = parasail.nw_banded(seqlist[i], seqlist[j], 8, 1, band_size, mymat).score
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
    
def find_centroids_from_file_vsearch (fastafile=None, output=None, identity=0.95, nthreads=0, scratch=None):
    if fastafile is None:
        logger.error("No fasta file provided")
        return
    if output is None:
        output = "centroids." + '%012x' % random.randrange(16**12) + ".fasta.gz"
        logger.warning (f"No output file specified, writing to file {output}")
    if scratch is None:
        scratch = "."
        logger.warning (f"No scratch directory provided, writing to current directory")
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


#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
from phylobarcode.pb_kmer import * 
import pandas as pd, numpy as np
import itertools, pathlib, shutil, gzip, parasail
from Bio import pairwise2
from sklearn import cluster

## TODO: fix task parameters

# legacy logging code, no need to create a distinct stream
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
    
def cluster_primers_from_tsv (tsv = None, output = None, min_samples = 10, subsample=100, kmer_length = 5, 
        threshold = 0.7, scratch = None, nthreads=1):
    if tsv is None: 
        logger.error("No tsv file provided")
        return
    if output is None:
        output = "clusters." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, writing to file {output}.tsv.xz")
    if scratch is None:
        scratch = "scratch." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No scratch directory specified, writing to directory {scratch}")
    if min_samples < 3: min_samples = 3
    if subsample < 1e-3: subsample = 1e-3
    if subsample > 100: subsample = 100
    if kmer_length < 3: kmer_length = 3
    if kmer_length > 18: kmer_length = 18
    if threshold < 0.01: threshold = 0.01
    if threshold > 1: threshold = 1
    if nthreads < 1: nthreads = 1

    df = pd.read_csv (tsv, compression="infer", sep="\t", dtype='unicode')
    logger.info(f"Read {len(df)} primers from file {tsv}")
    if "cluster" in df.columns:
        logger.warning (f"Column 'cluster' already exists in file {tsv}, overwriting")
        df.drop(columns=["cluster"], inplace=True)
    
    df = subsample_primers (df, subsample=subsample)
    df = cluster_primers_with_vsearch (df, nthreads = nthreads, scratch=scratch)
    df = cluster_profiles_with_vsearch (df, nthreads = nthreads, simple_cluster_id = True, scratch=scratch)

    logger.info (f"Found {len(df['v_cluster'].unique())} clusters with vsearch; Will now cluster profiles at similarity {threshold}")
    df = merge_vsearch_profiles (df, identity = threshold, nthreads = nthreads)
    df = reoder_dataframe_by_clusters (df)
    logger.info (f"Found {len(df['cluster'].unique())} clusters of vsearch profiles, writing to file {output}.tsv.xz")
    df.to_csv (f"{output}.tsv.xz", sep="\t", index=False)
    return

def reorder_dataframe_by_clusters (df):
    # python>=3.7 ensures order of dict keys is preserved; 
    # some samples don't have GTDB info (thus no genus) but all have taxon, which comes from GFF file
    sort_dic = {"genus_diversity":False,"taxon_diversity":False,"penalty":True,"frequency":False,"max_distance":True}
    col_sort = [i for i in sort_dic.keys() if i in df.columns]
    ord_sort = [sort_dic[i] for i in col_sort]
    df = df.sort_values(by=col_sort, ascending=ord_sort)
    df["clust_idx"] = df.groupby("cluster").cumcount() # idx = 1 best, idx = 2 second best, etc.
    col_sort.insert(0,"clust_idx")
    ord_sort.insert(0,True)
    df = df.sort_values(by=col_sort, ascending=ord_sort)
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

def cluster_primers_with_vsearch (df, maxdiffs=8, maxgaps=10, identity=0.6, scratch = None, nthreads=1):
    orig_cols = df.columns.tolist()
    seqs = df["primer"].tolist()
    logger.info(f"Clustering {len(seqs)} primers with vsearch")

    vclus = vsearch_clustering_sequences (seqs, seqs, maxdiffs=maxdiffs, maxgaps=maxgaps, identity=identity, scratch=scratch, 
            drop_profile = False, nthreads=nthreads)
    vclus.rename(columns={"seqname":"primer"}, inplace=True) # "key" 

    for col in ["v_cluster", "profile", "consensus"]: # these columns should be new from vclus
        if col in orig_cols: df.drop(columns=[col], inplace=True)
    df = df.merge(vclus, on="primer", how="left")
    return df[orig_cols + ["v_cluster", "profile"]] ## maintain original order of columns, excluding "consensus"

def cluster_profiles_with_vsearch (df, maxdiffs=8, maxgaps=10, identity=0.6, scratch = None, simple_cluster_id = True, nthreads=1):
    if identity < 0.1: identity = 0.1
    if identity > 0.9999: identity = 0.9999
    orig_cols = df.columns.tolist()
    profiles = df["profile"].unique().tolist()
    logger.info(f"Clustering {len(profiles)} profiles with vsearch")
    sequences = [re.sub("[a-z]", "N", str(x)) for x in profiles] # convert lowercase to N (ambiguous)
    # vclus will have 4 columns: profile, consensus, seqname, v_cluster; we only need seqname (which is original profile) and v_cluster
    vclus = vsearch_clustering_sequences (sequences, profiles, maxdiffs=maxdiffs, maxgaps=maxgaps, identity=identity, scratch=scratch, 
            drop_profile = True, nthreads=nthreads)
    vclus = vclus[["seqname", "v_cluster"]] # drop profile and consensus columns
    vclus.rename(columns={"seqname":"profile"}, inplace=True) # now we have "profile" and "v_cluster"
    # now we merge vclust with original df
    if "v_cluster" in orig_cols: df = df.rename(columns={"v_cluster":"v_cluster_prev"})
    df = df.merge(vclus, on="profile", how="left")
    if "v_cluster" in orig_cols:
        if simple_cluster_id is False: # cluster id will be newclust_oldclust; o.w. just newclust
            df["v_cluster"] = df[["v_cluster", "v_cluster_prev"]].agg("_".join, axis=1) 
        else:
            df["v_cluster"] = df["v_cluster"].fillna(df["v_cluster_prev"] + "_") # NaNs should not happen?
        df.drop("v_cluster_prev", axis=1, inplace=True)
    else:
        orig_cols += ["v_cluster"]
    orig_cols = [x for x in orig_cols if x not in ["consensus"]] # remove consensus but keep new profiles 
    return df[orig_cols] ## maintain original order of columns

def vsearch_clustering_sequences (seqs, seqnames, maxdiffs=8, maxgaps=10, identity=0.8, scratch = None, 
        drop_profile=False, nthreads=1):
    '''
    still tends to oversplit (from 241k to 65k clusters with maxgaps=maxdiffs=30 and id=0.2)
    '''
    def df_from_consensus_fasta (fasfile):
        falist = read_fasta_as_list (fasfile)
        logger.info(f"Read {len(falist)} sequences from consensus file {fasfile}")
        #seqnames are  ">centroid=CGGCCTTCTCGACCGC;seqs=1"
        falist = [[re.search('centroid=(.*);seqs', s.id).group(1), str(s.seq)] for s in falist]
        df = pd.DataFrame(falist, columns=['seqname', 'consensus'])
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
            prolist.append([str(seqname), str("".join(seq))])
        df = pd.DataFrame(prolist, columns=['seqname', 'profile'])
        logger.info(f"Finished dataframe with {len(df)} longer profiles (other rows will receive consensus sequence)")
        return df

    hash_name = '%012x' % random.randrange(16**12)
    if scratch is None:
        scratch = f"scratch_vsearch.{hash_name}"
    if not os.path.exists(scratch):
        pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # create scratch subdirectory
        scratch_exists = False
    else: scratch_exists = True

    fastafile = os.path.join(scratch, f"vseqs.fasta")
    ucfile = os.path.join(scratch, f"vseqs.uc")
    consfile = os.path.join(scratch, f"v_consensus.fasta")
    profile = os.path.join(scratch, f"v_profile.txt")
    if nthreads < 1: nthreads = 1
    
    with open_anyformat(fastafile, "w") as f:
        for seqname, seq in zip(seqnames, seqs):
            f.write(str(f">{seqname}\n{seq}\n").encode())

    # run vsearch and store centroids into unzipped tmpfile
    runstr = f"vsearch --cluster_fast {fastafile} --minseqlength 8 --maxdiffs {maxdiffs} --maxgaps {maxgaps} " + \
             f"--id {identity} --threads {nthreads} --uc {ucfile} --consout {consfile} --profile {profile}"
    proc_run = subprocess.check_output(runstr, shell=(sys.platform!="win32"), universal_newlines=True)
    # https://manpages.debian.org/stretch/vsearch/vsearch.1 : no headers, and first column is a hits (H), centroid (S), or cluster (C) info
    vclus = pd.read_csv(ucfile, sep="\t", names=["rectype","v_cluster","seqname"], usecols=[0,1,8]) 
    vclus = vclus[vclus["rectype"].isin(["S","H"])] # keep only centroids and non-centroids (Hits), excluding summary "C"
    vclus.drop("rectype", axis=1, inplace=True)
    vclus["v_cluster"] = vclus["v_cluster"].astype("string")
    ## add (shorter) consensus info (fast to calculate and needed even in absence of profile)
    vcons = df_from_consensus_fasta (consfile)
    vclus = pd.merge(vclus, vcons, on="seqname", how="left")
    # export consensus to all from same cluster
    vclus["consensus"] = vclus.groupby("v_cluster")["consensus"].transform(lambda x: x.mode()[0] if x.notna().any() else np.nan) 

    ## add (longer consensus) profile info
    if drop_profile is False:
        vprof = df_from_profile (profile)
        vclus = pd.merge(vclus, vprof, on="seqname", how="left")
        # export profile to all from same cluster (some clusters do not have a profile)
        vclus["profile"] = vclus.groupby("v_cluster")["profile"].transform(lambda x: x.mode()[0] if x.notna().any() else np.nan)
        vclus["profile"] = vclus["profile"].fillna(vclus["consensus"])
        vclus.drop("consensus", axis=1, inplace=True)
    else:
        vclus.rename(columns={"consensus": "profile"}, inplace=True)

    # delete tmp files
    for file in [fastafile, ucfile, consfile, profile]:
        pathlib.Path(file).unlink()
    if not scratch_exists:
        shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory
    return vclus ## will have 3 columns: seqname, v_cluster, profile

def pairwise_identity_matches (s1, s2, mode = "cdhit"):
    '''
    identity based on number of matches (https://manpages.debian.org/stretch/vsearch/vsearch.1)
    '''
    sim = parasail.sg_stats_striped_16(s1, s2, 10, 1, parasail.pam10)
    if mode == "cdhit":  x = sim.matches / min(sim.len_query, sim.len_ref) # CD-HIT similarity (always larger than others)
    elif mode == "edit": x = sim.matches / sim.length # edit distance (larger than MBL but lower than CD-HIT)
    else:                x = (sim.length - sim.matches) / max (sim.len_query, sim.len_ref) # MBL: all gaps are counted (more strict)
    return x

def merge_vsearch_profiles (df, identity=0.8, nthreads=1):
    profiles = df["profile"].unique().tolist()
    logger.info(f"Starting to merge {len(profiles)} profiles")
    sequences = [re.sub("[a-z]", "N", str(x)) for x in profiles] # convert lowercase to N (ambiguous)

    idx = cluster_profiles_parallel (sequences, identity=identity, nthreads=nthreads)
    df1 = pd.DataFrame({"cluster": idx, "profile": profiles}, dtype="string")
    logger.info(f"Finished merging {len(set(idx))} clusters")
    o_cols = df.columns.tolist()
    df = df.merge(df1, on="profile", how="left") ## all samples should have profile even if repeated
    df["cluster"] = df[["cluster", "v_cluster"]].agg("_".join, axis=1) # add vsearch cluster
    df = df[o_cols + ["cluster"]] ## force order
    df.drop("v_cluster", axis=1, inplace=True)
    return df

def cluster_profiles (sequences, identity = 0.8):
    n_seqs = len (sequences)
    clusters = []
    for i in range (n_seqs):
        for cl in clusters:
            if pairwise_identity_matches (sequences[i], sequences[cl[0]]) >= identity:
                # first position will be the longest sequence
                if len(sequence[i]) > len(sequences[cl[0]]): cl.insert(0, i) 
                else: cl.append (i)
                break
        else:
            clusters.append ([i])
    return map_clusters_to_indices (clusters, n_elements = n_seqs)
    #return pd.DataFrame ({'cluster': idx, 'profile': sequences})

def cluster_pair_of_profile_clusters (cluster_pair, sequences, identity = 0.8):
    for c1 in cluster_pair[1]:
        for c0 in cluster_pair[0]:
            if pairwise_identity_matches (sequences[c1[0]], sequences[c0[0]]) >= identity:
                if len(sequences[c1[0]]) > len(sequences[c0[0]]):
                    swap = c1[0]; c1[0] = c0[0]; c0[0] = swap
                c0.extend (c1) ## now c1[0] is the longest sequence
                break
        else:
            cluster_pair[0].append (c1)
    return cluster_pair[0]

def cluster_profiles_parallel (sequences, identity = 0.8, nthreads = 1):
    if nthreads < 1: nthreads = 1
    if nthreads == 1:
        return cluster_profiles_as_df (sequences, identity = identity)

    from multiprocessing import Pool
    from functools import partial
    n_seqs = len (sequences)
    cluster_chunks = [ [[j] for j in range (i, n_seqs, 2 * nthreads)] for i in range (2 * nthreads) ]
    with Pool (nthreads) as p:
        while (len (cluster_chunks) > 1):
            nt = len(cluster_chunks)//2 + len(cluster_chunks)%2
            n_clusters = sum ([len(i) for i in cluster_chunks])
            logger.info (f"Clustering {n_clusters} clusters using {nt} threads (pool size = {len(cluster_chunks)})")
            results = p.map (partial (cluster_pair_of_profile_clusters, sequences = sequences, identity = identity),
                zip(cluster_chunks, cluster_chunks[nt:]))
            if len(cluster_chunks) % 2 == 1: results.append (cluster_chunks[nt-1])
            cluster_chunks = results
    return map_clusters_to_indices (cluster_chunks[0], n_elements = n_seqs)


###    for cluster_flanks_from_fasta

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

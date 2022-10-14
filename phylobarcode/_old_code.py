### from task_cluster_primers.py

def cluster_primers_from_tsv_old (tsv = None, output = None, min_samples = 2, subsample=100, nthreads = 1):
    if tsv is None: 
        logger.error("No tsv file provided")
        return
    if output is None:
        output = "clusters." + '%012x' % random.randrange(16**12) 
        logger.warning (f"No output file specified, writing to file {output}.tsv")
    if subsample < 1e-5: subsample = 1e-5
    if subsample > 100: subsample = 100

    df = pd.read_csv (tsv, compression="infer", sep="\t", dtype='unicode')
    #df.set_index("primer", drop=False, inplace=True) # keep column with primer sequences
    logger.info(f"Read {len(df)} primers from file {tsv}")
    
    df = kmer_clustering_dataframe (df, nthreads=nthreads)
    logger.info(f"Clustering done, writing to file {output}.tsv")
    df.to_csv (f"{output}.tsv", sep="\t", index=False)

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

def kmer_clustering_dataframe_old (df):
    primers = df["primer"].tolist()
    df["frequency"] = df["frequency"].astype(int)
    df["max_distance"] = df["max_distance"].astype(int)
    df["penalty"] = df["penalty"].astype(float)
    df = df.sort_values(by=["frequency","max_distance","penalty"], ascending=[False,True,True])

    primers = [str(i) for i in primers]
    kmers = [short_kmer(i,length = [4,7]) for i in primers]
    logger.debug (f"Number of kmers: {len(kmers)}")
    cluster_1 = cluster_short_kmers (kmers, threshold=0.75, element="min", use_centroid=False, jaccard=True)
    n_clusters = len(set(cluster_1))
    logger.debug (f"Number of clusters: {n_clusters}")
    if n_clusters < len(primers)/100:
        cluster_2 = cluster_short_kmers (kmers, threshold=0.9, element="max", use_centroid=False, jaccard=False)
        n_clusters2 = len(set(cluster_2))
        logger.debug (f"Number of overlap clusters: {n_clusters2} for too few clusters")
        cluster_1 = consensus_clustering (cluster_1, cluster_2) # oversplits clusters
        n_clusters2 = len(set(cluster_1))
        logger.debug (f"Number of consensus clusters: {n_clusters2}")
    elif n_clusters > len(primers)/2:
        cluster_1 = cluster_short_kmers (kmers, threshold=0.6, use_centroid=True, element="min", jaccard=False)
        n_clusters = len(set(cluster_1))
        logger.debug (f"Number of overlap clusters: {n_clusters} for too many clusters")
    
    df["kmer_cluster"] = cluster_1
    df = df.groupby("kmer_cluster").head(5)
    return df

### old version of kmer clustering, has several kmer sizes at once

import copy

class short_kmer:
    kmers = {}
    seq = None 
    def __init__(self, seq=None, length=None):
        self.seq = str(seq)
        self.kmers = {}
        if length is None: length = [4]
        if isinstance(length, int):
            length = [length]
        for l in length:
            self.set_kmers (l)

    def set_kmers (self, length=None):
        if length is None: length = 4
        if length < 2: 
            self.kmers[length] = 1
            return
        kmer_list = []
        for i in range(len(self.seq)-length+1):
            kmer_list.append ( xxhash.xxh64_intdigest(self.seq[i:i+length]) )
#            print ("dbg0: ", self.seq[i:i+length], " :: ", xxhash.xxh64_intdigest(self.seq[i:i+length]))
        self.kmers[length] = set (kmer_list)
#        print (sorted(self.kmers[length]), "dbg1: ", length)
        return

    def get_kmers (self, length=None):
        if length is None: length = 4
#        print ("dbg2: ", self.kmers.keys(), "requested length: ", length, " :: ", self.seq, " :: ", sorted(self.kmers[length]))
        if length not in self.kmers.keys(): self.set_kmers (length)
        return self.kmers[length]

    def append (self, other): 
        if len(other.seq) > len(self.seq): self.seq = other.seq
        for k,v in other.kmers.items():
            if k in self.kmers.keys(): self.kmers[k].update(v)
            else: self.kmers[k] = v

    def get_seq (self): return str(self.seq)
    def get_k (self): return list(self.kmers.keys())

    def jaccard_similarity (self, other, length = None):
        if length is None: length = sorted(set(self.get_k() + other.get_k())) # sorted() returns a list
        if isinstance(length, int): length = [length]
        jac = [len(self.get_kmers(l).intersection(other.get_kmers(l))) / len(self.get_kmers(l).union(other.get_kmers(l))) for l in length]
#        print ("dbg3: ", self.get_seq(), " - ", other.get_seq(), " :")
#        for i,l in enumerate(length):
#            print ("\t", l, " [", jac[i], "]\n", sorted(self.get_kmers(l)), "\n", sorted(other.get_kmers(l)))
        return min(jac), sum(jac)/len(jac), max(jac)

    def set_overlap (self, other, length = None): # aka overlap coefficient or Szymkiewicz–Simpson coefficient
        if length is None: length = sorted(set(self.get_k() + other.get_k())) # sorted() returns a list
        if isinstance(length, int): length = [length]
        ove = [len(self.get_kmers(l).intersection(other.get_kmers(l))) / min(len(self.get_kmers(l)), len(other.get_kmers(l))) for l in length]
        return min(ove), sum(ove)/len(ove), max(ove)

def cluster_short_kmers (sequences, length=None, threshold=0.5, element="min", use_centroid=True, jaccard=True):
    if not isinstance (sequences, list): sequences = [sequences]
    if length is None:
        if isinstance (sequences[0], short_kmer):length = sequences[0].get_k() 
        else: length = [5] 
    if isinstance(length, int): length = [length]
    kmers = [short_kmer(seq, length) if isinstance(seq, str) else seq for seq in sequences]
#    for kmer in kmers: print (sorted(kmer.get_kmers(length[0])), "==" , kmer.get_seq(), length[0], " dbg4")
    if (element == "min"): element = 0
    elif (element == "max"): element = 2
    else: element = 1
    clusters = []
    centroids = []
    for i in range(len(kmers)):
        for cl, ce in zip(clusters, centroids):
            if (jaccard):
                if use_centroid: similarity = kmers[i].jaccard_similarity(ce)[element]
                else:            similarity = kmers[i].jaccard_similarity(kmers[cl[0]])[element]
            else:
                if use_centroid: similarity = kmers[i].set_overlap(ce)[element]
                else:            similarity = kmers[i].set_overlap(kmers[cl[0]])[element]
            if similarity > threshold:
                cl.append(i)
                ce.append(kmers[i]) 
                break
        else: 
            clusters.append([i])
            centroids.append(copy.deepcopy(kmers[i])) 
    idx = [None] * len(kmers)
    for i, c in enumerate(clusters):
        for j in c: idx[j] = i
    return idx

def consensus_clustering (cluster_1, cluster_2):
    clusters = [] 
    for i in range(len(cluster_1)):
        for c in clusters:
            if cluster_1[i] == cluster_1[c[0]] and cluster_2[i] == cluster_2[c[0]]:
                c.append(i)
                break
        else: clusters.append([i])
    idx = [None] * len(cluster_1)
    for i, c in enumerate(clusters):
        for j in c: idx[j] = i
    return idx


### kmer preclustering + optics (before using vsearch)

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
    if nthreads < 1: nthreads = 1

    df = pd.read_csv (tsv, compression="infer", sep="\t", dtype='unicode')
    logger.info(f"Read {len(df)} primers from file {tsv}")
    if "cluster" in df.columns:
        logger.warning (f"Column 'cluster' already exists in file {tsv}, overwriting")
        df.drop(columns=["cluster"], inplace=True)
    
    df = subsample_primers (df, subsample=subsample)

    df = cluster_primers_with_vsearch (df, nthreads = nthreads)
    logger.info (f"Found {len(df['v_cluster'].unique())} clusters with vsearch; Will now cluster at similarity {threshold}")
    df = merge_vsearch_profiles (df, identity = threshold, nthreads = nthreads)
    logger.info (f"Found {len(df['cluster'].unique())} clusters of vsearch profiles, writing to file {output}.tsv.xz")
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


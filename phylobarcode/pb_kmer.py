#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import copy

logger = logging.getLogger("phylobarcode_global_logger")

class single_kmer:
    kmers = set()
    seq = None 
    k = None
    size = 0
    def __init__(self, seq=None, k=None):
        if isinstance(seq, str): 
            self.seq = str(seq)
            self.kmers = set()
            if k is None: self.k = 4
            else:         self.k = k
            self.set_kmers()
            self.size = len(self.kmers)
        elif isinstance(seq, single_kmer):
            self.seq = copy.deepcopy(seq.seq)
            self.kmers = copy.deepcopy(seq.kmers)
            self.k = seq.k
            self.size = seq.size

    def set_kmers (self):
        kmer_list = []
        for i in range(len(self.seq)-self.k + 1):
            kmer_list.append ( xxhash.xxh64_intdigest(self.seq[i:i+self.k]) )
        self.kmers = set (kmer_list)
        return

    def get_kmers (self): return self.kmers

    def get_seq (self): return self.seq

    def get_k (self): return self.k 

    def get_size (self): return self.size

    def append (self, other): 
        if len(self.seq) > len(other.get_seq()): 
            self.seq = other.get_seq()
        self.kmers.update(other.get_kmers())
        self.size = len(self.kmers)

    def similarity (self, other):  # Jaccard similarity, overlap coefficient (aka Szymkiewicz–Simpson coefficient)
        intrsc = len(self.get_kmers().intersection(other.get_kmers()))
        unin   = len(self.get_kmers().union(other.get_kmers()))
        minsize = min(self.get_size(), other.get_size())
        return intrsc/unin, intrsc/minsize

def cluster_single_kmers (sequences, length=None, threshold=0.5, use_centroid=True, jaccard=True):
    if not isinstance (sequences, list): sequences = [sequences]
    if length is None:
        if isinstance (sequences[0], single_kmer):length = sequences[0].get_k() 
        else: length = 5 
    kmers = [single_kmer(seq, k=length) if isinstance(seq, str) else seq for seq in sequences]
    if jaccard is True: element = 0 ## which element from `similarity` to use
    else:               element = 1
    clusters = []
    centroids = []
    for i in range(len(kmers)):
        for cl, ce in zip(clusters, centroids):
            if use_centroid: similarity = kmers[i].similarity(ce)
            else:            similarity = kmers[i].similarity(kmers[cl[0]])
            if similarity[element] > threshold:
                cl.append(i)
                ce.append(kmers[i]) #single_kmer.append() 
                break
        else: 
            clusters.append([i])
            centroids.append(single_kmer(kmers[i])) #list.append()

    idx = map_clusters_to_indices(clusters, n_elements=len(kmers))
    return idx, clusters, centroids

def cluster_centroids (centroids, clusters, cent_2, clu_2, threshold=0.5, jaccard=True):
    if jaccard is True: element = 0 ## which element from `similarity` to use
    else:               element = 1
    if len(cent_2) != len(clu_2):
        logger.error("Length of centroids and clusters do not match")
        return
    for i in range(len(cent_2)):
        for cl,ce in zip(clusters, centroids):
            if ce.similarity(cent_2[i])[element] > threshold:
                cl.extend(clu_2[i]) # each element of cluster[] is a list (of indices to original sequences)
                ce.append(cent_2[i])
                break
        else:
            clusters.append(clu_2[i])
            centroids.append(single_kmer(cent_2[i]))
    return clusters, centroids

# TODO: use cluster_centroids() to run in parallel

def map_clusters_to_indices (clusters, n_elements = None):
    if n_elements is None: n_elements = max([max(c) for c in clusters]) + 1 
    idx = [None] * n_elements
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

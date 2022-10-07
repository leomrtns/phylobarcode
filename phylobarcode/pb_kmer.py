#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?

logger = logging.getLogger("phylobarcode_global_logger")

class short_kmer:
    kmers = {}
    seq = None 
    def __init__(self, seq=None, length=None):
        self.seq = str(seq)
        if length is None: length = [5]
        if isinstance(length, int):
            length = [length]
        for l in length:
            self.set_kmers (l)

    def set_kmers (self, length=None):
        if length is None: length = 5
        if length < 2: 
            self.kmers[length] = 1
            return
        kmer_list = []
        for i in range(len(self.seq)-length+1):
            kmer_list.append ( xxhash.xxh64_intdigest(self.seq[i:i+length]) )
            print (self.seq[i:i+length], " :: ", xxhash.xxh64_intdigest(self.seq[i:i+length]))
        self.kmers[length] = set (kmer_list)
        print (sorted(self.kmers[length]), "dbg0: ", length)
        return

    def get_kmers (self, length=None):
        if length is None: length = 5
        print ("dbg1: ", self.kmers.keys(), "requested length: ", length, " :: ", self.seq)
        if length not in self.kmers.keys(): self.set_kmers (length)
        print ("after:", sorted(self.kmers[length]))
        return self.kmers[length]

    def get_seq (self):
        return str(self.seq)

    def jaccard_similarity (self, other, length = None):
        if length is None:
            length = list(set(list(self.kmers.keys()) + list(other.kmers.keys())))
        if isinstance(length, int):
            length = [length]
        jac = [len(self.get_kmers(l).intersection(other.get_kmers(l))) / len(self.get_kmers(l).union(other.get_kmers(l))) for l in length]
        print (self.get_seq(), " - ", other.get_seq(), " :")
        for i,l in enumerate(length):
            print ("\t", l, " [", jac[i], "]\n", sorted(self.get_kmers(l)), "\n", sorted(other.get_kmers(l)))
        return min(jac), sum(jac)/len(jac), max(jac)

    def set_overlap (self, other, length = None): # aka overlap coefficient or Szymkiewicz–Simpson coefficient
        if length is None:
            length = list(set(list(self.kmers.keys()) + list(other.kmers.keys())))
        if isinstance(length, int):
            length = [length]
        ove = [len(self.get_kmers(l).intersection(other.get_kmers(l))) / min(len(self.get_kmers(l)), len(other.get_kmers(l))) for l in length]
        return min(ove), sum(ove)/len(ove), max(ove)

def cluster_short_kmers_by_similarity (sequences, length=None, threshold=0.5, element="min"):
    if length is None: length = [5]
    if isinstance(length, int): length = [length]
    if not isinstance (sequences, list): sequences = [sequences]
    kmers = [short_kmer(seq, length) if isinstance(seq, str) else seq for seq in sequences]
    for kmer in kmers:
        print (sorted(kmer.get_kmers(length[0])), "==" , kmer.get_seq(), length[0])
    if (element == "min"): element = 0
    elif (element == "max"): element = 2
    else: element = 1
    clusters = []
    for i in range(len(kmers)):
        for c in clusters:
            if kmers[i].jaccard_similarity(kmers[c[0]])[element] > threshold:
                c.append(i)
                break
        else: clusters.append([i])
    idx = [None] * len(kmers)
    for i, c in enumerate(clusters):
        for j in c: idx[j] = i
    return idx

def cluster_short_kmers_by_overlap (sequences, length=5, threshold=0.5, element="min"):
    if isinstance(length, int): length = [length]
    if not isinstance (sequences, list): sequences = [sequences]
    kmers = [short_kmer(seq, length) if isinstance(seq, str) else seq for seq in sequences]
    if (element == "min"): element = 0
    elif (element == "max"): element = 2
    else: element = 1
    clusters = []
    for i in range(len(kmers)):
        for c in clusters:
            if kmers[i].set_overlap(kmers[c[0]])[element] > threshold:
                c.append(i)
                break
        else: clusters.append([i])
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

#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?

logger = logging.getLogger("phylobarcode_global_logger")

class short_kmer:
    kmers = {}
    seq = None 
    def __init__(self, seq=None, length=5):
        self.seq = str(seq)
        if isinstance(length, int):
            length = [length]
        for l in length:
            self.set_kmers (l)

    def set_kmers (self, length):
        if length < 2: 
            self.kmers[length] = 1
            return
        kmer_list = []
        for i in range(len(self.seq)-length+1):
            kmer_list.append ( xxhash.xxh128_intdigest(self.seq[i:i+length]) )
        self.kmers[length] = set (kmer_list)
        return

    def get_kmers (self, length):
        if length in self.kmers:
            return self.kmers[length]
        else:
            self.set_kmers (length)
            return self.kmers[length]

    def jaccard_similarity (self, other, length = None):
        if length is None:
            length = list(set(self.kmers.keys() + other.kmers.keys()))
        if isinstance(length, int):
            length = [length]
        jac = [self.get_kmers(l).intersection(other.get_kmers(l)) / len(self.get_kmers(l).union(other.get_kmers(l))) for l in length]
        return min(jac), sum(jac)/len(jac), max(jac)

    def set_overlap (self, other, length = None): # aka overlap coefficient or Szymkiewicz–Simpson coefficient
        if length is None:
            length = list(set(self.kmers.keys() + other.kmers.keys()))
        if isinstance(length, int):
            length = [length]
        ove = [self.get_kmers(l).intersection(other.get_kmers(l)) / min(len(self.get_kmers(l)), len(other.get_kmers(l))) for l in length]
        return min(ove), sum(ove)/len(ove), max(ove)

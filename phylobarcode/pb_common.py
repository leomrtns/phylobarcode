import os, logging, xxhash
from Bio import Seq, SeqIO
import random, datetime, sys, re, glob, collections, subprocess, itertools, pathlib, base64, string
import lzma, gzip, bz2

# legacy code, now every module shares the same parent logger
#log_format = logging.Formatter(fmt='phylobarcode_common %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

base62 = string.digits + string.ascii_letters + '_.-+~@'  # 66 elements actually (62 is alphanum only)
len_base62 = len (base62)

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

def read_fasta_as_list (filename, clean_sequence=True):
    unaligned = []
    with open_anyformat (filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if clean_sequence:
                record.seq  = Seq.Seq(str(record.seq.upper()).replace(".","N"))
            unaligned.append(record)
    logger.debug("Read %s sequences from file %s", str(len(unaligned)), filename)
    return unaligned

def read_fasta_headers_as_list (filename):
    seqnames = []
    with open_anyformat (filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqnames.append(record.id)
    logger.debug("Read %s sequence headers from file %s", str(len(seqnames)), filename)
    return seqnames

def mafft_align_seqs (sequences=None, infile = None, outfile = None, reference_file = None, prefix = "/tmp/"):    # list not dict
    if (sequences is None) and (infile is None):
        print ("ERROR: You must give me a fasta object or a file")
    if prefix is None: prefix = "./"
    if infile is None: ifl = f"{prefix}/mafft.fasta"
    else: ifl = infile # if both infile and sequences are present, it will save (overwrite) infile
    if outfile is None: ofl = f"{prefix}/mafft.aln"
    else: ofl = outfile # in this case it will not exclude_reference
    if sequences: SeqIO.write(sequences, ifl, "fasta") ## else it should be present in infile

    runstr = f"mafft --auto --thread -1 {ifl} > {ofl}"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    aligned = AlignIO.read(ofl, "fasta")

    if infile is None:  os.system("rm -f " + ifl)
    if outfile is None: os.system("rm -f " + ofl)
    return aligned

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


#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import pandas as pd, numpy as np
import itertools, io
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='phylobarcode__blast %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

database_fallback = "/media/teradisk/03.n114312_bigdata/blast/ref_prok_rep_genomes"

def blast_primers_from_csv (csv=None, output=None, database = None, evalue=0.01, task="blastn", max_target_seqs = 1000, nthreads=1):
    if csv is None: 
        logger.error("No csv file provided")
        return
    if len(csv) !=2:
        logger.error("I need two csv files, one with left primers and one for right primers")
        return
    if output is None or len(output) != 2: ## this should not happen if function called from main script
        prefix = "blastprimers." + '%012x' % random.randrange(16**12) 
        output = [prefix + "_left", prefix + "_right"]
        logger.warning (f"No output file specified, using {output[0]} and {output[1]}")
    if database is None:
        logger.error("I need a database (full path with prefix of DB in blast format; usually filename without '.nal' or '.nsq')")
        return
    if evalue > 1000: evalue = 1000
    if evalue < 1e-9: evalue = 1e-9
    if max_target_seqs > 1e9: max_target_seqs = 1e9
    if max_target_seqs < 1: max_target_seqs = 1
    if task not in ["blastn-short", "blastn"]:
        logger.error("I need a task = 'blastn-short' or 'blastn'; I'll use 'blastn'")

    df_l = pd.read_csv (csv[0], compression="infer", sep=",", dtype='unicode')
    df_r = pd.read_csv (csv[1], compression="infer", sep=",", dtype='unicode')
    primers_l = df_l["primer"].tolist()
    primers_r = df_r["primer"].tolist()
    logger.info(f"Read {len(primers_l)} primers from file {csv[0]} and {len(primers_r)} from file {csv[1]};")
    blast_l = query_primers_blastn_on_database (primers_l, database=database, evalue=evalue, task=task,
            max_target_seqs=max_target_seqs, nthreads=nthreads)
    print (blast_l)


def query_string_blastn_on_database (seqstring, database=None, evalue = 0.01, task="blastn", max_target_seqs=1000, nthreads=1):
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

def query_primers_blastn_on_database (primer_list, database=None, evalue = 0.01, task="blastn-short", max_target_seqs=1000, nthreads=1):
    """ 
    task = "blastn-short" (more precise) or task = "blastn" 
    """
    if database is None:
        database = database_fallback
        logger.error(f"I need a database (full path with prefix of DB in blast format; usually filename without '.nal' \
                or '.nsq'); I will try to use {database} but there are no guarantees it will work")
    cline = f"blastn -query - -db {database} -evalue {evalue} -task \"{task}\" -word_size 7 -out - \
              -outfmt \"6 std qseq sseq\"  -max_target_seqs {max_target_seqs} -num_threads {nthreads}" ## format 6 plus aligned querey plus aligned subject
    child = subprocess.Popen(cline, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                            universal_newlines=True, shell=(sys.platform!="win32"))
    records = [SeqRecord(Seq.Seq(seq), id=seq) for seq in primer_list] ## head is primer itself (like "> ACGTGT") 
    SeqIO.write(records, child.stdin, "fasta")
    child.stdin.close()
    # print(child.stdout.read()) ## can be used only once since child is closed
    df = pd.read_table (io.StringIO(child.stdout.read()), header=None)
    default_outfmt6_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq ssqeq'.strip().split(' ')
    df.columns = default_outfmt6_cols
    return df


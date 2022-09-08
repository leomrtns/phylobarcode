#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import pandas as pd, numpy as np
import itertools, io, multiprocessing
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

# legacy code (no need to create a distinct stream)
#log_format = logging.Formatter(fmt='phylobarcode__blast %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

database_fallback = "/media/teradisk/03.n114312_bigdata/blast/ref_prok_rep_genomes"

def blast_primers_from_tsv (tsv=None, output=None, database = None, evalue=1, task="blastn", max_target_seqs = 1000, nthreads=1):
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
        logger.error("I need a blast database (full path with prefix of DB in blast format; usually filename without '.nal' or '.nsq')")
        return
    if evalue > 10000: evalue = 10000
    if evalue < 1e-3: evalue = 1e-3
    if max_target_seqs > 1e9: max_target_seqs = 1e9
    if max_target_seqs < 1: max_target_seqs = 1
    if task not in ["blastn-short", "blastn"]:
        logger.error("I need a task = 'blastn-short' or 'blastn'; I'll use 'blastn'")

    df_l = pd.read_csv (tsv[0], compression="infer", sep="\t", dtype='unicode')
    df_r = pd.read_csv (tsv[1], compression="infer", sep="\t", dtype='unicode')
    primers_l = df_l["primer"].tolist()
    primers_r = df_r["primer"].tolist()
    logger.info(f"Read {len(primers_l)} primers from file {tsv[0]} and {len(primers_r)} from file {tsv[1]}; Running blast now")

    if nthreads < 2: # single BLAST job, still this one job will use all available CPUs
        ncpus = multiprocessing.cpu_count()
        if ncpus < 1: ncpus = 1
        logger.info(f"Using single thread with {ncpus} CPUs per blast job")
        logger.info(f"Running blast on left primers from file {tsv[0]}")
        blast_l = query_primers_blastn_on_database (primers_l, database=database, evalue=evalue, task=task,
            max_target_seqs=max_target_seqs, ncpus=ncpus)
        logger.info(f"Running blast on right primers from file {tsv[1]}")
        blast_r = query_primers_blastn_on_database (primers_r, database=database, evalue=evalue, task=task,
            max_target_seqs=max_target_seqs, ncpus=ncpus)
    else: # user wants several BLAST jobs
        try:
            from multiprocessing import Pool
            from functools import partial
        except ImportError:
            ncpus = multiprocessing.cpu_count()
            logger.error(f"Multiprocessing not available for python, will run one blast process with {ncpus} CPUs")
            logger.info(f"Running blast on left primers from file {tsv[0]}")
            blast_l = query_primers_blastn_on_database (primers_l, database=database, evalue=evalue, task=task,
                    max_target_seqs=max_target_seqs, ncpus=ncpus)
            logger.info(f"Running blast on right primers from file {tsv[1]}")
            blast_r = query_primers_blastn_on_database (primers_r, database=database, evalue=evalue, task=task,
                    max_target_seqs=max_target_seqs, ncpus=ncpus)
        else:
            logger.info(f"Running blast on left primers from file {tsv[0]}")
            blast_l = query_primers_blastn_on_database_parallel (primers_l, database=database, evalue=evalue, task=task,
                    max_target_seqs=max_target_seqs, nthreads=nthreads)
            logger.info(f"Running blast on right primers from file {tsv[1]}")
            blast_r = query_primers_blastn_on_database_parallel (primers_r, database=database, evalue=evalue, task=task,
                    max_target_seqs=max_target_seqs, nthreads=nthreads)

    print (blast_l) # FIXME: need to save as tsv at least; STOPHERE

def query_primers_blastn_on_database_parallel (primer_list, database=None, evalue = 1, task="blastn-short", max_target_seqs=1000, nthreads=1):
    n_primers = len(primer_list)
    if nthreads > n_primers: nthreads = len(primer_list)
    chunk_size = n_primers//nthreads + 1 # no garantee that all threads will be used, specially if chumk_size is small like 2 or 3
    ## see extract_riboproterins, ids below might not be needed since python handles properly indices exceeding list length
    chunk_ids = [i for i in range(0, n_primers, chunk_size)] + [n_primers] # from [i] to [i+1] excusive
    nthreads = len(chunk_ids) - 1 
    primers_per_thread = [primer_list[chunk_ids[i]:chunk_ids[i+1]] for i in range(nthreads)]
    ncpus = multiprocessing.cpu_count() // nthreads
    if ncpus < 1: ncpus = 1
    logger.info (f"Using {nthreads} threads from those available for blast, and {ncpus} CPUs per blast job")

    from multiprocessing import Pool
    from functools import partial
    with Pool(nthreads) as p:
        results = p.map(
                partial(query_primers_blastn_on_database, database=database, evalue=evalue, task=task, max_target_seqs=max_target_seqs, ncpus=ncpus), 
                primers_per_thread)
    results = [i for i in results if i is not None]
    if len(results): return pd.concat(results)
    else:            return None

def query_primers_blastn_on_database (primer_list, database=None, evalue = 1, task="blastn-short", max_target_seqs=1000, ncpus=1):
    """ 
    task = "blastn-short" (more precise) or task = "blastn" 
    nthreads = number of threads for blast query
    """
    if database is None:
        database = database_fallback
        logger.error(f"I need a database (full path with prefix of DB in blast format; usually filename without '.nal' \
                or '.nsq'); I will try to use {database} but there are no guarantees it will work")
    cline = f"blastn -query - -db {database} -evalue {evalue} -task \"{task}\" -word_size 7 -out - \
              -outfmt \"6 std qseq sseq\"  -max_target_seqs {max_target_seqs} -num_threads {ncpus}" ## format 6 plus aligned querey plus aligned subject
    child = subprocess.Popen(cline, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                            universal_newlines=True, shell=(sys.platform!="win32"))
    records = [SeqRecord(Seq.Seq(seq), id=seq) for seq in primer_list] ## head is primer itself (like "> ACGTGT") 
    SeqIO.write(records, child.stdin, "fasta")
    child.stdin.close()
    # print(child.stdout.read()) ## can be used only once since child is closed
    try:
        df = pd.read_table (io.StringIO(child.stdout.read()), header=None)
    except pd.errors.EmptyDataError:
        logger.error("No data returned from blast; try increasing evalue and not using `accurate` task")
        return None
    default_outfmt6_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq ssqeq'.strip().split(' ')
    df.columns = default_outfmt6_cols
    return df

# legacy code
def query_string_blastn_on_database (seqstring, database=None, evalue = 1, task="blastn", max_target_seqs=1000, nthreads=1):
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

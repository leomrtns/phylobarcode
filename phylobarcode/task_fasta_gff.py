#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import pandas as pd, numpy as np
import itertools, io, multiprocessing, shutil
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

# legacy code, no need to create a separate logger
#log_format = logging.Formatter(fmt='phylobarcode_fasgff %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

#logger.propagate = False
#stream_log = logging.StreamHandler() 
#stream_log.setFormatter(log_format)
#stream_log.setLevel(logging.INFO)
#logger.addHandler(stream_log)

def merge_fasta_gff (fastadir=None, gffdir=None, scratch=None, output=None):
    hash_name = '%012x' % random.randrange(16**12)  # use same file random file name for all files (notice that main script should have taken care of these)
    if fastadir is None: 
        logger.error("No fasta directory provided"); return
    if not os.path.isdir (fastadir):
        logger.error(f"Fasta directory provided {fastadir} does not exist or is not a proper directory"); return
    if gffdir is None: 
        logger.error("No GFF3 directory provided"); return
    if not os.path.isdir (gffdir):
        logger.error(f"GFF3 directory provided {gffdir} does not exist or is not a proper directory"); return
    if output is None: ## this should not happen if function called from main script
        prefix = f"fastagff.{hash_name}"
        logger.warning (f"No output file specified, using {prefix} as prefix")
    if scratch is None: ## this should not happen if function called from main script; use current directory 
        scratch = f"scratch.{hash_name}"
    # create scratch directory (usually it's a subdirectory of the user-given scratch directory)
#    pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist
    fasta_files = \
            glob.glob(f"{fastadir}/*.fasta") + glob.glob(f"{fastadir}/*.fasta.*") + \
            glob.glob(f"{fastadir}/*.fa")    + glob.glob(f"{fastadir}/*.fa.*") + \
            glob.glob(f"{fastadir}/*.fna")   + glob.glob(f"{fastadir}/*.fna.*") + \
            glob.glob(f"{fastadir}/*.fas")   + glob.glob(f"{fastadir}/*.fas.*")
    print (f"Found {len(fasta_files)} fasta files in {fastadir}")
    #read_fasta_headers_as_list (filename);

    # delete scratch subdirectory and all its contents
#    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory



        


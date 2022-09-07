#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json? imports itertools, pathlib
import pandas as pd, numpy as np
import io, multiprocessing, shutil, gffutils
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

# legacy code, no need to create a separate logger
#log_format = logging.Formatter(fmt='phylobarcode_fasgff %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

def extract_riboproteins (tsvfile = None, gffdir=None, output=None, nthreads=1, scratch=None):
    hash_name = '%012x' % random.randrange(16**12) 
    if gffdir is None: 
        logger.error("No GFF3 directory provided"); return
    if not os.path.isdir (gffdir):
        logger.error(f"GFF3 directory provided {gffdir} does not exist or is not a proper directory"); return
    if output is None: ## this should not happen if function called from main script
        prefix = f"riboprot.{hash_name}"
        logger.warning (f"No output file specified, using {prefix} as prefix")
    if scratch is None: ## this should not happen if function called from main script; use current directory 
        scratch = f"scratch.{hash_name}"
    # create scratch directory (usually it's a subdirectory of the user-given scratch directory)
    pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

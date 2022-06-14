#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='phylobarcode_PRIM %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

def find_primers (defaults, fasta = None, output = None, entry_timestamp = None):
    print ("madadayo!")

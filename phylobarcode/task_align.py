#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json? imports itertools, pathlib
import pandas as pd, numpy as np
import io, multiprocessing, shutil, gffutils, json, collections
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
logger = logging.getLogger("phylobarcode_global_logger")

def cluster_align_gene_files (genefiles = None, outfile = None, nthreads = 1, csvfile = None, scratch = None):
    if genefiles is None:
        logger.error("No gene files provided")
        sys.exit(1)
    hash_name = '%012x' % random.randrange(16**12)
    if scratch is None:
        scratch = f"scratch.{hash_name}"
        logger.warning(f"No scratch directory provided, using {scratch}")
    if outfile is None:
        outfile = f"phyloaln.{hash_name}"
        logger.warning(f"No output file (prefix) provided, using {outfile}")

    if not os.path.exists(scratch):
        pathlib.Path(scratch).mkdir(parents=True, exist_ok=True)
        scratch_created = True

    if csvfile is not None:
        taxon_df = pd.read_csv(csvfile, sep='\t', header=0)
        taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, replace="unknown")
    else:
        taxon_df = None

    shortname = remove_prefix_suffix (genefiles)
    for short, long in zip(shortname, genefiles):
        cluster_align_each_gene (short, long, outfile, scratch, taxon_df, nthreads)

    if scratch_created:
        shutil.rmtree(pathlib.Path(scratch))

def cluster_align_each_gene (shortname, genefile, outfile, scratch, taxon_df, nthreads):
    fas = read_fasta_as_list (genefile)
    if len(fas) < 4:
        logger.warning (f"Too few sequences in {genefile}, skipping")
        return
    logger.info(f"Read {shortname} gene (file {genefile}) with {len(fas)} sequences")
    taxonomy = {}
    for x in fas:
        x.id = x.id.split("|")[0] # remove gene name e.g. ">NZ_CP028136.1|S31"
        taxonomy[x.id] = x.description.split(" ", 1)[1].split("|")[1:5] # |order|family|genus|species|

    print (taxonomy) # STOPHERE



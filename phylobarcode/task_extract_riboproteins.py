#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json? imports itertools, pathlib
import pandas as pd, numpy as np
import io, multiprocessing, shutil, gffutils, json
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

# legacy code, no need to create a separate logger
#log_format = logging.Formatter(fmt='phylobarcode_fasgff %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

def extract_riboproteins_from_gff (tsvfile = None, gffdir=None, output=None, jsonfile=None, nthreads=1, scratch=None):
    hash_name = '%012x' % random.randrange(16**12) 
    if gffdir is None: 
        logger.error("No GFF3 directory provided"); return
    if not os.path.isdir (gffdir):
        logger.error(f"GFF3 directory provided {gffdir} does not exist or is not a proper directory"); return
    if output is None: ## this should not happen if function called from main script
        prefix = f"riboprot.{hash_name}"
        logger.warning (f"No output file specified, using {prefix} as prefix")
    if jsonfile is None: 
        jsonfile = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/riboprotein_names.json")
    if scratch is None: ## this should not happen if function called from main script; use current directory 
        scratch = f"scratch.{hash_name}"
    # create scratch directory (usually it's a subdirectory of the user-given scratch directory)
    pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    df = pd.read_csv (tsvfile, sep="\t", dtype=str)
    df.dropna(subset=["gtdb_accession"], inplace=True) # same as df = df[~df["gtdb_accession"].isnull()]
    try:
        ribonames = json.load(open(jsonfile, "r"))
    except ValueError as e:
        logger.error(f"Error loading riboprotein names from {jsonfile}: {e}")
        ribonames = None

    gfiles = []
    for gf in df["gff_file"].unique():
        fullgf = os.path.join (gffdir, gf)
        if os.path.isfile (fullgf): gfiles.append (gf)
        else: logger.error (f"GFF file {fullgf} does not exist, skipping")

    if (nthreads > 1): # main() already checked that modules are available (o.w. nthreads=1)
        logger.info (f"Extracting ribosomal proteins from {len(gfiles)} GFF3 files using up to {nthreads} threads")
        from multiprocessing import Pool
        from functools import partial
        n_files = len (gfiles)
        if nthreads > n_files: nthreads = n_files
        chunk_size = n_files // nthreads + 1 
        gfile_chunks = [gfiles[i:i+chunk_size] for i in range(0, n_files, chunk_size)]
        with Pool (len(gfile_chunks)) as p:
            results = p.map (partial(get_features_from_gff, gff_dir=gffdir, scratch_dir=scratch, ribonames=ribonames), gfile_chunks)
        tbl = [row for chunk in results if chunk is not None for row in chunk] # [[[1,2],[4,5]], [[7,8],[10,11]]] -> [[1,2],[4,5],[7,8],[10,11]]

    else: ## one thread
        logger.info (f"Extracting ribosomal proteins from {len(gfiles)} GFF3 files using a single thread")
        tbl = get_features_from_gff (gff_file = gf, gff_dir = gffdir, scratch_dir = scratch, ribonames = ribonames)

    logger.info (f"Extracted information about {len(tbl)} ribosomal proteins")
    
    tbl = list(map(list, zip(*tbl))) # transpose s.t. each row is a feature
    tbl = {k:v for k,v in zip (["seqid","start", "end", "strand", "product"], tbl)}
    df = pd.DataFrame.from_dict (tbl, orient="columns")
    tsvfile = f"{output}.tsv.gz"
    df.to_csv (tsvfile, sep="\t", index=False)
    logger.info (f"Saved information about ribosomal proteins to {tsvfile}")

    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

def get_features_from_gff (gff_file_list, gff_dir, scratch_dir, ribonames):
    database = os.path.join (scratch_dir, os.path.basename(gff_file_list[0]) + ".db") ## unique name for the database, below scratch dir
    a = []
    n_files = len (gff_file_list)
    for i, gf in enumerate(gff_file_list):
        if i and i % (n_files//10) == 0: logger.info (f"{round((i*100)/n_files,1)}% of files processed in thread {gff_file_list[0]}")
        gff_file = os.path.join (gff_dir, gf) ## full path to GFF3 file
        db = gffutils.create_db(gff_file, dbfn=database, force=True, keep_order=False, merge_strategy="merge", sort_attribute_values=False)
        for i in db.features_of_type('CDS'):
            if any ("ribosomal protein" in str.lower(x) for x in i.attributes["product"]):
                prod = str.upper(i.attributes["product"][0])
                prod = prod[prod.find('RIBOSOMAL PROTEIN')+17:] # remove "ribosomal protein" from beginning; find() returns -1 if not found or idx of first match
                if prod in ribonames.keys():  # to inspect all possible names, remove this if statement and store all `prod`
                    prod = ribonames[prod] # replace with standard name
                    a.append ([i.seqid, i.start, i.end, i.strand, prod]) # gff_file doesnt know which seq from fasta file, seqid does

    #pathlib.Path(database).unlink() # delete database file (delete whole tree later)
    return a

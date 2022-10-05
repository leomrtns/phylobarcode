#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json? imports itertools, pathlib
import pandas as pd, numpy as np
import io, multiprocessing, shutil, gffutils, json, collections
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

# legacy code, no need to create a separate logger
#log_format = logging.Formatter(fmt='phylobarcode_fasgff %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
logger = logging.getLogger("phylobarcode_global_logger")

# generate pandas table with info from all GFF3 files

def extract_coordinates_from_gff (tsvfile=None, gffdir=None, output=None, jsonfiles=None, coord_tsvfile = None, nthreads=1, scratch=None):
    hash_name = '%012x' % random.randrange(16**12) 
    if tsvfile is None:
        logger.error ("No TSV file with matches between fasta and GFF3 files given, exiting"); sys.exit(1)
    if gffdir is None: 
        logger.error("No GFF3 directory provided"); sys.exit(1)
    if not os.path.isdir (gffdir):
        logger.error(f"GFF3 directory provided {gffdir} does not exist or is not a proper directory"); sys.exit(1)
    if output is None: ## this should not happen if function called from main script
        output = f"coordinates.{hash_name}"
        logger.warning (f"No output file specified, using {output} as prefix")
    if jsonfiles is None: 
        jsonfiles = {
                "riboproteins": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/riboproteins_names.json"),
                "ribogenes": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/ribogenes_names.json"),
                "extragenes": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/extragenes_names.json")
                }

    if scratch is None: ## this should not happen if function called from main script; use current directory 
        scratch = f"scratch.{hash_name}"

    df = pd.read_csv (tsvfile, sep="\t", dtype=str)
    # currently we work only with genomes included in GTDB (i.e. QC passed)
    df.dropna(subset=["gtdb_accession"], inplace=True) # same as df = df[~df["gtdb_accession"].isnull()]
    n_files = len(df["gff_file"].unique())
    logger.info (f"Found {len(df)} genomes from {n_files} GFF3 files and GTDB taxonomic info from file {tsvfile}")

    jmap = read_json_files (jsonfiles)

    if coord_tsvfile is not None:
        coord_df = pd.read_csv (coord_tsvfile, sep="\t", dtype=str)
        existing_seqids = coord_df["seqid"].unique()
        existing_gff_files = df.loc[df["seqid"].isin(existing_seqids), "gff_file"].unique() # files with already extracted coordinates
        df = df[~df["gff_file"].isin(existing_gff_files)]
        logger.info (f"It is assumed that coordinates from all genomes (seqid) from each file in coordinate table"
                f"{coord_tsvfile} have been extracted.\n If this is not the case, please run the command again without"
                f"the --coord_tsvfile option (if, for instance, a GFF3 or a fasta file was updated).")
        if len(df) == 0:
            logger.info(f"File {coord_tsvfile} has info about all {len(existing_seqids)} genomes")
            return 
        else:
            logger.info (f"File {coord_tsvfile} has info about {len(existing_seqids)} genomes; "
                f"{len(df)} genomes left to process from file {tsvfile}")

    # create scratch directory (usually it's a subdirectory of the user-given scratch directory)
    pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    gfiles = []
    for gf in df["gff_file"].unique():
        fullgf = os.path.join (gffdir, gf)
        if os.path.isfile (fullgf): gfiles.append (gf)
        else: logger.error (f"GFF file {fullgf} does not exist, skipping")
    if n_files < 1: 
        logger.error ("No GFF3 files available, exiting"); sys.exit(1)

    if (nthreads > 1): # main() already checked that modules are available (o.w. nthreads=1)
        logger.info (f"Extracting ribosomal proteins from {len(gfiles)} GFF3 files using up to {nthreads} threads")
        logger.info (f"Threads are named after first file in pool (i.e. names are arbitrary and do not relate to file itself)")
        from multiprocessing import Pool
        from functools import partial
        n_files = len (gfiles) # see below for alternative oneliner using slice
        if nthreads > n_files: nthreads = n_files
        chunk_size = n_files // nthreads + 1 
        gfile_chunks = [gfiles[i:i+chunk_size] for i in range(0, n_files, chunk_size)]
        with Pool (len(gfile_chunks)) as p:
            results = p.map (partial(
                get_features_from_gff, gff_dir=gffdir, scratch_dir=scratch, jmap=jmap), 
                gfile_chunks)
        tbl = [row for chunk in results if chunk is not None for row in chunk] # [[[1,2],[4,5]], [[7,8],[10,11]]] -> [[1,2],[4,5],[7,8],[10,11]]

    else: ## one thread
        logger.info (f"Extracting ribosomal proteins from {len(gfiles)} GFF3 files using a single thread")
        logger.info (f"Thread is named after first file in pool (i.e. name is arbitrary and does not relate to file itself)")
        tbl = get_features_from_gff (gff_file = gf, gff_dir = gffdir, scratch_dir = scratch, jmap = jmap)

    logger.info (f"Extracted information about {len(tbl)} ribosomal proteins")
    
    tbl = list(map(list, zip(*tbl))) # transpose s.t. each row is a feature
    tbl = {k:v for k,v in zip (["seqid","start", "end", "strand", "product"], tbl)}
    df = pd.DataFrame.from_dict (tbl, orient="columns")
    tsvfile = f"{output}.tsv.xz"
    df.to_csv (tsvfile, sep="\t", index=False)
    logger.info (f"Saved information about ribosomal proteins to {tsvfile}")

    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

def read_json_files (jsonfiles):
    """ Read json files with mappings from GFF gene and product names to standardised names; return a dictionary with 3 dicts"""
    logger.info (f"Reading json files with mappings from GFF gene and product names to standardised names")
    classnames = ["riboproteins", "ribogenes", "extragenes"]
    jmap = {}
    for k in classnames:
        if k not in jsonfiles:
            logger.warning (f"JSON file with {k} not provided, skipping")
            jmap[k] = None
        else:
            try:
                jmap[k] = json.load (open(jsonfiles[k], "r"))
            except ValueError as e:
                logger.error(f"Error loading {k} from {jsonfiles[k]}: {e}")
                logger.error("Coordinates will be extracted but gene/product names will not be corrected; " 
                        "information will be missing downstream")
                jmap[k] = None
            else:
                logger.info (f"Loaded {k} from {jsonfiles[k]}")
    return jmap

def get_features_from_gff (gff_file_list, gff_dir, scratch_dir, jmap):
    database = os.path.join (scratch_dir, os.path.basename(gff_file_list[0]) + ".db") ## unique name for the database, below scratch dir
    a = []
    n_files = len (gff_file_list)
    for i, gf in enumerate(gff_file_list):
        if i and i % (n_files//10) == 0: 
            logger.info (f"{round((i*100)/n_files,1)}% of files processed, {len(a)} riboprotein genes found so far from thread {gff_file_list[0]}")
        gff_file = os.path.join (gff_dir, gf) ## full path to GFF3 file
        db = gffutils.create_db(gff_file, dbfn=database, force=True, keep_order=False, merge_strategy="merge", sort_attribute_values=False)
        for i in db.features_of_type('CDS'):
            gene = i.attributes["gene"][0] if "gene" in i.attributes else None
            i.start = int(i.start - 1) # gffutils uses 1-based coordinates, we want 0-based
            i.end = int(i.end -1)
            if jmap["ribogenes"] and gene and gene in jmap["ribogenes"].keys(): # ribosomal genes have names like rpsJ, rplK, etc.
                gene = jmap["ribogenes"][gene]
                a.append ([i.seqid, i.start, i.end, i.strand, gene]) # gff_file doesnt know which seq from fasta file, seqid does
            elif jmap["extragenes"] and gene and gene in jmap["extragenes"].keys(): ## not RIBOSOMAL PROTEIN, but still close to operons
                gene = jmap["extragenes"][gene] + "_" # underscores marks non-riboprotein genes
                a.append ([i.seqid, i.start, i.end, i.strand, gene])
            elif any ("ribosomal protein" in str.lower(x) for x in i.attributes["product"]):
                prod = str.upper(i.attributes["product"][0])
                prod = prod[prod.find('RIBOSOMAL PROTEIN')+17:].lstrip() # remove "ribosomal protein" from beginning; find() returns -1 if not found or idx of first match
                if jmap["riboproteins"] and prod in jmap["riboproteins"].keys():  # to inspect all possible names, remove this if statement and store all `prod`
                    prod = jmap["riboproteins"][prod] # replace with standard name
                    a.append ([i.seqid, i.start, i.end, i.strand, prod])

    #pathlib.Path(database).unlink() # delete database file (delete whole tree later)
    return a

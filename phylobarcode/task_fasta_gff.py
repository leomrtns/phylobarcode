#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import pandas as pd, numpy as np
import itertools, io, multiprocessing, shutil, gffutils
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
    pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    # check if directories contain fasta and gff files first, before doing anything else
    fasta_files = list_of_files_by_extension (fastadir, ['fasta', 'fa', 'fna', 'faa', 'ffn', 'faa', 'fas'])
    logger.info(f"Found {len(fasta_files)} fasta files in {fastadir}")
    gff_files = list_of_files_by_extension (gffdir, ['gff', 'gff3'])
    logger.info(f"Found {len(gff_files)} gff files in {gffdir}")
    if (len(fasta_files) == 0):
        logger.error(f"No fasta files found in {fastadir}")
        return
    if (len(gff_files) == 0):
        logger.error(f"No gff files found in {gffdir}")
        return

    # get sequence names as dataframe
    df_fasta = pd.DataFrame()
    for fasfile in fasta_files[::100]:
        a = split_headers_in_fasta(fasfile)
        df_fasta = pd.concat ([df_fasta, a]) if not df_fasta.empty else a
    logger.info(f"Found {len(df_fasta)} sequences in {len(fasta_files)} fasta files")
    # get GFF chromosomes (excludes plasmids), using scratch dir to store the sqlite db, as dataframe
    dbfile = f"{scratch}/gff.db"
    df_gff = pd.DataFrame()
    for gffile in gff_files[::100]:
        a = split_region_elements_in_gff (gffile, dbfile)
        df_gff = pd.concat ([df_gff, a]) if not df_gff.empty else a
    logger.info(f"Found {len(df_gff)} chromosomes in {len(gff_files)} gff files; now finding the intersection between fasta and gff files")

    # merge dataframes using seqid as key, keeping only rows found in both
    df = pd.merge(df_fasta, df_gff, on='seqid', how='inner')
    print (df.shape)
    print (df.head())

    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory


def list_of_files_by_extension (dirname, extension):
    files = []
    for ext in extension:
        files += glob.glob(f"{dirname}/*.{ext}") + glob.glob(f"{dirname}/*.{ext}.*")
    return files

def split_headers_in_fasta (fasta_file):
    a = [x.split(",")[0] for x in read_fasta_headers_as_list (fasta_file)] # read_fasta_headers is defined in
    a = [[os.path.basename(fasta_file)] + x.split(" ",1) for x in a] # filename +  split on first space
    if (len(a) == 0): 
        return pd.DataFrame()
    else:
        a = list(map(list, zip(*a))) # transpose list of lists (https://stackoverflow.com/questions/6473679/python-transpose-list-of-lists)
    a = {"fasta_file": a[0], "seqid": a[1], "fasta_description": a[2]} # dictionary of lists (one row is chromosome and others are plasmids usually
    return pd.DataFrame.from_dict(a, orient='columns') # orient: each key is a column

def split_region_elements_in_gff (gff_file, database_file):
    db = gffutils.create_db (gff_file, database_file, merge_strategy='create_unique', keep_order=True, force=True) # force to overwrite existing db
    a = []
    for ft in db.features_of_type('region', order_by='start'):
        if ft.attributes['genome'][0] == 'chromosome': # skip plasmids
            longname = ""
            if ("old-name" in ft.attributes): 
                longname = ft.attributes["old-name"][0] + ";" ## alternative to ft["old-name"]
            if ("type-material" in ft.attributes): 
                longname = ft.attributes["type-material"][0] + ";" 
            if ("strain" in ft.attributes):
                longname = ft.attributes["strain"][0] + ";"
            a.append ([os.path.basename(gff_file), ft.seqid, longname, ft["Dbxref"][0].replace("taxon:","")]) # filename + seqid + longname + Dbxref
    a = list(map(list, zip(*a))) # transpose list of lists (https://stackoverflow.com/questions/6473679/python-transpose-list-of-lists)
    a = {"gff_file": a[0], "seqid": a[1], "gff_description": a[2], "gff_taxonid": a[3]} # dictionary of lists (usually one row only since chromosome)
    return pd.DataFrame.from_dict(a, orient='columns') # orient: each key is a column




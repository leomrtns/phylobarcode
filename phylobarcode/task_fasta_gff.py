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

#logger.propagate = False
#stream_log = logging.StreamHandler() 
#stream_log.setFormatter(log_format)
#stream_log.setLevel(logging.INFO)
#logger.addHandler(stream_log)


gtdb_columns_keep = ['accession', 'gtdb_genome_representative', 'gtdb_taxonomy', 'ssu_query_id', 
        'ncbi_assembly_name', 'ncbi_genbank_assembly_accession', 'ncbi_strain_identifiers', 'ncbi_taxid',
        'ncbi_taxonomy', 'ncbi_taxonomy_unfiltered']
gtdb_columns_rename = {
        'accession': 'gtdb_accession', 
        'ssu_query_id': 'seqid'}


def merge_fasta_gff (fastadir=None, gffdir=None, fasta_tsvfile = None, gff_tsvfile = None, gtdb = None, scratch=None, output=None):
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

    # check if directories contain fasta files first
    fasta_files = list_of_files_by_extension (fastadir, ['fasta', 'fa', 'fna', 'faa', 'ffn', 'faa', 'fas'])
    logger.info(f"Found {len(fasta_files)} fasta files in {fastadir}")
    if (len(fasta_files) < 1):
        logger.error(f"Not enough fasta files found in {fastadir}")
        return

    # get sequence names as dataframe
    fasta_files, fasta_tsv = update_tsv_from_filenames (fasta_files, fasta_tsvfile, "fasta_file")
    if len(fasta_files):
        a = [] # list of lists (samples=rows, features=columns)
        for fasfile in fasta_files:
            a.extend(split_headers_in_fasta (fasfile)) # append() would create 3x list; extend() "flattens" each element like "plus"
        a = list(map(list, zip(*a))) # transpose list of lists (https://stackoverflow.com/questions/6473679/python-transpose-list-of-lists)
        a = {"fasta_file": a[0], "seqid": a[1], "fasta_description": a[2]} # dictionary of lists (one row is chromosome and others are plasmids usually)
        df_fasta = pd.DataFrame.from_dict(a, orient='columns')
        if fasta_tsv is not None:
            logger.info(f"FASTA: {len(fasta_tsv)} sequences already in tsv file {fasta_tsvfile}")
            logger.info(f"FASTA: {len(df_fasta)} new sequences found in directory {fastadir}")
            df_fasta = pd.concat([df_fasta, fasta_tsv], ignore_index=True)
        else:
            logger.info(f"FASTA: found {len(df_fasta)} sequences in directory {fastadir}")
        tsvfilename = f"{output}_fasta.tsv.gz"
        df_fasta.to_csv (tsvfilename, sep="\t", index=False)
        logger.info(f"All fasta entries wrote to {tsvfilename}")
    else:
        logger.info (f"All fasta files already found in tsv file {fasta_tsvfile}, with {len(fasta_tsv)} entries")
        df_fasta = fasta_tsv

    # check if directories contain gff files first
    gff_files = list_of_files_by_extension (gffdir, ['gff', 'gff3'])
    logger.info(f"Found {len(gff_files)} gff files in {gffdir}")
    if (len(gff_files) < 1):
        logger.error(f"Not enough gff files found in {gffdir}")
        return

    # get GFF chromosomes (excludes plasmids), using scratch dir to store the sqlite db, as dataframe
    gff_files, gff_tsv = update_tsv_from_filenames (gff_files, gff_tsvfile, "gff_file")
    if len(gff_files):
        dbfile = f"{scratch}/gff.db"
        a = []
        for gffile in gff_files[::10]:
            #dbfile = f"{scratch}/{pathlib.Path(gffile).stem}.db" # Path = os.path.basename but "stem" removes extension
            a.extend(split_region_elements_in_gff (gffile, dbfile))
        a = list(map(list, zip(*a))) # transpose list of lists so that each row is one feature
        a = {"gff_file": a[0], "seqid": a[1], "gff_description": a[2], "gff_taxonid": a[3]} # dictionary of lists (usually one row only since chromosome)
        df_gff = pd.DataFrame.from_dict(a, orient='columns')
        if gff_tsv is not None:
            logger.info(f"GFF: {len(gff_tsv)} sequences already in tsv file {gff_tsvfile}")
            logger.info(f"GFF: {len(df_gff)} new sequences found in directory {gffdir}")
            df_gff = pd.concat([df_gff, gff_tsv], ignore_index=True)
        else:
            logger.info(f"GFF: found {len(df_gff)} sequences in directory {gffdir}")
        tsvfilename = f"{output}_gff.tsv.gz"
        df_gff.to_csv (tsvfilename, sep="\t", index=False)
        logger.info(f"All GFF entries wrote to {tsvfilename}")
    else:
        logger.info (f"All GFF files already found in tsv file {gff_tsvfile}, with {len(gff_tsv)} entries")
        df_gff = gff_tsv

    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

    # merge dataframes using seqid as key, keeping only rows found in both

    if (gtdb is None):
        logger.warning("No GTDB taxonomy file provided, the table with matching genomes will _not_ be created. Use the "
                "generated tsv files together with the GTDB file in order to produce the final table.")
        return

    df = pd.merge(df_fasta, df_gff, on='seqid', how='inner')
    if (len(df) < 1):
        logger.error(f"No common entries found in fasta and gff files. Cannot produce the final table with merged info.")
        return
    logger.info(f"Found {len(df)} common entries in fasta and gff files; Will now read GTDB file {gtdb} and merge")

    # get GTDB taxonomy as dataframe
    df_gtdb = pd.read_csv(gtdb, sep="\t", dtype=str)
    df_gtdb = df_gtdb[gtdb_columns_keep] # remove most columns
    df_gtdb.rename(columns=gtdb_columns_rename, inplace=True) # rename columns to match other tables
    df = pd.merge(df, df_gtdb, on='seqid', how='left')
    full_dlen = len(df) - df["gtdb_accession"].isnull().sum() # sum=count null values

    tsvfilename = f"{output}_merged.tsv.gz"
    logger.info(f"Found {full_dlen} samples with complete information; writing all (including incomplete) to {tsvfilename}")
    print (df.head())
    df.to_csv (tsvfilename, sep="\t", index=False)


def list_of_files_by_extension (dirname, extension):
    files = []
    for ext in extension:
        files += glob.glob(f"{dirname}/*.{ext}") + glob.glob(f"{dirname}/*.{ext}.*")
    return files

def update_tsv_from_filenames (files, tsvfile, columnname):
    if tsvfile is None:
        return files, None
    if not os.path.isfile (tsvfile):
        logger.error(f"tsv file {tsvfile} does not exist or is not a proper file"); return files, None
    df = pd.read_csv (tsvfile, sep="\t", dtype=str)
    if columnname not in df.columns:
        logger.error(f"tsv file {tsvfile} does not have column {columnname}"); return files, None
    new_files = [f for f in files if os.path.basename(f) not in df[columnname].unique()] 
    old_files = [os.path.basename(f) for f in list(set(files) - set(new_files))]
    df = df[df[columnname].isin(old_files)] # keep only rows with filenames found in files
    return new_files, df

def split_headers_in_fasta (fasta_file):
    a = [x.split(",")[0] for x in read_fasta_headers_as_list (fasta_file)] # read_fasta_headers is defined in
    a = [[os.path.basename(fasta_file), x.split(" ",1)[0], x.split(" ",1)[1]] for x in a] # filename +  split on first space
    return a

def split_region_elements_in_gff (gff_file, database_file):
    db = gffutils.create_db (gff_file, database_file, merge_strategy='create_unique', keep_order=True, force=True) # force to overwrite existing db
    a = []
    for ft in db.features_of_type('region', order_by='start'):
        if "genome" in ft.attributes and ft.attributes['genome'][0] == 'chromosome': # skip plasmids
            longname = ""
            if ("old-name" in ft.attributes): 
                longname = ft.attributes["old-name"][0] + ";" ## alternative to ft["old-name"]
            if ("type-material" in ft.attributes): 
                longname = ft.attributes["type-material"][0] + ";" 
            if ("strain" in ft.attributes):
                longname = ft.attributes["strain"][0] + ";"
            a.append ([os.path.basename(gff_file), ft.seqid, longname, ft["Dbxref"][0].replace("taxon:","")]) # filename + seqid + longname + Dbxref
    return a


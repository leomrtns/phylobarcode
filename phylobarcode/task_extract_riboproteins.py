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

# TODO: extra genes (rpoB, rpoC) are also close to riboproteins
# first command: generate pandas table with info from all GFF3 files

def extract_coordinates_from_gff (tsvfile=None, gffdir=None, output=None, jsonfile=None, nthreads=1, scratch=None):
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
    if jsonfile is None: 
        jsonfile = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/riboprotein_names.json")
    if scratch is None: ## this should not happen if function called from main script; use current directory 
        scratch = f"scratch.{hash_name}"
    # create scratch directory (usually it's a subdirectory of the user-given scratch directory)
    pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    df = pd.read_csv (tsvfile, sep="\t", dtype=str)
    # currently we work only with genomes included in GTDB (i.e. QC passed)
    df.dropna(subset=["gtdb_accession"], inplace=True) # same as df = df[~df["gtdb_accession"].isnull()]
    try:
        ribonames = json.load(open(jsonfile, "r"))
    except ValueError as e:
        logger.error(f"Error loading riboprotein names from {jsonfile}: {e}")
        logger.error("Coordinates will be extracted but gene names will not be corrected; some operons will be skipped downstream")
        ribonames = None

    gfiles = []
    for gf in df["gff_file"].unique():
        fullgf = os.path.join (gffdir, gf)
        if os.path.isfile (fullgf): gfiles.append (gf)
        else: logger.error (f"GFF file {fullgf} does not exist, skipping")

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
            results = p.map (partial(get_features_from_gff, gff_dir=gffdir, scratch_dir=scratch, ribonames=ribonames), gfile_chunks)
        tbl = [row for chunk in results if chunk is not None for row in chunk] # [[[1,2],[4,5]], [[7,8],[10,11]]] -> [[1,2],[4,5],[7,8],[10,11]]

    else: ## one thread
        logger.info (f"Extracting ribosomal proteins from {len(gfiles)} GFF3 files using a single thread")
        logger.info (f"Thread is named after first file in pool (i.e. name is arbitrary and does not relate to file itself)")
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
        if i and i % (n_files//10) == 0: 
            logger.info (f"{round((i*100)/n_files,1)}% of files processed, {len(a)} riboprotein genes found so far from thread {gff_file_list[0]}")
        gff_file = os.path.join (gff_dir, gf) ## full path to GFF3 file
        db = gffutils.create_db(gff_file, dbfn=database, force=True, keep_order=False, merge_strategy="merge", sort_attribute_values=False)
        for i in db.features_of_type('CDS'):
            if any ("ribosomal protein" in str.lower(x) for x in i.attributes["product"]):
                prod = str.upper(i.attributes["product"][0])
                prod = prod[prod.find('RIBOSOMAL PROTEIN')+17:].lstrip() # remove "ribosomal protein" from beginning; find() returns -1 if not found or idx of first match
                if ribonames and prod in ribonames.keys():  # to inspect all possible names, remove this if statement and store all `prod`
                    prod = ribonames[prod] # replace with standard name
                    a.append ([i.seqid, i.start, i.end, i.strand, prod]) # gff_file doesnt know which seq from fasta file, seqid does

    #pathlib.Path(database).unlink() # delete database file (delete whole tree later)
    return a

# second command: extract DNA sequences using pandas table with riboprotein info

def extract_operons_from_fasta (coord_tsvfile=None, merge_tsvfile=None, fastadir=None, output=None, nthreads=1, scratch=None):
    hash_name = '%012x' % random.randrange(16**12) 
    if coord_tsvfile is None:
        logger.error ("No TSV file with riboprot coordinates from GFF3 files given, exiting"); sys.exit(1)
    if merge_tsvfile is None:
        logger.error ("No TSV file with merged info from FASTA and GFF3 files given, exiting"); sys.exit(1)
    if fastadir is None: 
        logger.error ("No FASTA directory given, exiting"); sys.exit(1)
    if not os.path.isdir (fastadir):
        logger.error (f"FASTA directory {fastadir} does not exist or is not a proper directory, exiting"); sys.exit(1)
    if output is None: 
        output = f"operon.{hash_name}"
        logger.warning (f"No output file specified, using {output} as prefix")
    if scratch is None:
        scratch = f"scratch.{hash_name}"
        logger.warning (f"No scratch directory specified, using {scratch} as prefix")

    coord_df = pd.read_csv (coord_tsvfile, sep="\t", dtype = str)
    if coord_df.empty:
        logger.error (f"Coordinates file {coord_tsvfile} is empty, exiting"); sys.exit(1)
    merge_df = pd.read_csv (merge_tsvfile, sep="\t", dtype = str)
    if merge_df.empty:
        logger.error (f"Merged file {merge_tsvfile} with fasta x GFF3 info is empty, exiting"); sys.exit(1)

    genome_list = coord_df["seqid"].unique().tolist()
    # create scratch subdirectory
    pathlib.Path(scratch).mkdir(parents=True, exist_ok=True) # create scratch subdirectory

    if (nthreads > 1): ## multiple threads
        logger.info (f"Extracting operons from {len(genome_list)} genomes using {nthreads} threads")
        from multiprocessing import Pool
        from functools import partial
        genome_chunks = [genome_list[i::nthreads] for i in range(nthreads)]
        g_pool = []
        for g in genome_chunks:
            cdf = coord_df[coord_df["seqid"].isin(g)]
            mdf = merge_df[merge_df["seqid"].isin(g)]
            fname = f"{scratch}/coord.{g[0]}"
            g_pool.append ([cdf, mdf, fname])
        with Pool(len(genome_chunks)) as p:
            results = p.map(partial(extract_and_save_operons, fastadir=fastadir), g_pool)
    else: ## single thread
        logger.info (f"Extracting operons from {len(genome_list)} genomes using 1 thread")
        results = extract_and_save_operons ([coord_df, merge_df, f"{scratch}/coord"], fastadir)
    
    operon_seqs = []
    for g in g_pool:
        operon_seqs.extend (read_fasta_as_list (g[2]))
    ofile = f"{output}.fasta.xz"
    with open_anyformat (ofile, "w") as f:
        for rec in operon_seqs:
            f.write (f"> {rec.description}\n{rec.seq}\n").encode()
    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

def extract_and_save_operons (pool_info, fastadir):
    coord_df, merge_df, fname = pool_info
    genome_list = coord_df["seqid"].unique().tolist()
    fw = open_anyformat (fname, "a")

    for g in genome_list:
        cdf = coord_df[coord_df["seqid"] == g]
        mdf = merge_df[merge_df["seqid"] == g]

        genome_sequence = read_fasta_as_list (os.path.join (fastadir, mdf["fasta_file"].iloc[0]))
        genome_sequence = genome_sequence[0].seq
        operons = operon_from_coords (genome_sequence, cdf)
        for opr,seq in operons.items():
            name = f"> {g} {opr} " + mdf["fasta_description"].iloc[0]
            fw.write (f"{name}\n{seq}\n").encode()
    fw.close()

    def operon_from_coords (genome_sequence, coord_df, intergenic_space=100):
        hug_genes = ["S10","L3","L4","L2","S19","L22","S3","L16","S17","L14","L24","L5","S8","L6","L18","L15"] # 16 riboprots from Hug
        core_genes = hug_genes + ["L23","L29","S14","S5","L30"] 
        left_1_genes = ["L11","L1","L10","L7"] # L7 includes L7/L12 
        left_2_genes = ["S12","S7"]    
        right_genes  = ["L36","S13","S11","S4","L17"] # core between S7 and L36
        # remember that GFF is one-based
        cdf.sort_values(["start"], ascending=True, inplace=True)
        cdf["start"] = cdf["start"].astype(int) - 1 # convert to zero-based
        cdf["end"] = cdf["end"].astype(int) - 1
        minioperons = minioperon_merge_from_coords (cdf, intergenic_space)
        

    def minioperon_merge_from_coords (df, intergenic_space=100):
        # iterate over cdf rows
        for i, row in df.iterrows():
            ## STOPHERE



    
        






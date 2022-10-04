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

# TODO: extra genes (rpoB, rpoC) are also close to riboproteins
# first command: generate pandas table with info from all GFF3 files

def extract_coordinates_from_gff (tsvfile=None, gffdir=None, output=None, jsonfile=None, coord_tsvfile = None, nthreads=1, scratch=None):
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

    df = pd.read_csv (tsvfile, sep="\t", dtype=str)
    # currently we work only with genomes included in GTDB (i.e. QC passed)
    df.dropna(subset=["gtdb_accession"], inplace=True) # same as df = df[~df["gtdb_accession"].isnull()]
    n_files = len(df["gff_file"].unique())
    logger.info (f"Found {len(df)} genomes from {n_files} GFF3 files and GTDB taxonomic info from file {tsvfile}")
    try:
        ribonames = json.load(open(jsonfile, "r"))
    except ValueError as e:
        logger.error(f"Error loading riboprotein (product) names from {jsonfile}: {e}")
        logger.error("Coordinates will be extracted but gene names will not be corrected; some operons will be skipped downstream")
        ribonames = None

    ribogenes_jsonfile = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/ribogenes_names.json")
    try:
        ribogenes = json.load(open(ribogenes_jsonfile, "r"))
    except ValueError as e:
        logger.error(f"Error loading riboprotein (genes) names from {ribogenes_jsonfile}: {e}")
        logger.error("Coordinates will be extracted but gene names will not be corrected; some operons will be skipped downstream")
        ribogenes = None

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
                get_features_from_gff, gff_dir=gffdir, scratch_dir=scratch, ribonames=ribonames,ribogenes=ribogenes), 
                gfile_chunks)
        tbl = [row for chunk in results if chunk is not None for row in chunk] # [[[1,2],[4,5]], [[7,8],[10,11]]] -> [[1,2],[4,5],[7,8],[10,11]]

    else: ## one thread
        logger.info (f"Extracting ribosomal proteins from {len(gfiles)} GFF3 files using a single thread")
        logger.info (f"Thread is named after first file in pool (i.e. name is arbitrary and does not relate to file itself)")
        tbl = get_features_from_gff (gff_file = gf, gff_dir = gffdir, scratch_dir = scratch, ribonames = ribonames, ribogenes = ribogenes)

    logger.info (f"Extracted information about {len(tbl)} ribosomal proteins")
    
    tbl = list(map(list, zip(*tbl))) # transpose s.t. each row is a feature
    tbl = {k:v for k,v in zip (["seqid","start", "end", "strand", "product"], tbl)}
    df = pd.DataFrame.from_dict (tbl, orient="columns")
    tsvfile = f"{output}.tsv.gz"
    df.to_csv (tsvfile, sep="\t", index=False)
    logger.info (f"Saved information about ribosomal proteins to {tsvfile}")

    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

def get_features_from_gff (gff_file_list, gff_dir, scratch_dir, ribonames, ribogenes):
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
            if any ("ribosomal protein" in str.lower(x) for x in i.attributes["product"]):
                prod = str.upper(i.attributes["product"][0])
                prod = prod[prod.find('RIBOSOMAL PROTEIN')+17:].lstrip() # remove "ribosomal protein" from beginning; find() returns -1 if not found or idx of first match
            
                if ribonames and prod in ribonames.keys():  # to inspect all possible names, remove this if statement and store all `prod`
                    prod = ribonames[prod] # replace with standard name
                    a.append ([i.seqid, i.start, i.end, i.strand, prod]) # gff_file doesnt know which seq from fasta file, seqid does
            elif ribogenes and gene and gene in ribogenes.keys(): ## not RIBOSOMAL PROTEIN, but still in list of ribosomal genes
                gene = ribogenes[gene]
                a.append ([i.seqid, i.start, i.end, i.strand, gene]) # gff_file doesnt know which seq from fasta file, seqid does

    #pathlib.Path(database).unlink() # delete database file (delete whole tree later)
    return a

# second command: extract DNA sequences using pandas table with riboprotein info

def extract_operons_from_fasta (coord_tsvfile=None, merge_tsvfile=None, fastadir=None, output=None, 
        intergenic_space = 1000, short_operon = 1000, most_common_mosaics = 50, border = 50, nthreads=1, scratch=None):
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
    if intergenic_space < 1:
        logger.warning (f"Intergenic space is too small, setting to one")
        intergenic_space = 1
    if short_operon < 1:
        logger.warning (f"Short operon is too small, setting to one (i.e. single genes are removed")
        short_operon = 1
    if most_common_mosaics < 2:
        logger.warning (f"Most common mosaics {most_common_mosaics} is too small, setting to 2")
        most_common_mosaics = 2
    if border < 1:
        logger.warning (f"Border {border} is too small, setting to 1")
        border = 1

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
        logger.info (f"Thread is named after first file in pool (i.e. name is arbitrary and does not relate to file itself)")
        from multiprocessing import Pool
        from functools import partial
        genome_chunks = [genome_list[i::nthreads] for i in range(nthreads)]
        g_pool = []
        for g in genome_chunks:
            cdf = coord_df[coord_df["seqid"].isin(g)]
            mdf = merge_df[merge_df["seqid"].isin(g)]
            fname = f"{scratch}/coord.{g[0]}.gz"
            g_pool.append ([cdf, mdf, fname])
        with Pool(len(genome_chunks)) as p:
            results = p.map( partial(
                        extract_and_save_operons, 
                        fastadir=fastadir, 
                        intergenic_space=intergenic_space, 
                        short_operon=short_operon,
                        border=border),
                    g_pool)
        results = [elem for x in results for elem in x] # flatten list of lists [[1,2],[3,4]] -> [1,2,3,4]
    else: ## single thread
        logger.info (f"Extracting operons from {len(genome_list)} genomes using one thread")
        logger.info (f"Thread is named arbitrarily")
        g_pool = [[coord_df, merge_df, f"{scratch}/coord.fa.gz"]] # list of lists to be compatible with multithreaded
        results = extract_and_save_operons (g_pool[0], fastadir=fastadir, intergenic_space=intergenic_space,
            short_operon=short_operon, border=border)
    
    # results is a list of mosaics, we'll save the most common ones
    moscounter = collections.Counter(results)
    mosaics = [x[0] for x in moscounter.most_common(most_common_mosaics)] # get the 10 most common mosaics
    tolog = "\n".join([f"Found in {x[1]} genomes:\t{x[0]}" for x in moscounter.most_common(10)])
    logger.info (f"Finished scanning genomes, the {most_common_mosaics} most common mosaics will be saved, amongst them:\n{tolog}")

    operon_seqs = []
    for g in g_pool:
        operon_seqs.extend (read_fasta_as_list (g[2], substring=mosaics))

    save_mosaics_as_fasta (operon_seqs, output, mosaics)
    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

def save_mosaics_as_fasta (operon_seqs, output, mosaics):
    for i, m in enumerate(mosaics):
        ofile = f"{output}.{m}.fasta.xz"
        counter = 0
        with open_anyformat (ofile, "w") as f:
            for rec in operon_seqs:
                print (rec.description.split())
                if rec.description.split()[1] == m: # header has format "> genomeID mosaic description"
                    f.write (str(f"> {rec.description}\n{rec.seq}\n").encode())
                    counter += 1
        if i < 5:
            logger.info (f"Succesfully saved {counter} operons to {ofile}")
        elif i == 6:
            logger.info (f"etc.")

def extract_and_save_operons (pool_info, fastadir, intergenic_space=1000, short_operon=1000, border=50):
    coord_df, merge_df, fname = pool_info
    genome_list = coord_df["seqid"].unique().tolist()
    fw = open_anyformat (fname, "w")

    def operon_from_coords (genome_sequence, coord_df):
        hug_genes = ["S10","L3","L4","L2","S19","L22","S3","L16","S17","L14","L24","L5","S8","L6","L18","L15"] # 16 riboprots from Hug
        core_genes = hug_genes + ["L23","L29","S14","S5","L30"] 
        left_1_genes = ["L11","L1","L10","L7"] # L7 includes L7/L12 
        left_2_genes = ["S12","S7"]    
        right_genes  = ["L36","S13","S11","S4","L17"] # core between S7 and L36
        coord_df = coord_df.sort_values(by=["start"], ascending=True)
        coord_df = coord_df.reset_index(drop=True) # unlike iterate(), iterrows() returns index (before sorting)
        # remember that GFF is one-based
        coord_df["start"] = coord_df["start"].astype(int) - 1 # convert to zero-based
        coord_df["end"] = coord_df["end"].astype(int) - 1
        genome_length = len(genome_sequence)
        minioperons, extra_space = minioperon_merge_from_coords (coord_df, genome_length)
        minioperons = remove_short_operons (minioperons)
        minioperons = add_borders (minioperons, genome_length, border)
        return dict_of_operons (minioperons, genome_sequence, extra_space)
        
    def minioperon_merge_from_coords (df, genome_length=0):
        # iterate over cdf rows; minioperon will have [[gene1, ..., geneN], [start, end], strand]
        minioperons = []
        for i, row in df.iterrows():
            if (i == 0) or (row["start"] - df.iloc[i-1]["end"] > intergenic_space) or (row["strand"] != df.iloc[i-1]["strand"]):
                minioperons.append ([[row["product"]], [row["start"], row["end"]], row["strand"]])
            elif (row["start"] - df.iloc[i-1]["end"] <= intergenic_space) and (row["strand"] == df.iloc[i-1]["strand"]):
                minioperons[-1][0].append (row["product"])
                minioperons[-1][1][1] = row["end"]
        # if gene crosses zero, then GFF has two entries (first and last); genbank has one entry with 4 locations BTW
        extra_space = 0
        if ((minioperons[0][1][0] < 1) and  # gene crosses zero
            (minioperons[-1][1][1] == genome_length-1) and # gene crosses and of genome
            (minioperons[0][0][0] == minioperons[-1][0][-1]) and # same name for first and last gene
            (minioperons[0][2] == minioperons[-1][2])): # and same orientation --> _same_ _gene_
                # operon1 operon2 .... operonN-1 operonN --> operon2 .... operonN-1 operonN+operon1
                extra_space = minioperons[0][1][1] # space beyond genome length used by first minioperon
                minioperons[-1][1][1] = minioperons[0][1][1] + genome_length # new coordinate goes beyond genome length
                minioperons[-1][0].extend (minioperons[0][0]) # add first minioperon gene names to last minioperon
                minioperons.pop(0) # removes first minioperon which is now redundant with last minioperon
        return minioperons, extra_space

    def remove_short_operons (minioperons):
        # remove minioperons that are too short
        if (short_operon > 100): # remove short operons
            minioperons = [m for m in minioperons if (m[1][1] - m[1][0]) > short_operon]
        else: ## remove single-gene operons
            minioperons = [m for m in minioperons if len(m[0]) > 1]
        return minioperons

    def add_borders (minioperons, genome_length, border):
        for m in minioperons:
            m[1][0] = max (0, m[1][0] - border)
            m[1][1] = min (genome_length-1, m[1][1] + border)
        return minioperons

    def dict_of_operons (minioperons, genome_sequence, extra_space=0):
        operons = {}
        if extra_space > 0: genome_sequence = genome_sequence + genome_sequence[:extra_space]
        # minioperon is a list of [[gene1, ..., geneN], [start, end], strand]
        for m in minioperons:
            if (m[2] == "-"):
                seq = genome_sequence[m[1][0]:m[1][1]+1].reverse_complement()
                name = "".join (m[0][::-1]) # merge all gene names, like "L1L2L3"
            else:
                seq = genome_sequence[m[1][0]:m[1][1]+1]
                name = "".join (m[0]) # merge all gene names, like "L1L2L3"
            operons[name] = seq
        return operons

    # save all operon mosaics into same fasta file
    mosaics = []
    n_genomes = len(genome_list)
    for i, g in enumerate(genome_list):
        cdf = coord_df[coord_df["seqid"] == g]
        mdf = merge_df[merge_df["seqid"] == g]
        if len(cdf) < 2 or len(mdf) < 1: continue # some genomes, specially multi-chromosomal, have one gene 
        # e.g. Burkholderia multivorans strain P1Bm2011b has 3 chromosomes

        if i and i % (n_genomes//10) == 0: 
            logger.info (f"{round((i*100)/n_genomes,1)}% of files processed, {len(mosaics)} mosaics found so far from thread {fname[-32:]}")

        genome_sequence = read_fasta_as_list (os.path.join (fastadir, mdf["fasta_file"].iloc[0]))
        genome_sequence = [x for x in genome_sequence if x.id == g] # one fasta file can have multiple genomes, each
        genome_sequence = genome_sequence[0].seq # biopython SeqRecord object s.t. we can reverse_complement() if needed
        operons = operon_from_coords (genome_sequence, cdf)
        for opr,seq in operons.items():
            name = f"> {g} {opr} " + mdf["fasta_description"].iloc[0]
            fw.write (str(f"{name}\n{seq}\n").encode())
            mosaics.append (opr)
    fw.close()
    return mosaics

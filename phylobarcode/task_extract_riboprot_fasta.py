#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json? imports itertools, pathlib
import pandas as pd, numpy as np
import io, multiprocessing, shutil, gffutils, json, collections
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("phylobarcode_global_logger")

genesets = { # not used; noted here for  legacy purposes
        "hug"     : ["S10","L3","L4","L2","S19","L22","S3","L16","S17","L14","L24","L5","S8","L6","L18","L15"], # 16 riboprots from Hug
        "core"    : ["L23","L29","S14","S5","L30"], 
        "left"    : ["S12","S7", "fus_", "tuf_"],    
        "leftleft": ["secE", "nusG", "L11","L1","L10","L7", "rpoB", "rpoC"], # L7 includes L7/L12 
        "right"   : ["secY", "map", "infA", "L36","S13","S11","S4","L17"] # core between S7 and L36
    }
genesets["main"] = [x for y in genesets.keys() if y != "extended" for x in genesets[y]]
genesets["core"]     += genesets["hug"]
genesets["left"]     += genesets["core"]
genesets["leftleft"] += genesets["left"]
genesets["right"]    += genesets["core"]

# extract DNA sequences using pandas table with riboprotein info
# TODO: create structure storing original genewise info (and save it as tsv like "coords" table but pointing to
# generated fasta files)
# TODO: store only zero-based coordinates (when reading GFF)

def extract_operons_from_fasta (coord_tsvfile=None, merge_tsvfile=None, fastadir=None, output=None, 
        intergenic_space = 1000, short_operon = 1000, most_common_mosaics = 50, border = 50, riboprot_subset = None, 
        nthreads=1, scratch=None):
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
    coord_df = coord_df.drop_duplicates () # some sequences appear twice in the table
    if coord_df.empty:
        logger.error (f"Coordinates file {coord_tsvfile} is empty, exiting"); sys.exit(1)
    merge_df = pd.read_csv (merge_tsvfile, sep="\t", dtype = str)
    if merge_df.empty:
        logger.error (f"Merged file {merge_tsvfile} with fasta x GFF3 info is empty, exiting"); sys.exit(1)

    if isinstance (riboprot_subset, str) and riboprot_subset in genesets:
        coord_df = coord_df[coord_df["product"].isin(genesets[riboprot_subset])]
        logger.info (f"Using only {riboprot_subset} riboproteins: {genesets[riboprot_subset]}")
    elif isinstance (riboprot_subset, str): # unkonwn set; will just remove  non-riboproteins (end with "_")
        coord_df = coord_df[~coord_df["product"].str.endswith("_")]
        logger.info (f"Using only riboproteins (removing other genes)")
    else:
        logger.info (f"Using all genes from coordinates file")

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
    save_mosaic_frequency (moscounter, output)
    # delete scratch subdirectory and all its contents
    shutil.rmtree(pathlib.Path(scratch)) # delete scratch subdirectory

def save_mosaics_as_fasta (operon_seqs, output, mosaics):
    for i, m in enumerate(mosaics):
        ofile = f"{output}.seq-{m}.fasta.xz"
        counter = 0
        with open_anyformat (ofile, "w") as f:
            for rec in operon_seqs:
                if rec.description.split()[1] == m: # header has format "> genomeID mosaic description"
                    f.write (str(f">{rec.description}\n{rec.seq}\n").encode())
                    counter += 1
        if i < 5:
            logger.info (f"Succesfully saved {counter} operons to {ofile}")
        elif i == 6:
            logger.info (f"etc... (this may take a while)")

def save_mosaic_frequency (moscounter, output):
    ofile = f"{output}.mosaics.tsv"
    with open_anyformat (ofile, "w") as f:
        f.write (str(f"mosaic\tfrequency\n").encode())
        for k, v in moscounter.most_common():
            f.write (str(f"{k}\t{v}\n").encode())

def extract_and_save_operons (pool_info, fastadir, intergenic_space=1000, short_operon=1000, border=50):
    coord_df, merge_df, fname = pool_info
    genome_list = coord_df["seqid"].unique().tolist()
    fw = open_anyformat (fname, "w")

    def operon_from_coords (genome_sequence, coord_df):
        coord_df = coord_df.sort_values(by=["start"], ascending=True)
        coord_df = coord_df.reset_index(drop=True) # unlike iterate(), iterrows() returns index (before sorting)
        # remember that GFF is one-based, but recent version of phylobarcode stores coordinates as zero-based
        coord_df["start"] = coord_df["start"].astype(int)
        coord_df["end"] = coord_df["end"].astype(int)
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
            name = f">{g} {opr} " + mdf["fasta_description"].iloc[0]
            fw.write (str(f"{name}\n{seq}\n").encode())
            mosaics.append (opr)
    fw.close()
    return mosaics

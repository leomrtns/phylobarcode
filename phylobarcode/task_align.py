#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json? imports itertools, pathlib
import pandas as pd, numpy as np
import io, multiprocessing, shutil, gffutils, json, collections
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
logger = logging.getLogger("phylobarcode_global_logger")

def cluster_align_gene_files (genefiles = None, output = None, nthreads = 1, threshold = None, csvfile = None, scratch = None):
    if genefiles is None:
        logger.error("No gene files provided")
        sys.exit(1)
    hash_name = '%012x' % random.randrange(16**12)
    if scratch is None:
        scratch = f"scratch.{hash_name}"
        logger.warning(f"No scratch directory provided, using {scratch}")
    if output is None:
        output = f"phyloaln.{hash_name}"
        logger.warning(f"No output file (prefix) provided, using {output}")
    if threshold is None: threshold = 0.99
    if threshold < 0.0001: threshold = 0.0001
    if threshold > 1: threshold = 1

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
        cluster_align_each_gene (short, long, output, scratch, taxon_df, threshold, nthreads)

    if scratch_created:
        shutil.rmtree(pathlib.Path(scratch))

def cluster_align_each_gene (shortname, genefile, outfile, scratch, taxon_df, threshold, nthreads):
    fas = read_fasta_as_list (genefile)
    if len(fas) < 4:
        logger.warning (f"{genefile} has fewer than 4 sequences, skipping")
        return
    logger.info(f"Read {shortname} gene (file {genefile}) with {len(fas)} sequences")
    seqinfo = {}
    for x in fas:
        seqid = x.id.split("|")[0] # remove gene name e.g. ">NZ_CP028136.1|S31"
        if taxon_df is not None:
            tx = taxon_df.loc[taxon_df['seqid'] == seqid].iloc[0]
            seqinfo[x.id] = {
                    "seqid": seqid, 
                    "species": tx['species'], 
                    "genus": tx['genus'],
                    "family": tx['family'], 
                    "order": tx['order']}
        else:
            tx = x.description.split(" ", 1)[1].split("|")[1:5] # |order|family|genus|species|
            seqinfo[x.id] = {
                    "seqid": seqid,
                    "species": tx[3],
                    "genus": tx[2],
                    "family": tx[1],
                    "order": tx[0]}
    n_species = len(set([x['species'] for x in seqinfo.values()]))
    if n_species < 4:
        logger.warning (f"{genefile} represents fewer than 4 species, skipping")
        return
    logger.info(f"Found {n_species} species in {shortname} gene")
    cdfile = f"{scratch}/{shortname}.cdhit"
    alnfile = f"{outfile}.{shortname}.aln"
    # cdhit_cluster_seqs() will recreate (uncompressed) fasta file since infile may be compressed
    rep_seqs, clusters = cdhit_cluster_seqs (sequences = fas, outfile=cdfile, nthreads=nthreads, id = threshold)
    logger.info(f"CD-HIT clustered {shortname} gene into {len(clusters)} clusters")
    rep_aln = mafft_align_seqs (infile=cdfile, outfile=alnfile, nthreads=nthreads)
    logger.info(f"Finished MAFFT alignment of {shortname}")
    os.remove(cdfile)
    n_species = [len(set([seqinfo[x]["species"] for x in cl])) for cl in clusters if len(cl) > 1]
    n_species = len([x for x in n_species if x > 1])
    n_genus   = [len(set([seqinfo[x]["genus"]   for x in cl])) for cl in clusters if len(cl) > 1]
    n_genus   = len([x for x in n_genus if x > 1])
    logger.info(f"Clustered {shortname} gene into {len(clusters)} clusters, with {n_species} with >1 species, and {n_genus} with >1 genus")

#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json? imports itertools, pathlib
import pandas as pd, numpy as np
import io, multiprocessing, shutil, gffutils, json, collections
from Bio.Blast import NCBIXML
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
logger = logging.getLogger("phylobarcode_global_logger")

def cluster_align_gene_files (genefiles = None, output = None, nthreads = 1, threshold = None, tsvfile = None, scratch = None):
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

    if tsvfile is not None:
        taxon_df = pd.read_csv(tsvfile, sep='\t', header=0)
        taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, replace="unknown")
    else:
        taxon_df = None

    shortname = remove_prefix_suffix (genefiles)
    tbl = []
    for short, long in zip(shortname, genefiles):
        tbl_row = cluster_align_each_gene (short, long, output, scratch, taxon_df, threshold, nthreads)
        tbl.append(tbl_row)

    if scratch_created:
        shutil.rmtree(pathlib.Path(scratch))

    df = pd.DataFrame(tbl, dtype=object)
    ofilename = f"{output}.cluster.tsv"
    df.to_csv(f"{ofilename}", sep='\t', index=False)
    logger.info(f"statistics save to {ofilename}")

def get_seqinfo_from_sequence_header (xid, xdescription, taxon_df = None):
    seqid = xid.split("|")[0] # remove gene name e.g. ">NZ_CP028136.1|S31"
    if taxon_df is not None:
        tx = taxon_df.loc[taxon_df['seqid'] == seqid].iloc[0]
        seqinfo = {
                "seqid": seqid, 
                "species": tx['species'], 
                "genus": tx['genus'],
                "family": tx['family'], 
                "order": tx['order']}
    else:
        tx = xdescription.split(" ", 1)[1].split("|")[1:5] # |order|family|genus|species|
        seqinfo = {
                "seqid": seqid,
                "species": tx[3],
                "genus": tx[2],
                "family": tx[1],
                "order": tx[0]}
    return seqinfo

def cluster_align_each_gene (shortname, genefile, outfile, scratch, taxon_df, threshold, nthreads):
    fas = read_fasta_as_list (genefile)
    if len(fas) < 4:
        logger.warning (f"{genefile} has fewer than 4 sequences, skipping")
        return stats
    logger.info(f"Read {shortname} gene (file {genefile}) with {len(fas)} sequences")
    seqinfo = {}
    for x in fas: seqinfo[x.id] = get_seqinfo_from_sequence_header (x.id, x.description, taxon_df)

    stats = {"gene": shortname, 
            "n_sequences": len(fas),
            "n_species": len(set([x['species'] for x in seqinfo.values()])),
            "n_genus": len(set([x['genus'] for x in seqinfo.values()]))}
    if stats["n_species"] < 4:
        logger.warning (f"{genefile} represents fewer than 4 species, skipping")
        return
    logger.info(f"Found {stats['n_species']} species in {shortname} gene")
    cdfile = f"{scratch}/{shortname}.cdhit"
    alnfile = f"{outfile}.{shortname}.aln"
    # cdhit_cluster_seqs() will recreate (uncompressed) fasta file since infile may be compressed
    rep_seqs, clusters = cdhit_cluster_seqs (sequences = fas, outfile=cdfile, nthreads=nthreads, id = threshold)
    logger.info(f"CD-HIT clustered {shortname} gene into {len(clusters)} clusters")
    n_species = [len(set([seqinfo[x]["species"] for x in cl])) for cl in clusters if len(cl) > 1]
    max_species = max(n_species)
    n_species = len([x for x in n_species if x > 1])
    n_genus   = [len(set([seqinfo[x]["genus"] for x in cl])) for cl in clusters if len(cl) > 1]
    max_genus = max(n_genus)
    n_genus   = len([x for x in n_genus if x > 1])
    stats = {**stats, **{
        "max_paraphyl_species": max_species,
        "max_paraphyl_genus": max_genus,
        "n_paraphyl_species": n_species, 
        "n_paraphyl_genus": n_genus, 
        "n_representatives": len(clusters)}}
    logger.info(f"Clustered {shortname} gene into {len(clusters)} clusters, with {n_species} with >1 species, and {n_genus} with >1 genus")
    seqlen = [len(x.seq) for x in rep_seqs]
    stats = {**stats, **{
        "min_seqlen": min(seqlen),
        "first_quart_seqlen": round(np.quantile(seqlen, 0.25),0),
        "mean_seqlen": round(np.mean(seqlen),0),
        "third_quart_seqlen": round(np.quantile(seqlen, 0.75),0),
        "max_seqlen": max(seqlen),
        }}

    # align representative sequences
    rep_aln = mafft_align_seqs (infile=cdfile, outfile=alnfile, nthreads=nthreads)
    if rep_aln:
       logger.info(f"Finished MAFFT alignment of {shortname}")
       stats["alignment_length"] = len(rep_aln[0].seq)
    else:
         logger.error(f"Failed to align {shortname}")
    os.remove(cdfile)

    return stats

### tree and silhouette 

def estimate_compare_trees (alnfiles = None, output = None, scratch = None, tsvfile = None, 
        gtdb_tree = None, prev_tsv = None, nthreads = 1):
    hash_name = '%012x' % random.randrange(16**12)
    if alnfiles is None:
        logger.error("No alignment files specified")
        sys.exit(1)
    if scratch is None:
        scratch = f"scratch.{hash_name}"
        logger.warning(f"No scratch directory provided, using {scratch}")
    if output is None:
        output = f"phylotree.{hash_name}"
        logger.warning(f"No output file (prefix) provided, using {output}")
    if not os.path.exists(scratch):
        pathlib.Path(scratch).mkdir(parents=True, exist_ok=True)
        scratch_created = True
    if tsvfile is not None:
        taxon_df = pd.read_csv(tsvfile, sep='\t', header=0)
        taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, replace="unknown")
    else:
        taxon_df = None
    if gtdb_tree is not None and taxon_df is not None:
        ref_tree = read_translate_gtdb_tree_dendropy (gtdb_tree, taxon_df)
    else:
        ref_tree = None
        logger.warning("No reference tree provided or taxon table, skipping tree comparison")

    shortname = remove_prefix_suffix (alnfiles)
    tbl = []
    for short, long in zip(shortname, alnfiles):
        tbl_row = generate_tree (short, long, output, scratch, ref_tree, taxon_df, nthreads)
        tbl.append(tbl_row)

    if scratch_created:
        shutil.rmtree(pathlib.Path(scratch))

    df = pd.DataFrame(tbl, dtype=object)
    if prev_tsv is not None:
        logger.info(f"Mergin with other (clustering) results from {prev_tsv}")
        prev_df = pd.read_csv(prev_tsv, sep='\t')
        df = df.merge(prev_df, how='outer', on='gene')

    ofilename = f"{output}.treestats.tsv"
    df.to_csv(f"{ofilename}", sep='\t', index=False)
    logger.info(f"statistics save to {ofilename}")

def read_translate_gtdb_tree_dendropy (treefile, taxon_df): # both have to be present (i.e. not None)
    tree = dendropy.Tree.get_from_path(treefile, schema="newick", preserve_underscores=True)
    tree.retain_taxa_with_labels(taxon_df["gtdb_accession"].unique())
    for node in tree.postorder_node_iter():
        node.label = None # gtdb has internal node annotations (which in dendropy are distinct from node.taxon.label)
    for tx in tree.taxon_namespace:
        tx.label = taxon_df[taxon_df["gtdb_accession"] == tx.label]["seqid"].values[0]
    return tree

def generate_tree (shortname, alnfile, output, scratch, reference_tree, taxon_df, nthreads):
    treefile = f"{output}.{shortname}.tre"

    seqinfo = read_fasta_headers_as_list (alnfile)
    seqinfo = [x[1:].split(" ", 1) for x in seqinfo] # remove leading ">" and split id and description
    seqinfo = {x[0]:get_seqinfo_from_sequence_header (x[0], x[1], taxon_df) for x in seqinfo}
    logger.info(f"Read seqinfo from {len(seqinfo)} sequences in {shortname} alignment")
    
    if not os.path.exists(treefile):
        logger.info(f"Generating tree for {shortname}")
        treestring = newick_string_from_alignment (infile=alnfile, outfile=treefile, simple_names = True, nthreads=1)
    else:
        logger.info(f"Tree file {treefile} already exists, no estimation needed")
        treestring = open(treefile).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
  
    logger.info(f"Now comparing gene and reference trees")
    ## remember that tree labels have gene name like ">NZ_CP032229.1|L7"
    commontaxa = dendropy.TaxonNamespace([x for x in seqinfo.keys()])
    tree = dendropy.Tree.get_from_string (treestring, schema="newick", 
            preserve_underscores=True, taxon_namespace=commontaxa)
    if reference_tree is not None:
        # reduce reftree to genetree
        reftree = copy.deepcopy(reference_tree) 
        reftree.retain_taxa (commontaxa)
        # reduce genetree to reftree
        tree.retain_taxa (reftree.taxon_namespace)
        commontaxa = tree.taxon_namespace
        reftree = dendropy.Tree.get (data=reftree.as_string(schema="newick"), schema="newick",
                preserve_underscores=True, taxon_namespace=commontaxa) # to make sure taxon_namespace is the same
        tree = dendropy.Tree.get (data=tree.as_string(schema="newick"), schema="newick",
                preserve_underscores=True, taxon_namespace=commontaxa)
        ref_only, gene_only = dendropy.calculate.treecompare.false_positives_and_negatives(reftree, tree, is_bipartitions_updated=False)
        # store unnormalised tree as string
        rtre_str = reftree.as_string(schema="newick", suppress_rooting=True, suppress_edge_lengths=False)
        r_blens = [x.edge_length for x in reftree.postorder_node_iter() if x.edge_length is not None]
        stats = {
            "gene": shortname,
            "ref_only": ref_only,
            "gene_only": gene_only,
            "common_taxa": len(commontaxa), ## this may be smaller than the number of seqs if are missing from gtdb
            "ref_phylodiversity": sum(r_blens),
            }
    else:
        # brlens stats and normalise trees (first copy strings, which are used in silhouette)
        gtre_str = tree.as_string(schema="newick", suppress_rooting=True, suppress_edge_lengths=False)
        g_blens = [x.edge_length for x in tree.postorder_node_iter() if x.edge_length is not None]
        stats = {
            "gene": shortname,
            "gene_phylodiversity": sum(g_blens),
            }

    ## calc silhouette (use dict for both genetree and reftree)
    class_dict_sp = {x:y["species"] for x,y in seqinfo.items()}
    class_dict_genus = {x:y["genus"] for x,y in seqinfo.items()}
    stats["gene_sscore_species"] = silhouette_score_from_newick_swift (gtre_str, class_dict_sp)
    stats["gene_sscore_genus"] = silhouette_score_from_newick_swift (gtre_str, class_dict_genus)

    if reference_tree is not None:
        for x in reftree.postorder_node_iter():
            if x.edge_length is not None:
                x.edge_length /= stats["ref_phylodiversity"]
        for x in tree.postorder_node_iter():
            if x.edge_length is not None:
                x.edge_length /= stats["gene_phylodiversity"]
        d = dendropy.calculate.euclidean_distance_matrix(reftree, tree)
        stats["blen_distance"] = d
        stats["ref_sscore_species"]  = silhouette_score_from_newick_swift (rtre_str, class_dict_sp)
        stats["ref_sscore_genus"]  = silhouette_score_from_newick_swift (rtre_str, class_dict_genus)

    return stats

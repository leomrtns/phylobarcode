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
        if tbl_row is not None: tbl.append(tbl_row)

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
        return None
    logger.info(f"Read {shortname} gene (file {genefile}) with {len(fas)} sequences")
    seqinfo = {}
    for x in fas: seqinfo[x.id] = get_seqinfo_from_sequence_header (x.id, x.description, taxon_df)

    stats = {"gene": shortname, 
            "n_sequences": len(fas),
            "n_species": len(set([x['species'] for x in seqinfo.values()])),
            "n_genus": len(set([x['genus'] for x in seqinfo.values()]))}
    if stats["n_species"] < 4:
        logger.warning (f"{genefile} represents fewer than 4 species, skipping")
        return None
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
        gtdb_tree = None, prev_tsv = None, rapidnj = None, nthreads = 1):
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
    if rapidnj is None: rapidnj = False
    if rapidnj is not False: rapidnj = True
    if not os.path.exists(scratch):
        pathlib.Path(scratch).mkdir(parents=True, exist_ok=True)
        scratch_created = True
    if tsvfile is not None:
        logger.info(f"Reading taxonomic information from {tsvfile}")
        taxon_df = pd.read_csv(tsvfile, sep='\t', header=0)
        taxon_df = split_gtdb_taxonomy_from_dataframe (taxon_df, replace="unknown")
    else:
        taxon_df = None
    if gtdb_tree is not None and taxon_df is not None:
        logger.info(f"Reading GTDB tree from {gtdb_tree}")
        ref_tree = read_translate_gtdb_tree_dendropy (gtdb_tree, taxon_df)
        ref_tree.write_tree_newick (f"gtdb.tree") ## treeswift 
        #ref_tree.write(path=f"gtdb.tree", schema="newick") ## dendropy
    else:
        ref_tree = None
        logger.warning("No reference tree provided or taxon table, skipping tree comparison")

    shortname = remove_prefix_suffix (alnfiles)
    tbl = []
    for short, long in zip(shortname, alnfiles):
        tbl_row = generate_tree (short, long, output, scratch, ref_tree, taxon_df, rapidnj, nthreads)
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
    uniq_taxa = taxon_df["gtdb_genome_representative"].unique() # taxon_df["gtdb_accession"]
    logger.info(f"Read tree with {len(tree.taxon_namespace)} leaves; will now remove internal node annotations")
    for node in tree.postorder_node_iter():
        node.label = None # gtdb has internal node annotations (which in dendropy are distinct from node.taxon.label)
    ## treeswift is much faster than dendropy
    swtree = treeswift.read_tree_dendropy(tree) # faster than dendropy
    logger.info(f"Will now reduce tree to {len(uniq_taxa)} known genomes")
    swtree = swtree.extract_tree_with (uniq_taxa)

    logger.info(f"Reduced tree to {swtree.num_nodes(internal=False)} taxa; will now map representative leaves to genome names")
    leaf_map = swtree.label_to_node(selection="leaves") ## cannot use traverse_leaves() since new leaves confuse the iterator
    for lab, node in leaf_map.items():
        if lab == "": continue ## treeswift sometimes thinks an internal node is a leaf
        newlabel = list(set(taxon_df[taxon_df["gtdb_genome_representative"] == lab]["seqid"].tolist()))
        if len(newlabel) == 1:
            node.label = newlabel[0]
        else:
            node.set_label("") # becomes internal node
            for newlab in newlabel:
                node.add_child(treeswift.Node(label=newlab, edge_length=0))

    logger.info(f"Expanded tree has {swtree.num_nodes(internal=False)} leaves (genomes)")
    #return dendropy.Tree.get_from_string(swtree.newick(), schema="newick", preserve_underscores=True)
    return swtree

def generate_tree (shortname, alnfile, output, scratch, reference_tree, taxon_df, rapidnj, nthreads):
    treefile = f"{output}.{shortname}.tre"
    seqinfo = read_fasta_headers_as_list (alnfile)
    seqinfo = [x.split(" ", 1) for x in seqinfo] #  split id and description
    seqinfo = {x[0]:get_seqinfo_from_sequence_header (x[0], x[1], taxon_df) for x in seqinfo}
    logger.info(f"Read seqinfo from {len(seqinfo)} sequences in {shortname} alignment")

    def stats_silhouette (treestring, class_dict):
        try:
            labdic = silhouette_score_from_newick_swift (treestring, class_dict)
            vals = [x for x in labdic.values()]
            return [np.quantile(vals, 0.01), np.quantile(vals, 0.05), np.quantile(vals, 0.5)]
        except:
            return [None, None, None]
    
    if not os.path.exists(treefile):
        logger.info(f"Generating tree for {shortname}")
        gtre_str = newick_string_from_alignment (infile=alnfile, outfile=treefile, 
                rapidnj = rapidnj, simple_names = True, nthreads=nthreads)
    else:
        logger.info(f"Tree file {treefile} already exists, no estimation needed")
        gtre_str = open(treefile).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")

    ## remember that tree labels have gene name like ">NZ_CP032229.1|L7"
    for k,v in seqinfo.items():
        gtre_str = gtre_str.replace(k, v["seqid"])

    if reference_tree is None:
        logger.info(f"Calculating silhouette scores (reference tree not provided)")
    else:
        logger.info(f"Comparing gene and reference trees, and calculating silhouette scores")

    # brlens stats and normalise trees (keep copy gtre_str, which is used in silhouette)
    tree = treeswift.read_tree_newick (gtre_str)
    g_blens = [x.edge_length for x in tree.traverse_postorder() if x.edge_length is not None]
    stats = {
        "gene": shortname,
        "gene_phylodiv_full": sum(g_blens),
        }

    if reference_tree is not None:
        c1 = set([x.label for x in reference_tree.traverse_leaves()])
        c2 = set([x.label for x in tree.traverse_leaves()])
        commontaxa = list(c1.intersection(c2)) # debug note: offending ref taxa not here (NC_011891.1)
        # reduce reftree to genetree
        reftree = reference_tree.extract_tree(None, False, False) ## deep copy
        reftree = reftree.extract_tree_with (commontaxa)
        # reduce genetree to reftree
        tree = tree.extract_tree_with (commontaxa)
        # normalise tree lengths
        g_blens = [x.edge_length for x in tree.traverse_postorder() if x.edge_length is not None]
        r_blens = [x.edge_length for x in reftree.traverse_postorder() if x.edge_length is not None]
        # store unnormalised tree as string and save to file
        rtre_str = reftree.newick()
        outreffile = f"{output}.{shortname}.ref.tre"
        reftree.write_tree_newick (outreffile)
        # use dendropy to calc false pos and negs (RF distance)
        ct = dendropy.TaxonNamespace(commontaxa) # same namespace must be used for both trees
        rdendro = dendropy.Tree.get_from_string(reftree.newick(), schema="newick", taxon_namespace=ct, preserve_underscores=True)
        gdendro = dendropy.Tree.get_from_string(tree.newick(), schema="newick", taxon_namespace=ct, preserve_underscores=True)
        r_only, g_only = dendropy.calculate.treecompare.false_positives_and_negatives(rdendro, gdendro, is_bipartitions_updated=False)
        stats = {**stats, **{
            "ref_only": r_only,
            "gene_only": g_only,
            "common_taxa": len(commontaxa), ## this may be smaller than the number of seqs if are missing from gtdb
            "ref_phylodiversity": sum(r_blens),
            "gene_phylodiversity": sum(g_blens),
            }}
        logger.info(f"Reference tree has {len(commontaxa)} leaves in common with gene tree, saved to {outreffile}")

    ## calc silhouette using treestrings (use dict for both genetree and reftree)
    class_dict_sp = {y["seqid"]:y["species"] for y in seqinfo.values()}
    class_dict_genus = {y["seqid"]:y["genus"] for y in seqinfo.values()}
    sp_stats = stats_silhouette (gtre_str, class_dict_sp)
    ge_stats = stats_silhouette (gtre_str, class_dict_genus)
    stats = {**stats, **{
        "gene_sscore_species_1pct": sp_stats[0],
        "gene_sscore_species_5pct": sp_stats[1],
        "gene_sscore_species_median": sp_stats[2],
        "gene_sscore_genus_1pct": ge_stats[0],
        "gene_sscore_genus_5pct": ge_stats[1],
        "gene_sscore_genus_median": ge_stats[2],
        }}
    logger.info(f"Silhouette scores for {shortname} calculated")

    if reference_tree is not None: # gdendro and tdendro are dendropy trees, used in Euclidean distance
        for x in rdendro.postorder_node_iter():
            if x.edge_length is not None:
                x.edge_length /= stats["ref_phylodiversity"]
        for x in gdendro.postorder_node_iter():
            if x.edge_length is not None:
                x.edge_length /= stats["gene_phylodiversity"]
        stats["blen_distance"] =  dendropy.calculate.treecompare.euclidean_distance(rdendro, gdendro)
        sp_stats = stats_silhouette (rtre_str, class_dict_sp)
        ge_stats = stats_silhouette (rtre_str, class_dict_genus)
        stats = {**stats, **{
            "ref_sscore_species_1pct": sp_stats[0],
            "ref_sscore_species_5pct": sp_stats[1],
            "ref_sscore_species_median": sp_stats[2],
            "ref_sscore_genus_1pct": ge_stats[0],
            "ref_sscore_genus_5pct": ge_stats[1],
            "ref_sscore_genus_median": ge_stats[2],
            }}

    return stats


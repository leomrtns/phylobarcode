import os, logging, argparse, sys, pathlib, multiprocessing, datetime, itertools, pathlib, random
from phylobarcode.__version__ import __version__

stream_log = logging.StreamHandler()
#log_format = logging.Formatter(fmt='phylobarcode___main %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d%H:%M") # now it's shared 
log_format = logging.Formatter(fmt='phylobarcode %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)

#logger = logging.getLogger(__name__) # this creates one logger per module, but I want one logger for the whole program
logger = logging.getLogger("phylobarcode_global_logger") # creates a named global logger
logger.addHandler(stream_log)
logger.propagate = False

defaults = {
    "current_dir": os.getcwd() + "/",
    "timestamp": datetime.datetime.now().strftime("%y%m%d_%H%M%S"),
    "nthreads": multiprocessing.cpu_count(),
    "scratch": 'phylobar.%016x' % random.randrange(16**16),
    "json_riboproteins": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/riboproteins_names.json"),
    "json_ribogenes": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/ribogenes_names.json"),
    "json_extragenes": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/extragenes_names.json")
    }

long_description = """
Phylobarcode is a tool to search and analyse long operons with phylogenetic signal.

It can search for potential long segments ('operons' or 'barcodes') which can potentially serve as phylogenetic markers. 
It can then search for potential primers for these markers, which can be used in the lab to amplify the operon.
And it can be used to reconstruct the evolutionary history of the amplicons (not implemented yet).
"""
epilogue = """
SPDX-License-Identifier: GPL-3.0-or-later; Copyright (C) 2022-today Leonardo de Oliveira Martins
https://github.com/quadram-institute-bioscience/phylobarcode\n
"""

def generate_prefix_for_task (args, taskname):
    if args.prefix:
        args.prefix = os.path.join(defaults["current_dir"], args.prefix) # "current_dir" = args.outdir, buth updated with user choice
    else:
        args.prefix = os.path.join(defaults["current_dir"], f"pb.{defaults['timestamp']}_{taskname}")

def run_merge_fasta_gff (args):
    from phylobarcode import task_fasta_gff
    generate_prefix_for_task (args, "fastagff")
    if not args.nthreads: args.nthreads = defaults["nthreads"]
    task_fasta_gff.merge_fasta_gff (fastadir=args.fasta, gffdir=args.gff, fasta_tsvfile = args.tsv_fasta, 
            gff_tsvfile = args.tsv_gff, scratch=args.scratch, gtdb = args.gtdb, output=args.prefix, nthreads = args.nthreads)

def run_extract_coordinates_from_gff (args):
    from phylobarcode import task_extract_riboprot_gff
    generate_prefix_for_task (args, "coordinates")
    if not args.nthreads: args.nthreads = defaults["nthreads"]
    jsonfiles = {"riboproteins": defaults["json_riboproteins"], "ribogenes": defaults["json_ribogenes"], "extragenes": defaults["json_extragenes"]}
    task_extract_riboprot_gff.extract_coordinates_from_gff (tsvfile=args.tsv, gffdir = args.gff, output=args.prefix, 
            jsonfiles = jsonfiles, coord_tsvfile = args.coords, nthreads=args.nthreads, scratch=args.scratch)

def run_extract_operons_from_fasta (args):
    from phylobarcode import task_extract_riboprot_fasta
    generate_prefix_for_task (args, "operons")
    if not args.nthreads: args.nthreads = defaults["nthreads"]

    task_extract_riboprot_fasta.extract_operons_from_fasta (coord_tsvfile = args.coords, merge_tsvfile = args.tsv, 
            fastadir=args.fasta, output=args.prefix, intergenic_space = args.intergenic, short_operon = args.short,
            most_common_mosaics = args.most_common, border = args.border, riboprot_subset = args.subset,
            nthreads=args.nthreads, scratch=args.scratch)

def run_extract_genes_from_fasta (args):
    from phylobarcode import task_extract_riboprot_fasta
    generate_prefix_for_task (args, "genes")
    if not args.nthreads: args.nthreads = defaults["nthreads"]
    task_extract_riboprot_fasta.extract_genes_from_fasta (coord_tsvfile = args.coords, merge_tsvfile = args.tsv,
            fastadir=args.fasta, output=args.prefix, scratch=args.scratch, keep_paralogs = args.paralogs, 
            nthreads=args.nthreads)

def run_find_primers (args):
    from phylobarcode import task_find_primers
    generate_prefix_for_task (args, "primers")
    if args.nthreads and args.nthreads < 2:
        logger.info("Single-threaded mode requested by user")
        task_find_primers.find_primers (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, 
                num_return=args.n_primers, output=args.prefix)
        return
    if defaults["nthreads"] < 2:
        logger.warning("Multiprocessing not available, falling back to single-threaded mode")
        task_find_primers.find_primers (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, 
                num_return=args.n_primers, output=args.prefix)
        return

    if not args.nthreads: 
        args.nthreads = defaults["nthreads"]
        logger.info(f"{args.nthreads} threads are available (actual pool may be smaller)")
    else:
        logger.info(f"{args.nthreads} threads were requested by user (actual pool may be smaller)")
    task_find_primers.find_primers_parallel (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, 
            num_return=args.n_primers, taxon=args.taxon, output=args.prefix, nthreads=args.nthreads)
    return

def run_cluster_flanks (args):
    from phylobarcode import task_cluster
    generate_prefix_for_task (args, "flank")
    if not args.nthreads: args.nthreads = defaults["nthreads"]
    task_cluster.cluster_flanks_from_fasta (fastafile=args.fasta, output=args.prefix, border=args.border,
            identity=args.id, nthreads=args.nthreads, min_samples = args.min_samples, scratch=args.scratch)
    return

def run_cluster_primers (args):
    from phylobarcode import task_cluster
    generate_prefix_for_task (args, "cluster")
    if not args.nthreads: args.nthreads = defaults["nthreads"]
    if len(args.tsv) < 2: ## "nargs='+'" always returns a list of at least one element (or None, but here it's positional)
        task_cluster.cluster_primers_from_tsv (tsv=args.tsv[0], output=args.prefix, 
                subsample=args.subsample, threshold = args.threshold, scratch=args.scratch,
                nthreads=args.nthreads)
        return
    uniq = remove_prefix_suffix (args.tsv)
    for infile, outfile in zip (args.tsv, uniq):
        task_cluster.cluster_primers_from_tsv (tsv=infile, output=f"{args.prefix}_{outfile}", 
                subsample=args.subsample, threshold = args.threshold, scratch=args.scratch,
                nthreads=args.nthreads)
    return

def run_subsample_primers (args):
    from phylobarcode import task_cluster
    generate_prefix_for_task (args, "subsample")
    if len(args.tsv) < 2:
        task_cluster.subsample_primers_from_tsv (tsv=args.tsv[0], output=args.prefix, 
                subsample=args.subsample, n_elements = args.n_primers)
        return
    uniq = remove_prefix_suffix (args.tsv)
    for infile, outfile in zip (args.tsv, uniq):
        task_cluster.subsample_primers_from_tsv (tsv=infile, output=f"{args.prefix}_{outfile}", 
                subsample=args.subsample, n_elements = args.n_primers)
    return

def run_blast_primers (args):
    from phylobarcode import task_blast_primers
    generate_prefix_for_task (args, "blast")
    if not args.nthreads: 
        args.nthreads = defaults["nthreads"]
        logger.info(f"{args.nthreads} threads are available (actual pool may be smaller)")
    else:
        logger.info(f"{args.nthreads} threads were requested by user (actual pool may be smaller)")
    if args.accurate: task = "blastn-short"
    else: task = "blastn"
    if args.max_target_seqs > 1e9: args.max_target_seqs = 1e9
    if args.max_target_seqs < 1: args.max_target_seqs = 1
    if len(args.tsv) < 2:
        task_blast_primers.blast_primers_from_tsv (tsv=args.tsv[0], output=args.prefix, 
                database=args.database, evalue=args.evalue, task=task, 
                max_target_seqs=args.max_target_seqs, taxon=args.taxon, nthreads=args.nthreads)
        return
    uniq = remove_prefix_suffix (args.tsv)
    for infile, outfile in zip (args.tsv, uniq):
        task_blast_primers.blast_primers_from_tsv (tsv=infile, output=f"{args.prefix}_{outfile}", 
                database=args.database, evalue=args.evalue, task=task, 
                max_target_seqs=args.max_target_seqs, taxon=args.taxon, nthreads=args.nthreads)
    return

def pair_blast_primers (args):
    from phylobarcode import task_pair_primers
    generate_prefix_for_task (args, "paired")
    if len(args.tsv) != 2:
        logger.error("Exactly two tsv files must be provided (left and right, in order)")
        return
    uniq = remove_prefix_suffix (args.tsv)
    output = [f"{args.prefix}_{outfile}" for outfile in uniq]
    if not args.nthreads: 
        args.nthreads = defaults["nthreads"]
        logger.info(f"{args.nthreads} threads are available (actual pool may be smaller)")
    else:
        logger.info(f"{args.nthreads} threads were requested by user (actual pool may be smaller)")
    task_pair_primers.pair_primers_from_raw_tsv (tsv=args.tsv, output=output, taxon=args.taxon, nthreads=args.nthreads)
    return

def run_select_primers (args):
    from phylobarcode import task_blast_primers
    generate_prefix_for_task (args, "selected")
    if len(args.tsv) < 2:
        task_blast_primers.select_primers_from_tsv (tsv=args.tsv[0], subsample=args.subsample, 
                n_elements = args.n_primers, output=args.prefix)
        return
    uniq = remove_prefix_suffix (args.tsv)
    for infile, outfile in zip (args.tsv, uniq):
        task_blast_primers.select_primers_from_tsv (tsv=infile, subsample=args.subsample, 
                n_elements = args.n_primers, output=f"{args.prefix}_{outfile}")
    return

def remove_prefix_suffix (strlist):
    def all_same(x):  # https://stackoverflow.com/a/6719272/204903
        return all(x[0] == y for y in x)
    char_tuples = zip(*strlist)
    fix_tuples  = itertools.takewhile(all_same, char_tuples)
    prefix = ''.join(x[0] for x in fix_tuples)
    inverse = [x[::-1] for x in strlist]
    char_tuples = zip(*inverse)
    fix_tuples  = itertools.takewhile(all_same, char_tuples)
    suffix = ''.join(x[0] for x in fix_tuples)
    suffix = suffix[::-1]

    l_pre = len(prefix) ## we could skip "prefix" and store lenght of fix_tuples but this is more readable
    l_suf = len(suffix)
    return [x[l_pre:len(x)-l_suf] for x in strlist] # does not work well for 'lefT' and 'righT' 

class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help(sys.stderr)
        sys.exit(2)

def main():
    parent_parser = ParserWithErrorHelp(description=long_description, add_help=False)
    parent_group = parent_parser.add_argument_group('Options common to all commands (some commands may not use them)')
    parent_group.add_argument('--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.INFO, 
            help="Print debugging statements (most verbose)")
    parent_group.add_argument('--silent', action="store_const", dest="loglevel", const=logging.WARNING, 
            help="Less verbose output (only warnings and errors)")
    parent_group.add_argument('--nthreads', metavar='int', type=int, 
            help="Number of threads requested (default = maximum available)")
    parent_group.add_argument('--outdir', metavar="<dir/>", action="store", 
            help="Output directory. Default: working directory")
    parent_group.add_argument('--scratch', action="store", 
            help="Existing scratch directory (i.e. path must exist). Default: working directory")
    parent_group.add_argument('--prefix', metavar="string", 
            help="optional file name --- just the prefix, since suffixes (a.k.a. extensions) are added by the program")
    parent_group.add_argument('--version', action='version', version=f"%(prog)s {__version__}") ## called with subcommands

    # alternative to subp= parent_parser.add_subparsers(dest='command', description=None, title="Commands")
    main_parser = ParserWithErrorHelp(description=long_description, formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    main_parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}") ## called without  subcommands (grouped together with --help)
    subp= main_parser.add_subparsers(dest='Commands', description=None, title="Commands", required=True)

    # TODO: deduplicate (using genus+sourmash)
    this_help = "Given one folder with fasta, one with GFF files of reference genomes, and a GTDB metadata file, creates a table with matches. "
    extra_help= '''\n
    The fasta and GFF files are recognised by their extensions (fasta, fna, faa, gff, gff3) with optional compression.
    You can add a TSV file with GFF and fasta info from a previous run (e.g. to add more genomes or to include the GTDB
    taxonomy and thus generate a TSV table with the matches).
    This program works without a GTDB metadata file, in which case it will only generate only the individual fasta and GFF3
    tables. However the information from the GTDB metadata file is required in downstream analyses, and thus the final
    table with matches is only generated if a GTDB metadata file is provided.
    Notice that we use TSV instead of CSV because the latter is not compatible with NCBI taxonomy shenanigans (commas in names).
    '''
    up_findp = subp.add_parser('merge_fasta_gff', help=this_help, description=this_help + extra_help, parents=[parent_parser], 
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('-a', '--fasta', metavar="<dir>", required=True, 
            help="directory where fasta genomic files can be found (required)") # technically not needed if gff3 contains fasta
    up_findp.add_argument('-g', '--gff',   metavar="<dir>", required=True, help="directory with GFF3 files (required)")
    up_findp.add_argument('-A', '--tsv_fasta', metavar="tsv", help="tsv file with fasta entries from previous run (optional)")
    up_findp.add_argument('-G', '--tsv_gff',   metavar="tsv", help="tsv file with GFF entries from previous run (optional)")
    up_findp.add_argument('-d', '--gtdb', metavar="tsv[.xz,.gz]", help="GTDB metadata tsv file (optional but needed downstream)")
    up_findp.set_defaults(func = run_merge_fasta_gff)

    this_help = "Given a table with matches, extracts the coordinates of all riboproteins from GFF files"
    extra_help= '''\n
    This program extracts the annotations (coordinates info) from the GFF files with a match and GTDB taxonomy. 
    The table with matches should have been generated by the "merge_fasta_gff" command with the exact same GFF
    directory, which should contain all files mentioned in the table.
    It does not need the fasta files, and it will store all needed info from the GFF files in a single TSV file.
    '''
    up_findp = subp.add_parser('extract_coordinates', help=this_help, description=this_help + extra_help, parents=[parent_parser],
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('tsv', help="tsv file with file matches between fasta and GFF3 (required)")
    up_findp.add_argument('-g', '--gff',   metavar="<dir>", required=True, help="directory with GFF3 files (required)")
    up_findp.add_argument('-c', '--coords', metavar="tsv", help="tsv file with coordinates from previous run (optional)")
    up_findp.set_defaults(func = run_extract_coordinates_from_gff)

    this_help = "Given riboprotein coordinates and table with fasta x GFF matches, extracts the riboprotein operons"
    extra_help= '''\n
    This program extracts the riboprotein operons from the fasta files given the coordinates of the riboproteins.
    The table with matches should have been generated by the "merge_fasta_gff" command with the exact same fasta
    directory, which should contain all files mentioned in the table.
    The coordinates of the riboproteins should have been generated by the "extract_riboprots_from_gff" command.
    You can tell to use just a subset of the riboprotein (and closeby) genes: "main", "hug", "core", "left", "leftleft", or "right".
    Any other name (e.g. "only") will exclude the non-riboprotein genes which are usually present in the operon.
    '''
    up_findp = subp.add_parser('extract_operons', help=this_help, description=this_help + extra_help, parents=[parent_parser],
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('tsv', help="tsv file with file matches between fasta and GFF3 (required)")
    up_findp.add_argument('-a', '--fasta', metavar="<dir>", required=True, 
            help="directory where fasta genomic files can be found (required)")
    up_findp.add_argument('-c', '--coords', metavar="tsv", required=True, 
            help="tsv file with riboprotein coordinates for all genomes (required)")
    up_findp.add_argument('-i', '--intergenic', metavar="int", default=1000, type=int,
            help="maximum distance between riboproteins to be considered part of the same operon (default: 1000)")
    up_findp.add_argument('-s', '--short', metavar="int", default=1500, type=int, 
            help="minimum number of riboproteins in an operon if <= 100 or minimum operon length if > 100 (default: 1500)")
    up_findp.add_argument('-m', '--most_common', metavar="int", default=100, type=int,
            help="number of most common operon mosaic (gene patterns) to be saved to fasta files (default: 100)")
    up_findp.add_argument('-S', '--subset', metavar="str", type=str, help="subset of riboproteins to be extracted (default: all)")
    up_findp.add_argument('-b', '--border', metavar="int", default=50, type=int,
            help="number of nucleotides to be added to the start and end of the operon (default: 50)")
    up_findp.set_defaults(func = run_extract_operons_from_fasta)

    this_help = "Given riboprotein coordinates and table with fasta x GFF matches, extracts individual genes"
    extra_help= '''\n
    This program extracts the riboprotein genes from the fasta files given the coordinates of the riboproteins.
    The table with matches should have been generated by the "merge_fasta_gff" command with the exact same fasta
    directory, which should contain all files mentioned in the table.
    If paralogs are included (i.e. all gene copies from each genome) then sequence IDs will include their number.
    The coordinates of the riboproteins should have been generated by the "extract_riboprots_from_gff" command.
    '''
    up_findp = subp.add_parser('extract_genes', help=this_help, description=this_help + extra_help, parents=[parent_parser],
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('tsv', help="tsv file with file matches between fasta and GFF3 (required)")
    up_findp.add_argument('-a', '--fasta', metavar="<dir>", required=True,
            help="directory where fasta genomic files can be found (required)")
    up_findp.add_argument('-c', '--coords', metavar="tsv", required=True,
            help="tsv file with riboprotein coordinates for all genomes (required)")
    up_findp.add_argument('-p', '--paralogs',  action="store_true", default = False,
            help="include all copies (paralogs) in fasta files (default=keep only one copy)")
    up_findp.set_defaults(func = run_extract_genes_from_fasta)

    this_help = "Find primers given a fasta file." # help is shown in "prog -h", description is shown in "prog this -h"
    extra_help= '''\n
    Runs `primer3_core` on each sequence in the alignment file, generating two tables, of left (l) and right (r) primers. Can use multiple threads.
    If file with taxonomic information is provided (from `merge_fasta_gff`), it will be used to calculate
    representativity of primers (called "taxon_diversity")
    '''
    up_findp = subp.add_parser('find_primers', help=this_help, description=this_help + extra_help, parents=[parent_parser], 
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('fasta', help="unaligned sequences")
    up_findp.add_argument('-l', '--length', metavar='int', type=int, 
            help="optimal primer size, in bp (default=20)")
    up_findp.add_argument('-b', '--border', metavar='int', type=int, default=400, 
            help="how far from sequence borders, in bp, we should look for primers (default=400)")
    up_findp.add_argument('-n', '--n_primers', metavar='int', type=int, default=100, 
            help="how many primers, per sequence, per end, should be returned (default=100)")
    up_findp.add_argument('-x', '--taxon', metavar='tsv',
            help="tsv file with taxonomy information (output from 'merge_fasta_gff')")
    up_findp.set_defaults(func = run_find_primers)

    this_help = "Extract and cluster flanking regions of operons where primers may be found."
    extra_help= '''\n
    The left and right flanking regions are short sequences (default=400 bp) which are extracted from the beginning (l)
    or end (r) of each sequence, and clustered using vsearch or OPTICS. Both original and clustered sequences are returned.
    OPTICS is used by default, but if you set an identity threshold (option `-i` or `--id`), vsearch is used instead.
    '''
    up_findp = subp.add_parser('get_flanks', help=this_help, description=this_help + extra_help, parents=[parent_parser], 
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('fasta', help="unaligned sequences")
    up_findp.add_argument('-b', '--border', metavar='int', type=int, default=400, 
            help="how far from sequence borders, in bp, we should look for primers (default=400)")
    up_findp.add_argument('-m', '--min_samples', type=int, default=5, 
            help="in OPTICS, minimum number of neighbours for sequence to be a core point (default=5 ;  should be larger than 2)")
    up_findp.add_argument('-i', '--id', type=float, default=None, 
            help="If defined, identity threshold for vsearch (default is to use OPTICS clustering instead of vsearch; suggested value is 0.5 ~ 0.9)")
    up_findp.set_defaults(func = run_cluster_flanks)

    this_help = "Cluster primers described in tsv file, adding a column to table with cluster number."
    extra_help= '''\n
    If several tsv files are given, default output files will keep their unique names (i.e. without common prefix or
    suffix). 

    If subsample percentage is given, then only "best" primers are considered from the beginning, based on genus/taxon
    frequency where primer was found, longest distance of primer to border, and penalty score as given by primer3.
    The subsample percentage should be close to 100 (percent) here. If you want to subsample _after_ clustering, then
    use `subsample_primers` afterwards. It is strongly suggested to cluster primers before subsampling (thus do not
    change the default value of `-s` or `--subsample` here unless you have too many primers to cluster).

    The cluster is based on vsearch, which is very fast but overestimates the number of clusters. Since it provides a
    profile for each cluster, we can use this to further cluster the profiles (as an extended consensus sequence) and
    try to reduce the number of clusters. Finally the profiles clusters are merged through single-linkage, and the final
    cluster IDs will be "id(single-likage)_id(profile)". Thus we may have `0_1` and `0_2` clusters, meaning that their
    profiles are `1` and `2`, but they are connected through single-linkage. The threshold controls the single linkage
    step, since it has little effect on the vsearch-based clusterings. This threshold can be a low value since it just
    marks similar clusters --- the total number of clusters will be given by the vsearch profile clustering.

    The input tsv files here should be the output of `find_primers`.
    '''
    up_findp = subp.add_parser('cluster_primers', help=this_help, description=this_help + extra_help, parents=[parent_parser], 
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('tsv', nargs="+", help="tsv files with primers (each ouput file from 'find_primers')")
#    up_findp.add_argument('-k', '--kmer', type=int, default=5, help="kmer size for primer clustering (default=5)")
    up_findp.add_argument('-t', '--threshold', type=float, default=0.8, help="threshold for primer clustering (default=0.8)")
    up_findp.add_argument('-s', '--subsample', metavar='float', type=float, default=100,
            help="subsample percentage of best primers to be clustered (default=100, i.e. all primers)")
    up_findp.set_defaults(func = run_cluster_primers)

    this_help = "Subsample primers from tsv file, keeping only the best ones based on clustering and taxonomic frequency."
    extra_help= '''\n
    The input tsv files here should be the output of `cluster_primers`. If possible (i.e. you have a manageable number
    of primers), use this command to subsample only _after_ clustering (since it can use the clustering information to
    maximilse the diversity of primers).
    This command can use a subsample percentage (with `-s` or `--subsample`) or a number of primers to keep (with `-n` or `--n_primers`).
    When using the percentage option, the univariate quantile of each variable is used, before finding the best primers
    according to _all_ variables. Thus the total number of primers will likely be (quite) smaller than the percentage given.
    If you just want the best N primers (using the same variables), then use the `-n` or `--n_primers` option.
    The variables used are the same as in `cluster_primers`, however the cluster ID is taken into account such that
    distinct clusters are preferred.  The other variables include the number of genera/taxa/sequences where primer was
    found, longest distance of primer to border, primer length, and penalty score as given by primer3.

    If both `-s` and `-n` are given, then first the quantiles are used to exclude primers, and from the remaining the best N are selected.

    '''
    up_findp = subp.add_parser('subsample_primers', help=this_help, description=this_help + extra_help, parents=[parent_parser],
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('tsv', nargs="+", help="tsv files with primers (each ouput file from 'cluster_primers')")
    up_findp.add_argument('-s', '--subsample', metavar='float', type=float, default=100,
            help="quantile percentage of best primers to be clustered (default=100, i.e. all primers)")
    up_findp.add_argument('-n', '--n_primers', metavar='int/float', type=float, default=None,
            help="number (or proportions if < 1) of primers to keep (default is to ignore this and use subsample percentage)")
    up_findp.set_defaults(func = run_subsample_primers)

    this_help = "Blast primers against database (in the future will check for left-right pairs)."
    extra_help= '''\n
    Needs exactly two tsv files, with left and right primers respectively. These tsv files should be the output of
    `find_primers`, `cluster_primers`, or `subsample_primers`. 
    Optionally you can give the merged tsv file from `merge_fasta_gff` to use taxonomy information of hits
    Output files will keep their unique names (i.e. without common prefix or suffix)
    Can use multiple threads within each BLAST job, or to create several concurrent BLAST jobs.
    '''
    up_findp = subp.add_parser('blast_primers', help=this_help, description=this_help + extra_help, parents=[parent_parser], 
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('tsv', nargs=2, 
            help="tsv files with left and right primers, respectively (out from 'find_primers', 'cluster_primers', or 'subsample_primers')")
    up_findp.add_argument('-d', '--database', required=True, 
            help="full path to database prefix")
    up_findp.add_argument('-e', '--evalue', type=float, default=1., 
            help="E-value threshold for blast (default=1, beware of low values)")
    up_findp.add_argument('-m', '--max_target_seqs', type=int, default=1000, 
            help="max number of hist per primer (default=1000; recommended to be close to number of genomes in blast DB)")
    up_findp.add_argument('-a', '--accurate', action="store_true", 
            help="use accurate blastn-short (default=regular blastn)")
    up_findp.add_argument('-x', '--taxon', metavar='tsv',
            help="tsv file with taxonomy information (output from 'merge_fasta_gff')")
    up_findp.set_defaults(func = run_blast_primers)

    this_help = "Select best primers based on blast results."
    extra_help= '''\n
    Can use a single tsv file or a pair (with left and right primers respectively). First it uses subsampling (`-s` or
    `--subsample`) on each of the variables it tracks (like number of unique families or genera hit by primer), to
    select only the "best" according to each variable. Then it sorts the resulting primers according to all variables
    and selectes the best `-n` or `--n_primers`.
    It also removes unused columns calculated by `blast_primers` (those related to "_all" hits, not just the best ones).
    The input files will be the `stats.tsv` files from `blast_primers`.
    '''
    up_findp = subp.add_parser('select_primers', help=this_help, description=this_help + extra_help, parents=[parent_parser],
            formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('tsv', nargs="+", help="tsv files with primers (each ouput file from 'blast_primers')")
    up_findp.add_argument('-s', '--subsample', metavar='float', type=float, default = 90,
            help="quantile percentage of best primers to be selected (default=90 percent)")
    up_findp.add_argument('-n', '--n_primers', metavar='int/float', type=float, default = None,
            help="number (or proportions if < 1) of primers to keep (default is to ignore this and use subsample percentage)")
    up_findp.set_defaults(func = run_select_primers)

    args = main_parser.parse_args()
    logger.setLevel(args.loglevel)
    for handler in logger.handlers: handler.setLevel(args.loglevel)

    if args.outdir: 
        defaults["current_dir"] = args.outdir = os.path.join(defaults["current_dir"], args.outdir)
    else:
        args.outdir = defaults["current_dir"]
    pathlib.Path(defaults["current_dir"]).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    if args.scratch:
        tmpdir = os.path.join(defaults["current_dir"], args.scratch)
        if not os.path.isdir(tmpdir):
            logger.error(f"Scratch directory {tmpdir} does not exist or it's not a directory")
            sys.exit(1)
        defaults["scratch"] = args.scratch = os.path.join(defaults["current_dir"], args.scratch, defaults["scratch"])
    else: # in any case, scratch will be a subdirectory of existing scratch area (which we can delete it with shutil.rmtree())
        defaults["scratch"] = args.scratch = os.path.join(defaults["current_dir"], defaults["scratch"])

    try:
        from multiprocessing import Pool
        from functools import partial
    except ImportError:
        logger.warning ("Multiprocessing not available, setting single thread as default")
        defaults["threads"] = 1
    if args.nthreads:
        if args.nthreads < 1: args.nthreads = 1
        if args.nthreads > defaults["nthreads"]:
            logger.warning(f"Requested {args.nthreads} threads, but only {defaults['nthreads']} are available. Using {defaults['nthreads']} instead.")
            args.nthreads = defaults["nthreads"]

    args.func(args) # calls task 

if __name__ == '__main__':
    main()

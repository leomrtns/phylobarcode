import os, logging, argparse, sys, pathlib, multiprocessing, datetime, itertools, pathlib, random
from phylobarcode.__version__ import __version__

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='phylobarcode___main %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

defaults = {
    "current_dir": os.getcwd() + "/",
    "timestamp": datetime.datetime.now().strftime("%y%m%d_%H%M%S"),
    "nthreads": multiprocessing.cpu_count(),
    "scratch": '%012x' % random.randrange(16**12),
#   "reference": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas") 
    }

long_description = """
Phylobarcode is a tool to search and analyse long operons with phylogenetic signal.

It can search for potential long segments ('operons' or 'barcodes') which can serve as phylogenetic markers. 
It can then search for potential primers which can be used to amplify the operon.
And it can be used to reconstruct the evolutionary history of the amplicons.
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
    task_fasta_gff.merge_fasta_gff (fastadir=args.fasta, gffdir=args.gff, scratch=args.scratch, output=args.prefix)

def run_find_primers (args):
    from phylobarcode import task_find_primers
    generate_prefix_for_task (args, "primers")
    if args.nthreads and args.nthreads < 2:
        logger.info("Single-threaded mode requested by user")
        task_find_primers.find_primers (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, num_return=args.n_primers, output=args.prefix)
        return
    try: # I could have used the same function with nthreads=1 but I already had this function to benchmark single threaded
        from multiprocessing import Pool
        from functools import partial
    except ImportError:
        logger.warning("Multiprocessing not available, falling back to single-threaded mode")
        task_find_primers.find_primers (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, num_return=args.n_primers, output=args.prefix)
        return

    if not args.nthreads: 
        args.nthreads = defaults["nthreads"]
        logger.info(f"{args.nthreads} threads are available (actual pool may be smaller)")
    else:
        logger.info(f"{args.nthreads} threads were requested by user (actual pool may be smaller)")
    task_find_primers.find_primers_parallel (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, num_return=args.n_primers, output=args.prefix, nthreads=args.nthreads)
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
    if len(args.csv) < 2: ## "nargs='+'" always returns a list of at least one element (or None, but here it's positional)
        task_cluster.cluster_primers_from_csv (csv=args.csv[0], output=args.prefix, nthreads=args.nthreads)
        return
    if not args.nthreads: args.nthreads = defaults["nthreads"]
    uniq = remove_prefix_suffix (args.csv)
    for infile, outfile in zip (args.csv, uniq):
        task_cluster.cluster_primers_from_csv (csv=infile, output=f"{args.prefix}_{outfile}", 
                min_samples = args.min_samples, nthreads=args.nthreads)
    return

def run_blast_primers (args):
    from phylobarcode import task_blast_primers
    generate_prefix_for_task (args, "blast")
    if len(args.csv) != 2:
        logger.error("Exactly two CSV files must be provided (left and right, in order)")
        return
    uniq = remove_prefix_suffix (args.csv)
    output = [f"{args.prefix}_{outfile}" for outfile in uniq]
    if args.accurate: task = "blastn-short"
    else: task = "blastn"
    if args.max_target_seqs > 1e9: args.max_target_seqs = 1e9
    if args.max_target_seqs < 1: args.max_target_seqs = 1
    if not args.nthreads: 
        args.nthreads = defaults["nthreads"]
        logger.info(f"{args.nthreads} threads are available (actual pool may be smaller)")
    else:
        logger.info(f"{args.nthreads} threads were requested by user (actual pool may be smaller)")
    task_blast_primers.blast_primers_from_csv (csv=args.csv, output=output, database = args.database, evalue=args.evalue, 
            task=task, max_target_seqs = args.max_target_seqs, nthreads=args.nthreads)
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
    parent_group.add_argument('--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING, 
            help="Print debugging statements (most verbose)")
    parent_group.add_argument('--verbose', action="store_const", dest="loglevel", const=logging.INFO, 
            help="Add verbosity")
    parent_group.add_argument('--nthreads', metavar='int', type=int, 
            help="Number of threads requested (default = maximum available)")
    parent_group.add_argument('--outdir', metavar="</path/dir/>", action="store", 
            help="Output directory. Default: working directory")
    parent_group.add_argument('--scratch', action="store", 
            help="Existing scratch directory (i.e. path must exist). Default: working directory")
    parent_group.add_argument('--prefix', metavar="<no_extension>", 
            help="optional file name --- just the prefix, since suffixes (a.k.a. extensions) are added by the program")
    parent_group.add_argument('--version', action='version', version=f"%(prog)s {__version__}") ## called with subcommands

    # alternative to subp= parent_parser.add_subparsers(dest='command', description=None, title="Commands")
    main_parser = ParserWithErrorHelp(description=long_description, formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    main_parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}") ## called without  subcommands (grouped together with --help)
    subp= main_parser.add_subparsers(dest='Commands', description=None, title="Commands", required=True)

    # TODO: deduplicate (using genus+sourmash) and add GTDB info
    this_help = "Given one folder with fasta files and one with GFF files of reference genomes database, creates a table with matches"
    up_findp = subp.add_parser('merge_fasta_gff', help=this_help, description=this_help, parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('-a', '--fasta', required=True, help="directory where fasta genomic files can be found") # technically not needed if gff3 contains fasta
    up_findp.add_argument('-g', '--gff', required=True, help="directory with GFF3 files")
    up_findp.set_defaults(func = run_merge_fasta_gff)

    this_help = "Find primers given an alignment file" # help is shown in "prog -h", description is shown in "prog this -h"
    up_findp = subp.add_parser('find_primers', help=this_help, description=this_help, parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('fasta', help="unaligned sequences")
    up_findp.add_argument('-l', '--length', metavar='int', type=int, 
            help="optimal primer size, in bp (default=20)")
    up_findp.add_argument('-b', '--border', metavar='int', type=int, default=400, 
            help="how far from sequence borders, in bp, we should look for primers (default=400)")
    up_findp.add_argument('-n', '--n_primers', metavar='int', type=int, default=100, 
            help="how many primers, per sequence, per end, should be returned (default=100)")
    up_findp.set_defaults(func = run_find_primers)

    this_help = "Extract and cluster flanking regions of operons where primers may be found"
    up_findp = subp.add_parser('get_flanks', help=this_help, description=this_help, parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('fasta', help="unaligned sequences")
    up_findp.add_argument('-b', '--border', metavar='int', type=int, default=400, 
            help="how far from sequence borders, in bp, we should look for primers (default=400)")
    up_findp.add_argument('-m', '--min_samples', type=int, default=5, 
            help="in OPTICS, minimum number of neighbours for sequence to be a core point (default=5 ;  should be larger than 2)")
    up_findp.add_argument('-i', '--id', type=float, default=None, 
            help="If defined, identity threshold for vsearch (default is to use OPTICS clustering instead of vsearch; suggested value is 0.5 ~ 0.9)")
    up_findp.set_defaults(func = run_cluster_flanks)

    this_help = "Cluster primers described in csv file;\nIf several csv files are given, default output files will keep their unique names (i.e. without common prefix or suffix)"
    up_findp = subp.add_parser('cluster_primers', help=this_help, description=this_help, parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('csv', nargs="+", help="csv files with primers (each ouput file from 'find_primers')")
    up_findp.add_argument('-m', '--min_samples', type=int, default=5, 
            help="in OPTICS, minimum number of neighbours for sequence to be a core point (default=5 ; should be larger than 2)")
    up_findp.set_defaults(func = run_cluster_primers)

    this_help = "Blast primers against database, checking for left-right pairs;\nNeeds exactly two CSV files, with " \
    "left and right primers respect. Output files will keep their unique names (i.e. without common prefix or suffix)"
    up_findp = subp.add_parser('blast_primers', help=this_help, description=this_help, parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('csv', nargs=2, 
            help="csv files with left and right primers, respectively (out from 'find_primers' or 'cluster_primers')")
    up_findp.add_argument('-d', '--database', required=True, 
            help="full path to database prefix")
    up_findp.add_argument('-e', '--evalue', type=float, default=1., 
            help="E-value threshold for blast (default=1, beware of low values)")
    up_findp.add_argument('-m', '--max_target_seqs', type=int, default=100, 
            help="max number of hist per primer (default=100)")
    up_findp.add_argument('-a', '--accurate', action="store_true", 
            help="use accurate blastn-short (default=regular blastn)")
    up_findp.set_defaults(func = run_blast_primers)

    args = main_parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    if args.outdir: 
        defaults["current_dir"] = args.outdir = os.path.join(defaults["current_dir"], args.outdir)
    else:
        args.outdir = defaults["current_dir"]
    pathlib.Path(defaults["current_dir"]).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    if args.scratch:
        tmpdir = os.path.join(defaults["current_dir"], args.scratch)
        if not os.path.isdir(tmpdir):
            logger.error(f"{tmpdir} does not exist or it's not a directory")
            sys.exit(1)
        defaults["scratch"] = args.scratch = os.path.join(defaults["current_dir"], args.scratch, defaults["scratch"])
    else: # in any case, scratch will be a subdirectory of existing scratch area (we can delete it with shutil.rmtree())
        defaults["scratch"] = args.scratch = os.path.join(defaults["current_dir"], defaults["scratch"])

    if args.nthreads:
        if args.nthreads < 1: args.nthreads = 1
        if args.nthreads > defaults["nthreads"]:
            logger.warning(f"Requested {args.nthreads} threads, but only {defaults['nthreads']} are available. Using {defaults['nthreads']} instead.")
            args.nthreads = defaults["nthreads"]

    args.func(args) # calls task 

if __name__ == '__main__':
    main()

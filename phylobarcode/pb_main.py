import os, logging, argparse, sys, pathlib, multiprocessing, datetime

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='phylobarcode_MAIN %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

defaults = {
    "current_dir": os.getcwd() + "/",
    "timestamp": datetime.datetime.now().strftime("%y%m%d_%H%M%S"),
    "n_threads": multiprocessing.cpu_count(),
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


def run_find_primers (args):
    from phylobarcode import task_find_primers
    if args.prefix:
        args.prefix = os.path.join(defaults["current_dir"], args.prefix)
    else:
        args.prefix = os.path.join(defaults["current_dir"], f"{defaults['timestamp']}_primers")
    if args.nthreads and args.nthreads < 2:
        logger.warning("Single-threaded mode requested by user")
        task_find_primers.find_primers (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, num_return=args.n_primers, output=args.prefix)
        return
    try:
        from multiprocessing import Pool
        from functools import partial
    except ImportError:
        logger.warning("Multiprocessing not available, falling back to single-threaded mode")
        task_find_primers.find_primers (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, num_return=args.n_primers, output=args.prefix)
        return

    if not args.nthreads or args.nthreads > defaults["n_threads"]: args.nthreads = defaults["n_threads"]
    logger.info(f"{args.nthreads} threads are available (actual pool may be smaller)")
    task_find_primers.find_primers_parallel (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, num_return=args.n_primers, output=args.prefix, nthreads=args.nthreads)
    return

class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help(sys.stderr)
        sys.exit(2)

def main():
    parent_parser = ParserWithErrorHelp(description=long_description, add_help=False)
    parent_group = parent_parser.add_argument_group('Options common to all commands')
    parent_group.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING, help="Print debugging statements (most verbose)")
    parent_group.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parent_group.add_argument('--nthreads', metavar='int', type=int, help="Number of threads requested (default = maximum available)")
    parent_group.add_argument('--outdir', metavar="</path/dir/>", action="store", help="Output directory. Default: working directory")
    parent_group.add_argument('--prefix', metavar="<no_extension>", help="optional file name --- just the prefix, since suffixes (a.k.a. extensions) are added by the program")

    # alternative to subp= parent_parser.add_subparsers(dest='command', description=None, title="Commands")
    main_parser = ParserWithErrorHelp(description=long_description, formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    subp= main_parser.add_subparsers(dest='Commands', description=None, title="Commands", required=True)

    this_help = "Find primers given an alignment file" # help is shown in "prog -h", description is shown in "prog this -h"
    up_findp = subp.add_parser('find_primers', help=this_help, description=this_help, parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter, epilog=epilogue)
    up_findp.add_argument('fasta', help="unaligned sequences")
    up_findp.add_argument('-l', '--length', metavar='int', type=int, help="optimal primer size, in bp (default=20)")
    up_findp.add_argument('-b', '--border', metavar='int', type=int, help="how far from sequence borders, in bp, we should look for primers (default=400)")
    up_findp.add_argument('-n', '--n_primers', metavar='int', type=int, help="how many primers, per sequence, per end, should be returned (default=100)")
    up_findp.set_defaults(func = run_find_primers)

    args = main_parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    if args.outdir: 
        defaults["current_dir"] = args.outdir = os.path.join(defaults["current_dir"], args.outdir)
        common.pathlib.Path(defaults["current_dir"]).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist

    if args.nthreads:
        if args.nthreads < 1: args.nthreads = 1
        max_threads = multiprocessing.cpu_count()
        if args.nthreads > max_threads:
            logger.warning(f"Requested {args.nthreads} threads, but only {max_threads} are available. Using {max_threads} instead.")
            args.nthreads = max_threads
        defaults["n_threads"] = args.nthreads

    args.func(args) # calls task 

if __name__ == '__main__':
    main()

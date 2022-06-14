import os, logging, argparse, sys, pathlib, multiprocessing, datetime

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='phylobarcode_MAIN  %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M%S")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

defaults = {
    "current_dir": os.getcwd() + "/",
    "timestamp": datetime.datetime.now().strftime("%y%m%d_%H%M%S"),
    "n_threads": multiprocessing.cpu_count(),
#   "reference": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas") 
    }

def run_find_primers (args):
    from phylobarcode import task_find_primers
#    if args.timestamp is not None and not args.timestamp.isdigit():
#        logger.warning(f"provided timestamp '{args.timestamp}' is not numeric, may cause problems downstream. Ideally it would be a YearMonthDay")
    if args.output:
        args.output = os.path.join(defaults["current_dir"], args.output)
    else:
        args.output = os.path.join(defaults["current_dir"], f"{defaults['timestamp']}_primers")
    task_find_primers.find_primers (fastafile=args.fasta, primer_opt_size=args.length, border=args.border, num_return=args.n_primers, output=args.output)

class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp()
    #description="""
    #peroba top level
    #""") 

    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING, help="Print debugging statements (most verbose)")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('--nthreads', metavar='int', help="Number of threads requested (default = maximum available)")
#    parser.add_argument('--timestamp', help="optional timestamp for naming new files. You can safely ignore it, otherwise use format YYMMDD")
    parser.add_argument('--outdir', action="store", help="Output database directory. Default: working directory")
    subp= parser.add_subparsers(dest='command')

    up_findp = subp.add_parser('find_primers', help="find primers given an alignment")
    up_findp.add_argument('fasta', help="unaligned sequences")
    up_findp.add_argument('-o', '--output', metavar='prefix', help="optional file name prefix of table with primers")
    up_findp.add_argument('-l', '--length', metavar='int', help="optimal primer size, in bp (default=20)")
    up_findp.add_argument('-b', '--border', metavar='int', help="how far from sequence borders, in bp, we should look for primers (default=400)")
    up_findp.add_argument('-n', '--n_primers', metavar='int', help="how many primers, per sequence, per end, should be returned (default=100)")
    up_findp.set_defaults(func = run_find_primers)

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    if args.outdir: 
        defaults["current_dir"] = args.outdir = os.path.join(defaults["current_dir"], args.outdir)
        common.pathlib.Path(defaults["current_dir"]).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist
    if args.nthreads:
        if args.nthreads < 1: args.nthreads = 1
        defaults["n_threads"] = args.nthreads

    args.func(args) # calls task 

if __name__ == '__main__':
    main()

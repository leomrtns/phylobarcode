#!/usr/bin/env python
from phylobarcode.pb_common import *  ## better to have it in json?
import re

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='phylobarcode_PRIM %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

def find_primers (defaults, fastafile = None, output = None, entry_timestamp = None):
    fas = read_fasta_as_list (fastafile)
    for i in fas[:4]:
        left, right = get_primers (i.seq, i.id, num_return=5)
        print (i.id, "\n", left, "\n", right)
    

def get_primers (sequence, seqname, primer_opt_size = 21, border = 400, num_return=100):
    seqlen = len(sequence)
    primer_task = "pick_primer_list" # "pick_sequencing_primers" "generic" "pick_primer_list"
    primer_min_size = 14
    primer_max_size = 28
    if primer_opt_size < primer_min_size:
        primer_opt_size = primer_min_size
    if primer_opt_size > primer_max_size:
        primer_opt_size = primer_max_size
    max_n_accepted = 1
    if border < primer_opt_size:
        border = primer_opt_size
    if num_return < 1:
        num_return = 1
    if seqlen - 2 * border - 2 * primer_opt_size < seqlen * 0.1:
        border = int(seqlen * 0.45 - primer_opt_size)
        if border < primer_opt_size:
            logger.info (f"Border is too small for sequence {seqname}")
            return None, None;
    min_product_size = seqlen - 2 * border - 2 * primer_opt_size;
    if (min_product_size > 0.8 * seqlen):
        min_product_size = int(0.8 * seqlen)

    arguments = (     # sequence must be in one line, and there must be no extra newlines or spaces
        "PRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_INTERNAL_OLIGO=0\nPRIMER_PICK_RIGHT_PRIMER=1\n"
        f"SEQUENCE_ID={seqname}\nPRIMER_TASK={primer_task}\n"
        f"PRIMER_MIN_SIZE={primer_min_size}\nPRIMER_MAX_SIZE={primer_max_size}\n"
        f"PRIMER_MAX_NS_ACCEPTED={max_n_accepted}\n"
        f"PRIMER_PRODUCT_SIZE_RANGE={min_product_size}-{seqlen}\n" # 5000-9500
        f"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,{border},{seqlen-border},{border}\n"  # left, left_len, right, right_len : 0,500,9000,554
        f"PRIMER_NUM_RETURN={num_return}\nPRIMER_OPT_SIZE={primer_opt_size}\n"
        f"SEQUENCE_TEMPLATE={sequence}\n=\n"
        )
    child = subprocess.Popen(" primer3_core " , stdin=subprocess.PIPE, stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
    output = child.communicate(input=arguments)[0].split("\n")

    return extract_primer_from_output (output)

def extract_primer_from_output (output):
    '''
    Extracts the primers from the output of primer3_core. Returns lists of left and right primers, with
    [sequence,penalty, position] per primer.
    '''
    def xtract_left_or_right (out, side):
        y=[re.match(f"PRIMER_{side}_(\d+)_SEQUENCE=(.*)",x) for x in out] # most will be "None" since won't match
        pseq = [[int(i.group(1)),i.group(2)] for i in y if i is not None]
        y=[re.match(f"PRIMER_{side}_(\d+)_PENALTY=(.*)",x) for x in out] # most will be "None" since won't match
        ppen = [[int(i.group(1)),float(i.group(2))] for i in y if i is not None]
        y=[re.match(f"PRIMER_{side}_(\d+)=(\d+),(\d+)",x) for x in out] ## =start,length (but we dont use length)
        pstart = [[int(i.group(1)),int(i.group(2))] for i in y if i is not None]
        # pseq=[idx, sequence]; ppen=[idx, penalty]; pstart=[idx, start]; we want [sequence,penalty, start]

        idx = [i[0] for i in pseq] + [i[0] for i in ppen] + [i[0] for i in pstart] # flattened version of two index lists (e.g. [0,1,2...,0,1,2...])
        min_idx = min(idx) ## make sure it starts from 0 (which is already the case in current version of primer3...)
        spl = [[None,None,None] for i in range (max(idx) + 1 - min_idx)] ## length = max between two lists flattened (base_zero -> max +1 is length)
        for i in pseq:
            spl[i[0]-min_idx][0] = i[1]
        for i in ppen:
            spl[i[0]-min_idx][1] = i[1]
        for i in pstart:
            spl[i[0]-min_idx][2] = i[1]
        return [i for i in spl if i[0] is not None]
    left = xtract_left_or_right (output, "LEFT")
    right = xtract_left_or_right (output, "RIGHT")
    return left, right

def extract_primer_from_output_legacy (output):
    pseq_l = [[int(x.split("_")[2]), x.split("=")[1]] for x in output if "_SEQUENCE=" in x and "PRIMER_LEFT_" in x]  # only lines with PRIMER_LEFT_xxx_SEQUENCE
    pseq_r = [[int(x.split("_")[2]), x.split("=")[1]] for x in output if "_SEQUENCE=" in x and "PRIMER_RIGHT_" in x] # only lines with PRIMER_RIGHT_xxx_SEQUENCE
    penalty_l = [[int(x.split("_")[2]), x.split("=")[1]] for x in output if "_PENALTY=" in x and "PRIMER_LEFT_" in x]
    penalty_r = [[int(x.split("_")[2]), x.split("=")[1]] for x in output if "_PENALTY=" in x and "PRIMER_RIGHT_" in x]
    def seq_penalty_list(pseq, penalty):
        idx = [i[0] for i in pseq] + [i[0] for i in penalty] # flattened version of two index lists (e.g. [0,1,2...,0,1,2...])
        min_idx = min(idx) ## make sure it starts from 0 (which is already the case in current version of primer3...)
        spl = [[None,None] for i in range (max(idx) - min_idx)] ## length = max between two lists flattened
        for i in pseq:
            spl[i[0]-min_idx][0] = i[1]
        for i in penalty:
            spl[i[0]-min_idx][1] = i[1]
        return spl
    pseq_l = seq_penalty_list(pseq_l, penalty_l)
    pseq_r = seq_penalty_list(pseq_r, penalty_r)
    return pseq_l, pseq_r




#! /usr/bin/env python3
import sys
import os
import numpy as np
import argparse
"""
Convert output from genmap (mappability program) to the format created by
umap (mappability program). This can then be used together with a BAM file
of ChiP input data by the ENCODE blacklist program to create a blacklist
BED file.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", "-o", help="The output directory. \
Output files will be named <chromosome>.uint8.unique", required=True)
    parser.add_argument("--infiles", "-i", help="The input file(s). \
If more than one is given, they will be consolidated", required=True, \
    nargs="+")
    return parser.parse_args()

def nextseq(f):
    seqid = None
    seq = []
    for line in f:
        if line[0] == '>':
            if seqid is not None:
                yield (seqid, np.array(seq, dtype=float))
            seqid = line.strip()[1:].split(' ')[0]
        else:
            seq += line.strip().split(' ')
    if seqid is not None:
        yield (seqid, np.array(seq, dtype=float))

def main(args):
    options = parse_args()
    
    ptrs = []
    for filename in options.infiles:
        ptrs.append(open(filename, 'r'))
    iters = []
    for ptr in ptrs:
        iters.append(nextseq(ptr))

    for seqid, mapvals in iters[0]:
        # Convert map to 1 (uniquely mappable) or 0 (not)
        maps_bool = np.array(mapvals == 1)
        for itr in iters[1:]:
            seqid2, maps2 = next(itr)
            if seqid2 != seqid:
                print("ERROR: files not in same order; seqs = {} {}"\
                    .format(seqid, seqid2), file=sys.stderr)
                exit(1)
            maps_bool = np.logical_and(maps_bool, np.array(maps2 == 1))
        maps_bool = maps_bool.astype(np.uint8)
        outf = open('{}/{}.uint8.unique'.format(options.outdir, seqid), 'wb')
        outf.write(maps_bool.tobytes())

if __name__ == "__main__":
    sys.exit(main(sys.argv))

#! /usr/bin/env python3
import sys
import os
import numpy as np
import argparse
"""
Convert output from genmap (mappability program) to a BED file listing
regions not uniquely mappable.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infile", "-i", help="The input file.", required=True)
    parser.add_argument("--minlen", "-m", help="Minimum length run to print", type=int, default=100)
    parser.add_argument("--fai", "-f", help="FASTA index for reference genome", required=True)
    return parser.parse_args()

def nextseq(f):
    seqid = None
    seq = []
    for line in f:
        if len(line) == 0:
            continue
        if line[0] == '>':
            if seqid is not None:
                yield (seqid, np.array(seq, dtype=float))
            seqid = line.strip()[1:].split(' ')[0]
        else:
            seq = line.strip().split(' ')
    if seqid is not None:
        yield (seqid, np.array(seq, dtype=float))

def main(args):
    options = parse_args()
    
    chromlens = {}
    f = open(options.fai, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        chromlens[dat[0]] = int(dat[1])
    f.close()

    ptr = open(options.infile, 'r')
    
    for seqid, mapvals in nextseq(ptr):
        if len(mapvals) > chromlens[seqid]:
            print("ERROR: {} length > fasta index length ({})".format(seqid, chromlens[seqid]), file=sys.stderr)
            exit(1)

        runstart = None
        runend = None
        # Convert map to 1 (uniquely mappable) or 0 (not)
        maps_bool = np.array(mapvals == 1)
        for idx, elt in enumerate(maps_bool):
            if not elt:
                # Not mappable
                if runstart is None:
                    runstart = idx
                    runend = idx
                else:
                    runend = idx
            else:
                if runstart is not None and runend-runstart+1 >= options.minlen:
                    print("{}\t{}\t{}".format(seqid, runstart, runend+1))
                runstart = None
                runend = None
        if runstart is not None and runend-runstart+1 >= options.minlen:
            print("{}\t{}\t{}".format(seqid, runstart, runend+1))

if __name__ == "__main__":
    sys.exit(main(sys.argv))

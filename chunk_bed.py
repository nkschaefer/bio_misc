#! /usr/bin/env python3
import sys
import os
import gzip
import math
import subprocess

def main(args):
    if len(args) < 4:
        print("USAGE: chunk_bed.py BED numchunks outbase", file=sys.stderr)
        exit(1)

    vcf = args[1]
    numchunks = int(args[2])
    outbase = args[3]

    # Learn how many lines
    nlines = 0
    
    f = None
    gzipped_input = False
    if vcf[-3:] == '.gz':
        gzipped_input = True
        f = gzip.open(vcf, 'r')
    else:
        f = open(vcf, 'r')

    for line in f:
        if gzipped_input:
            line = line.decode()
        if line[0] == '#':
            continue
        else:
            nlines += 1
    f.close()

    print("Input has {} lines".format(nlines), file=sys.stderr)

    vars_per_chunk = math.ceil(nlines / numchunks)

    print("Writing {} lines per chunk".format(vars_per_chunk), file=sys.stderr)

    f = None

    if vcf[-3:] == '.gz':
        f = gzip.open(vcf, 'r')
    else:
        f = open(vcf, 'r')

    outs = []
    outnames = []
    for i in range(0, numchunks):
        fname = "{}.{}.bed".format(outbase, i)
        outnames.append(fname)
        fh = open(fname, 'w')
        outs.append(fh)

    chunk_idx = 0
    n_written = 0

    for line in f:
        if gzipped_input:
            line = line.decode()
        line = line.rstrip()
        if line[0] == '#':
            for out in outs:
                print(line, file=out)
        else:
            print(line, file=outs[chunk_idx])
            n_written += 1
            if n_written >= vars_per_chunk:
                chunk_idx += 1
                n_written = 0
                if chunk_idx > numchunks-1:
                    chunk_idx = 0

    for i in range(0, len(outs)):
        outs[i].close()

if __name__ == '__main__':
    sys.exit(main(sys.argv))

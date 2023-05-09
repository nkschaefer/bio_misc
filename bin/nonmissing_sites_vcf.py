#! /usr/bin/env python3
import sys
import os
import gzip
from collections import Counter

def main(args):
    if len(args) < 2:
        print("USAGE: nonmissingsites_vcf.py FILE.vcf.gz", file=sys.stderr)
        print("Prints the number of sites where each individual has a non-missing genotype.", file=sys.stderr)
        exit(1)
    f = gzip.open(args[1], 'r')
    hdr = []
    counts = Counter()
    for line in f:
        line = line.decode().rstrip()
        if line[0:2] == "##":
            continue
        elif line[0] == '#':
            hdr = line.split('\t')[9:]
        else:
            dat = line.split('\t')
            for idx in range(9, len(dat)):
                gt = dat[idx].split(':')[0]
                if gt != '.' and gt != './.':
                    counts[hdr[idx-9]] += 1
    f.close()

    for item in sorted(counts.items(), key=lambda x: -x[1]):
        print("{}\t{}".format(item[0], item[1]))

if __name__ == '__main__':
    sys.exit(main(sys.argv))


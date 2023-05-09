#! /usr/bin/env python3
import sys
import os
import random
import argparse
import gzip
"""
Generates random 10X cell barcodes.
Ensures that no barcodes conflict.
Optionally can select from a whitelist instead of generating random ones.
"""

def parse_whitelist(filename):
    f = None
    if filename[-3:] == '.gz':
        f = gzip.open(filename, 'r')
    else:
        f = open(filename, 'r')
    bcs = []
    for line in f:
        if filename[-3:] == '.gz':
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        bcs.append(line)
    f.close()
    return bcs

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--num", "-n", type=int, help="Number of desired barcodes", required=True)
    parser.add_argument("--whitelist", "-w", help="Barcode whitelist to select from (OPTIONAL)", required=False)
    return parser.parse_args()

def main(args):
    options = parse_args()
    if options.whitelist is not None:
        wl = parse_whitelist(options.whitelist)
        if options.num > len(wl):
            print("ERROR: requested number ({}) greater than whitelist size ({})".format(\
                options.num, len(wl)), file=sys.stderr)
            exit(1)
        samp = random.sample(wl, options.num)
        for s in samp:
            print("{}-1".format(s))
    else:
        # Randomly generate.
        chosen = set([])
        while len(chosen) < options.num:
            bc = [None] * 16
            for idx in range(0, 16):
                bc[idx] = random.choice(['A', 'C', 'G', 'T'])
            bc = "".join(bc)
            chosen.add(bc)
        for bc in chosen:
            print("{}-1".format(bc))


if __name__ == "__main__":
    sys.exit(main(sys.argv))

#! /usr/bin/env python3
import sys
import os
import argparse
import pysam
import random
from numpy.random import poisson
"""
Generates a chosen number of random cell barcodes and inserts into a BAM file for bulk
(ATAC/RNA-seq) data. Assumes the whole BAM came from a single individual.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", "-b", help="Input BAM file", required=True)
    parser.add_argument("--ncell", "-n", help="Number of cells to generate", required=True, type=int)
    parser.add_argument("--barcodes_out", "-bo", \
        help="Optional file to store generated barcodes", required=False, default=None)
    parser.add_argument("--barcode_blacklist", "-bb", \
        help="Optional file containing barcodes not to generate", required=False, default=None)
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    barcode_out = None
    if options.barcodes_out is not None:
        barcode_out = open(options.barcodes_out, "w")
    
    barcode_bl = set([])
    if options.barcode_blacklist is not None:
        f = open(options.barcode_blacklist, 'r')
        for line in f:
            line = line.rstrip()
            barcode_bl.add(line)
        f.close()

    # Generate barcodes
    print("Generating {} random cell barcodes".format(options.ncell), file=sys.stderr)
    chosen = set([])
    while len(chosen) < options.ncell:
        bc = [None] * 16
        for idx in range(0, 16):
            bc[idx] = random.choice(['A', 'C', 'G', 'T'])
        bc = "".join(bc) + "-1"
        if bc not in barcode_bl:
            if barcode_out is not None:
                print(bc, file=barcode_out)
            chosen.add(bc)
        
    if barcode_out is not None:
        barcode_out.close()

    chosen = list(chosen)
    
    # Count total reads
    print("Counting reads", file=sys.stderr)
    totreads = int(pysam.view("-c", options.bam))
    print("Found {} total reads".format(totreads), file=sys.stderr)

    meancov = int(round(totreads / options.ncell))
    
    print("Determining cell barcode probabilities", file=sys.stderr)
    # Create CDF-like structure
    covs = list(poisson(meancov, len(chosen)))
    covstot = 0
    cdf = []
    for idx, bc in enumerate(chosen):
        cdf.append((covs[idx], bc))
        covstot += covs[idx]
    
    cdf.sort(key=lambda x: -x[0])

    # Parse input BAM
    bam = pysam.AlignmentFile(options.bam, 'r')
    
    # Create file to write to
    outfile = pysam.AlignmentFile("-", "wb", template=bam)
    
    print("Traversing input BAM file", file=sys.stderr)

    for rec in bam:

        # Choose a cell barcode using weighted probabilities from Poisson
        r = random.random()
        cumulative = 0
        bc = None
        for elt in cdf:
            cumulative += (elt[0]/covstot)
            if cumulative >= r:
                # This is the bin.
                bc = elt[1]
                break
    
        if bc is not None:
            # Insert it.
            rec.set_tag("CB", bc)

        outfile.write(rec)




if __name__ == '__main__':
    sys.exit(main(sys.argv))

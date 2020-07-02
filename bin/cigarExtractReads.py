#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import re
"""
Given a BAM alignment -- especially if it's Omni-C data,
extracts just the mapped portion of reads, excluding the junction
sequence. Outputs FASTQ format.

cigarExtractReads.py
Created by Nathan Schaefer on 01/13/20 at 14:33"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()

def main(args):
    """Main method"""
    cigparse = re.compile(r'([0-9]+)([MIDNSHP=X])')
    for line in sys.stdin:
        line = line.rstrip()
        if line[0] != "#" and line[0] != "@":
            dat = line.split('\t')
            cigstr = dat[5]
            read_id = dat[0]
            # Build up sequence and quality strings as we go
            seq = ""
            qual = ""
            
            if False:
                seq_index = 0
                cigops = cigparse.findall(cigstr)
                for cigop in cigops:
                    if cigop[1] == 'M' or cigop[1] == '=' or cigop[1] == 'X':
                        seq += dat[9][seq_index:seq_index + int(cigop[0])]
                        qual += dat[10][seq_index:seq_index + int(cigop[0])]
                    if cigop[1] == 'M' or cigop[1] == 'I' or cigop[1] == 'S' or \
                        cigop[1] == '=' or cigop[1] == 'X':
                        seq_index += int(cigop[0])
            else:
                seq = dat[9]
                qual = dat[10]
                    
            print("@{}\n{}\n+\n{}".format(read_id, seq, qual))
            
if __name__ == "__main__" :
    sys.exit(main(sys.argv))


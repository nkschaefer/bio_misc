#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import re
"""
Pipe in SAM format. Parses the CIGAR string to figure out the end
coordinate of the alignment. Prints out BED format for each 
alignment.

cigar2bed.py
Created by Nathan Schaefer on 01/13/20 at 14:06
"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()

def main(args):
    """Main method"""
    cigparse = re.compile(r'([0-9]+)([MIDNSHP=X])')
    for line in sys.stdin:
        line = line.rstrip()
        if line[0] != "@":
            dat = line.split('\t')
            cigstr = dat[5]
            map_start = int(dat[3]) - 1
            map_end = map_start
            cigops = cigparse.findall(cigstr)
            for cigop in cigops:
                if cigop[1] == 'M' or cigop[1] == 'D' or cigop[1] == 'N' or \
                    cigop[1] == '=' or cigop[1] == 'X':
                    map_end += int(cigop[0])
            print("{}\t{}\t{}\t{}".format(dat[2], map_start, map_end+1, dat[0]))
            
if __name__ == "__main__" :
    sys.exit(main(sys.argv))


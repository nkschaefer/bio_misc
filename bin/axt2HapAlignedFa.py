#! /usr/bin/env python2.7
from __future__ import division, print_function
import sys
import argparse
import re
import gzip
"""
axt2HapAlignedFa.py
Created by Nathan Schaefer on 09/04/14 at 14:40

Converts an AXT formatted alignment to a haploid FASTA file in the coordinates
of the reference genome. Ignores gaps in the reference sequence and replaces
all gaps in the query sequence with ambiguous "N" bases.

"""

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()

def read_axt(filename):
    if filename[-3:] == ".gz":
        fp = gzip.open(filename, "r")
    else:
        fp = open(filename, "r")
    
    readingAlignment = False
    nextAlnNumber = 0
    
    refStart = -1
    refEnd = -1
    chrom = None
    
    refSeqBlock = None
    qSeqBlock = None
    
    for line in fp:
        line = line.rstrip()
        if len(line) == 0 or line[0] == "#":
            continue
        elif not readingAlignment and "{}".format(nextAlnNumber) == line.split()[0]:
            data = line.split()
            alnNumber = int(data[0])
            if alnNumber == nextAlnNumber:
                nextAlnNumber += 1
                readingAlignment = True
                refStart = int(data[2])
                refEnd = int(data[3])
                chrom = data[1]
            else:
                print("ERROR: alignment block indices do not match: {} {}".format(alnNumber, nextAlnNumber), file=sys.stderr)
                print(line, file=sys.stderr)
                exit(1)
        else:
            # Read sequences
            if refSeqBlock is None:
                refSeqBlock = line
            elif qSeqBlock is None:
                qSeqBlock = line
                
                # Get ready to return data
                
                # Transform gaps in query seq into "N"
                #qSeqBlock = qSeqBlock.replace('-', 'N')
                
                # Remove gaps in reference seq from query seq
                refGaps = re.finditer('\-+', refSeqBlock)
                
                qseqFixed = ''
                qseqIndex = 0
                
                for refGap in refGaps:
                    #gaplen = len(refGap.group(0))
                    gapStart = refGap.start(0)
                    gapEnd = refGap.end(0)
                    for i in range(qseqIndex, gapStart):
                        qseqFixed += qSeqBlock[i]
                    qseqIndex = gapEnd  
                
                for i in range(qseqIndex, len(qSeqBlock)):
                    qseqFixed += qSeqBlock[i]
                
                # Yield data
                yield(chrom, refStart, refEnd, qseqFixed)
                
                # Reset data
                readingAlignment = False
                refSeqBlock = None
                qSeqBlock = None
    fp.close()
    
def main(args):
    """Main method"""
    # Arguments: input AXT files
    for axtfilename in args[1:]:
        lastRefIndex = 1
        lastChromName = None
        needLineBreak = False
        for chromName, refStart, refEnd, qSeq in read_axt(axtfilename):
            if chromName != lastChromName:
                if needLineBreak:
                    print("\n", end="")
                print(">{}".format(chromName))
                lastChromName = chromName
                
            # Add any sequence since the end of the previous sequence
            if (refStart - lastRefIndex) > 0:
                fillerSeq = "-" * (refStart - lastRefIndex)
                print(fillerSeq, end="")
            
            lastRefIndex = refEnd + 1
            
            # Print this sequence block
            print(qSeq.upper(), end="")
            
            needLineBreak = True
        
if __name__ == "__main__" :
    sys.exit(main(sys.argv))


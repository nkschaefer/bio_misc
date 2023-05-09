#! /usr/bin/bash

if [ $# -lt 2 ]; then
    >&2 echo "USAGE: merge_many_bams.sh bamlist.txt outfilename"
    exit 1
fi
bamlist=$1
outfile=$2

outbase="${outfile%.bam}"

# Merge in groups of 50
cat $bamlist | awk '{printf("%s\t%d\n", $line, NR % 50); }' | sort -k2,2n | cut -f2 | uniq | while read idx; do
    files=$( cat $bamlist | awk '{printf("%s\t%d\n", $line, NR % 50);}' | grep -P "\t${idx}$" | cut -f1 | tr '\n' ' ' )
    samtools merge "${outbase}.${idx}.bam" $files
done

out_merge=$( cat $bamlist | awk '{printf("%s\t%d\n", $line, NR % 50); }' | sort -k2,2n | cut -f2 | uniq | sed "s/^/${outbase}./" | sed "s/$/.bam/" | tr '\n' ' ' )

samtools merge $outfile $out_merge
rm $out_merge


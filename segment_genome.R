#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
    cat("USAGE: segment_genome.R fai nblocks outprefix\n")
    cat("Given a fai file for a genome FASTA (indexed using samtools faidx),\n")
    cat("Creates BED files for genomic regions segmenting the genome into nblocks chunks.\n")
    cat("Cannot subdivide chromosomes - this is intended to break up a genome with many\n")
    cat("contigs/scaffolds into a manageable number of chunks for parallel processing.\n")
    q()
}

fai <- read.table(args[1])
nblocks <- as.numeric(args[2])
outprefix <- args[3]

fai <- fai[,c(1,2)]
# Preserve original sort order
fai_orig <- fai 
fai <- fai[order(fai$V2, decreasing=T),]
fai$start <- 0 

totsize <- sum(as.numeric(fai$V2))

blocksize <- round(totsize / nblocks)

block_start <- 1
block_idx <- 1
blocksum <- 0

for (i in seq(1, length(rownames(fai)))){
    blocksum <- blocksum + fai[i,2]
    if (blocksum >= blocksize){
        chunk <- fai[block_start:i,c(1,3,2)]
        #chunk <- chunk[order(chunk$V1),]
        chunk_orig <- fai_orig[which(fai_orig$V1 %in% chunk$V1),]
        chunk <- chunk[match(chunk_orig$V1, chunk$V1),]
        write.table(chunk, file=paste(outprefix, block_idx, 'bed', sep='.'), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)        
        block_idx <- block_idx + 1
        blocksum <- 0
        block_start <- i + 1
    }
}


#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    cat("Please provide stats file (stdout from vcf_depth_filter)\n")
    q()
}

statsf <- args[1]
stats <- read.table(statsf)
colnames(stats) <- c("type", "id", "num1", "num2")

stats <- stats[which(stats$type=="DP"),]
maxes <- aggregate(stats$num2, by=list(id=stats$id), FUN=max)
colnames(maxes)[2] <- "max"

stats <- merge(stats, maxes, by="id")
stats <- stats[which(stats$num2==stats$max),]

stats <- stats[,which(colnames(stats) %in% c("id", "num1"))]
colnames(stats)[2] <- "depth"
stats <- stats[order(stats$depth, decreasing=T),]

write.table(stats, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

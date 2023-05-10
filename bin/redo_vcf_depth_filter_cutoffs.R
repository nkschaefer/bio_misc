#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
    cat("Please provide stats file (stdout from vcf_depth_filter), and a new lower and upper cutoff percentile.\n")
    q()
}

statsf <- args[1]
stats <- read.table(statsf)
colnames(stats) <- c("type", "id", "num1", "num2")

lower <- as.numeric(args[2])
upper <- as.numeric(args[3])

stats <- stats[which(stats$type=="DP"),]

res <- data.frame(id=c(), lower=c(), upper=c())

for (id in unique(stats$id)){
    subs <- stats[which(stats$id==id),]
    subs <- subs[which(subs$num1 != 0),]
    subs$area <- subs$num1*subs$num2
    tot <- sum(subs$area)
    subs$cumulative <- cumsum(subs$area / tot)
    dp_lower <- head(subs[which(subs$cum >= lower),], 1)$num1
    dp_upper <- head(subs[which(subs$cum >= upper),], 1)$num1
    row <- data.frame(id=c(id), lower=c(dp_lower), upper=c(dp_upper))
    res <- rbind(res, row)
}

write.table(res, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


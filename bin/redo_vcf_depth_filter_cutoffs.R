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

modes <- stats[which(stats$type=="MODE_DP_GQ"),]
modes <- modes[,c(2,3)]

stats <- stats[which(stats$type=="DP"),]

res <- data.frame(id=c(), peak=c(), lower=c(), upper=c(), lower_exp=c(), upper_exp=c())

for (id in unique(stats$id)){
    subs <- stats[which(stats$id==id),]
    subs <- subs[which(subs$num1 != 0),]
    subs$area <- subs$num1*subs$num2
    tot <- sum(subs$area)
    subs$cumulative <- cumsum(subs$area / tot)
    dp_lower <- head(subs[which(subs$cum >= lower),], 1)$num1
    dp_upper <- head(subs[which(subs$cum >= upper),], 1)$num1
    expmean <- 1.0 / modes[which(modes$id==id),]$num1
    dp_lower_alt <- round(-log(1-lower)/expmean)
    dp_upper_alt <- round(-log(1-upper)/expmean)
    row <- data.frame(id=c(id), peak=c(modes[which(modes$id==id),]$num1), 
        lower=c(dp_lower), upper=c(dp_upper), lower_exp=c(dp_lower_alt), upper_exp=c(dp_upper_alt))
    res <- rbind(res, row)
}

write.table(res, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


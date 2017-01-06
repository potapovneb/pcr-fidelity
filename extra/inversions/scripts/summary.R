args <- commandArgs(TRUE)

ifile <- args[1]
rname <- args[2]

## read data
x <- read.table( ifile, sep=",", header=TRUE )

## compute precent coverage
qcov <- x$QueryCoverage/x$RefLength

## define filtering criteria
ix <- x$Reference == rname & qcov >= 0.95 & x$NP >= 3 & x$Mismatches == 0

## count number of reads
num <- sum(ix)
cat(num)

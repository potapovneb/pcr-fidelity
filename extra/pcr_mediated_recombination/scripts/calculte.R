args <- commandArgs(TRUE)

ifile <- args[1]
amplicon.length <- as.numeric(args[2])
input <- as.numeric(args[3])
yield <- as.numeric(args[4])

table <- NULL

x <- read.table(ifile, sep=",", header=TRUE)

## number of recombination events
table$N.re <- sum(apply(x,1,function(z) z["Recombinations"] * z["Reads"]))

## total number of sequenced bases
table$N.total <- sum(x["Reads"]*amplicon.length)

## recombination rate
doubling <- log(yield/input,2)
table$Recombination.rate <- (table$N.re / table$N.total) * (2/doubling) * 2

## strand w/ at least 1 recombination event
table$N.strands.with.re <- sum(x$Reads[x$Recombination > 0]) / sum(x$Reads)

write.table(table,"",quote=TRUE,sep=",",col.names=TRUE,row.names=FALSE)

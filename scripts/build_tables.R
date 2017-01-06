args <- commandArgs(TRUE)
ifile <- args[1]

## load data
df <- read.table( ifile, sep=",", header=TRUE, stringsAsFactors=FALSE )

## define column name lists
colsDesc <- c("Enzyme","Amplicon","Input","Yield")
colsData <- c("AA","AC","AT","AG","CA","CC","CT","CG","TA","TC","TT","TG","GA","GC","GT","GG","Deletion","Insertion")

colsAT <- c("AA","AC","AT","AG","TA","TC","TT","TG")
colsGC <- c("CA","CC","CT","CG","GA","GC","GT","GG")

subsAT <- c("AC","AT","AG","TA","TC","TG")
subsGC <- c("CA","CT","CG","GA","GC","GT")


## ----- Aggregate PacBio data -------------------------------------------------

## group by GroupID
desc <- aggregate(df[colsDesc], df["GroupID"], unique)
data <- aggregate(df[colsData], df["GroupID"], sum)
x <- merge(desc,data,by="GroupID")

## PCR & plasmid sample definitions (used later)
PcrSamplesIx <- x$GroupID >= 1 & x$GroupID <= 31
PlasmidSamplesIx <- x$GroupID >= 34 & x$GroupID <= 46

## total sequenced bases
x$TotalAT  <- apply(x[colsAT],1,sum)
x$TotalGC  <- apply(x[colsGC],1,sum)
x$Total    <- x$TotalAT + x$TotalGC + x$Insertion

## raw error rates
x$ErrorAT  <- apply(x[subsAT],1,sum) / x$TotalAT
x$ErrorGC  <- apply(x[subsGC],1,sum) / x$TotalGC
x$ErrorSub <- (x$ErrorAT + x$ErrorGC) / 2
x$ErrorDel <- x$Deletion / x$Total
x$ErrorIns <- x$Insertion / x$Total

## normalized error rates
x$Doubling <- log2(x$Yield/x$Input)

x$RateSub <- x$ErrorSub * 2 / x$Doubling
x$RateDel <- x$ErrorDel * 2 / x$Doubling
x$RateIns <- x$ErrorIns * 2 / x$Doubling

## Fidelity relative to Taq
taqfind <- function(y) {
    ix <- PcrSamplesIx & x$Enzyme == "Taq" & x$Amplicon == y["Amplicon"];
    if( sum(ix) == 1 )
        return(x$RateSub[ix])
    else
        return(NA)
}
TaqRateSub <- apply(x, 1, function(x) taqfind(x))
x$Fidelity <- TaqRateSub / x$RateSub



## ----- Table 1 ---------------------------------------------------------------
## Error rate of Taq DNA polymerase
## -----------------------------------------------------------------------------

## find Taq samples
TaqSampleIx <- PcrSamplesIx & x$Enzyme == "Taq"

table <- x[TaqSampleIx,c("Amplicon","RateSub","RateDel","RateIns","Total")]

table$RateTotal  <- table$RateSub + table$RateDel + table$RateIns
table$RateSubPct <- table$RateSub / table$RateTotal
table$RateDelPct <- table$RateDel / table$RateTotal
table$RateInsPct <- table$RateIns / table$RateTotal

write.table( table, "table1.csv", quote=TRUE, sep=",", row.names=FALSE,
            col.names=c("Amplicon",
                "Substitution rate",
                "Deletion rate",
                "Insertion rate",
                "Total bases",
                "Total error rate",
                "Sub, %",
                "Del, %",
                "Ins, %") )



## ----- Table 2 ---------------------------------------------------------------
## Distribution of individaul error types for Taq DNA polymerase
## -----------------------------------------------------------------------------

table <- NULL

## Find Taq LacZ-1 sample
TaqSampleIx <- PcrSamplesIx & x$Enzyme == "Taq" & x$Amplicon == "LacZ-1"

taq <- x[TaqSampleIx,]

## Total number of substitutions
total <- sum(taq[c(subsAT,subsGC)])

## Compute error types
table["AG,TC"] <- sum(taq[c("AG","TC")]) / total
table["GA,CT"] <- sum(taq[c("GA","CT")]) / total
table["AT,TA"] <- sum(taq[c("AT","TA")]) / total
table["AC,TG"] <- sum(taq[c("AC","TG")]) / total
table["GC,CG"] <- sum(taq[c("GC","CG")]) / total
table["GT,CA"] <- sum(taq[c("GT","CA")]) / total

write.table( t(table), "table2.csv", quote=TRUE, sep=",", row.names=FALSE, col.names=TRUE )



## ----- Table 3 ---------------------------------------------------------------
## Substitution error rates measured by PacBio single-molecule sequencing
## -----------------------------------------------------------------------------

table <- NULL

pcr <- x[PcrSamplesIx,]

## Substitution rate, standard deviation
rate  <- aggregate(pcr["RateSub"], pcr["Enzyme"], mean)
sd    <- aggregate(pcr["RateSub"], pcr["Enzyme"], sd)
total <- aggregate(pcr["Total"],   pcr["Enzyme"], sum)
fidel <- aggregate(pcr["Fidelity"],pcr["Enzyme"], mean)
table <- merge( rate, sd, by="Enzyme" )

## Accuracy
table["Accuracy"]  <- 1 / table$RateSub.x

## Average fidelity
table <- merge( table, fidel, by="Enzyme" )

## Total bases sequenced
table <- merge( table, total, by="Enzyme" )
table <- table[order(table["RateSub.x"]),]

## Save table
write.table( table, "table3.csv", quote=TRUE, sep=",", row.names=FALSE,
            col.names=c("DNA Polymerase","Substitution Rate","SD","Accuracy","Fidelity","Total Bases") )



## ----- Table S1 --------------------------------------------------------------
## DNA Polymerase base substitution rates for individual amplicons
## -----------------------------------------------------------------------------

table <- x[PcrSamplesIx,c("Enzyme","Amplicon","RateSub","Doubling","Total","Fidelity")]

write.table( table, "tableS1.csv", quote=TRUE, sep=",", row.names=FALSE,
            col.names=c("DNA Polymerase","Amplicon","Substitution rate","Doubling","Total bases","Fidelity") )



## ----- Table S2 --------------------------------------------------------------
## Distribution of individual error types for DNA polymerases
## -----------------------------------------------------------------------------

## get PCR samples
pcr <- x[PcrSamplesIx,c("Enzyme",subsAT,subsGC)]

## group by enzyme
pcr <- aggregate( pcr[c(subsAT,subsGC)], pcr["Enzyme"], sum )

## compute error types
total <- apply( pcr[c(subsAT,subsGC)], 1, sum)

table <- pcr["Enzyme"]
table["AG,TC"] <- apply( pcr[c("AG","TC")], 1, sum) / total
table["GA,CT"] <- apply( pcr[c("GA","CT")], 1, sum) / total
table["AT,TA"] <- apply( pcr[c("AT","TA")], 1, sum) / total
table["AC,TG"] <- apply( pcr[c("AC","TG")], 1, sum) / total
table["GC,CG"] <- apply( pcr[c("GC","CG")], 1, sum) / total
table["GT,CA"] <- apply( pcr[c("GT","CA")], 1, sum) / total

write.table( table, "tableS2.csv", quote=TRUE, sep=",", row.names=FALSE, col.names=TRUE )



## ----- Table S3 --------------------------------------------------------------
## Number of individual error types in plasmid sequencing
## -----------------------------------------------------------------------------

plasmid <- x[PlasmidSamplesIx,c("Enzyme",subsAT,subsGC,"Total")]

plasmid["Subs"] <- apply( plasmid[c(subsAT,subsGC)], 1, sum)
plasmid["Rate"] <- plasmid["Subs"] / plasmid["Total"]

write.table( plasmid, "tableS3.csv", quote=TRUE, sep=",", row.names=FALSE, col.names=TRUE )



## ----- Table S4 --------------------------------------------------------------
## Polymerase substitution rates normalized by number of PCR cycles
## -----------------------------------------------------------------------------

## get PCR samples
pcr <- x[PcrSamplesIx,c("Enzyme","ErrorSub")]

## normalize by number of PCR cycles
pcr["ErrorSub"] <- pcr["ErrorSub"] / 16

## group by enzyme
mean <- aggregate( pcr["ErrorSub"], pcr["Enzyme"], mean )
sd   <- aggregate( pcr["ErrorSub"], pcr["Enzyme"], sd )

table <- merge(mean,sd,by="Enzyme")
table <- table[order(table["ErrorSub.x"]),]

## save table
write.table( table, "tableS4.csv", quote=TRUE, sep=",", row.names=FALSE,
            col.names=c("DNA Polymerase","Substitution rate","SD") )

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

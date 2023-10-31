#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#Extract variables from bash submission
INFILE=args[1]
OUTFILE=args[2]
NO_SITES=args[3]

library(dplyr)
VCF_positions <- read.delim(gzfile(INFILE), header = FALSE)
n <- NO_SITES
randomSubset <- VCF_positions[sample(nrow(VCF_positions), n), ]
randomSubset <- randomSubset[with(randomSubset, order(V1,V2)),]
write.table(randomSubset, OUTFILE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



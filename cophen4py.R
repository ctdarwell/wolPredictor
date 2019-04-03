library(ape)
rm(list=ls())

myArg <- commandArgs(trailingOnly = TRUE)

tree <- read.nexus(myArg)
phydist <- cophenetic(tree)
phydist <- phydist[order(rownames(phydist)), order(colnames(phydist))]

write.csv(phydist, file=paste(myArg, "_cophen.csv", sep=""), quote = F, row.names = F, col.names = F)

#NEED TO OUTPUT DIST CSV!!!!!!!!!!!!!!!!!!!!!!
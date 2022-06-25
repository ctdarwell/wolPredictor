if (!require("ape")) install.packages("ape")
library(ape)
rm(list=ls())

myArg <- commandArgs(trailingOnly = TRUE)

tree <- read.nexus(myArg)
phydist <- cophenetic(tree)
phydist <- phydist[order(rownames(phydist)), order(colnames(phydist))]

write.csv(as.matrix(phydist), stdout(), quote = F, row.names = F)

##################################
# R source code file used to analyse TDT results
# Created by Sahir, June 17, 2014
# Updated June 17, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

setwd("~/git_repositories/GAW19/data")

# TDT results -------------------------------------------------------------
DT.tdt <- as.data.table(read.table("everyone_affected.tdt", header=TRUE))
setkey(DT.tdt,SNP)

str(DT.tdt)

manhattan(DT.tdt)

DT.tdt[which.min(DT.tdt$P), ,][,"SNP",with=FALSE]
manhattan(DT.tdt, highlight="21_10742678")

qq(DT.tdt$P, main="Q-Q plot of TDT p-values")


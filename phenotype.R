##################################
# R source code file used to make the phenotype file for analysis
# Created by Sahir, June 16, 2014
# Updated June 16, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################


# Phenotype information ---------------------------------------------------

DT.phen <- fread("PHEN.csv")
phen.names <- colnames(DT.phen)
#extract family ID from individual ID
set(DT.phen, i=NULL, j="fam_id", value = gsub("(?<![0-9])0+","",substr(DT.phen[["t2dgID"]],5,6),perl=TRUE))
setcolorder(DT.phen, c("fam_id",phen.names[-30]))

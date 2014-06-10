##################################
# R source code file used to create variables used for GAW analysis
# Created by Sahir, March 27, 2014
# Updated March 27, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

setwd("~/git_repositories/GAW19/data/")
library(data.table)
library(bit64)
library(plotrix)


# Genotype Files ----------------------------------------------------------

filename <- "chr21-geno.csv"
filename <- "chr21t.csv"
DT.geno <- fread(filename)
DT.geno.t <- as.data.table(t(DT.geno[,2:ncol(DT.geno), with=FALSE]), keep.rownames=TRUE)
setnames(DT.geno.t,"rn","ID")
setkey(DT.geno.t,ID)


# PED file -------------------------------------------------------------

DT.ped <- fread("PED.csv")
setkey(DT.ped,ID)
tables()

DT.all <- DT.ped[DT.geno.t]

setcolorder(DT.all,c("PEDNUM","ID",colnames(DT.all)[-1:-2]))


source("http://bioconductor.org/biocLite.R")
biocLite("GWASTools")

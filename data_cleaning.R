##################################
# R source code file used to create variables used for GAW analysis
# Created by Sahir, March 27, 2014
# Updated June 10, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

#setwd("~/git_repositories/GAW19/data/")
library(data.table)
library(bit64)
library(plotrix)
library(stringr)
library(doParallel)
registerDoParallel(cores = 4)
library(foreach)



# Genotype Files ----------------------------------------------------------
#needs to be in the following 4 column formate
# chromosome number | snp identifier | genetic distance (cM) | base pair
filename <- "chr9-geno.csv"
filename <- "chr21-geno.csv"
DT.geno <- fread(filename, select=1) #, nrows=200000)
setkey(DT.geno, snp)

#extract chromosome and base pair
DT.geno$chr <- as.numeric(sub(" .*", "", gsub("\\_", " ", as.character(DT.geno$snp))))
DT.geno$base <- as.numeric(sub(" ", "", sub(DT.geno$chr[1],"",gsub("\\_", " ", as.character(DT.geno$snp)),perl=TRUE)))
DT.geno$cM <- 0
DT.geno[1:5,1:4,with=FALSE]
setcolorder(DT.geno, c("chr","snp","cM","base"))
write.table(DT.geno, file=paste("chr",DT.geno$chr[1],".map",sep=""),col.names=FALSE,
            row.names=FALSE, quote=FALSE)


#bring in SNP data
DT <- fread(filename, select=2)#, nrows=10) #, nrows=200000)
people <- colnames(DT)[-1]
setkey(DT,snp)

DT.m <- DT[DT.geno]
#DT.m[1:10,1:5,with=FALSE]
setcolorder(DT.m, c("chr","snp","cM","base",people))
rm(DT.geno,DT)
write.table(DT.m, file=paste("trans",DT.m$chr[1],".tped",sep=""),col.names=FALSE,
            row.names=FALSE, quote=FALSE)

read.fwf

# Break up genotypes ------------------------------------------------------
# filename <- "chr21-geno.csv"
# DT <- fread(filename, select=2:5, nrows=10) #, nrows=200000)
# DT[,snp:=NULL]
#DT[,c("v1","v2"):=list(substr(DT$T2DG0200001, 1, 1),substr(DT$T2DG0200001, 2, 2))]

#id <- colnames(DT)
#ids <- unlist(lapply(id, function(i) rep(i,2)))

#f <- function(x) {list(substr(x, 1, 1),substr(x, 2, 2))} 

# ids <- colnames(DT)
# gene.split <- function(i) 
#   {
#   as.data.table(do.call(rbind, strsplit(as.vector(eval(parse(text=paste("DT$",ids[i])))), split = "")))
# }
# 
# DT
# all.gene <- foreach(i=1:length(ids)) %dopar% gene.split(i)
# do.call(cbind,all.gene)


# Break up genotypes NEW - This works --------------------------------------------------

rm(DT)
filename <- "chr21-geno.csv.gz"
DT <- fread(filename)#, select=2:300)#, nrows=10) #, nrows=200000)
out_names <- paste("V", 1:(2*ncol(DT)), sep="_")
invar1 <- names(DT)
alloc.col(DT,3000)
for (i in seq_along(invar1)) {
  set(DT, i=NULL, j=out_names[2*i-1], value=substr(DT[[invar1[i]]],1,1))
  set(DT, i=NULL, j=out_names[2*i], value=substr(DT[[invar1[i]]],2,2))
}


system(paste("gunzip ",filename))

# DT[1:25,c(1,2,299:305), with=FALSE]
# for (i in seq_along(invar1)) {
#   set(DT, i=NULL, j=out_names[2*i-1], value=do.call(rbind, strsplit(DT[[invar1[i]]], split = ""))[,1])
#   set(DT, i=NULL, j=out_names[2*i], value=do.call(rbind, strsplit(DT[[invar1[i]]], split = ""))[,2])
# }

# DT[1:5,c(18,19,56:57),with=FALSE]


# method 1 from SO
# pairs <- c(outer(c('C', 'T', 'G', 'A'), c('C', 'T', 'G', 'A'), 'paste0'))
# l <- replicate(1000, sample(pairs, 260000, replace=TRUE), simplify=FALSE)



# PED file -------------------------------------------------------------

# leave this one as is, this is in the proper format as per
# http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr 
DT.ped <- fread("PED.csv")
DT.ped
write.table(DT.ped, file="trans.tfam",col.names=FALSE,row.names=FALSE, quote=FALSE)


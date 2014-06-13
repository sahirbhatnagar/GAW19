##################################
# R source code file used to split genotypes into two columns
# Created by Sahir, June 12, 2014
# Updated June 13, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

LIB_LOC = "~/share/greenwood.group/Rlibs/sahir.bhatnagar"

#LIB_LOC = "/home/sahir/R/x86_64-pc-linux-gnu-library/3.1"
#setwd("~/git_repositories/GAW19/data")

library(data.table, lib.loc=LIB_LOC)
library(bit, lib.loc=LIB_LOC)
library(bit64, lib.loc=LIB_LOC)
library(plotrix, lib.loc=LIB_LOC)
library(stringr, lib.loc=LIB_LOC)
library(doParallel, lib.loc=LIB_LOC)
registerDoParallel(cores = 4)
library(foreach, lib.loc=LIB_LOC)

filename <- "~/share/sy/GAW19/data/chr21-geno.csv"
#filename = "chr21-geno.csv"

#file.names <- list.files(pattern = "*geno.csv")

# function to create the data files in proper format
# this function will create the .tped file in the format:
# chromosome number | snp identifier | genetic distance (cM) | base pair | Subject Genotypes (each allele in individual column)
# and the map file in the format: chromosome number | snp identifier | genetic distance (cM) | base pair
data.clean <- function(filename, select.all=FALSE, nrows=10, ncols=10) {
    #filename: path and filename or just the filename of genotype file
    #select: use all the data or a subset of it, if "all" then nrows and ncols is ignored
    #nrows: # of rows to select
    #ncols: # of columns to select
  
    #read in SNPs
    DT <- if (select.all) fread(filename) else fread(filename, select=1:ncols, nrows=nrows) 
    
    #remove column of SNP Id's because we only want to split the genotypes 
    DT[,snp:=NULL]
    
    #used for deleting columns of grouped SNPs after
    n.col.DT <- ncol(DT)
    
    #read in SNP ID's
    DT.snp <- if (select.all) fread(filename, select=1) else fread(filename, select=1, nrows=nrows)
    #set the keys 
    setkey(DT.snp, snp)
    
    out_names <- paste("V", 1:(2*ncol(DT)), sep="_")
    invar1 <- names(DT)
    if (!select.all & ncols>50 | select.all) alloc.col(DT,3000)
    
    for (i in seq_along(invar1)) {
          set(DT, i=NULL, j=out_names[2*i-1], value=substr(DT[[invar1[i]]],1,1))
          set(DT, i=NULL, j=out_names[2*i], value=substr(DT[[invar1[i]]],2,2))
    }

    #remove grouped SNPs
    DT[,seq.int(1,n.col.DT,1):=NULL, with=FALSE]

    #extract chromosome and base pair
    DT.snp$chr <- as.numeric(sub(" .*", "", gsub("\\_", " ", as.character(DT.snp$snp))))
    DT.snp$base <- as.numeric(sub(" ", "", sub(DT.snp$chr[1],"",gsub("\\_", " ", as.character(DT.snp$snp)),perl=TRUE)))
    DT.snp$cM <- 0
    
    setcolorder(DT.snp, c("chr","snp","cM","base"))
    write.table(DT.snp, file=paste("chr",DT.snp$chr[1],".map",sep=""),col.names=FALSE,
                row.names=FALSE, quote=FALSE)

    DT.final <- cbind(DT.snp,DT)

    write.table(DT.final, file=paste("chr",DT.snp$chr[1],".tped",sep=""), col.names=FALSE,
                row.names=FALSE, quote=FALSE)

}

data.clean(filename, select.all=TRUE)

# PED file -------------------------------------------------------------

# leave this one as is, this is in the proper format as per
# http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr 
#DT.ped <- fread("PED.csv")
#write.table(DT.ped, file="trans.tfam",col.names=FALSE,row.names=FALSE, quote=FALSE)
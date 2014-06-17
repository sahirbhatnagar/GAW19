##################################
# R source code file used to split genotypes into two columns
# Created by Sahir, June 12, 2014
# Updated June 14, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

#LIB_LOC = "~/share/greenwood.group/Rlibs/sahir.bhatnagar"

LIB_LOC = "/home/sahir/R/x86_64-pc-linux-gnu-library/3.1"
setwd("~/git_repositories/GAW19/data")

library(data.table, lib.loc=LIB_LOC)
library(bit, lib.loc=LIB_LOC)
library(bit64, lib.loc=LIB_LOC)
library(plotrix, lib.loc=LIB_LOC)
library(stringr, lib.loc=LIB_LOC)
library(doParallel, lib.loc=LIB_LOC)
registerDoParallel(cores = 4)
library(foreach, lib.loc=LIB_LOC)

#filename <- "~/share/sy/GAW19/data/chr21-geno.csv"
filename = "chr21-geno.csv"

#file.names <- list.files(pattern = "*geno.csv")

# function to create the data files in proper format
# this function will create the .tped file in the format:
# chromosome number | snp identifier | genetic distance (cM) | base pair | Subject Genotypes (each allele in individual column)
# and the map file in the format: chromosome number | snp identifier | genetic distance (cM) | base pair
# use this for creating the file used in PLINK (the output of this function will include
# the every individual in the families, including those 430 for which there is no genotype data)
data.clean <- function(filename, select.all=FALSE, nrows=10, ncols=10, col.alloc=100) {
    #filename: path and filename or just the filename of genotype file
    #select: use all the data or a subset of it, if "all" then nrows and ncols is ignored
    #nrows: # of rows to select
    #ncols: # of columns to select
    #col.alloc: # # of columns to allocate to data.table::alloc.col 
    #select.all=FALSE; nrows=5; ncols=960
  
    #read in SNPs
    DT <- if (select.all) fread(filename) else fread(filename, select=1:ncols, nrows=nrows) 
    
    #remove column of SNP Id's because we only want to split the genotypes 
    DT[,snp:=NULL]
    
    #get subject ID's for which we have genotype information
    id.geno <- colnames(DT)
    
    # Differences between genotype and family file 
    # read in family file
    DT.tfam <- fread("PED.csv")
    fam.id <- DT.tfam$ID
    
    #ID's with missing genotypes
    missing.geno <- setdiff(fam.id, id.geno)
    
    alloc.col(DT,col.alloc)
    for (i in seq_along(missing.geno)) {
      set(DT, i=NULL, j=missing.geno[i], value="XX")
    }
    
    setcolorder(DT, fam.id)
    
    #used for deleting columns of grouped SNPs after
    n.col.DT <- ncol(DT)
    
    #out_names <- paste("V", 1:(2*ncol(DT)), sep="_")
    out_names=NULL
    for (i in 1:length(fam.id)){
      out_names <- c(out_names, paste(fam.id[i],1,sep="_"), paste(fam.id[i],2,sep="_"))
    }
        
    invar1 <- names(DT)
    alloc.col(DT, col.alloc)
    for (i in seq_along(invar1)) {
          set(DT, i=NULL, j=out_names[2*i-1], value=substr(DT[[invar1[i]]],1,1))
          set(DT, i=NULL, j=out_names[2*i], value=substr(DT[[invar1[i]]],2,2))
    }

    #read in SNP ID's
    DT.snp <- if (select.all) fread(filename, select=1) else fread(filename, select=1, nrows=nrows)
    #set the keys 
    setkey(DT.snp, snp)
    
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
    
    return(DT.final)

}


# use this function for local calculations only (this will not include missing genotype data)
data.clean.incomplete <- function(filename, select.all=FALSE, select.all.cols=FALSE, nrows, ncols, col.alloc=4000) {
  #filename: path and filename or just the filename of genotype file
  #select: use all the data or a subset of it, if "all" then nrows and ncols is ignored
  #nrows: # of rows to select
  #ncols: # of columns to select
  #col.alloc: # # of columns to allocate to data.table::alloc.col 
  #select.all=FALSE; nrows=5; ncols=960
  
  #read in SNPs
  DT <- if (select.all) fread(filename) else if (select.all.cols) fread(filename, nrows=nrows) else fread(filename, select=1:ncols, nrows=nrows) 
  
  #remove column of SNP Id's because we only want to split the genotypes 
  DT[,snp:=NULL]
  
  #get subject ID's for which we have genotype information
  id.geno <- colnames(DT)
  
#   # Differences between genotype and family file 
#   # read in family file
#   DT.tfam <- fread("PED.csv")
#   fam.id <- DT.tfam$ID
#   
#   #ID's with missing genotypes
#   missing.geno <- setdiff(fam.id, id.geno)
#   
#   alloc.col(DT,col.alloc)
#   for (i in seq_along(missing.geno)) {
#     set(DT, i=NULL, j=missing.geno[i], value="XX")
#   }
#   
#   setcolorder(DT, fam.id)
  
  #used for deleting columns of grouped SNPs after
  n.col.DT <- ncol(DT)
  
  #out_names <- paste("V", 1:(2*ncol(DT)), sep="_")
  out_names=NULL
  for (i in 1:length(id.geno)){
    out_names <- c(out_names, paste(id.geno[i],1,sep="_"), paste(id.geno[i],2,sep="_"))
  }
  
  invar1 <- names(DT)
  alloc.col(DT, col.alloc)
  for (i in seq_along(invar1)) {
    set(DT, i=NULL, j=out_names[2*i-1], value=substr(DT[[invar1[i]]],1,1))
    set(DT, i=NULL, j=out_names[2*i], value=substr(DT[[invar1[i]]],2,2))
  }
  
  #read in SNP ID's
  #DT.snp <- if (select.all) fread(filename, select=1) else fread(filename, select=1, nrows=nrows)
  DT.snp <- if (select.all) fread(filename, select=1) else fread(filename, select=1, nrows=nrows) 

  #set the keys 
  setkey(DT.snp, snp)
  
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
  
  return(DT.final)
  
}



#DT <- fread(filename, nrows=5) 

#dat <- data.clean(filename, select.all=FALSE, nrows=239352, ncols=2, col.alloc=3000)


DT <- data.clean.incomplete(filename, select.all=FALSE, select.all.cols=TRUE, nrows=6, col.alloc=3025)
warnings()
f=1
#################################################################3
DT[snp=="21_9411318",1:7,with=FALSE]



do.call(rbind,apply(DT[,-1:-4, with=FALSE],1,table))

library(plyr)

ddply(DT, "snp", function(x) {
  apply(DT[f,5:12, with=FALSE],1,function(i) sum(i=="T"))
  apply(DT[f,-1:-4, with=FALSE],1,function(i) sum(i=="T"))
  apply(DT[f,-1:-4, with=FALSE],1,function(i) sum(i=="G"))
  apply(DT[f,-1:-4, with=FALSE],1,function(i) sum(i=="C"))
  apply(DT[f,-1:-4, with=FALSE],1,function(i) sum(i=="X"))
  data.frame(cv.count = cv)
})

example(data.table)
# PED file -------------------------------------------------------------

# leave this one as is, this is in the proper format as per
# http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr 
DT.tfam <- fread("PED.csv")
write.table(DT.tfam, file="trans.tfam",col.names=FALSE,row.names=FALSE, quote=FALSE)




DT<-fread("chr21.tped", )



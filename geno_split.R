##################################
# R source code file used to split genotypes into two columns
# Created by Sahir, June 12, 2014
# Updated June 20, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

# LIB_LOC = "~/share/greenwood.group/Rlibs/sahir.bhatnagar"
# 
# #LIB_LOC = "/home/sahir/R/x86_64-pc-linux-gnu-library/3.1"
# #setwd("~/git_repositories/GAW19/data")
# 
# library(data.table, lib.loc=LIB_LOC)
# library(bit, lib.loc=LIB_LOC)
# library(bit64, lib.loc=LIB_LOC)
# library(plotrix, lib.loc=LIB_LOC)
# library(stringr, lib.loc=LIB_LOC)
# library(doParallel, lib.loc=LIB_LOC)
# registerDoParallel(cores = 10)
# library(foreach, lib.loc=LIB_LOC)


#filename <- "~/share/sy/GAW19/data/chr21-geno.csv"
#filename = "chr21-geno.csv"

#file.names <- list.files(pattern = "*geno.csv")

# function to create the data files in proper format
# this function will create the .tped file in the format:
# chromosome number | snp identifier | genetic distance (cM) | base pair | Subject Genotypes (each allele in individual column)
# and the map file in the format: chromosome number | snp identifier | genetic distance (cM) | base pair
# use this for creating the file used in PLINK (the output of this function will include
# the every individual in the families, including those 430 for which there is no genotype data)
data.clean <- function(filename, select.all=FALSE, nrows=10, ncols=10, col.alloc=100, chr.digits=2) {
    #filename: path and filename or just the filename of genotype file
    #select: use all the data or a subset of it, if "all" then nrows and ncols is ignored
    #nrows: # of rows to select
    #ncols: # of columns to select
    #col.alloc: # # of columns to allocate to data.table::alloc.col 
    #chr.digits: # of digits in chromosome number
    #select.all=FALSE; nrows=400; ncols=960; filename = "chr21-geno.csv";col.alloc=4500;chr.digits=2
  
    #read in SNPs
    DT <- if (select.all) fread(filename) else fread(filename, select=1:ncols, nrows=nrows) 
    setkey(DT,snp)
    
    #get subject ID's for which we have genotype information
    id.geno <- colnames(DT)[-1]
    
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
    
    setcolorder(DT, c("snp",fam.id))
    
    #used for deleting columns of grouped SNPs after
    n.col.DT <- ncol(DT)-1
    
    out_names <- paste("V", 1:(2*n.col.DT), sep="_")
    invar1 <- names(DT)[-1]

    cat("Starting to split genotypes\n")

    for (i in seq_along(invar1)) {
          set(DT, i=NULL, j=out_names[2*i-1], value=substr(DT[[invar1[i]]],1,1))
          set(DT, i=NULL, j=out_names[2*i], value=substr(DT[[invar1[i]]],2,2))
    }
    
    cat("Splitting genotypes complete\n")

    #remove grouped SNPs
    DT[,seq.int(2,n.col.DT+1,1):=NULL, with=FALSE]
  
    cat("Starting to extract chromosome position\n")

    set(DT, i=NULL, j="chr", value= substr(DT[["snp"]],1,chr.digits))
    set(DT, i=NULL, j="base", value = if (chr.digits==2) sub("[0-9][0-9]_","",DT[["snp"]]) else sub("[0-9]_","",DT[["snp"]]) )
    set(DT, i=NULL, j="cM", value=0) 
    
    cat("Finished extracting chromosome position")

    setcolorder(DT, c("chr","snp","cM","base",out_names))  
    
#    setwd(paste(getwd(),"/chr",DT.snp$chr[1],sep=""))

#     write.table(as.matrix(DT), file=paste("chr",DT[["chr"]][1],".tped",sep=""), col.names=FALSE,
#                              row.names=FALSE, quote=FALSE)

    #rm(DT.final)

#   write.table(DT.tfam, file="trans.tfam",col.names=FALSE,row.names=FALSE, quote=FALSE)
    
    return(DT)
#     system(paste("mkdir ","chr",DT[["chr"]][1], sep=""))
#     system(paste("mv ","chr",DT[["chr"]][1],".tped ",getwd(),"/chr",DT[["chr"]][1],sep=""))
#     system(paste("mv ","trans.tfam ",getwd(),"/chr",DT[["chr"]][1],sep=""))
#     system(paste("cp mztwins.txt ",getwd(),"/chr",DT[["chr"]][1],sep=""))
#     system(paste("cp plink.sh ",getwd(),"/chr",DT[["chr"]][1],sep=""))
#     system(paste("./plink.sh ",DT.snp$chr[1],sep=""))

}


#raw_data <- list.files(pattern="*-geno.csv")
#data.clean("chr9-geno.csv", select.all=FALSE, nrows=5,ncols=950, col.alloc=4500)
#l <- data.clean("chr21-geno.csv", select.all=FALSE, nrows=50, ncols=960, col.alloc=4500, chr.digits=2)

#command line args should be --filename chr.digits 
args <- commandArgs(trailingOnly=TRUE)

DT <- data.clean(args[1], select.all=TRUE, col.alloc=4500, chr.digits=args[2]) 

write.table(DT, file=tempfile(pattern=paste("chr",DT[["chr"]][1],".tped",sep=""), tmpdir=getwd()) ,
            col.names=FALSE, row.names=FALSE, quote=FALSE)


# rm(DT.final)
# 
#   write.table(DT.tfam, file="trans.tfam",col.names=FALSE,row.names=FALSE, quote=FALSE)
# 
# 
# DT <- data.clean(commandArgs(trailingOnly=TRUE), select.all=TRUE, col.alloc=4500, chr.digits=1)
# 



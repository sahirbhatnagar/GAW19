##################################
# R source code file used to split genotypes into two columns
# Created by Sahir, June 12, 2014
# Updated 
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

LIB_LOC = "~/share/greenwood.group/Rlibs/sahir.bhatnagar"

library(bit, lib.loc=LIB_LOC)
library(bit64, lib.loc=LIB_LOC)
library(plotrix, lib.loc=LIB_LOC)
library(stringr, lib.loc=LIB_LOC)
library(doParallel, lib.loc=LIB_LOC)
registerDoParallel(cores = 4)
library(foreach, lib.loc=LIB_LOC)
library(data.table, lib.loc=LIB_LOC)

filename <- "~/share/sy/GAW19/data/chr21-geno.csv"
DT <- fread(filename)#, select=2:300)#, nrows=10) #, nrows=200000)
out_names <- paste("V", 1:(2*ncol(DT)), sep="_")
invar1 <- names(DT)
alloc.col(DT,3000)
for (i in seq_along(invar1)) {
  set(DT, i=NULL, j=out_names[2*i-1], value=substr(DT[[invar1[i]]],1,1))
  set(DT, i=NULL, j=out_names[2*i], value=substr(DT[[invar1[i]]],2,2))
}
write.table(DT, file="chr21.tped")
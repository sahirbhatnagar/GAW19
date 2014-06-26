##################################
# R source code file that will read in HWE and MAF file, and output HWE manhattan and 
# HWE vs MAF plots
# Created by Sahir, June 26, 2014
# Updated June 26, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

#command line args should be frq_file hwe_file 
args <- commandArgs(trailingOnly=TRUE)


# -log10(HWE pvalue) vs MAF -----------------------------------------------

DT.frq <- fread(args[1])
DT.hwe <- fread(args[2]) 
head(DT.frq);head(DT.hwe)
DT.hwe <- DT.hwe[V3=="ALL"]
setkey(DT.hwe,V2);setkey(DT.frq, SNP)
DT.all <- DT.hwe[DT.frq]
DT.all[order(DT.all$V9)[1:3000]]
#DT.all2 <- DT.hwe2[DT.frq]


xyplot(-log10(V9) ~ MAF, DT.all, grid = TRUE, main=paste("chromosome ",chromosome, sep=" "))

##################################
# R source code file for manhattan plots
# Created by Sahir, June 5, 2014
# Updated July 6, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: this is different from geno_split.R in that it no longer needs to split genotypes
# since PLINK 1.90 can handle grouped alleles for regular ped files only
##################################
# PDT ---------------------------------------------------------------------

#command line args should be chromosome number 
args <- commandArgs(trailingOnly=TRUE)

#all
DT.pdt.all <- fread(paste("chr",args[1],"/pdt",args[1],".csv", sep=""), sep=";")
set(DT.pdt.all, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.pdt.all[["V1"]])))
set(DT.pdt.all, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.pdt.all[["V1"]])) )
DT.pdt.all <- DT.pdt.all[V5!="bad_data"]
setnames(DT.pdt.all, c("V1", "V5"), c("SNP","P"))
set(DT.pdt.all, i=NULL, j="Pval", value= as.numeric(DT.pdt.all[["P"]])) 
png(paste("chr", args[1],"/PDT_chr",args[1],"all.png",sep=""));
manhattan(DT.pdt.all, p="Pval", main=paste("PDT, chr ",args[1],", all", sep=""));
dev.off()

#sequenced
DT.pdt.seq <- fread(paste("chr",args[1],"/pdt",args[1],"seq.csv", sep=""), sep=";")
set(DT.pdt.seq, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.pdt.seq[["V1"]])))
set(DT.pdt.seq, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.pdt.seq[["V1"]])) )
DT.pdt.seq <- DT.pdt.seq[V5!="bad_data"]
setnames(DT.pdt.seq, c("V1", "V5"), c("SNP","P"))
set(DT.pdt.seq, i=NULL, j="Pval", value= as.numeric(DT.pdt.seq[["P"]])) 
png(paste("chr", args[1],"/PDT_chr",args[1],"seq.png",sep=""));
manhattan(DT.pdt.seq, p="Pval", main=paste("PDT, chr ",args[1],", seq only", sep=""));
dev.off()


#nuclear families
DT.pdt.nuc <- fread(paste("chr",args[1],"/pdt",args[1],"nuc.csv", sep=""), sep=";")
set(DT.pdt.nuc, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.pdt.nuc[["V1"]])))
set(DT.pdt.nuc, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.pdt.nuc[["V1"]])) )
DT.pdt.nuc <- DT.pdt.nuc[V5!="bad_data"]
setnames(DT.pdt.nuc, c("V1", "V5"), c("SNP","P"))
set(DT.pdt.nuc, i=NULL, j="Pval", value= as.numeric(DT.pdt.nuc[["P"]])) 
png(paste("chr", args[1],"/PDT_chr",args[1],"nuc.png",sep=""));
manhattan(DT.pdt.nuc, p="Pval", main=paste("PDT, chr ",args[1],", nuclear families only", sep=""));dev.off()
dev.off()

rm(list=ls())

# # Manhattan all -----------------------------------------------------------
# 
# setwd("~/share/sy/GAW19/PDT/")
# DT <- fread("all.pdt", sep=";")
# head(DT)
# 
# set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
# set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
# DT <- DT[V5!="bad_data"]
# setnames(DT, c("V1", "V5"), c("SNP","P"))
# set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
# 
# 
# png("PDT_all.png");manhattan(DT, p="Pval", main="PDT, All chr, seq only");dev.off()
# 
# DT[order(DT$Pval)[1:5]]

# TDT all ---------------------------------------------------------------------
DT <- fread(paste("chr",args[1],"/chr",args[1],"all.csv", sep=""), sep=";")
png(paste("chr", args[1],"/TDT_chr",args[1],"all.png",sep=""));
manhattan(DT, main=paste("TDT, chr ",args[1],", all", sep=""));dev.off()
rm(DT)

# TDT seq ---------------------------------------------------------------------
DT <- fread(paste("chr",args[1],"/chr",args[1],"seq.csv", sep=""), sep=";")
png(paste("chr", args[1],"/TDT_chr",args[1],"seq.png",sep=""));
manhattan(DT, main=paste("TDT, chr ",args[1],", seq only", sep=""));dev.off()

# TDT nuc ---------------------------------------------------------------------
DT <- fread(paste("chr",args[1],"/chr",args[1],"nuc.csv", sep=""), sep=";")
png(paste("chr", args[1],"/TDT_chr",args[1],"nuc.png",sep=""));
manhattan(DT, main=paste("TDT, chr ",args[1],", nuclear families only", sep=""));dev.off()



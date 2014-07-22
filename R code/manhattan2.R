##################################
# R source code file used to create manhattan plots
# Created by Sahir, July 6, 2014
# Updated 
# 
# NOTE: run this from GAW19 folder
##################################

#command line args should be PDT/TDT and ALL/NUC/SEQ 
#args <- commandArgs(trailingOnly=TRUE)

#setwd("~/share/sy/GAW19/REAL")
#source("~/share/sy/GAW19/scripts/manhattan2.R")


#args<-c("PDT","ALL")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
if (args[1]=="PDT") {
  set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
  set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
  DT <- DT[V5!="bad_data"]
  setnames(DT, c("V1", "V5"), c("SNP","P"))
  set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
  png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
  manhattan(DT, p="Pval", main=paste(args[1],"_",args[2],"_REAL_Pheno.png",sep=""));
  dev.off() } else if (args[1]=="TDT") {
    png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
    manhattan(DT, main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
    dev.off()
  } 



head(DT)
DT[is.na(P)]

# rm(DT)
# 
# 
# args<-c("PDT","SEQ")
# DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
# if (args[1]=="PDT") {
#   set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
#   set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
#   DT <- DT[V5!="bad_data"]
#   setnames(DT, c("V1", "V5"), c("SNP","P"))
#   set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
#   png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#   manhattan(DT, p="Pval", main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#   dev.off() } else if (args[1]=="TDT") {
#     png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#     manhattan(DT, main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#     dev.off()
#   } 
# 
# rm(DT)
# 
# args<-c("PDT","NUC")
# DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
# if (args[1]=="PDT") {
#   set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
#   set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
#   DT <- DT[V5!="bad_data"]
#   setnames(DT, c("V1", "V5"), c("SNP","P"))
#   set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
#   png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#   manhattan(DT, p="Pval", main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#   dev.off() } else if (args[1]=="TDT") {
#     png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#     manhattan(DT, main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#     dev.off()
#   } 
# 
# rm(DT)
# 
# args<-c("TDT","ALL")
# DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
# if (args[1]=="PDT") {
#   set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
#   set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
#   DT <- DT[V5!="bad_data"]
#   setnames(DT, c("V1", "V5"), c("SNP","P"))
#   set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
#   png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#   manhattan(DT, p="Pval", main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#   dev.off() } else if (args[1]=="TDT") {
#     png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#     manhattan(DT, main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#     dev.off()
#   } 
# 
# rm(DT)
# 
# args<-c("TDT","SEQ")
# DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
# if (args[1]=="PDT") {
#   set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
#   set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
#   DT <- DT[V5!="bad_data"]
#   setnames(DT, c("V1", "V5"), c("SNP","P"))
#   set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
#   png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#   manhattan(DT, p="Pval", main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#   dev.off() } else if (args[1]=="TDT") {
#     png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#     manhattan(DT, main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#     dev.off()
#   } 
# 
# rm(DT)
# 
# args<-c("TDT","NUC")
# DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
# if (args[1]=="PDT") {
#   set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
#   set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
#   DT <- DT[V5!="bad_data"]
#   setnames(DT, c("V1", "V5"), c("SNP","P"))
#   set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
#   png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#   manhattan(DT, p="Pval", main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#   dev.off() } else if (args[1]=="TDT") {
#     png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan_REAL.png",sep=""));
#     manhattan(DT, main=paste(args[1],"_",args[2],"_REAL_Pheno",sep=""));
#     dev.off()
#   } 
# 
# rm(DT)


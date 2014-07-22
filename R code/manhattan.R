##################################
# R source code file used to create manhattan plots
# Created by Sahir, July 6, 2014
# Updated 
# 
# NOTE: run this from GAW19 folder
##################################

#command line args should be PDT/TDT and ALL/NUC/SEQ 
#args <- commandArgs(trailingOnly=TRUE)


##################################
#             PDT                #
##################################


setwd("~/share/sy/GAW19/")

args<-c("PDT","ALL")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
DT <- DT[V5!="bad_data"]
setnames(DT, c("V1", "V5"), c("SNP","P"))
set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
manhattan(DT, p="Pval", main="PDT: All Individuals", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off() 

rm(DT)


args<-c("PDT","SEQ")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
  set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
  set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
  DT <- DT[V5!="bad_data"]
  setnames(DT, c("V1", "V5"), c("SNP","P"))
  set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
  png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
  manhattan(DT, p="Pval", main="PDT: Sequenced Only", suggestiveline=FALSE, genomewideline=-log10(1e-03));
  dev.off() 

rm(DT)

args<-c("PDT","NUC")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
  set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["V1"]])))
  set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["V1"]])) )
  DT <- DT[V5!="bad_data"]
  setnames(DT, c("V1", "V5"), c("SNP","P"))
  set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
  png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
  manhattan(DT, p="Pval", main="PDT: 1 Nuclear Family per Pedigree", suggestiveline=FALSE, genomewideline=-log10(1e-03));
  dev.off() 

rm(DT)

args<-c("PDT","BEST")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["LOCUS_NAME"]])))
set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["LOCUS_NAME"]])) )
#DT <- DT[V5!="bad_data"]
setnames(DT, c("LOCUS_NAME", "P-VALUE"), c("SNP","P"))
set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
manhattan(DT, p="Pval", main="PDT: Best Imputed Within 1%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off() 

rm(DT)

args<-c("PDT","BEST")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["LOCUS_NAME"]])))
set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["LOCUS_NAME"]])) )
#DT <- DT[V5!="bad_data"]
setnames(DT, c("LOCUS_NAME", "P-VALUE"), c("SNP","P"))
set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
manhattan(DT, p="Pval", main="PDT: Best Imputed Within 1%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off() 

rm(DT)

args<-c("PDT","BEST5")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["LOCUS_NAME"]])))
set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["LOCUS_NAME"]])) )
#DT <- DT[V5!="bad_data"]
setnames(DT, c("LOCUS_NAME", "P-VALUE"), c("SNP","P"))
set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
manhattan(DT, p="Pval", main="PDT: Best Imputed Within 5%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off() 

rm(DT)

args<-c("PDT","BEST10")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["LOCUS_NAME"]])))
set(DT, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT[["LOCUS_NAME"]])) )
#DT <- DT[V5!="bad_data"]
setnames(DT, c("LOCUS_NAME", "P-VALUE"), c("SNP","P"))
set(DT, i=NULL, j="Pval", value= as.numeric(DT[["P"]])) 
png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
manhattan(DT, p="Pval", main="PDT: Best Imputed Within 10%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off() 

rm(DT)

##################################
#             TDT                #
##################################

args<-c("TDT","ALL")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
    png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
    manhattan(DT, main="TDT: All Individuals", suggestiveline=FALSE, genomewideline=-log10(1e-03));
    dev.off()

rm(DT)

args<-c("TDT","SEQ")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
    png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
    manhattan(DT, main="TDT: Sequenced Only", suggestiveline=FALSE, genomewideline=-log10(1e-03));
    dev.off()

rm(DT)

args<-c("TDT","NUC")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
    png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
    manhattan(DT, main="TDT: 1 Nuclear Family per Pedigree", suggestiveline=FALSE, genomewideline=-log10(1e-03));
    dev.off()

rm(DT)

args<-c("TDT","BEST5")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
manhattan(DT, main="TDT: Best Imputed Within 5%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()

rm(DT)

args<-c("TDT","BEST10")
DT <- fread(paste(args[1], "/",args[2],"/",args[2],".",args[1],sep=""), sep=";")
png(paste(args[1], "/",args[2],"/",args[1],"_",args[2],"_Manhattan.png",sep=""));
manhattan(DT, main="TDT: Best Imputed Within 10%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()

rm(DT)


##################################
#             FBAT              #
##################################

setwd("~/share/sy/GAW19/FBAT/ALL")

DT.fbat <- fread("ALL.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_ALL_Manhattan.png");
manhattan(DT.fbat, snp="Marker", main="FBAT: All Individuals", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()

rm(DT.fbat)
setwd("~/share/sy/GAW19/FBAT/SEQ")

DT.fbat <- fread("SEQ.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_SEQ_Manhattan.png");
manhattan(DT.fbat, snp="Marker", main="FBAT: Sequenced Only", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()

rm(DT.fbat)
setwd("~/share/sy/GAW19/FBAT/NUC")

DT.fbat <- fread("NUC.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_NUC_Manhattan.png");
manhattan(DT.fbat, snp="Marker", main="FBAT: 1 Nuclear Family per Pedigree", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()

DT.fbat <- fread("BEST.FBAT", sep=";")
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_BEST_Manhattan.png");
manhattan(DT.fbat, snp="Marker", main="FBAT: Best Imputed Within 1%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()

DT.fbat <- fread("BEST5.FBAT", sep=";")
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_BEST5_Manhattan.png");
manhattan(DT.fbat, snp="Marker", main="FBAT: Best Imputed Within 5%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()

DT.fbat <- fread("BEST10.FBAT", sep=";")
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_BEST10_Manhattan.png");
manhattan(DT.fbat, snp="Marker", main="FBAT: Best Imputed Within 10%", suggestiveline=FALSE, genomewideline=-log10(1e-03));
dev.off()


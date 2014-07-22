##################################
# R source code file that will read in HWE and MAF file, and output HWE manhattan and 
# HWE vs MAF plots
# Created by Sahir, June 26, 2014
# Updated July 16, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

#command line args should be frq_file hwe_file 
#args <- commandArgs(trailingOnly=TRUE)

# -log10(HWE pvalue) vs MAF -----------------------------------------------

# freqname <- "~/share/sy/GAW19/MAF/all/ALL.MAF"
 hwename <- "~/share/sy/GAW19/HWE/ALL/ALL.HWE"
hwename <- "~/share/sy/GAW19/HWE/SEQ/SEQ.HWE"



#mafhweplots <- function (freqname, hwename){
#DT.frq <- fread(freqname)
DT.hwe <- fread(hwename, sep=";") 
head(DT.hwe)
DT.hwe <- DT.hwe[V3=="ALL"]
#head(DT.frq)
setnames(DT.hwe, c("V1","V2","V9"), c("CHR", "SNP", "P"))
#set(DT.hwe, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.hwe[["SNP"]])))
ll <- unlist(strsplit(DT.hwe$V6, "/", fixed=TRUE))
idx <- seq(1, length(ll), by=3)
DT.hwe[,`:=`(aa=ll[idx], Aa=ll[idx+1], AA=ll[idx+2])]
setkey(DT.hwe,SNP);
set(DT.hwe, i=NULL, j="a.a",value=as.numeric(DT.hwe[["aa"]]))
set(DT.hwe, i=NULL, j="A.a",value=as.numeric(DT.hwe[["Aa"]]))
set(DT.hwe, i=NULL, j="A.A",value=as.numeric(DT.hwe[["AA"]]))
str(DT.hwe)
head(DT.hwe)

#setkey(DT.frq, SNP)
#DT.all <- DT.hwe[DT.frq]
#return(DT.all)
#DT.all[order(DT.all$V9)[1:3000]]
#DT.all2 <- DT.hwe2[DT.frq]
#xyplot(-log10(V9) ~ MAF, DT.all, grid = TRUE, main=paste("chromosome ",chromosome, sep=" "))
}

# head(DT.all)
# getwd()
# setwd("/home/data1/share/sy/GAW19/MAF_HWE")
# png("MAF_HWE_ALL.png")
# DT.all[,plot(MAF,-log10(P))];dev.off()
# str(DT.all)


a<-ggplot(DT.hwe, aes(A.A)) + geom_bar()
b<-ggplot(DT.hwe, aes(A.a)) + geom_bar()
c<-ggplot(DT.hwe, aes(a.a)) + geom_bar()

source("~/share/sy/GAW19/scripts/multiplot.R")
png("geno_hist_seq.png");multiplot(a,b,c,cols=3);dev.off()

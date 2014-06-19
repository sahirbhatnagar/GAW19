##################################
# R source code file used to make plots
# Created by Sahir, June 16, 2014
# Updated June 16, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

setwd("~/git_repositories/GAW19/data/")
rm(list=ls())
DT.hwe <- DT.hwe[TEST=="ALL"]
#system(paste("mkdir", DT.frq[["CHR"]][1], sep=" "))

DT.frq <- as.data.table(read.table("chr21.frq", header=TRUE))
DT.hwe <- as.data.table(read.table("chr21.hwe", header=TRUE))

DT.hwe2 <- as.data.table(read.table("chr21.hwe2.hwe", header=TRUE))
DT.hwe2 <- DT.hwe2[TEST=="ALL"]


setkey(DT.hwe2,SNP)
#setkey(DT.hwe,NULL)
setkey(DT.frq, SNP)
DT.all <- DT.hwe2[DT.frq]
DT.all[MAF!=0][,plot(MAF,-log10(P))]
#rm(DT.all)

DT.all <- cbind(DT.frq,DT.hwe["ALL"])
DT.all[,log.expect.p:=-log10((1:nrow(DT.all)-0.5)/nrow(DT.all))]
DT.all[,log.sorted.p:=-log10(sort(P))]

#data for founders only
DT.frq.found <- as.data.table(read.table("founders.frq", header=TRUE))
DT.hwe.found <- as.data.table(read.table("founders.hwe", header=TRUE))
setkey(DT.hwe.found,SNP)
setkey(DT.frq.found,SNP)
DT.all.found <- DT.hwe.found[TEST=="ALL"][DT.frq.found]
DT.all.found[MAF!=0][,plot(MAF,-log10(P))]

DT.all.found[,log10p:=-log10(P+1e-10)]
DT.all.found[,plot(MAF,log10p)]


plot(DT.all.found[,"log10p", with=FALSE])


DT.all.found[,log.expect.p:=-log10((1:nrow(DT.all.found)-0.5)/nrow(DT.all.found))]
DT.all.found[,log.sorted.p:=-log10(sort(P))]


DT.tdt.every.affect <- as.data.table(read.table("everyone_affected.tdt", header=TRUE))
setkey(DT.tdt.every.affect,SNP)
DT.tdt.switch <- as.data.table(read.table("switch_pheno.tdt", header=TRUE))
setkey(DT.tdt.switch,SNP)



DT.frq[,hist(MAF)]
DT.frq.found[,hist(MAF)]



DT.all[,plot(MAF,-log10(P))]
DT.all.found[,plot(MAF,-log10(P))]



DT.all[,plot(log.expect.p,log.sorted.p)]
abline(a=0,b=1)
DT.all.found[,plot(log.expect.p,log.sorted.p)]
abline(a=0,b=1)



manhattan(DT.tdt.every.affect)
manhattan(DT.tdt.switch)



qq(DT.tdt.every.affect$P, main="Q-Q plot of TDT p-values")
qq(DT.tdt.switch$P, main="Q-Q plot of TDT p-values")



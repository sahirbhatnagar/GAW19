##################################
# R source code file used to make plots
# Created by Sahir, June 16, 2014
# Updated June 16, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

setwd("~/git_repositories/GAW19/data/")

DT.frq <- as.data.table(read.table("chr21.frq", header=TRUE))
DT.frq[,hist(MAF)]

DT.hwe <- as.data.table(read.table("chr21.hwe", header=TRUE))
DT.hwe[1:5,1:9,with=FALSE]

setkey(DT.hwe,TEST)

#p-values
DT.hwe["ALL"][,"P",with=FALSE]

DT.all <- cbind(DT.frq,DT.hwe["ALL"] )

DT.all[,plot(MAF,-log10(P))]



DT.all[,log.expect.p:=-log10((1:nrow(DT.all)-0.5)/nrow(DT.all))]
DT.all[,sorted.p:=-log10(sort(P))]
DT.all[,plot(log.expect.p,sorted.p)]
abline(a=0,b=1)




# DT.all[,sorted.p:=NULL]
# DT.all[,expect.p:=NULL]
# DT.all[,log.expect.p:=NULL]

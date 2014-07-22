##################################
# R source code file for FBAT Manhattan plots
# Created by Sahir, June 7, 2014
# Updated July 7, 2014
# 
# NOTE: 
##################################



# FBAT all affected -------------------------------------------------------

setwd("~/share/sy/GAW19/FBAT/ALL")

DT.fbat <- fread("ALL.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_ALL_Manhattan.png");manhattan(DT.fbat, snp="Marker", main="FBAT_ALL");dev.off()


setwd("~/share/sy/GAW19/FBAT/SEQ")

DT.fbat <- fread("SEQ.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_SEQ_Manhattan.png");manhattan(DT.fbat, snp="Marker", main="FBAT_SEQ");dev.off()



setwd("~/share/sy/GAW19/FBAT/NUC")

DT.fbat <- fread("NUC.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_NUC_Manhattan.png");manhattan(DT.fbat, snp="Marker", main="FBAT_NUC");dev.off()






# FBAT for REAL PHENO -----------------------------------------------------

setwd("~/share/sy/GAW19/REAL/FBAT/ALL")

DT.fbat <- fread("ALL.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_ALL_REAL_PHENO_Manhattan.png");manhattan(DT.fbat, snp="Marker", main="FBAT_ALL_REAL_PHENO");dev.off()


setwd("~/share/sy/GAW19/REAL/FBAT/SEQ")

DT.fbat <- fread("SEQ.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_SEQ_REAL_PHENO_Manhattan.png");manhattan(DT.fbat, snp="Marker", main="FBAT_SEQ_REAL_PHENO");dev.off()



setwd("~/share/sy/GAW19/REAL/FBAT/NUC")

DT.fbat <- fread("NUC.FBAT", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_NUC_REAL_PHENO_Manhattan.png");manhattan(DT.fbat, snp="Marker", main="FBAT_NUC_REAL_PHENO");dev.off()




# FBAT all affected without -e option -------------------------------------------------------

setwd("~/share/sy/GAW19/FBAT/ALL")

DT.fbat <- fread("fbat21all.csv", sep=";")
head(DT.fbat)
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_ALL_we_Manhattan.png");manhattan(DT.fbat, snp="Marker", main="FBAT_ALL without e option");dev.off()


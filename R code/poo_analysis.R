##################################
# R source code file for Parent of Origin analysis
# Created by Sahir, June 10, 2014
# Updated July 10, 2014
# 
# NOTE: 
##################################

setwd("~/share/sy/GAW19/data/cleangeno/chr21")
DT <- fread("chr21poo.csv", sep=";")
set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["SNP"]])))

png("mat21_poo.png")
manhattan(DT,p="P_MAT", main="Maternal P-val, chr 21, seq");dev.off()

png("pat21_poo.png")
manhattan(DT,p="P_PAT", main="Paternal P-val, chr 21, seq");dev.off()

png("21_poo.png")
manhattan(DT,p="P_POO", main="POO P-val, chr 21, seq");dev.off()


#genome wide
setwd("/home/data1/share/sy/GAW19/POO")
DT <- fread("all.chr", sep=";")
set(DT, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT[["SNP"]])))

png("matall_poo.png")
manhattan(DT,p="P_MAT", main="Maternal P-val,TDT, seq");dev.off()

png("patall_poo.png")
manhattan(DT,p="P_PAT", main="Paternal P-val,TDT, seq");dev.off()

png("all_poo.png")
manhattan(DT,p="P_POO", main="POO P-val, TDT, seq");dev.off()

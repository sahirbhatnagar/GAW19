setwd("~/Documents/GAW19/")
library(data.table)
library(bit64)
library(plotrix)

filename <- "chr21-geno.csv"
DT <- fread(filename, nrows=10, select=c("snp","T2DG0200001"))
DT <- fread(filename, select=c("snp"))
DT

str(gene.exp)
colnames(gene.exp)

read.table
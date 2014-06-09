##################################
# R source code file used to create variables used for GAW analysis
# Created by Sahir, March 27, 2014
# Updated March 27, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

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
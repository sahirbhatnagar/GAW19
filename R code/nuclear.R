##################################
# R source code file used to create nuclear families
# Created by Sahir, July 13, 2014
# Updated 
# 
# NOTE: 
##################################

# Differences between genotype and family file 
# read in family file
DT.tfam <- fread("~/share/sy/GAW19/data/PED.csv")
fam.id <- DT.tfam$ID

#geno.id
DT.geno <- fread("~/share/sy/GAW19/data/rawgeno/chr21-geno.csv", nrows=2)
#get subject ID's for which we have genotype information
id.geno <- colnames(DT.geno)[-1]

#sequenced id
DT.seq <- fread("~/share/sy/GAW19/data/sequenced.id")
seq.id <- DT.seq$V2 

#nuclear families
DT.nuc <- fread("~/share/sy/GAW19/data/nuclear.fam")

#real phenotype
DT.real <- fread("/home/data1/share/sy/GAW19/data/real.pheno")

set(DT.tfam, i=NULL, j="hyper", value=DT.real[["hyper"]])

set(DT.tfam, i=NULL, j="genotyped", value=DT.tfam[["ID"]] %in% id.geno)
sum(DT.tfam$genotyped)
set(DT.tfam, i=NULL, j="sequenced", value=DT.tfam[["ID"]] %in% seq.id)
sum(DT.tfam$sequenced)

set(DT.tfam, i=NULL, j="nuclear", value=DT.tfam[["ID"]] %in% DT.nuc[["V2"]])

DT.tfam[nuclear==TRUE][genotyped==TRUE]

sum(DT.tfam[nuclear==TRUE]$sequenced)

table(DT.tfam[nuclear==TRUE]$PEDNUM)


DT.tfam[sequenced==TRUE][hyper==1]


#number of founders who are genotyped.. This is where the 117 comes from
sum(DT.tfam[FA==0]$ID %in% id.geno)

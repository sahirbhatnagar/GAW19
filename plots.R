##################################
# R source code file used to make plots
# Created by Sahir, June 16, 2014
# Updated June 16, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

setwd("~/git_repositories/GAW19/data/")
library(lattice)
chromosome <- 1
#rm(list=ls())
#system(paste("mkdir", DT.frq[["CHR"]][1], sep=" "))


# -log10(HWE pvalue) vs MAF -----------------------------------------------

DT.frq <- fread(paste("chr",chromosome,"frq.csv",sep=""))
DT.hwe <- fread(paste("chr",chromosome,"hwe.csv",sep=""), nrows=1000) 
head(DT.frq);head(DT.hwe)
setnames(DT.hwe, c("V1","V2","V9"), c("CHR", "SNP", "P"))
set(DT.hwe, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.hwe[["SNP"]])))
#DT.hwe <- DT.hwe[V3=="ALL"]
setkey(DT.hwe,SNP);setkey(DT.frq, SNP)
#setkey(DT.hwe2, V2)
str(DT.hwe)
DT.all <- DT.hwe[DT.frq]
DT.all[order(DT.all$P)[1:100]]
#DT.all2 <- DT.hwe2[DT.frq]
xyplot(-log10(V9) ~ MAF, DT.all, grid = TRUE, main=paste("chromosome ",chromosome, sep=" "))
xyplot(-log10(V9) ~ MAF, DT.all2, grid = TRUE, main=paste("chromosome ",chromosome," HWE standard", sep=" "))
#DT.all[P!=1][,plot(MAF,-log10(P))]
#DT.all[,plot(MAF,-log10(P))]
library(reshape2)

DT.hwe[, `:=` (aa=colsplit(DT.hwe[["V6"]],pattern="\\/", c("aa","Aa","aa"))[,1],
               Aa=colsplit(DT.hwe[["V6"]],pattern="\\/", c("aa","Aa","aa"))[,2],
               AA=colsplit(DT.hwe[["V6"]],pattern="\\/", c("aa","Aa","aa"))[,3])]

ll <- unlist(strsplit(DT.hwe$V6, "/", fixed=TRUE))
idx <- seq(1, length(ll), by=3)
DT.hwe[,`:=`(aa=ll[idx], Aa=ll[idx+1], AA=ll[idx+2])]

#DT.hwe2 <- as.data.table(read.table("chr21.hwe2.hwe", header=TRUE))
#DT.hwe2 <- DT.hwe2[TEST=="ALL"]
#setkey(DT.hwe,NULL)

# # Vince data to check my code
# DT <- fread("21.snp-stats")
# DT[HWE!=0][,plot(MAF, HWE)]
# 
# a=-log10(DT[HWE!=0][1:10, ,][["HWE"]])
# b=DT[HWE!=0][1:10, ,][["MAF"]]
# plot(b,a, main="plot ab")
xyplot(HWE ~ MAF, DT, grid = TRUE)


# Possible explanation for unexpected HWE vs MAF plot ---------------------
setwd("~/git_repositories/GAW19/data/")
DT.tfam <- fread("PED.csv")

#get subject ID's for which we have genotype information
filename = "chr21-geno.csv"
DT <- fread(filename, nrows=1) 
id.geno <- colnames(DT)[-1]
fam.id <- DT.tfam$ID

#ID's with missing genotypes
missing.geno <- setdiff(fam.id, id.geno)

# founder ID's 
founder.id <- DT.tfam[FA==0][["ID"]]

# 296 founders are missing genotype information
# so out of the 413 founders only 413-296=117 (30% have genotype information)
sum(!(founder.id %in% id.geno))
missing.geno.founder <- setdiff(founder.id, id.geno)

seqid <- fread("~/git_repositories/GAW19/data/sequenced.id", header=F)

#of the 413 founders, 91 have been actually sequenced
sum((founder.id %in% seqid$V2))

#id's who are founders and have been sequenced
id_seq_found <- intersect(founder.id,seqid$V2)

# id's who are not founders of the 959 individuals. so there should be 959-91=868
id_not_seq_not_found <- setdiff(id.geno,id_seq_found)

write.csv(id_not_found, file="not_seq_not_found.id", row.names=F, quote=F, col.names=F)



# Manhattan Plot ----------------------------------------------------------

DT.tdt.every.affect <- as.data.table(read.table("everyone_affected.tdt", header=TRUE))
setkey(DT.tdt.every.affect,SNP)
DT.tdt.switch <- as.data.table(read.table("switch_pheno.tdt", header=TRUE))
setkey(DT.tdt.switch,SNP)

manhattan(DT.tdt.every.affect)
manhattan(DT.tdt.switch)

qq(DT.tdt.every.affect$P, main="Q-Q plot of TDT p-values")
qq(DT.tdt.switch$P, main="Q-Q plot of TDT p-values")




# Q-Q plot ----------------------------------------------------------------

DT.all[,log.expect.p:=-log10((1:nrow(DT.all)-0.5)/nrow(DT.all))]
DT.all[,log.sorted.p:=-log10(sort(P))]

DT.all[,plot(log.expect.p,log.sorted.p)]
abline(a=0,b=1)
DT.all.found[,plot(log.expect.p,log.sorted.p)]
abline(a=0,b=1)




# Checking HWE calculations -----------------------------------------------
library(lattice)
setwd("~/git_repositories/GAW19/data/chr21")
chromosome <- 21;chr.digits=nchar(chromosome)

DT.frq <- fread(paste("chr",chromosome,"frq.csv",sep=""))
DT.hwe <- fread(paste("chr",chromosome,"hwe.csv",sep=""), header=FALSE) 
colnames(DT.hwe) <- c("CHR","SNP","TEST","A1","A2","aa","Aa","AA","O(HET)","E(HET)","pvalue")
head(DT.hwe);head(DT.frq)
setkey(DT.hwe,SNP);setkey(DT.frq, SNP)

DT.all <- DT.hwe[DT.frq]
head(DT.all)
xyplot(-log10(pvalue) ~ MAF, DT.all, grid = TRUE, main=paste("chromosome ",chromosome, sep=" "))
set(DT.all, i=NULL, j="my.maf", value=(2*DT.all[["aa"]]+DT.all[["Aa"]])/(2*(DT.all[["aa"]]+DT.all[["Aa"]]+DT.all[["AA"]])))
set(DT.all, i=NULL, j="NOBS", value=(DT.all[["aa"]]+DT.all[["Aa"]]+DT.all[["AA"]]))
set(DT.all, i=NULL, j="my.chisq", value= (DT.all[["aa"]] - DT.all[["NOBS"]]*(DT.all[["my.maf"]])^2)^2/(DT.all[["NOBS"]]*(DT.all[["my.maf"]])^2)+
(DT.all[["Aa"]] - 2*DT.all[["NOBS"]]*(DT.all[["my.maf"]])*(1-DT.all[["my.maf"]]))^2/(2*DT.all[["NOBS"]]*(DT.all[["my.maf"]])*(1-DT.all[["my.maf"]]))+
(DT.all[["AA"]] - DT.all[["NOBS"]]*(1-DT.all[["my.maf"]])^2)^2/(DT.all[["NOBS"]]*(1-DT.all[["my.maf"]])^2))
set(DT.all, i=NULL, j="P", value=pchisq(DT.all[["my.chisq"]], 1, lower.tail=F))
set(DT.all, i=NULL, j="BP", value = as.numeric(if (chr.digits==2) sub("[0-9][0-9]_","",DT.all[["SNP"]]) else sub("[0-9]_","",DT.all[["SNP"]])) )
str(DT.all)

#my pvalue
manhattan(DT.all)

#p-value from plink
manhattan(DT.all, p="pvalue")




p <- (2*81+10)/(2*(10+81))
q <-1-p

chisq <- (p^2*91-81)^2/(p^2*91)+(2*p*q*91-10)^2/(2*p*q*91)+(q^2*91-0)^2/(q^2*91)
pchisq(chisq,1, lower.tail=F)

head(DT.all[P<=0.1], n=10)


table(DT.all[["NOBS"]])

seqid <- fread("~/git_repositories/GAW19/data/seq.id", header=F)
DT.tped <- fread("chr21.tped", header=F)

hist(DT.all[["my.chisq"]])



# Manhattan HWE -----------------------------------------------------------

setwd("~/git_repositories/GAW19/data/")
DT.hwe.seq <- fread("allhwe.csv")
head(DT.hwe.seq)
set
nchar(DT.hwe.seq[["V1"]][1])
set(DT.hwe.seq[1:5, , ], i=NULL, j="BP", value = sub('\\d*_','',DT.hwe.seq[["V2"]]))

k <- DT.hwe.seq$V2[1:10]

grep(k,".{1,2}_", perl=TRUE)
gsub("[0-9]_.", "", k, perl=TRUE)
str_split(k[1], "_")[[1]][2]
library(stringr)


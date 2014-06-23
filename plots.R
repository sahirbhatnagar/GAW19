##################################
# R source code file used to make plots
# Created by Sahir, June 16, 2014
# Updated June 16, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

setwd("~/git_repositories/GAW19/data/")
library(lattice)
chromosome <- 7
#rm(list=ls())
#system(paste("mkdir", DT.frq[["CHR"]][1], sep=" "))


# -log10(HWE pvalue) vs MAF -----------------------------------------------

DT.frq <- as.data.table(read.table(paste("chr",chromosome,".frq",sep=""), header=TRUE))
DT.hwe <- as.data.table(read.table(paste("chr",chromosome,".hwe",sep=""), header=TRUE))[TEST=="ALL"]
DT.hwe2 <- as.data.table(read.table(paste("chr",chromosome,".hwe2.hwe",sep=""), header=TRUE))[TEST=="ALL"]

setkey(DT.hwe,SNP)
setkey(DT.frq, SNP)
setkey(DT.hwe2, SNP)

DT.all <- DT.hwe[DT.frq]
DT.all2 <- DT.hwe2[DT.frq]
xyplot(-log10(P) ~ MAF, DT.all, grid = TRUE, main=paste("chromosome ",chromosome, sep=" "))
xyplot(-log10(P) ~ MAF, DT.all2[MAF!=0], grid = TRUE, main=paste("chromosome ",chromosome," HWE standard", sep=" "))
#DT.all[P!=1][,plot(MAF,-log10(P))]
#DT.all[,plot(MAF,-log10(P))]



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







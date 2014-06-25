##################################
# R source code file used to analyse TDT results
# Created by Sahir, June 17, 2014
# Updated June 17, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

setwd("~/git_repositories/GAW19/data")
chromosome <- 1
# TDT results -------------------------------------------------------------
DT.tdt <- fread(paste("chr",chromosome,"tdt.csv",sep=""))
setkey(DT.tdt,SNP)

#check their chi-square statistic
set(DT.tdt, i=NULL, j="my.chi", value=((DT.tdt[["T"]]-DT.tdt[["U"]])^2)/(DT.tdt[["T"]]+DT.tdt[["U"]]))
set(DT.tdt, i=NULL, j="my.p", value=pchisq(DT.tdt[["my.chi"]],df=1,lower.tail=FALSE))
DT.tdt[order(DT.tdt$P)[1:5]]


str(DT.tdt)

manhattan(DT.tdt)

DT.tdt[which.min(DT.tdt$P), ,][,"SNP",with=FALSE]

as.character(DT.tdt[order(DT.tdt$P)[1:5]][["SNP"]])

manhattan(DT.tdt, highlight=as.character(DT.tdt[order(DT.tdt$P)[1:5]][["SNP"]]))

qq(as.numeric(DT.tdt[!is.na(P)][["P"]]), main="Q-Q plot of TDT p-values")
str(DT.tdt)
which.min
max



# Use kinship package to break pedigrees ----------------------------------
setwd("~/git_repositories/GAW19/data")
DT.tfam <- fread("PED.csv")

library(kinship2)

set(DT.tfam, i=NULL, j="id", value=as.integer(substr(DT.tfam[["ID"]],5,11)))
set(DT.tfam, i=NULL, j="dadid", value=as.integer(substr(DT.tfam[["FA"]],5,11)))
set(DT.tfam, i=NULL, j="momid", value=as.integer(substr(DT.tfam[["MO"]],5,11)))
#replace NA with 0
f_dowle3 = function(DT) {
  # either of the following for loops
  
  # by name :
  #for (j in names(DT))
  # set(DT,which(is.na(DT[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,value=0)
}

f_dowle3(DT.tfam)
str(DT.tfam)

ped <- with(DT.tfam,pedigree(id,dadid,momid,SEX,affected=rep(1,nrow(DT.tfam)),famid=PEDNUM))
plot(ped['3'])

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


#create geno availability column
set(DT.tfam, i=NULL,j="avail", value=DT.tfam[["ID"]] %in% id.geno)

shrink1.B30 <- pedigree.shrink(ped=ped['2'], avail=DT.tfam[PEDNUM==2]$avail, maxBits=10)
as.data.frame(shrink1.B30$pedObj)

plot(shrink1.B30$pedObj)

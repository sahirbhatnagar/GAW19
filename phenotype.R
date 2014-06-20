##################################
# R source code file used to make the phenotype file for analysis
# Created by Sahir, June 16, 2014
# Updated June 16, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################


# Phenotype information ---------------------------------------------------
setwd("~/git_repositories/GAW19/data")
#rm(DT.phen)
DT.phen <- fread("PHEN.csv")
phen.names <- colnames(DT.phen)
#extract family ID from individual ID
set(DT.phen, i=NULL, j="FID", value = gsub("(?<![0-9])0+","",substr(DT.phen[["t2dgID"]],5,6),perl=TRUE))
setcolorder(DT.phen, c("FID",phen.names[-30]))
setnames(DT.phen,"t2dgID","IID")

f_dowle3 = function(DT) {
  # either of the following for loops
  
  # by name :
  #for (j in names(DT))
  # set(DT,which(is.na(DT[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,value=-9)
}

#replace NA's with -9 which is default missing for plink phenotype files
f_dowle3(DT.phen)

# add affected status for everyone
set(DT.phen, i=NULL, j="EVERYONE_AFFECTED", value=1)
#set(DT.phen, which(DT.phen[["HTN_1"]]==-9), "EVERYONE_AFFECTED", value=-9)
#set(DT.phen, which(DT.phen[["HTN_1"]]!=-9), "EVERYONE_AFFECTED", value=1)

# add phenotype where Cases --> Missing and Controls --> Cases
# table(DT.phen$HTN_1);table(DT.phen$HTN_2)
# table(DT.phen$HTN_3);table(DT.phen$HTN_4)
set(DT.phen, which(DT.phen[["HTN_1"]]==1),"SWITCH_PHENO",value=-9)
set(DT.phen, which(DT.phen[["HTN_1"]]==0),"SWITCH_PHENO",value=1)
set(DT.phen, which(DT.phen[["HTN_1"]]==-9),"SWITCH_PHENO",value=-9)

#write phenotype file
write.table(DT.phen, file="pheno.txt", col.names=TRUE,
            row.names=FALSE, quote=FALSE)






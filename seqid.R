##################################
# R source code file used to create ID's who were actually sequenced
# Created by Sahir, June 24, 2014
# Updated June 24, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
# I used pseq (shell script program) to extract ID of sequenced people from vcf files 
##################################

setwd("~/git_repositories/GAW19/data")
DT <- fread("sequenced.id", header=F)
str(DT)
set(DT, i=NULL, j="FID", value=as.numeric(substr(DT[["V1"]],5,6)))
setcolorder(DT, c("FID","V1"))

write.table(DT, file="sequenced.id",
            col.names=FALSE, row.names=FALSE, quote=FALSE)

##################################
# R source code file for creating hypertension phenotype written to real.pheno
# Created by Sahir, June 11, 2014
# Updated July 12, 2014
# 
# NOTE: 
##################################
`%nin%` <- Negate(`%in%`) 

setwd("~/share/sy/GAW19/data")

rm(DT)
DT <- fread("PHEN.csv")
DT <- DT[, c("t2dgID","HTN_1","HTN_2","HTN_3","HTN_4"), with=F]
set(DT, i=NULL, j="IID", value=sub("T2DG", "", DT[["t2dgID"]]))
set(DT, i=NULL, j="tempFID", value=substr(DT[["IID"]],1,2))
set(DT, i=which(DT$tempFID %in% c("02","03","04","05","06","07","08","09")), j="FID", value=sub("0","", substr(DT[tempFID %in% c("02","03","04","05","06","07","08","09")][["IID"]],1,2)))
set(DT, i=which(DT$tempFID %nin% c("02","03","04","05","06","07","08","09")), j="FID", value=DT[tempFID %nin% c("02","03","04","05","06","07","08","09")][["tempFID"]])

set(DT, i=NULL, j="na1", value = is.na(DT[["HTN_1"]]));set(DT, i=NULL, j="na2", value = is.na(DT[["HTN_2"]]))
set(DT, i=NULL, j="na3", value = is.na(DT[["HTN_3"]]));set(DT, i=NULL, j="na4", value = is.na(DT[["HTN_4"]]))
set(DT, i=NULL, j="sums", value=rowSums(DT[,c("na1","na2","na3","na4"),with=F]))
set(DT, i=NULL, j="hyper", value = pmin(1,rowSums(DT[,c("HTN_1","HTN_2","HTN_3","HTN_4"),with=F], na.rm=T)) )
set(DT, i=which(DT$sums==4), j="hyper", value = -9 )

# there are 455 missing phenotypes
# 564 controls, 370 cases
table(DT$hyper)


write.table(DT[,c("FID","IID","hyper"), with=F], file="/home/data1/share/sy/GAW19/data/real.pheno", 
            row.names=FALSE, quote=FALSE)






j<- read.table("real.pheno")
j


d <- fread("~/share/sy/GAW19/data/cleangeno/realpheno/chr21/FBAT21.ped", select=1:6)
table(d$V6)
setkey(d,V6)
setkey(DT, IID)



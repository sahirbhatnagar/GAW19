##################################
# R source code file for HWE vs p-values of test statistics
# Created by Sahir, July 14
# Updated July 14, 2014
# 
# NOTE: 
##################################



# PDT SEQ  ----------------------------------------------------------------
rm(DT.hwe)
setwd("~/share/sy/GAW19/HWE/SEQ_ONLY")
DT.hwe <- fread("allhwe.csv", sep=";") 
head(DT.hwe)
setnames(DT.hwe, c("V1","V2","V9"), c("CHR", "SNP", "HWE"))
#set(DT.hwe, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.hwe[["SNP"]])))
setkey(DT.hwe,SNP)

#hypertension pdt seq
setwd("~/share/sy/GAW19/REAL/PDT/SEQ")
DT.hyper <- fread("SEQ.PDT",sep=";")
head(DT.hyper)
setnames(DT.hyper, c("LOCUS_NAME","P-VALUE"), c("SNP","P_HYPERTENSION"))
setkey(DT.hyper,SNP)

#all pdt seq
setwd("~/share/sy/GAW19/PDT/SEQ")
DT.all <- fread("SEQ.PDT",sep=";")
head(DT.all)
setnames(DT.all,c("V1","V5"), c("SNP","P_ALL"))
setkey(DT.all, SNP)

DT <- DT.hyper[,c("SNP","P_HYPERTENSION"),with=F][ DT.all[,c("SNP","P_ALL"),with=F][P_ALL!="bad_data"]][DT.hwe[,c("SNP","HWE"),with=F]]
head(DT)

summary(DT[,lm(P_ALL~HWE)])
summary(DT[,lm(P_HYPERTENSION~HWE)])

summary(DT[!is.na(P_ALL)][P_ALL<=0.001][,lm(-log10(as.numeric(P_ALL))~I(-log10(as.numeric(HWE))))])
summary(DT[!is.na(P_HYPERTENSION)][P_HYPERTENSION<=0.001][,lm(P_HYPERTENSION~HWE)])

DT[!is.na(P_ALL)][P_ALL<=0.001][,plot(HWE,-log10(as.numeric(P_ALL)))]
DT[!is.na(P_ALL)][P_HYPERTENSION<=0.001][,plot(HWE,-log10(as.numeric(P_HYPERTENSION)))]




# PDT ALL -----------------------------------------------------------------

setwd("~/share/sy/GAW19/HWE/")
DT.hwe <- fread("all.hwe") 
head(DT.hwe)
setnames(DT.hwe, c("V1","V2","V9"), c("CHR", "SNP", "HWE"))
#set(DT.hwe, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.hwe[["SNP"]])))
setkey(DT.hwe,SNP)


#hypertension pdt all
setwd("~/share/sy/GAW19/REAL/PDT/ALL")
DT.hyper <- fread("ALL.PDT",sep=";")
head(DT.hyper)
setnames(DT.hyper,c("V1","V5"), c("SNP","P_HYPERTENSION"))
setkey(DT.hyper,SNP)

#everyone affected pdt all
setwd("~/share/sy/GAW19/PDT/ALL")
DT.all <- fread("ALL.PDT",sep=";")
head(DT.all)
setnames(DT.all,c("V1","V5"), c("SNP","P_ALL"))
setkey(DT.all,SNP)

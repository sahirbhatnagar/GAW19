##################################
# R source code file for all affected vs. hypertension plots
# Created by Sahir, July 13
# Updated July 13, 2014
# 
# NOTE: 
##################################


# PDT SEQ, All affected vs Hypertension -----------------------------------

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

DT <- DT.hyper[,c("SNP","P_HYPERTENSION"),with=F][ DT.all[,c("SNP","P_ALL"),with=F][P_ALL!="bad_data"]]
head(DT)

#none of the markers that both have small p-values
DT[order(DT$P_HYPERTENSION)[1:50]]



setwd("~/share/sy/GAW19/REAL_VS_ALL")
png("PDT_SEQ_REAL_VS_ALL.png");
DT[,plot(P_HYPERTENSION,P_ALL, main="PDT SEQ, All affected vs Hypertension")];dev.off()

# PDT All, All affected vs Hypertension -----------------------------------

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


DT <- DT.hyper[,c("SNP","P_HYPERTENSION"),with=F][P_HYPERTENSION!="bad_data"][ DT.all[,c("SNP","P_ALL"),with=F][P_ALL!="bad_data"]]
head(DT)

DT[order(DT$P_HYPERTENSION)[1:50]]



setwd("~/share/sy/GAW19/REAL_VS_ALL")
png("PDT_ALL_REAL_VS_ALL.png");
DT[,plot(P_HYPERTENSION,P_ALL, main="PDT ALL, All affected vs Hypertension")];dev.off()



# PDT NUC, All affected vs Hypertension -----------------------------------

#hypertension pdt all
setwd("~/share/sy/GAW19/REAL/PDT/NUC")
DT.hyper <- fread("NUC.PDT",sep=";")
head(DT.hyper)
setnames(DT.hyper,c("V1","V5"), c("SNP","P_HYPERTENSION"))
setkey(DT.hyper,SNP)

#everyone affected pdt all
setwd("~/share/sy/GAW19/PDT/NUC")
DT.all <- fread("NUC.PDT",sep=";")
head(DT.all)
setnames(DT.all,c("V1","V5"), c("SNP","P_ALL"))
setkey(DT.all,SNP)

DT[order(DT$P_HYPERTENSION)[1:50]]


DT <- DT.hyper[,c("SNP","P_HYPERTENSION"),with=F][P_HYPERTENSION!="bad_data"][ DT.all[,c("SNP","P_ALL"),with=F][P_ALL!="bad_data"]]
head(DT)
setwd("~/share/sy/GAW19/REAL_VS_ALL")
png("PDT_NUC_REAL_VS_ALL.png");
DT[,plot(P_HYPERTENSION,P_ALL, main="PDT NUC, All affected vs Hypertension")];dev.off()





# TDT SEQ, All affected vs Hypertension ---------------------------------------------------------------------

#hypertension tdt seq
setwd("~/share/sy/GAW19/REAL/TDT/SEQ")
DT.hyper <- fread("SEQ.TDT",sep=";")
head(DT.hyper)
setnames(DT.hyper, c("P"), c("P_HYPERTENSION"))
setkey(DT.hyper,SNP)

#all tdt seq
setwd("~/share/sy/GAW19/TDT/SEQ")
DT.all <- fread("SEQ.TDT",sep=";")
head(DT.all)
DT.all[order(DT.all$P)[1:10]]

setnames(DT.all,c("P"), c("P_ALL"))

DT <- DT.hyper[,c("SNP","P_HYPERTENSION"),with=F][ DT.all[,c("SNP","P_ALL"),with=F]]
head(DT)

DT[order(DT$P_HYPERTENSION)[1:50]]




setwd("~/share/sy/GAW19/REAL_VS_ALL")
png("TDT_SEQ_REAL_VS_ALL.png");
DT[!is.na(P_HYPERTENSION)][,plot(P_HYPERTENSION,P_ALL, main="TDT SEQ, All affected vs Hypertension")];dev.off()



# TDT ALL, All affected vs Hypertension ---------------------------------------------------------------------

#hypertension tdt all
setwd("~/share/sy/GAW19/REAL/TDT/ALL")
DT.hyper <- fread("ALL.TDT",sep=";")
head(DT.hyper)
setnames(DT.hyper, c("P"), c("P_HYPERTENSION"))
setkey(DT.hyper,SNP)

#all tdt all
setwd("~/share/sy/GAW19/TDT/ALL")
DT.all <- fread("ALL.TDT",sep=";")
head(DT.all)
setnames(DT.all,c("P"), c("P_ALL"))

DT.all[order(DT.all$P_ALL)[1:10]]
DT.all[P_ALL<=1e-10]


DT <- DT.hyper[,c("SNP","P_HYPERTENSION"),with=F][ DT.all[,c("SNP","P_ALL"),with=F]]
head(DT)

setwd("~/share/sy/GAW19/REAL_VS_ALL")
png("TDT_ALL_REAL_VS_ALL.png");
DT[!is.na(P_HYPERTENSION)][,plot(P_HYPERTENSION,P_ALL, main="TDT ALL, All affected vs Hypertension")];dev.off()

DT[order(DT$P_HYPERTENSION)[1:10]]


# FBAT SEQ, All affected vs Hypertension ---------------------------------------------------------------------

#hypertension fbat seq
setwd("~/share/sy/GAW19/REAL/FBAT/SEQ")
DT.hyper <- fread("SEQ.FBAT",sep=";")
head(DT.hyper)
setnames(DT.hyper, c("Marker","P"), c("SNP","P_HYPERTENSION"))
setkey(DT.hyper,SNP)
DT.hyper <- unique(DT.hyper)

#all fbat seq
setwd("~/share/sy/GAW19/FBAT")
DT.all <- fread("all.temp",sep=";")
head(DT.all)
setnames(DT.all,c("Marker","P"), c("SNP","P_ALL"))
setkey(DT.all,SNP)
DT.all <- unique(DT.all)


DT <- DT.all[,c("SNP","P_ALL"),with=F][ DT.hyper[,c("SNP","P_HYPERTENSION"),with=F]]
head(DT)

setwd("~/share/sy/GAW19/REAL_VS_ALL")
png("FBAT_SEQ_REAL_VS_ALL.png");
DT[,plot(P_HYPERTENSION,P_ALL, main="FBAT SEQ, All affected vs Hypertension")];dev.off()

DT[!is.na(P_ALL)][order(DT$P_HYPERTENSION)[1:20]]

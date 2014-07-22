##################################
# R source code file for regression analysis
# Created by Sahir, June 7, 2014
# Updated July 14, 2014
# 
# NOTE: The analysis here is for Everyone Affected
# we are performing 18 different analyses here:
# SEQUENCED/ALL/NUC then TDT/PDT/FBAT then All_SNPS/SNPS<0.001
##################################


# VCF file for covariates -------------------------------------------------
# Information fields provided in the vcf files included 
# NS: number of samples with fully called data 
# AF: allele frequency 
# DB: dbSNP membership
# RSID: dbSNP rs identifier
# STR: strand bias Pearson's correlation
# STZ: strand bias z-score 
# CBR: cycle bias Pearson's correlation
# CBZ: cycle bias z-score 
# CSR: cycle-strand Pearson's correlation
# IOZ: base-quality inflation z-score 
# IOR: ratio of base-quality inflation 
# AOZ: alternate allele quality z-score 
# AOI: alternate allele inflation score
# MQ0: fraction of bases with map quality of 0, 
# MQ10: less than 10, MQ20: less than 20, and MQ30: less than 30.

setwd("~/share/sy/GAW19/VCF")
#chromosome <- 21
rm(DT.cov)
#DT.cov <- fread(paste("covar",chromosome,".csv",sep=""), sep=";")
DT.cov <- fread(paste("all.covar",sep=""), sep=";")
head(DT.cov)

#non missing values
# nrow(DT.cov[NS!="."]);nrow(DT.cov[AF!="."]);nrow(DT.cov[STR!="."]);nrow(DT.cov[STZ!="."])
# nrow(DT.cov[CBR!="."]);nrow(DT.cov[CBZ!="."]);nrow(DT.cov[CSR!="."]);nrow(DT.cov[IOZ!="."])
# nrow(DT.cov[IOR!="."]);nrow(DT.cov[AOZ!="."]);nrow(DT.cov[AOI!="."]);nrow(DT.cov[MQ0!="."])
# nrow(DT.cov[MQ10!="."]);nrow(DT.cov[MQ20!="."]);nrow(DT.cov[MQ30!="."]);

set(DT.cov, i=NULL, j="SNP", value=paste(DT.cov[["CHR"]],"_",DT.cov[["BP"]], sep=""))
setkey(DT.cov, SNP)
DT.cov <- DT.cov[,c("NS","AF","SNP","STR","STZ","CBR","CBZ","IOZ","IOR","AOZ","AOI","MQ0","MQ10","MQ20","MQ30"),with=F]
#str(DT.cov)
DT.cov[,"STR":=as.numeric(STR)];DT.cov[,"STZ":=as.numeric(STZ)];DT.cov[,"CBR":=as.numeric(CBR)]
DT.cov[,"CBZ":=as.numeric(CBZ)];DT.cov[,"IOZ":=as.numeric(IOZ)];DT.cov[,"IOR":=as.numeric(IOR)]
DT.cov[,"AOZ":=as.numeric(AOZ)];DT.cov[,"AOI":=as.numeric(AOI)];DT.cov[,"MQ0":=as.numeric(MQ0)]
DT.cov[,"MQ10":=as.numeric(MQ10)];DT.cov[,"MQ20":=as.numeric(MQ20)]
DT.cov <- DT.cov[complete.cases(DT.cov),]


#########################################################
#                   SEQUENCED-------------------
#########################################################

# TDT ---------------------------------------------------
#########################################################

setwd("~/share/sy/GAW19/TDT/SEQ")
#DT.tdt <- fread(paste("chr",chromosome,"seq.csv",sep=""), sep=";")
DT.tdt <- fread("SEQ.TDT", sep=";")
head(DT.tdt,20)
setnames(DT.tdt, "P", "P_TDT")
setkey(DT.tdt, SNP)
str(DT.tdt)
DT.tdt <- DT.tdt[!is.na(P_TDT)]
#plot(density(-log10(DT.tdt[["P_TDT"]]), na.rm=TRUE))


# PDT ---------------------------------------------------
#########################################################
rm(DT.pdt)
setwd("~/share/sy/GAW19/PDT/SEQ")
#DT.pdt <- fread(paste("pdt",chromosome,"seq.csv",sep=""), sep=";")
DT.pdt <- fread(paste("SEQ.PDT",sep=""), sep=";")
head(DT.pdt)
setnames(DT.pdt, c("V1","V5"), c("SNP","P"))
DT.pdt <- DT.pdt[P!="bad_data"]
set(DT.pdt, i=NULL, j="P_PDT", value= as.numeric(DT.pdt[["P"]])) 
setkey(DT.pdt, SNP)
#DT.pdt[P_PDT==1]

#DT.pdt[,c("SNP","P_PDT"),with=F]

#plot(density(-log10(DT.pdt[["P_PDT"]])))

#DT.fbat["21_9923775"]

# k <- substr(DT.pdt[["V18"]],1,5)
# unique(k)
# 
# DT.pdt[DT.fbat]
# DT.frq[DT.pdt[P_PDT==1]][,plot("MAF","P_PDT")]

# FBAT ---------------------------------------------------
#########################################################

setwd("/home/data1/share/sy/GAW19/FBAT/SEQ")
DT.fbat <- fread(paste("SEQ.FBAT"), sep=";")
#set(DT.fbat, i=NULL, j="CHR", value=as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])))
#set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
head(DT.fbat)
setnames(DT.fbat, c("Marker","P"), c("SNP","P_FBAT"))
setkey(DT.fbat, SNP)
# these figures show that A1 and A2 Pvalues are the same
# png("FBAT_seq_A1.png");manhattan(DT.fbat[Allele==1], snp="Marker", main="FBAT,A1, sequenced only");dev.off()
# png("FBAT_seq_A2.png");manhattan(DT.fbat[Allele==2], snp="Marker", main="FBAT,A2, sequenced only");dev.off()

#DT.fbat <- DT.fbat[CHR==chromosome]
DT.fbat <- unique(DT.fbat)
#DT.fbat[,c("SNP","P_FBAT"),with=F]
head(DT.fbat, 10)

# DT.fbat[P_FBAT==1]
# 
# 
# DT.fbat[,plot(density(P_FBAT))]
# 
# rm(DT.fbat)
# 
# ls()

# MAF (SEQ ONLY) ---------------------------------------
#########################################################
setwd("/home/data1/share/sy/GAW19/MAF/SEQ")
DT.frq <- fread("SEQ.MAF", sep=";")
head(DT.frq)
#DT.frq <- DT.frq[CHR==chromosome]
setkey(DT.frq, SNP)
#DT.frq[,c("SNP","MAF"), with=F]
#table(DT.frqseq$NCHROBS)


# HWE (SEQ only) ---------------------------------------
#########################################################

setwd("/home/data1/share/sy/GAW19/HWE/SEQ")
DT.hwe <- fread("SEQ.HWE", sep=";")
#DT.hwe <- DT.hwe[V1==chromosome]
DT.hwe <- DT.hwe[V3=="ALL"]
head(DT.hwe)
setnames(DT.hwe, c("V2","V9"), c("SNP","HWE"))
setkey(DT.hwe, SNP)
#DT.hwe[,c("SNP","HWE"), with=F]


# TDT Regression Sequenced ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.cov[,c("NS","SNP","STR","STZ","CBR","CBZ","IOZ","IOR","AOZ","AOI","MQ10","MQ20","MQ30"),with=F][DT.frq[,c("SNP","MAF"), with=F]][DT.hwe[,c("SNP","HWE"), with=F]][DT.tdt[,c("SNP","P_TDT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_TDT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_TDT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#TDT sequenced all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"
gg2 <- "log10p ~ log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#TDT sequenced SNPs with TDT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"
gg4 <- "log10p ~ log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"

fit3 <- DT[P_TDT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_TDT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_TDT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]



# PDT Regression Sequenced ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.cov[,c("NS","SNP","STR","STZ","CBR","CBZ","IOZ","IOR","AOZ","AOI","MQ10","MQ20","MQ30"),with=F][DT.frq[,c("SNP","MAF"), with=F]][DT.hwe[,c("SNP","HWE"), with=F]][DT.pdt[,c("SNP","P_PDT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_PDT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_PDT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#PDT sequenced all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"
gg2 <- "log10p ~ log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#PDT sequenced SNPs with PDT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"
gg4 <- "log10p ~ log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"

fit3 <- DT[P_PDT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_PDT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_PDT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]






# FBAT Regression Sequenced ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.cov[,c("NS","SNP","STR","STZ","CBR","CBZ","IOZ","IOR","AOZ","AOI","MQ10","MQ20","MQ30"),with=F][DT.frq[,c("SNP","MAF"), with=F]][DT.hwe[,c("SNP","HWE"), with=F]][DT.fbat[,c("SNP","P_FBAT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_FBAT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_FBAT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#FBAT sequenced all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"
gg2 <- "log10p ~ log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#FBAT sequenced SNPs with FBAT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"
gg4 <- "log10p ~ log10MAF+NS+STR+STZ+CBR+CBZ+IOZ+IOR+AOZ+AOI+MQ10+MQ20+MQ30"

fit3 <- DT[P_FBAT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_FBAT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_FBAT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]




#########################################################
#                   ALL-------------------
#########################################################

# TDT ---------------------------------------------------
#########################################################

setwd("~/share/sy/GAW19/TDT/ALL")
#DT.tdt <- fread(paste("chr",chromosome,"seq.csv",sep=""), sep=";")
DT.tdt <- fread("ALL.TDT", sep=";")
head(DT.tdt,20)
setnames(DT.tdt, "P", "P_TDT")
setkey(DT.tdt, SNP)
str(DT.tdt)
DT.tdt <- DT.tdt[!is.na(P_TDT)]
#plot(density(-log10(DT.tdt[["P_TDT"]]), na.rm=TRUE))


# PDT ---------------------------------------------------
#########################################################
rm(DT.pdt)
setwd("~/share/sy/GAW19/PDT/ALL")
#DT.pdt <- fread(paste("pdt",chromosome,"seq.csv",sep=""), sep=";")
DT.pdt <- fread(paste("ALL.PDT",sep=""), sep=";")
head(DT.pdt)
setnames(DT.pdt, c("V1","V5"), c("SNP","P"))
DT.pdt <- DT.pdt[P!="bad_data"]
set(DT.pdt, i=NULL, j="P_PDT", value= as.numeric(DT.pdt[["P"]])) 
setkey(DT.pdt, SNP)
#DT.pdt[P_PDT==1]

#DT.pdt[,c("SNP","P_PDT"),with=F]

#plot(density(-log10(DT.pdt[["P_PDT"]])))

#DT.fbat["21_9923775"]

# k <- substr(DT.pdt[["V18"]],1,5)
# unique(k)
# 
# DT.pdt[DT.fbat]
# DT.frq[DT.pdt[P_PDT==1]][,plot("MAF","P_PDT")]

# FBAT ---------------------------------------------------
#########################################################

setwd("/home/data1/share/sy/GAW19/FBAT/ALL")
DT.fbat <- fread(paste("ALL.FBAT"), sep=";")
#set(DT.fbat, i=NULL, j="CHR", value=as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])))
#set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
head(DT.fbat)
setnames(DT.fbat, c("Marker","P"), c("SNP","P_FBAT"))
setkey(DT.fbat, SNP)
# these figures show that A1 and A2 Pvalues are the same
# png("FBAT_seq_A1.png");manhattan(DT.fbat[Allele==1], snp="Marker", main="FBAT,A1, sequenced only");dev.off()
# png("FBAT_seq_A2.png");manhattan(DT.fbat[Allele==2], snp="Marker", main="FBAT,A2, sequenced only");dev.off()

#DT.fbat <- DT.fbat[CHR==chromosome]
DT.fbat <- unique(DT.fbat)
#DT.fbat[,c("SNP","P_FBAT"),with=F]
head(DT.fbat, 10)

# DT.fbat[P_FBAT==1]
# 
# 
# DT.fbat[,plot(density(P_FBAT))]
# 
# rm(DT.fbat)
# 
# ls()

# MAF (ALL) ---------------------------------------
#########################################################
setwd("/home/data1/share/sy/GAW19/MAF/ALL")
DT.frq <- fread("ALL.MAF", sep=";")
head(DT.frq)
#DT.frq <- DT.frq[CHR==chromosome]
setkey(DT.frq, SNP)
#DT.frq[,c("SNP","MAF"), with=F]

table(DT.frq$NCHROBS)


# HWE (ALL) ---------------------------------------
#########################################################

setwd("/home/data1/share/sy/GAW19/HWE/ALL")
DT.hwe <- fread("ALL.HWE", sep=";")
#DT.hwe <- DT.hwe[V1==chromosome]
DT.hwe <- DT.hwe[V3=="ALL"]
head(DT.hwe)
setnames(DT.hwe, c("V2","V9"), c("SNP","HWE"))
setkey(DT.hwe, SNP)
#DT.hwe[,c("SNP","HWE"), with=F]





# TDT Regression all ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.frq[,c("SNP","MAF"), with=F][DT.hwe[,c("SNP","HWE"), with=F]][DT.tdt[,c("SNP","P_TDT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_TDT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_TDT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#TDT all all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF"
gg2 <- "log10p ~ log10MAF"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#TDT all SNPs with TDT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF"
gg4 <- "log10p ~ log10MAF"

fit3 <- DT[P_TDT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_TDT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_TDT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]


# PDT Regression All ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.frq[,c("SNP","MAF"), with=F][DT.hwe[,c("SNP","HWE"), with=F]][DT.pdt[,c("SNP","P_PDT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_PDT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_PDT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#PDT all, all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF"
gg2 <- "log10p ~ log10MAF"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#PDT all SNPs with PDT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF"
gg4 <- "log10p ~ log10MAF"

fit3 <- DT[P_PDT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_PDT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_PDT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]


# FBAT Regression all ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.frq[,c("SNP","MAF"), with=F][DT.hwe[,c("SNP","HWE"), with=F]][DT.fbat[,c("SNP","P_FBAT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_FBAT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_FBAT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#FBAT all all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF"
gg2 <- "log10p ~ log10MAF"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#FBAT all SNPs with FBAT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF"
gg4 <- "log10p ~ log10MAF"

fit3 <- DT[P_FBAT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_FBAT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_FBAT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]



#########################################################
#                   NUCLEAR-------------------
#########################################################

# TDT ---------------------------------------------------
#########################################################

setwd("~/share/sy/GAW19/TDT/NUC")
#DT.tdt <- fread(paste("chr",chromosome,"seq.csv",sep=""), sep=";")
DT.tdt <- fread("NUC.TDT", sep=";")
head(DT.tdt,20)
setnames(DT.tdt, "P", "P_TDT")
setkey(DT.tdt, SNP)
str(DT.tdt)
DT.tdt <- DT.tdt[!is.na(P_TDT)]
#plot(density(-log10(DT.tdt[["P_TDT"]]), na.rm=TRUE))


# PDT ---------------------------------------------------
#########################################################
rm(DT.pdt)
setwd("~/share/sy/GAW19/PDT/NUC")
#DT.pdt <- fread(paste("pdt",chromosome,"seq.csv",sep=""), sep=";")
DT.pdt <- fread(paste("NUC.PDT",sep=""), sep=";")
head(DT.pdt)
setnames(DT.pdt, c("V1","V5"), c("SNP","P"))
DT.pdt <- DT.pdt[P!="bad_data"]
set(DT.pdt, i=NULL, j="P_PDT", value= as.numeric(DT.pdt[["P"]])) 
setkey(DT.pdt, SNP)
#DT.pdt[P_PDT==1]

#DT.pdt[,c("SNP","P_PDT"),with=F]

#plot(density(-log10(DT.pdt[["P_PDT"]])))

#DT.fbat["21_9923775"]

# k <- substr(DT.pdt[["V18"]],1,5)
# unique(k)
# 
# DT.pdt[DT.fbat]
# DT.frq[DT.pdt[P_PDT==1]][,plot("MAF","P_PDT")]

# FBAT ---------------------------------------------------
#########################################################

setwd("/home/data1/share/sy/GAW19/FBAT/NUC")
DT.fbat <- fread(paste("NUC.FBAT"), sep=";")
#set(DT.fbat, i=NULL, j="CHR", value=as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])))
#set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
head(DT.fbat)
setnames(DT.fbat, c("Marker","P"), c("SNP","P_FBAT"))
setkey(DT.fbat, SNP)
# these figures show that A1 and A2 Pvalues are the same
# png("FBAT_seq_A1.png");manhattan(DT.fbat[Allele==1], snp="Marker", main="FBAT,A1, sequenced only");dev.off()
# png("FBAT_seq_A2.png");manhattan(DT.fbat[Allele==2], snp="Marker", main="FBAT,A2, sequenced only");dev.off()

#DT.fbat <- DT.fbat[CHR==chromosome]
DT.fbat <- unique(DT.fbat)
#DT.fbat[,c("SNP","P_FBAT"),with=F]
head(DT.fbat, 10)

# DT.fbat[P_FBAT==1]
# 
# 
# DT.fbat[,plot(density(P_FBAT))]
# 
# rm(DT.fbat)
# 
# ls()

# MAF (NUC) ---------------------------------------
#########################################################
setwd("/home/data1/share/sy/GAW19/MAF/NUC")
DT.frq <- fread("NUC.MAF", sep=";")
head(DT.frq)
#DT.frq <- DT.frq[CHR==chromosome]
setkey(DT.frq, SNP)
#DT.frq[,c("SNP","MAF"), with=F]

table(DT.frq$NCHROBS)


# HWE (NUC) ---------------------------------------
#########################################################

setwd("/home/data1/share/sy/GAW19/HWE/NUC")
DT.hwe <- fread("NUC.HWE", sep=";")
#DT.hwe <- DT.hwe[V1==chromosome]
DT.hwe <- DT.hwe[V3=="ALL"]
head(DT.hwe)
setnames(DT.hwe, c("V2","V9"), c("SNP","HWE"))
setkey(DT.hwe, SNP)
#DT.hwe[,c("SNP","HWE"), with=F]






# TDT Regression nuc ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.frq[,c("SNP","MAF"), with=F][DT.hwe[,c("SNP","HWE"), with=F]][DT.tdt[,c("SNP","P_TDT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_TDT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_TDT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#TDT all all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF"
gg2 <- "log10p ~ log10MAF"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#TDT all SNPs with TDT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF"
gg4 <- "log10p ~ log10MAF"

fit3 <- DT[P_TDT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_TDT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_TDT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]



# PDT Regression nuc ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.frq[,c("SNP","MAF"), with=F][DT.hwe[,c("SNP","HWE"), with=F]][DT.pdt[,c("SNP","P_PDT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_PDT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_PDT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#PDT all, all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF"
gg2 <- "log10p ~ log10MAF"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#PDT all SNPs with PDT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF"
gg4 <- "log10p ~ log10MAF"

fit3 <- DT[P_PDT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_PDT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_PDT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]



# FBAT Regression nuc ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.frq[,c("SNP","MAF"), with=F][DT.hwe[,c("SNP","HWE"), with=F]][DT.fbat[,c("SNP","P_FBAT"),with=F]]
str(DT)
DT <- DT[complete.cases(DT),]
#remove pvalues of 1
DT <- DT[P_FBAT!=1]
set(DT, i=NULL, j="log10p",value=-log10(DT[["P_FBAT"]]))
set(DT, i=NULL, j="log10HWE",value=-log10(DT[["HWE"]]))
set(DT, i=NULL, j="log10MAF",value=-log10(DT[["MAF"]]))
DT <- DT[!is.infinite(log10MAF)]
head(DT)

#FBAT all all SNPs-----
gg <- "log10p ~ log10HWE+log10MAF"
gg2 <- "log10p ~ log10MAF"

fit1 <- DT[,lm(as.formula(gg), .SD)]
fit2 <- DT[,lm(as.formula(gg2), .SD)]

#FBAT all SNPs with FBAT p-value < 0.001-----
gg3 <- "log10p ~ log10HWE+log10MAF"
gg4 <- "log10p ~ log10MAF"

fit3 <- DT[P_FBAT<=0.001][,lm(as.formula(gg3), .SD)];summary(fit3)
fit4 <- DT[P_FBAT<=0.001][,lm(as.formula(gg4), .SD)];summary(fit4)

#anova(reduced,full)
summary(fit1)$r.squared-summary(fit2)$r.squared
nrow(DT)
summary(fit3)$r.squared-summary(fit4)$r.squared
nrow(DT[P_FBAT<=0.001])
anova(fit2,fit1)[2,6]
anova(fit4,fit3)[2,6]

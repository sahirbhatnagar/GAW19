# =====================================================
#   FBAT
# =====================================================
#   
# Get FBAT p-values
# fbat # starts fbat
# load map FBAT21.map # load map before ped
# load ped FBAT21.ped
# log fbat21 # turns the logging on
# fbat # calculates p-values
# manually delete first few lines of log file, 
# then run ./makecsv.sh to import into R

setwd("~/share/sy/GAW19/data/transpose/chr21ped")

DT.fbat <- fread("fbat21seq.csv")
set(DT.fbat, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.fbat[["Marker"]])))
set(DT.fbat, i=NULL, j="CHR", value= as.numeric(sub("_\\d*","", DT.fbat[["Marker"]])) )
png("FBAT_chr21.png");manhattan(DT.fbat, snp="Marker", main="FBAT, chr21, seq only");dev.off()

DT.plink <- fread("chr21seqtdt.csv")
png("TDT_chr21.png");manhattan(DT.plink, main="TDT from PLINK, chr21, seq only");dev.off()

png("TRDcompare21.png")
par(mfrow=c(2,2));manhattan(DT.fbat, snp="Marker", main="FBAT, chr21, seq only");
manhattan(DT.plink, main="TDT from PLINK, chr21, seq only");
manhattan(DT, p="Pval", main="PDT, chr 21, seq only");dev.off()


# Info from VCF files -----------------------------------------------------

setwd("~/share/sy/GAW19/data")
DT.frq <- fread("chr1frq.csv")
set(DT.frq, i=NULL, j="A1", value=as.numeric(sub("[ATGC]:","",DT.frq[["V5"]])))
set(DT.frq, i=NULL, j="A2", value=as.numeric(sub("[ATGC]:","",DT.frq[["V6"]])))
set(DT.frq, i=NULL, j="MAF", value=pmin(DT.frq[["A1"]],DT.frq[["A2"]] )   )
DT.hwe <- fread("chr1hwe.csv")

setkey(DT.hwe,POS )
setkey(DT.frq,V2)
head(DT.hwe);head(DT.frq)

DT.all <- DT.frq[DT.hwe]
head(DT.all)
png("maf_hwe_vcf_chr1.png")
DT.all[,plot(MAF, -log10(P))];dev.off()


fam <- fread("PED.csv")
fam
table(fam$PEDNUM)






# VCF files covariate -----------------------------------------------------

setwd("~/share/sy/GAW19/data/vcf/")
DTc <- fread("covar21.csv", sep=";")


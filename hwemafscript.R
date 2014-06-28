##################################
# R source code file that will read in HWE and MAF file, and output HWE manhattan and 
# HWE vs MAF plots
# Created by Sahir, June 26, 2014
# Updated June 27, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
##################################

#command line args should be frq_file hwe_file 
#args <- commandArgs(trailingOnly=TRUE)

# -log10(HWE pvalue) vs MAF -----------------------------------------------

mafhweplots <- function (freqname, hwename){
DT.frq <- fread(freqname)
DT.hwe <- fread(hwename) 
setnames(DT.hwe, c("V1","V2","V9"), c("CHR", "SNP", "P"))
set(DT.hwe, i=NULL, j="BP", value=as.numeric(sub("\\d*_","", DT.hwe[["SNP"]])))
ll <- unlist(strsplit(DT.hwe$V6, "/", fixed=TRUE))
idx <- seq(1, length(ll), by=3)
DT.hwe[,`:=`(aa=ll[idx], Aa=ll[idx+1], AA=ll[idx+2])]
setkey(DT.hwe,SNP);setkey(DT.frq, SNP)
DT.all <- DT.hwe[DT.frq]
return(DT.all)
#DT.all[order(DT.all$V9)[1:3000]]
#DT.all2 <- DT.hwe2[DT.frq]
#xyplot(-log10(V9) ~ MAF, DT.all, grid = TRUE, main=paste("chromosome ",chromosome, sep=" "))
}







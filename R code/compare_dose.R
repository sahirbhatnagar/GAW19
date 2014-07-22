##################################
# R source code file used to check agreement between DOSE and Genotype file
# Created by Sahir, June 27, 2014
# Updated June 27, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: you need to change your working directory to where the data is
# I used pseq (shell script program) to extract ID of sequenced people from vcf files 
##################################


#command line args should be frq_file hwe_file chromosome 
args <- commandArgs(trailingOnly=TRUE)

DT.all <- mafhweplots(args[1],args[2])
DT.dose <- fread(paste("chr",args[3],"-dose.csv", sep=""))
DT.all <- DT.all[CHR==args[3]]

not_seq_not_found_id <- fread("~/share/sy/GAW19/data/not_seq_not_found.id")
DT.dose[, not_seq_not_found_id$x:=NULL]
#number of ids who are founders and sequenced
ids <- ncol(DT.dose)
#counts of Aa, and aa
DT.dose[ ,`:=`(countaa=rowSums(.SD==2), countAa=rowSums(.SD==1), countAA=rowSums(.SD==0)), .SDcol=2:ids]
setnames(DT.dose, 'snp', 'SNP')
setkey(DT.dose, SNP)
DT.geno.dose <- DT.all[DT.dose]

#seq_id <- fread("~/share/sy/GAW19/data/sequenced.id")
#DT.geno.dose[,seq_id$V2:=NULL]

# png(paste("compare_dose_geno_chr_",args[3],sep=""))
# par(mfrow=c(2,2))
# DT.geno.dose[,plot(aa,countaa,xlab='Dose file count of aa', ylab='Genotype file count of aa',main=paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep=""))]
# DT.geno.dose[,plot(Aa,countAa,xlab='Dose file count of Aa', ylab='Genotype file count of Aa',main=paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep="")]
# DT.geno.dose[,plot(AA,countAA,xlab='Dose file count of AA', ylab='Genotype file count of AA',main=paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep="")]
# dev.off()

# png(paste("compare_dose_geno_chr_",args[3],sep=""))
# par(mfrow=c(2,2))
# ggplot(DT.geno.dose, aes(x=countaa, y=as.numeric(aa))) + geom_point(aes(size=-log10(P))) + xlab("Dose file count of aa")+ylab("Genotype file count of aa")+ggtitle(paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep=""))
# ggplot(DT.geno.dose, aes(x=countAa, y=as.numeric(Aa))) + geom_point(aes(size=-log10(P))) + xlab("Dose file count of Aa")+ylab("Genotype file count of Aa")+ggtitle(paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep=""))
# ggplot(DT.geno.dose, aes(x=countAA, y=as.numeric(AA))) + geom_point(aes(size=-log10(P))) + xlab("Dose file count of AA")+ylab("Genotype file count of AA")+ggtitle(paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep=""))
# dev.off()

a <- ggplot(DT.geno.dose, aes(x=countaa, y=as.numeric(aa))) + geom_point(aes(size=-log10(P))) + xlab("Dose file count of aa")+ylab("Genotype file count of aa")+ggtitle(paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep=""))
b <- ggplot(DT.geno.dose, aes(x=countAa, y=as.numeric(Aa))) + geom_point(aes(size=-log10(P))) + xlab("Dose file count of Aa")+ylab("Genotype file count of Aa")+ggtitle(paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep=""))
c <- ggplot(DT.geno.dose, aes(x=countAA, y=as.numeric(AA))) + geom_point(aes(size=-log10(P))) + xlab("Dose file count of AA")+ylab("Genotype file count of AA")+ggtitle(paste('Chromosome ',args[3],", ",nrow(DT.geno.dose)," SNP's",sep=""))

png(paste("compare_dose_geno_chr_",args[3],sep=""))
multiplot(a,b,c,cols=2)
dev.off()

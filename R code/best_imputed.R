##################################
# R source code file used to get the 'Best' Imputed Individuals
# Created by Sahir, July 18, 2014
# Updated July 18, 2014
# 
# NOTE: We are trying to incorporate only the individuals who are best imputed
# 
##################################

#command line args should be chromosome-number
#args <- commandArgs(trailingOnly=TRUE)

#DT.all <- mafhweplots(args[1],args[2])
for (i in seq(1,21,2)){
args<-i
#setwd("~/share/sy/GAW19/DOSE")
DT.dose <- fread(paste("/home/data1/share/sy/GAW19/DOSE/chr",args[1],"-dose.csv", sep=""))
seq_id <- fread("/home/data1/share/sy/GAW19/data/sequenced.id")
#head(seq_id)
setnames(seq_id,c("V1","V2"),c("PEDNUM","ID"))
#remove sequenced, because you just want to see who were the best imputed
DT.dose[, seq_id$ID:=NULL]
#head(DT.dose)
setkey(DT.dose, snp)

near_int <- function(x,range){
  zero.lim <- c(0,0+range)
  one.lim <- c(1-range,1+range)*1
  two.lim <- c(1-range,1+range)*2
  (x>=zero.lim[1] & x <= zero.lim[2]) | (x>=one.lim[1] & x <= one.lim[2])  | (x>=two.lim[1] & x <= two.lim[2])
}
cols.to.sum = colnames(DT.dose)[-1]

#filter for counts within 1% of 0, 1, or 2
DT.dose[ , cols.to.sum := lapply(DT.dose[,cols.to.sum,with=FALSE], 
                    near_int,0.20), with=FALSE]

#head(DT.dose)

k <- data.frame(n.true=colSums(DT.dose[,cols.to.sum,with=FALSE]))
#ID's who are the best imputed
best.imp <- rownames(k)[which(k$n.true==nrow(DT.dose))]

result <- rbind(data.frame(PEDNUM=as.numeric(substr(best.imp,5,6)), ID=best.imp),seq_id)          

result$ID <- sub("T2DG","",result$ID)

#write to file with sequenced.id
cat("writing table to file\n")
write.table(result,file=paste("/home/data1/share/sy/GAW19/data/cleangeno/chr",
                              args[1],"/bestimputed20",args[1],sep=""), 
            row.names=F, col.names=F, quote=F)
cat("DONE writing table for chromosome",i)
}



# Plots for best imputed data ---------------------------------------------

chr <- paste0("chr",seq(1,21,2))
chr <- seq(1,21,2)
within.01<-c(12,16,16,14,23,20,37,34,37,32,119)
within.05<-c(17,27,21,23,31,31,45,46,49,39,143)
within.10<-c(20,36,24,30,40,38,60,57,55,51,157)
within.20<-c(32,47,52,45,67,60,89,80,79,73,195)

dat <- data.frame(chr,within.01,within.05,within.10,within.20)
dat.m <- melt(dat,"chr")
str(dat.m)

p <- ggplot(dat.m, aes(x=variable, y=value, group=chr))
p + geom_line(aes(colour=chr),size=2)

png("plots/best_imputed_summary.png")
ggplot(dat.m, aes(x=variable, y=value, fill=factor(chr))) + geom_bar(position="dodge", stat="identity")
dev.off()


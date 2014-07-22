##################################
# R source code file used to create pedigree plots
# USED also for creating dataset for PEDSTATS program. The PEDSTATS
# program is not working on the server so run locally (installed on ACER laptop)
# Created by Sahir, July 13, 2014
# Updated 
# 
# NOTE: 
##################################


# Bring in the data and modify -------------------------------------------

# Differences between genotype and family file 
# read in family file
DT.tfam <- fread("~/share/sy/GAW19/data/PED.csv")
fam.id <- DT.tfam$ID

#geno.id
DT.geno <- fread("~/share/sy/GAW19/data/rawgeno/chr21-geno.csv", nrows=2)
#get subject ID's for which we have genotype information
id.geno <- colnames(DT.geno)[-1]

#sequenced id
DT.seq <- fread("~/share/sy/GAW19/data/sequenced.id")
seq.id <- DT.seq$V2 

#nuclear families
DT.nuc <- fread("~/share/sy/GAW19/data/nuclear.fam")

#real phenotype
DT.real <- fread("/home/data1/share/sy/GAW19/data/real.pheno")

set(DT.tfam, i=NULL, j="hyper", value=DT.real[["hyper"]])
set(DT.tfam, i=NULL, j="genotyped", value=DT.tfam[["ID"]] %in% id.geno)
sum(DT.tfam$genotyped)
set(DT.tfam, i=NULL, j="sequenced", value=DT.tfam[["ID"]] %in% seq.id)
sum(DT.tfam$sequenced)
set(DT.tfam, i=NULL, j="nuclear", value=DT.tfam[["ID"]] %in% DT.nuc[["V2"]])

#Cranefoot only handles numeric ID's
set(DT.tfam, i=NULL, j="id",value=as.numeric(sub("T2DG","",DT.tfam[["ID"]])))
set(DT.tfam, i=NULL, j="dadid",value=as.numeric(sub("T2DG","",DT.tfam[["FA"]])))
set(DT.tfam, i=NULL, j="momid",value=as.numeric(sub("T2DG","",DT.tfam[["MO"]])))
#replace missing phenotype with 0(cranefoot doesnt handle missings)
#set(DT.tfam, i=which(DT.tfam[["hyper"]]==-9),j="hyper", value=0)

#replace missing phenotype with X for PEDSTATS program
set(DT.tfam, i=which(DT.tfam[["hyper"]]==-9), j="hyper2", value=0)
set(DT.tfam, i=which(DT.tfam[["hyper"]]==1), j="hyper2", value=2)
set(DT.tfam, i=which(DT.tfam[["hyper"]]==0), j="hyper2", value=1)

DT.tfam[which(DT.tfam[["SEX"]]==1), GENDER:="M" ,]
DT.tfam[which(DT.tfam[["SEX"]]==2), GENDER:="F" ,]
DT.tfam$SEQ <- as.numeric(DT.tfam[["sequenced"]])+1
DT.tfam$allaff <- 2




# Some stats on pedigrees -------------------------------------------------

DT.tfam[nuclear==TRUE][genotyped==TRUE]
sum(DT.tfam[nuclear==TRUE]$sequenced)
table(DT.tfam[nuclear==TRUE]$PEDNUM)
DT.tfam[sequenced==TRUE][hyper==1]
table(DT.tfam$hyper)
DT.tfam[ID=="T2DG0300174"|ID=="T2DG0400281"]

#number of founders who are genotyped.. This is where the 117 comes from
sum(DT.tfam[FA==0]$ID %in% id.geno)

#pedstats ped file
setwd("~/share/sy/GAW19/data")
write.table(DT.tfam[,c("PEDNUM","ID","FA","MO","SEX","hyper2","SEQ","allaff"),with=F], 
            file="pedstats.ped", row.names=F, quote=F, sep="\t", col.names=F)



# Founders ----------------------------------------------------------------

DT.tfam[FA==0][sequenced==TRUE]
DT.tfam[FA==0][genotyped==TRUE]
#sequenced
DT.tfam[sequenced==TRUE]
#not sequenced
DT.tfam[sequenced==FALSE]

#genotyped
DT.tfam[genotyped==TRUE]

# Plot Pedigrees ----------------------------------------------------------

head(DT.tfam)

#reduce values of id numbers
# DT.tfam$id2<-mapvalues(DT.tfam$id,DT.tfam$id,seq(1,nrow(DT.tfam),1))
# DT.tfam$dadid2<-mapvalues(DT.tfam$dadid,DT.tfam$id,seq(1,nrow(DT.tfam),1))
# DT.tfam$momid2<-mapvalues(DT.tfam$momid,DT.tfam$id,seq(1,nrow(DT.tfam),1))

library(kinship2)
ped <- with(DT.tfam, pedigree(id=id2, dadid=dadid2, momid=momid2,sex=SEX, affected=hyper, famid=PEDNUM))
png("~/share/sy/GAW19/Pedigrees/pedall.png");
plot(ped);
#pedigree.legend(ped['2'], location="topright")
dev.off()

png("~/share/sy/GAW19/Pedigrees/ped2.png");
ped2<-ped['2']
plot(ped2, col=ifelse(DT.tfam[PEDNUM==2][["sequenced"]],2,1));
#pedigree.legend(ped2, location="topright", radius=0.3)
dev.off()

str(DT.tfam)



setwd("~/share/sy/GAW19/Pedigrees")
write.table(DT.tfam, file="pedigree.txt", row.names=F, quote=F, sep="\t")
table(DT.tfam[["hyper"]])

#entire pedigrees not sequenced
DT.tfam[PEDNUM==14]
DT.tfam[PEDNUM==15]
DT.tfam[PEDNUM==23]
DT.tfam[PEDNUM==25]



# Pair Counts of Sequenced from PEDSTATS ----------------------------------
status <- c("sib_pairs","half_sibs", "cousins", "parent_child","grandparent_grandchild")
neither <- c(843,73,1442,740,1040)
discordant <- c(500,119,1364,953,626)
both <- c(106,62,494,259,156)

dat <- rbind(neither,discordant,both)

colnames(dat)<-status
dat
library(reshape2)
dat.m <- melt(dat)
png("paircounts.png");
ggplot(dat.m, aes(x=Var2, y=value, fill=Var1)) + geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values=c( "blue","red" ,"green")) + xlab("") + ylab("count")+
  guides(fill=guide_legend(title=NULL))+
  scale_x_discrete(labels=c("Sib Pairs", "Half Sibs", "Cousins", "Parent-\nChild", "Grandparent-\nGrandchild"))+
  theme(axis.text.x  = element_text(size=13), axis.title.y = element_text(size=13),
        legend.text = element_text(size = 13), legend.position="top")
dev.off()
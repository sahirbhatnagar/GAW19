##################################
# R source code file for finding genes or rs# associated with small p-value
# to see if any of the small p-values have been previously found for TRD
# Created by Sahir, July 15, 2014
# Updated July 15, 2014
# 
# NOTE: Finding the rs # for smallest p-value
# you need to get the datasets from the regression.R script first
# Based on Everyone affected, sequenced individuals
##################################

head(DT.cov, 50)

#########################################################
#                   SEQUENCED-------------------
#########################################################

# PDT Sequenced ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.cov[,c("NS","SNP","RSID", "DB"),with=F][DT.frq[,c("SNP","MAF"), with=F]][DT.pdt[,c("SNP","P_PDT"),with=F]]
str(DT)
head(DT)

#pdt_seq_rs <- unlist(str_extract_all(DT[order(DT$P_PDT)[1:50]][["RSID"]], "rs\\d*"))

pdt_seq_rs <- unlist(str_extract_all(DT[P_PDT<0.001][["RSID"]], "rs\\d*"))

setwd("~/share/sy/GAW19/Imprinting/rs")
write.table(pdt_seq_rs, file="rs.significant", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)


# TDT Regression Sequenced ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.cov[,c("NS","SNP","RSID", "DB"),with=F][DT.frq[,c("SNP","MAF"), with=F]][DT.tdt[,c("SNP","P_TDT"),with=F]]
str(DT)
head(DT)

tdt_seq_rs <- unlist(str_extract_all(DT[P_TDT<0.0001][MAF!=0][["RSID"]], "rs\\d*"))

write.table(tdt_seq_rs, file="rs.significant", append=TRUE,
            row.names=FALSE, col.names=FALSE, quote=FALSE)


# FBAT Regression Sequenced ---------------------------------------------------------------
rm(DT)
#dont include MQ0 as 9.27 million datapoints are 0
DT <- DT.cov[,c("NS","SNP","RSID", "DB"),with=F][DT.frq[,c("SNP","MAF"), with=F]][DT.fbat[,c("SNP","P_FBAT"),with=F]]
str(DT)
head(DT)

fbat_seq_rs <- unlist(str_extract_all(DT[P_FBAT<0.001][MAF!=0][["RSID"]], "rs\\d*"))

write.table(fbat_seq_rs, file="rs.significant", append=TRUE,
            row.names=FALSE, col.names=FALSE, quote=FALSE)




# Overlap of SNPs and imprinting genes ------------------------------------
#import significant SNPs
#rs.snp <- fread("rs/rs.significant", header=F)

setwd("/home/data1/share/sy/GAW19/Imprinting")
#import rs positions. these were extracteed from ucsc website by pasting in the
# rs/rs.significant text file in the browser
DT.rspos <- fread("rs_position")
setnames(DT.rspos, c("V2","V3"), c("start","stop"))

#write.table(rs.snp, file="", quote=F, col.names=F, row.names=F)

#rs.snp[["V1"]] %in% rspos[["V4"]]

library(GenomicRanges)
# mySNP <- GRanges(seqnames = rspos$V4,
#              ranges = IRanges(start=rspos$V2, end = rspos$V3, names = rspos$V1))#,
#             strand = Rle(strand(rspos$V6)))

mySNP <- GRanges(seqnames = rspos$V1,
             ranges = IRanges(start=rspos$V2, end = rspos$V3))#,
            strand = Rle(strand(rspos$V6)))

seqnames(mySNP)

# import GO genes for impriting
setwd("/home/data1/share/sy/GAW19/Imprinting")
gogenes <- read.table("GOgenes.txt", header=F)

#convert Gene symbol to ENSId
library(biomaRt)
#to use the hsapiens dataset of ensembl dataset
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
head(attributes)

#get gene positions
#attributes[grep("start",attributes$name),]
#attributes[grep("end",attributes$name),]
#attributes[grep("chr",attributes$name),]
gene.pos <- getBM(attributes=c("entrezgene", "ensembl_gene_id","hgnc_symbol","description",
                         "start_position","end_position", "chromosome_name" ), filters = "hgnc_symbol", 
            values = as.vector(gogenes), mart = ensembl)

DT.gene <- as.data.table(gene.pos)
set(DT.gene, i=NULL, j="V1", value=paste("chr",DT.gene[["chromosome_name"]], sep=""))



setkey(DT.gene,V1);setkey(DT.rspos,V1)
tables()
DT.all <- DT.rspos[DT.gene]

str(DT.gene)
str(DT.rspos)



library(IRanges)
rangesA <- split(IRanges(DT.rspos$start, DT.rspos$stop), DT.rspos$V1)
rangesB <- split(IRanges(DT.gene$start_position, DT.gene$end_position), DT.gene$V1)
#which rangesA wholly contain at least one rangesB?
ov <- countOverlaps(rangesB, rangesA, type="within")



#create grange for imprinting genes
# gene <- GRanges(seqnames = gene.pos$hgnc_symbol,
#                 ranges = IRanges(start=gene.pos$start_position, end = gene.pos$end_position))#, 
#                                  names = gene.pos$ensembl_gene_id),
#                 strand = Rle(strand(c("*"))))

gene <- GRanges(seqnames=paste("chr",gene.pos$chromosome_name,sep=""),
                ranges = IRanges(start=gene.pos$start_position, end = gene.pos$end_position))



mtch <- findOverlaps(gene,mySNP)

as.matrix(mtch)


seqnames(gene);ranges(gene)
seqnames(mySNP);ranges(mySNP)

# library(AnnotationHub)
# hub <- AnnotationHub()
# ## data exploration
# length(names(hub))                  # resources available
# md <- metadata(hub)                 # DataFrame
# md[grep("SNP",md$Description),]
# hub$goldenpath.hg19.encodeDCC  # tab completion
#   
#   ## retrieval
#   res <- hub$goldenpath.hg19.encodeDCC.wgEncodeUwTfbs.wgEncodeUwTfbsNhlfCtcfStdPkRep1.narrowPeak_0.0.1.RData
# class(res)         # GRanges representation of ENCODE bed file 

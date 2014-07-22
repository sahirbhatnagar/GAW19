##################################
# R source code file used to input Genotypes of missing ID's
# Created by Sahir, June 30, 2014
# Updated July 1, 2014
# hosted on Github repo 'sahirbhatnagar/gaw19'
# NOTE: this is different from geno_split.R in that it no longer needs to split genotypes
# since PLINK 1.90 can handle grouped alleles for regular ped files only
##################################
data.clean <- function(filename, pedfile, col.alloc=100){
    DT <- fread(filename)
    #get subject ID's for which we have genotype information
    id.geno <- colnames(DT)[-1]

    # Differences between genotype and family file 
    # read in family file
    DT.tfam <- fread(pedfile)
    fam.id <- DT.tfam$ID

    #ID's with missing genotypes
    missing.geno <- setdiff(fam.id, id.geno)
  
    alloc.col(DT,col.alloc)
    for (i in seq_along(missing.geno)) {
    set(DT, i=NULL, j=missing.geno[i], value="XX")
      }

    setcolorder(DT, c("snp",fam.id))

    # we need to add in map information after since we have to
    # transpose first
#     set(DT, i=NULL, j="chr", value= as.numeric(sub("_\\d*","", DT[["snp"]])) )
#     set(DT, i=NULL, j="base", value = as.numeric(sub("\\d*_","", DT[["snp"]])) )
#     set(DT, i=NULL, j="cM", value=0) 
# 
#     setcolorder(DT, c("chr","snp","cM","base",fam.id))  
    
    return(DT)

}

#command line args should be --filename chr.digits 
args <- commandArgs(trailingOnly=TRUE)

DT <- data.clean(filename=args[1], pedfile=args[2], col.alloc=4500) 
chr <- as.numeric(sub("_\\d*","", DT[["snp"]][1]))

write.table(DT, file=tempfile(pattern=paste("chr",chr,".tped",sep=""), tmpdir=args[3]) ,
            col.names=TRUE, row.names=FALSE, quote=FALSE)



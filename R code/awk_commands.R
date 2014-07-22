#write snp id's in 1 column to new file
awk 'NR==1 {for (i=2;i<=NF;i++) print $i}' chr00.csv.ped > id1

#count number of columns
awk '{print NF; exit}' merge.ped 

#print to terminal first few columns and lines -d is the delimiter, 
# change ' ' to ',' for .csv files
head chr21.temp3.ped | cut -d' ' -f1-7

#print specific records
awk 'NR>=5&&NR<=9' "file"



DT<-fread("chr21.temp.ped")

DT[1:10,1:10,with=F]

ped <- fread("PED.csv")
ped[, MZTWIN:=NULL];set(ped, i=NULL, j="pheno", value=1)
write.table(ped, file="PEDheader.fam",
            col.names=TRUE, row.names=FALSE, quote=FALSE)


DT2 <- fread("chr21.temp2.ped")
DT2[1:10,1:10,with=F]

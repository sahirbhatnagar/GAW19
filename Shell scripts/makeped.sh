#!/bin/bash
# file to create PED file. Use this script only after
# having run steps 1 to 4 from README.PED

# command line args: chromosome_number

# 5) Grab the SNP names to make the .MAP file
awk 'NR==1 {for (i=2;i<=NF;i++) print $i}' chr$1/$1.temp.ped > chr$1/snp.id$1

# 6) Create .MAP file (keep separate from snp.id in case we want to filter SNPs and create different MAP files):
awk '{$2=$1 ; $3=0 ; $4=$1; t=sub(/_[0-9]*/,"",$1);y=sub(/[0-9]*_/,"",$4) }1' chr$1/snp.id$1 > chr$1/chr$1.map

# 7) Bring in family data from PEDheader.fam. this file is the same as PED.fam, except it has a header:
join -12 -21 ~/share/sy/GAW19/data/PEDheader.fam chr$1/$1.temp.ped > chr$1/$1.temp2.ped

# 8) Remove T2DG from ID's as FBAT does not accept characters, and swap first two columns to get it in the write PLINK format
awk '{$1=substr($1,5)}1{if ($3==0) $3=0; else $3=substr($3,5)}1{if ($4==0) $4=0; else $4=substr($4,5)}1' chr$1/$1.temp2.ped | awk '{ t = $1; $1 = $2; $2 = t; print; }' > chr$1/$1.temp3.ped

# 9) Run PLINK to recode the alleles to 0/1/2, since FBAT only accepts this (note that this recodes the phenotypes to 2 which means affected. The --1 flag tells PLINK that control/case status is coded 0/1, and then --recode changes this to 1/2)
plink.1.90betadev --ped chr$1/$1.temp3.ped --map chr$1/chr$1.map --missing-genotype X --recode 12 --1 --make-bed --remove ~/share/sy/GAW19/data/mztwins.txt --out chr$1/chr$1;

# 10) Replace missing genotypes with 0 for FBAT
sed 's/X/0/g' chr$1/chr$1.ped > chr$1/FBAT$1.ped

# 11) Create FBAT map file which is in the format: 
# marker_name chr# genetic_pos physical_pos sex_link
awk '{$2=$1 ; $3=0 ; $4=$1; t=sub(/_[0-9]*/,"",$2);y=sub(/[0-9]*_/,"",$4);$5=0 }1' chr$1/snp.id$1 > chr$1/FBAT$1.map

# 12) Create PDT map file which is in the format : 
#Chromosome marker_name location
awk '{$2=$1 ; $3=$1 ; $4=2; t=sub(/_[0-9]*/,"",$1);y=sub(/[0-9]*_/,"",$3);}1' chr$1/snp.id$1 > chr$1/PDT$1.map
#or use the map file from PLINK (this is better option if you are filtering for HWE or MAF)
# awk '{ $3=$2 ; $4=2; y=sub(/[0-9]*_/,"",$3);}1' chr$1/chr$1seq_hwe.map > chr$1/PDT$1.map



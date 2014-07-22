#!/bin/bash
# Script for plink commands

# make binary data files using TPED files
#plink.1.90beta --noweb --tped chr$1.tped --tfam trans.tfam --missing-genotype X --make-bed --remove mztwins.txt --keep sequenced.id --no-pheno --out chr$1;

# make binary data files using PED files 
#plink.1.90betadev --file chr$1 --missing-genotype X --recode 12 --1 --make-bed --out chr$1;

# make binary file sequenced only
#plink.1.90betadev --file chr$1 --missing-genotype X --make-bed --recode --keep sequenced.id.fbat --out chr$1seq;

# make binary file bestimputed(within 1 percent of 0,1,2 for all markers) + sequenced
 plink.1.90betadev --file chr$1/chr$1 --missing-genotype X --make-bed --recode --keep chr$1/bestimputed$1 --out chr$1/chr$1best;

# load vcf file
#plink.1.90beta --noweb --vcf chr1-seq.vcf --vcf-half-call 'missing' --make-bed --out chr$1vcf

# MAF all
# plink.1.90betadev --bfile chr$1/chr$1 --freq --out chr$1/chr$1all ;

# MAF seq
# plink.1.90betadev --bfile chr$1/chr$1seq --freq --out chr$1/chr$1seq ;

# MAF nuc
# plink.1.90betadev --bfile chr$1/chr$1nuc --freq --out chr$1/chr$1nuc ;

# HWE all
# plink.1.90betadev --bfile chr$1/chr$1 --hardy midp --out chr$1/chr$1all 

# HWE seq
# plink.1.90betadev --bfile chr$1/chr$1seq --hardy midp --out chr$1/chr$1seq 

# HWE nuc
# plink.1.90betadev --bfile chr$1/chr$1nuc --hardy midp --out chr$1/chr$1nuc 

# HWE on sequenced people (vcf files only have sequenced people)
# plink --bfile chr$1 --hardy --out chr$1 

# check for mendelian errors
# plink --noweb --bfile chr21 --mendel

# bring in phenotype file, the --1 recodes for unaffected/affected to be 0/1 instead of 1/2
# plink --noweb --bfile chr21 --pheno pheno.txt --pheno-name HTN_1 --assoc --1

# extract founders and write binary file for founders
# plink --noweb --bfile mydata21 --filter-founders --make-bed --out founders --recode --remove mztwins.txt

# MAF on founders (by default, PLINK only calculates MAF using founders)
# plink.1.90beta --bfile chr$1 --freq --keep sequenced.id --out chr$1 ;

# HWE exact test on founders (by default, PLINK only calculates HWE using founders)
# plink.1.90beta --bfile chr$1 --hardy midp --keep sequenced.id --out chr$1 ;

# HWE standard asymptotic test on founders (by default, PLINK only calculates HWE using founders)
# plink --noweb --bfile chr$1 --hardy2 --maf 0.01 --geno 0.1 --mind 0.1 --out chr$1.hwe2 

# family based association tdt based on sequenced individuals, parent of origin analysis
#plink.1.90betadev --bfile chr$1/chr$1seq --tdt poo --out chr$1/chr$1poo

# family based association tdt based on 1 nuclear family per pedigree
# plink.1.90beta --bfile chr$1 --pheno pheno.txt --pheno-name EVERYONE_AFFECTED --tdt --1 --keep nuclear.fam --out chr$1

# family based association tdt based on sequenced individuals
# plink.1.90betadev --bfile chr$1seq --tdt --out chr$1seq

# family based association tdt with cases switched to missing, and controls to cases
# plink --noweb --bfile mydata21 --pheno pheno.txt --pheno-name SWITCH_PHENO --tdt --1 --out switch_pheno

# family based association tdt with real hypertension phenotype called hyper based on sequenced only
# plink.1.90betadev --bfile chr$1/chr$1seq --pheno ~/share/sy/GAW19/data/real.pheno --pheno-name hyper --tdt --1 --out switch_pheno

# Check for mendel errors (i am using to extract nuclear families)
# plink --noweb --bfile chr$1 --mendel --out chr$1

#!/bin/bash
# script to extract VCF information

# command line argument is chromosome number

bcftools query -f '%CHROM;%POS;%ID;%INFO/NS;%INFO/AF;%INFO/DB;%INFO/RSID;%INFO/STR;%INFO/STZ;%INFO/CBR;%INFO/CBZ;%INFO/CSR;%INFO/IOZ;%INFO/IOR;%INFO/AOZ;%INFO/AOI;%INFO/MQ0;%INFO/MQ10;%INFO/MQ20;%INFO/MQ30\n' chr$1-seq.vcf > temp$1

cat ~/share/sy/GAW19/data/vcf.header temp$1 > covar$1.csv

rm temp$1

mv -f covar$1.csv ~/share/sy/GAW19/VCF/

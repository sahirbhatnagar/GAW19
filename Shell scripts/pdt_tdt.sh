#!/bin/bash
# script to run PDT analysis and TDT analysis. This will output manhattan plot for 
# sequenced, all and nuclear fam for each chromosome. 
# Do the combined manhattan plot manually.
# note that the makeped.sh script creates the ped file FBAT$1.ped which contains everyone
# (i.e. sequenced + non sequenced people), so this script creates the non-sequenced. 
# Note also that FBAT.ped can be used for PDT software
# In this script, we will create another ped file which is based on the sequenced only

# 1) Run plink to create ped file for sequenced only 
#plink.1.90betadev --bfile chr$1/chr$1 --missing-genotype X --make-bed --recode 12 --1 --keep ~/share/sy/GAW19/data/sequenced.id.fbat --out chr$1/chr$1seq;

# 1.2) Run plink to create ped file for nuclear only 
#plink.1.90betadev --bfile chr$1/chr$1 --missing-genotype X --make-bed --recode 12 --1 --keep ~/share/sy/GAW19/data/nuclear.fam.fbat --out chr$1/chr$1nuc;

# 1.3) Run plink to create ped file for best imputed 
plink.1.90betadev --bfile chr$1/chr$1 --missing-genotype X --make-bed --recode 12 --1 --keep chr$1/bestimputed10$1 --out chr$1/chr$1best;

# 2) Replace missing genotypes with 0 for FBAT and PDT
# NOTE THAT PED FILE PRODUCED FOR FBAT CAN BE USED DIRECTLY FOR PDT
#sed 's/X/0/g' chr$1/chr$1seq.ped > chr$1/FBAT$1seq.ped
#sed 's/X/0/g' chr$1/chr$1nuc.ped > chr$1/FBAT$1nuc.ped
sed 's/X/0/g' chr$1/chr$1best.ped > chr$1/FBAT$1best.ped

# 3) Create pdt control text files for everyone. I have written a script for this:
#sh ~/share/sy/GAW19/scripts/makepdtctrl.sh $1 chr$1/pdtctrl$1.txt
#sh ~/share/sy/GAW19/scripts/makepdtctrlseq.sh $1 chr$1/pdtctrl$1seq.txt
#sh ~/share/sy/GAW19/scripts/makepdtctrlnuc.sh $1 chr$1/pdtctrl$1nuc.txt
sh ~/share/sy/GAW19/scripts/makepdtctrlbest.sh $1 chr$1/pdtctrl$1best.txt

# need to make it executable otherwise pdt2 doesnt work for some odd reason
#chmod a+x chr$1/pdtctrl$1.txt
#chmod a+x chr$1/pdtctrl$1seq.txt
#chmod a+x chr$1/pdtctrl$1nuc.txt
chmod a+x chr$1/pdtctrl$1best.txt

# 4) Run the pdt:
#pdt2 chr$1/pdtctrl$1.txt
#pdt2 chr$1/pdtctrl$1seq.txt
#pdt2 chr$1/pdtctrl$1nuc.txt
pdt2 chr$1/pdtctrl$1best.txt

# 5) Delete the second line from each output file:
#sed '2d' chr$1/pdt$1\_marker.out > chr$1/pdt$1.out
#sed '2d' chr$1/pdt$1seq\_marker.out > chr$1/pdt$1seq.out
#sed '2d' chr$1/pdt$1nuc\_marker.out > chr$1/pdt$1nuc.out
sed '2d' chr$1/pdt$1best\_marker.out > chr$1/pdt$1best.out


# 6) Run makecsv.sh on each of the output files for easy read into R
#sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/pdt$1.out chr$1/pdt$1.csv
#sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/pdt$1seq.out chr$1/pdt$1seq.csv
#sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/pdt$1nuc.out chr$1/pdt$1nuc.csv
sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/pdt$1best.out chr$1/pdt$1best.csv

# 8) Run TDT analysis from PLINK
# family based association tdt based on everyone
#plink.1.90betadev --bfile chr$1/chr$1 --tdt --out chr$1/chr$1all

# family based association tdt based on sequenced individuals
#plink.1.90betadev --bfile chr$1/chr$1seq --tdt --keep ~/share/sy/GAW19/data/sequenced.id.fbat --out chr$1/chr$1seq

# family based association tdt based on sequenced individuals
#plink.1.90betadev --bfile chr$1/chr$1nuc --tdt --keep ~/share/sy/GAW19/data/nuclear.fam.fbat --out chr$1/chr$1nuc

# family based association tdt based on sequenced and best imputed individuals
plink.1.90betadev --bfile chr$1/chr$1best --tdt --keep chr$1/bestimputed$1 --out chr$1/chr$1best

# 9) run makecsv on TDT files
#sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/chr$1all.tdt chr$1/chr$1all.csv
#sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/chr$1seq.tdt chr$1/chr$1seq.csv
#sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/chr$1nuc.tdt chr$1/chr$1nuc.csv
sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/chr$1best.tdt chr$1/chr$1best.csv

#10) Create PLots
#Rscript ~/share/sy/GAW19/scripts/PDT_analysis.R $1

# 11) Move files to PDT and TDT folders:
#mv -f chr$1/pdt$1.csv ~/share/sy/GAW19/PDT/ALL/
#mv -f chr$1/pdt$1seq.csv ~/share/sy/GAW19/PDT/SEQ/
#mv -f chr$1/pdt$1nuc.csv ~/share/sy/GAW19/PDT/NUC/
#mv -f chr$1/pdt$1best.csv ~/share/sy/GAW19/PDT/BEST/
#mv -f chr$1/pdt$1best.csv ~/share/sy/GAW19/PDT/BEST5/
mv -f chr$1/pdt$1best.csv ~/share/sy/GAW19/PDT/BEST10/

#mv -f chr$1/chr$1all.csv ~/share/sy/GAW19/TDT/ALL/
#mv -f chr$1/chr$1seq.csv ~/share/sy/GAW19/TDT/SEQ/
#mv -f chr$1/chr$1nuc.csv ~/share/sy/GAW19/TDT/NUC/
#mv -f chr$1/chr$1best.csv ~/share/sy/GAW19/TDT/BEST/
#mv -f chr$1/chr$1best.csv ~/share/sy/GAW19/TDT/BEST5/
mv -f chr$1/chr$1best.csv ~/share/sy/GAW19/TDT/BEST10/

#!/bin/bash

#script for post processing output of FBAT software
# command line inputs are chromosome number| all or seq or nuc or best | ALL or SEQ or NUC or BEST
# for i in $(seq 1 2 21); do echo "sh ~/share/sy/GAW19/scripts/fbatplot.sh $i all ALL" | qsub -cwd -N run$i -V; done
# and
# for i in $(seq 1 2 21); do echo "sh ~/share/sy/GAW19/scripts/fbatplot.sh $i seq SEQ" | qsub -cwd -N run$i -V; done



#1) Remove the unwanted rows
sed '1,4d;6d' chr$1/fbat$1$2 | sed '$d' | sed '$d' > chr$1/fbat$1$2temp

#2) make csv
sh ~/share/sy/GAW19/scripts/makecsv.sh chr$1/fbat$1$2temp chr$1/fbat$1$2.csv

#3) delete temp file
rm chr$1/fbat$1$2temp

#3) move to FBAT folder
mv -f chr$1/fbat$1$2.csv ~/share/sy/GAW19/FBAT/$3/

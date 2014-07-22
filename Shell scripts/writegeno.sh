#!/bin/bash

for file in chr*
do
       echo "Rscript geno_split.R $file 2" | qsub -cwd -N run1 -V
done

#!/bin/bash
#takes three user inputs
#run using ./datasplit.sh filename output_dir_to_create num_of_lines_to_split


# check if an input filename was passed as a command
# line argument:
if [ ! $# == 3 ]; then
  echo "Please specify 3 arguments: filename outputdir no.lines"
  exit
fi

# create a directory to store the output:
mkdir $2

# create a temporary file containing the header without
# the content:
head -n 1 $1 > $2/header.csv

# create a temporary file containing the content without
# the header:
tail -n +2 $1 > $2/content.csv

# split the content file into multiple files of 5 lines each:
split -l $3 -a 2 -d $2/content.csv $2/chr

# loop through the new split files, adding the header
# and a '.csv' extension:
for f in $2/chr*; do cat $2/header.csv $f > $f.csv; rm $f; done;

# remove the temporary files:
rm $2/header.csv
rm $2/content.csv

#copy necessary files to newly created directory
#cp pheno.txt mztwins.txt PED.csv plink.sh trans.tfam geno_split.R writegeno.sh transpose.sh $2/

#!/bin/bash
# commands to combine results
# the command line argument is the extension of the files returned by PLINK

#this command will keep the header from the first file and remove the header from subsequent files and combine them
awk -F ";" 'FNR==1 && NR!=1{next;}{print}' $1* > all.$1 ;

# this command will remove spaces and add ; to seperate each column, this will allow you to use the data.table::fread function
#sed 's/^[ ]*//g' all.$1 |sed -e 's/\([0-9a-zA-Z\.]*\)  */\1;/g' > all$1.csv

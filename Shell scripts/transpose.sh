#!/bin/bash
# this script will read in the split genotype csv files, 
# and transpose the data
# cmd line inputs are filename directory

NEW_UUID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $1 > $2/ped$NEW_UUID

#done

# sed 's/,/ /g' chr21-geno.csv | ./transpose.sh > chr21.ped
# this command removes first column and first row and puts a , between each allele
# head chr00.csv -n 5 | cut -d "," -f 2- | tail -n +2 | sed 's/.\{1\}/&,/g' | sed 's/,,,/,/g' | sed 's/,/ /g'

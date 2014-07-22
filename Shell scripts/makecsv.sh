#!/bin/bash


# this command will remove spaces and add ; to seperate each column, this will allow you to use the data.table::fread function
# first argument is input, second argument is output
sed 's/^[ ]*//g' $1 |sed -e 's/\([0-9a-zA-Z\.]*\)  */\1;/g' > $2

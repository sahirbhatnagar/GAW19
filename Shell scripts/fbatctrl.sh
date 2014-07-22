#!/bin/bash
# script for making control file for FBAT so that you can 
# run fbat ---> run fbat.script
# there is no other way of automating the FBAT for each chromosome
# so you need to manually enter fbat for each chromosome
# =========================================

echo "load FBAT$1.map
load FBAT$1best.ped
log fbat$1best10
fbat -e
log off
quit " > chr$1/fbat$1.ctrl


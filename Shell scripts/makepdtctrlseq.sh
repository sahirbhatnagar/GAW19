#!/bin/bash
# this is the same as makepdtctrl.sh except this is for sequenced only people

echo "#first line is a comment, it must begin with #
#file that contains genotypes
pedigree_file:                         chr$1/FBAT$1seq.ped
#<standard frequency file MERLIN uses. Will NOT print outfile_by_marker if this option is selected>
#freq_file:                            example.freq
map_file:                              chr$1/PDT$1.map
# Result file, sorted by p-value of markers
outfile_by_pvalue:                   chr$1/pdt$1seq_pvalue.out
# Result file, sorted by marker
outfile_by_marker:                     chr$1/pdt$1seq_marker.out
# The fields in the .ped file that are not alleles. Usually 6, when using SIMLA 8
non_maker_fields:                      6
# Always leave this at 1 unless you really know what you're doing
window:                                1
# Use geno-PDT, 0 = no (means regular PDT), 1 = use geno-PDT
geno_pdt:                              1
# Number of CPUs that PDT2 will try to use. 4 maximum
max_cpus:                              1
# 0=use all info; 1=triads only; 2=sibs only; 3=only triads if available, otherwise use discordant sibs
options:                               0
# size of largest pedigree in .ped file. Going under will cause crash, going over will waste some memory.
max_fam_size:                         107" > $2
# list for "special" PDT2 markers. ONLY those markers will be printed to the .out file that is "sorted by marker", all others will not be printed
#special_mrk_list:                     special_markers.txt
verbose:                              1" 



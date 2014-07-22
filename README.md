GAW19 
=====

Assessing Transmission Ratio Distortion in Extended Families: A Comparison of Analysis Methods
--------------------------------------------------------------------------------------------------------------------------
### Bhatnagar SR, Greenwood CMT, Labbe A
#### McGill University

#### Abstract:
A statistical departure from Mendel's law of segregation is known as transmission ratio distortion (TRD). Though well documented in many other organisms, the extent of TRD and its influence in the human genome remains incomplete. Using GAW19 whole genome sequence (WGS) data from 20 large Mexican American pedigrees, our goal was to identify potentially distorted regions in the genome using family based association methods such as the transmission disequilibrium test (TDT), the pedigree disequilibrium test (PDT) and the family based association test (FBAT). Preliminary results showed an unusually high number of TRD signals identified by the TDT, but this phenomenon could not be replicated by the PDT or FBAT. Applying these tests to different subsets of the data we found the TDT to be very sensitive to imputed genotypes. Regression analysis of TRD test p-values controlling for minor allele frequency and quality control checks showed that Hardy Weinberg p-values are associated with this inflation. While the TDT appears confounded by imputation of SNP's, the PDT and FBAT seem to offer a more robust alternative when searching for TRD loci in WGS data from extended families. 

Source code used in the analysis for [GAW19](http://www.gaworkshop.org/gaw19/index.html)

Refer to README files to see how `R` and `shell` scripts were used. 

Family Based Association Software Used:
* [FBAT - Harvard University](http://www.hsph.harvard.edu/fbat/fbat.htm)
* [PLINK 1.9](https://www.cog-genomics.org/plink2)
* [Pedstats](http://www.sph.umich.edu/csg/abecasis/Pedstats/index.html)
* [PDT](http://hihg.med.miami.edu/cgesg/statistical-programming-page)




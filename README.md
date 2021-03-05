# Assessment of significance of conditionally independent GWAS signals

#Open codes-simulation under the null.zip and follow the steps bellow:

#step1: 50 times randomly pick a set of 22 independent SNPs, one from each autosome with MAF>0.1.
run 01-50times_22SNPs_MAF.R

#step2: extract the imputed genotypes of pre-defined 22 SNPs to the separate files.
run 02-extract-imputed-GT.22.sh

#step3: Calculate the risk disease and make 3 binary phenotypes for binary association tests
run 03-simulation-nMAF.Binary.R

#step4: Make PLINK format fam file and run 3 binary association tests
run 04-GWAS.Binary.sh

#step5: Perform meta-analysis of the three GWAS results with METAL 

#step6: Identify all SNPs (Top22.txt) between pre-defined 22 SNPs that reach genomewide significance at established significant criterion (5*10^(-8)) and make 1-MB sourrouning regions as an input for --cojo-file parameter (GCTA analysis).
run 06-GCTA.R

#step7: Run condition GCTA analysis on all index SNPs (Top22.txt) in the respective 1-Mb surrounding regions.
run 07-GCTA.sh

#step8: For each index SNP ascertain 1-Mb window size surrounding LD region for INTERSNP tool
run 08-5Mb_LD_region.sh

#step9: run INTERSNP test (INTERSNP.zip file)

#step10: Determin the SNP-specific alpha-threshold for each candidate SNP and find secondary signals with "quasi adaptive method".
run 10-Final_G5-function.R

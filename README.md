# BovineQTL - TWAS study

## Main points of workflow

## Preprocessing
 - This script performs QC, alignment and quantification of RNA-seq data from all n = 123 animals
 - The script performs TPM normalisation of counts
 - The script remaps variants from UMD3.1 to ARS1.2 taking account of potential strand flips
 - The script phases target data using beagle and imputes variants up to WGS using Minimac4
 - The script filters out variants to leave ~ 3.8 million variants for the eQTL analysis
 
  

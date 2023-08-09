# BovineQTL - TWAS study

## Main points of workflow

 ## Reference files
 - Reference genome downloaded from https://ftp.ensembl.org/pub/release-110/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
 - Reference gtf file downloaded from https://ftp.ensembl.org/pub/release-110/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.110.gtf.gz
 - WGS reference file from Dutta et al., 2020 https://www.nature.com/articles/s41467-020-18550-1

## Preprocessing
 - This script performs QC, alignment and quantification of RNA-seq data from all n = 123 animals
 - The script performs TPM normalisation of counts
 - The script remaps variants from UMD3.1 to ARS1.2 taking account of potential strand flips
 - The script phases target data using beagle and imputes variants up to WGS using Minimac4
 - The script filters out variants to leave ~ 3.8 million variants for the eQTL analysis
 


  

import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_UPDATED.vcf.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.expr_tmm_inv.bed.gz", matrix = config["groups"])


rule split_groups:
    input:
        filtered_imputed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.FILTERED.dose.vcf.gz",
        sample_names = lambda wildcards: expand(f'/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/{config["groups"][wildcards.group]}_names.txt') 
    output:
        subsetted =  expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{{group}}_IMPUTED_CLEAN.vcf.gz")

    params:
        prefix = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/",
        cohort = "{group}"

    shell:
        '''
        # use VCF tools as it calculates allele frequencies

        vcftools --gzvcf {input.filtered_imputed} --keep {input.sample_names} --recode --recode-INFO-all --stdout | bgzip -c > {params.prefix}{params.cohort}_IMPUTED_CLEAN.vcf.gz
        '''

rule af_calculation:
    input: 
        subsetted = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_CLEAN.vcf.gz"

    output:
        updated_subset =   "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_UPDATED.vcf.gz"
    
    shell:
        '''
        bcftools +fill-tags {input.subsetted} -Oz -o {output.updated_subset}
        '''
 
rule prepare_4_eQTL:
    input:
        file_counts = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_clean.txt",
        file_tpm = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_TPM.txt", 
        annot_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_annotation_MF2.csv",
        vcf_fn = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}_IMPUTED_UPDATED.vcf.gz",
        known_cov = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt", 
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/prepare_data_tensorqtl.R"
    output:
        bed_4_tensor = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.expr_tmm_inv.bed.gz",
        scree = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}_qtl_scree_plot.pdf",
        geno_eigenvectors = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.PCA_eigenvect.txt",
        geno_eigenvalues = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.PCA_var.txt",
        covariates = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.covariates.txt"
    params:
        cohort = "{matrix}",
        temporary_bed =  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.expr_tmm_inv.bed",
        temporary_snprelate = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.ccm.gds"
    threads: 10
    shell:
        '''
        Rscript {input.script} {input.file_counts} {input.file_tpm} {input.annot_file} {input.vcf_fn} {params.cohort} {params.temporary_bed} {threads} \
        {output.scree} {params.temporary_snprelate} {output.geno_eigenvectors} {output.geno_eigenvalues} {input.known_cov} {output.covariates}

        bgzip {params.temporary_bed}

        '''
import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_UPDATED.vcf.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.expr_tmm_inv.bed.gz", matrix = config["groups"]),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_tss_file.gtf",
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor.fam", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl.txt.gz", group = config["groups"]),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_fdr0.05.txt', group = config["groups"]),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.{chr}.parquet', group = config["groups"], chr = autosomes),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.{chr}.txt.gz', group = config["groups"], chr = autosomes),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_independent_qtl.txt.gz", group = config["groups"]),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_eGene_upset.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Gtex_thresholds_determined.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_replication_correlation.pdf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_pairs.ALL.txt.gz",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/numbers_stats_cis.txt"


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
rule edit_annotation:
    input:
        annotation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"

    output:
        tss_annotation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_tss_file.gtf"
    
    shell:
        '''
        awk 'BEGIN{{OFS="\t"}}{{if($3=="gene"){{if($7=="+"){{$5=$4;$4=$5-1}}else{{$4=$5;$4=$4-1}};print $0}}}}' {input.annotation} > {output.tss_annotation}
        '''
    

rule prepare_4_eQTL:
    input:
        file_counts = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_clean.txt",
        file_tpm = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_TPM.txt", 
        annot_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_tss_file.gtf",
        vcf_fn = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}_IMPUTED_UPDATED.vcf.gz",
        known_cov = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt", 
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/prepare_data_tensorqtl.R"
    output:
        bed_4_tensor = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{matrix}.expr_tmm_inv", ".bed.gz", ".bed.gz.tbi"),
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

        tabix -p bed {output.bed_4_tensor[0]}
        '''
rule prepare_genotype_data:
    input:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_UPDATED.vcf.gz"
    output:
        plink_format = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor", ".bed", ".bim", ".fam")
    
    params:
        plink = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor"
    shell:
        '''
        plink --cow --vcf {input.vcf} --keep-allele-order --make-bed --out {params.plink}
        '''

rule qtl_mapping_permute:
    input:
        bed_4_tensor = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv", ".bed.gz", ".bed.gz.tbi"),# Phenotypes
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt", # Covariates
    output:
        outputdir="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl.txt.gz"
    resources:
        mem_mb = 6000,
        threads = 20
    params:
       plink_prefix_path="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor", # Genotypes
       group = "{group}"
    shell:
        '''
        mkdir -p /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results

        prefix="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{params.group}"

        python3 -m tensorqtl {params.plink_prefix_path} {input.bed_4_tensor[0]} ${{prefix}} \
        --mode cis \
        --window 1000000 \
        --seed 1856 \
        --permutations 10000 \
        --covariates {input.covariates_file} \
        --fdr 0.05
        '''

rule eGenedetection:
    input:
        eqtl_permute="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl.txt.gz", # Phenotypes
        script="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/eGene_detection.R"
    output:
        eGenes= '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_fdr0.05.txt'
    shell:
        '''
        Rscript {input.script} {input.eqtl_permute} {output.eGenes} 0.05 # FDR
        '''

rule nominal_mappinng:
    input:
        expression_bed= multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv", ".bed.gz", ".bed.gz.tbi"), # Phenotypes
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt", # Covariates
    output:
        parquet_files= '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.{chr}.parquet'
    resources:
        mem_mb = 6000,
        threads = 20
    params:
       plink_prefix_path="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor", # Genotypes
       group = "{group}",
    shell:
        '''
        mkdir -p /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results

        prefix="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{params.group}"
        python3 -m tensorqtl {params.plink_prefix_path} {input.expression_bed[0]} ${{prefix}} \
        --covariates {input.covariates_file} \
        --window 1000000 \
        --mode cis_nominal
        '''

rule parquet_2_txt:
    input:
        parquet_files= '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.{chr}.parquet',
        script = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/parquet2txt.py'
    output:
        txt_files = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.{chr}.txt.gz'
    shell:
        '''
        python3 {input.script} {input.parquet_files} {output.txt_files}
        '''

rule conditional_eQTL:
    input:
        bed_4_tensor = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv", ".bed.gz", ".bed.gz.tbi"),# Phenotypes
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt", # Covariates
    output:
        independent = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_independent_qtl.txt.gz"
    params:
        plink_prefix_path="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_genotypes_tensor", # Genotypes
        group = "{group}"
    shell:
        '''
        mkdir -p /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/result

        prefix="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{params.group}"


        python3 -m tensorqtl {params.plink_prefix_path} {input.bed_4_tensor[0]} ${{prefix}} \
        --covariates {input.covariates_file} \
        --cis_output ${{prefix}}.cis_qtl.txt.gz \
        --mode cis_independent
        '''
rule cis_conditional_upset_plotting:
    input:
        perm = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_fdr0.05.txt", group = config["groups"]),
        conditional = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_independent_qtl.txt.gz", group = config["groups"]),
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/eGene_Intersection_conditional_number.R"
    output:
        upset = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_eGene_upset.pdf",
        conditional = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/conditional_eGenes.pdf",
        conditional_zoom = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/conditional_zoomed.pdf",
        ridge_plot_distance = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/ggridges_TSS.pdf",
        cor_TSS_slope = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cor_TSS_slope.pdf"
    shell:
        '''
        Rscript {input.script} {input.perm} {output.upset} {input.conditional} {output.conditional} \
        {output.conditional_zoom} {output.ridge_plot_distance} {output.cor_TSS_slope}
        '''

rule Cattle_Gtex_correction:
    input:
        perm = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Blood.permutations.2rd.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/eGene_detection_Gtex.R"
    
    output:
        eQTL_output = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Gtex_thresholds_determined.txt"
    shell:
        '''
        Rscript {input.script} {input.perm} {output.eQTL_output} 0.05
        '''
rule cis_eQTL_replication:
    input:
        perm = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_fdr0.05.txt", group = config["groups"]),
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/eQTL_Replication.R",
        gtex_nominal = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Blood.nominals.2rd.txt",
        gtex_corrected = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Gtex_thresholds_determined.txt",
    
    output:
        pioest = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_replication_pi1_bootstrap.pdf",
        correlation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_replication_correlation.pdf"
    
    shell:
        '''
        Rscript {input.script} {input.perm} {input.gtex_nominal} {input.gtex_corrected} {output.pioest} {output.correlation}
        '''
rule merged_nominal:
    input:
        ALL_nominal = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_pairs.{autosome}.txt.gz", autosome = autosomes),
        CONTROL_nominal = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_pairs.{autosome}.txt.gz", autosome = autosomes),
        INFECTED_nominal = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_pairs.{autosome}.txt.gz", autosome = autosomes)

    output:
        ALL_nominal_merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_pairs.ALL.txt.gz",
        CONTROL_nominal_merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_pairs.ALL.txt.gz",
        INFECTED_nominal_merged = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_pairs.ALL.txt.gz"

    shell:
        '''
        cat {input.ALL_nominal} > {output.ALL_nominal_merged}
        cat {input.CONTROL_nominal} > {output.CONTROL_nominal_merged}
        cat {input.INFECTED_nominal} > {output.INFECTED_nominal_merged}
        '''

rule number_eQTLs:
    input:
        permute = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_fdr0.05.txt", group = config["groups"]),
        nominal = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.cis_qtl_pairs.ALL.txt.gz", group = config["groups"]),
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/cis-eQTL_numbers.R"
    output:
        stats = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/numbers_stats_cis.txt"
    shell:
        '''
        Rscript {input.script} {input.permute} {input.nominal} {output.stats}
        '''
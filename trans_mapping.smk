import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.nominal.best.txt.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.permute.best.txt.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_qtltools.bed.gz", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates_qtltools.txt", group = config["groups"])

rule edit_phenotype_bed:
    input:
        bed_4_tensor = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv", ".bed.gz", ".bed.gz.tbi")

    output:
        bed_4_qtltools = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_qtltools", ".bed.gz", ".bed.gz.tbi")

    shell:
        '''
        zcat {input.bed_4_tensor[0]} | awk '{{ $4=$4" . +"; print $0 }}' | tr " " "\t" | bgzip -c > {output.bed_4_qtltools[0]}
        tabix -p bed {output.bed_4_qtltools[0]}
        '''

rule edit_covariates:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/EQTL/Covariate_edit_qtltools.R",
        covs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates.txt"
    output:
        qtl_tools_cov = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates_qtltools.txt"
    shell:
        '''
        Rscript {input.script} {input.covs} {output.qtl_tools_cov}
        '''
rule trans_mapping:
    input:
        vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}_IMPUTED_UPDATED.vcf.gz",
        bed_4_qtltools = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_qtltools", ".bed.gz", ".bed.gz.tbi"),
        covariates_file="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.covariates_qtltools.txt", # Covariates
    output:
        nominal=multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.nominal",".best.txt.gz", ".bins.txt.gz", ".hits.txt.gz"),
        permute=multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.permute",".best.txt.gz", ".bins.txt.gz", ".hits.txt.gz")
    resources:
        mem_mb = 6000,
        threads = 20
    params:
       prefix_nominal = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.nominal",
       prefix_permute = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/{group}.trans.permute", 
       group = "{group}"
    singularity: "docker://jogrady/qtltools:1.3.1"
    shell:
        '''
        # nominal
        QTLtools trans --vcf {input.vcf} --bed {input.bed_4_qtltools[0]} --nominal --cov {input.covariates_file} --normal --threshold 1e-5 --out {params.prefix_nominal}


        # permutation
        QTLtools trans --vcf {input.vcf} --bed {input.bed_4_qtltools[0]} --threshold 1e-5 --cov {input.covariates_file} --normal --permute --out {params.prefix_permute} --seed 1894
        '''
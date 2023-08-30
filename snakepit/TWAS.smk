import pandas as pd
import os
import subprocess

# File to preprocess RNA-seq data and genotype data

autosomes = [str(i) for i in range(1,30)] # bovine autosomes
rule all:
    input:
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{gwas}_GWAS_final.txt", gwas = config["gwas_sets"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_correlations_FDR1.txt", group = config["groups"]),
        expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}.cis_qtl_pairs_ALL_chr.txt.gz", group = config["groups"])


rule GWAS_check:
    input:
        ARS_GWAS  = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/{gwas}_ARS.csv",
        UMD_GWAS = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/{gwas}_sires.txt",
        Run6 = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/1000Bull_UMD_GWAS_REF_ALT_SNPs.tab",
        ARS_alleles = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/1_GWAS_flipping.R"
    output:
        final_GWAS = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{gwas}_GWAS_final.txt"
    params:
        set = "{gwas}"      

    shell:
        '''
        mkdir -p /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs
        # Note these are massive files, please have the computational resources necessary
        # Otherwise you may need to split by chromosome

        Rscript {input.script} {input.ARS_GWAS} {input.UMD_GWAS} {input.Run6} {input.ARS_alleles} {params.set} {output.final_GWAS}
        '''
rule MeQTL_correlation:
    input:
        exp = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/{group}.expr_tmm_inv.bed.gz",
        TF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_TF.txt",
        coTF = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_Cof.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/TWAS/2_MediatorQTLAnalysis.R"

    output:
        med_loc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_Mediators_location.txt",
        med_intensity = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_Mediators_intensity.txt",
        med_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/{group}_correlations_FDR1.txt"

    shell:
        '''
        Rscript {input.script} {input.exp} {input.TF} {input.coTF} {output.med_loc} {output.med_intensity} {output.med_results}
        '''

rule Get_nominal_associations:
    input:
        all = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_pairs.{autosome}.txt.gz", autosome = autosomes),
        con = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_pairs.{autosome}.txt.gz", autosome = autosomes),
        infec = expand("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_pairs.{autosome}.txt.gz", autosome = autosomes),
    output:
        nominal_eqtl = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/ALL.cis_qtl_pairs_ALL_chr.txt.gz",
        con_eqtl = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/CONTROL.cis_qtl_pairs_ALL_chr.txt.gz",
        infec_eqtl = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/inputs/INFECTED.cis_qtl_pairs_ALL_chr.txt.gz"
    shell:
        '''
        cat {input.all} > {output.nominal_eqtl}
        cat {input.con} > {output.con_eqtl}
        cat {input.infec} > {output.infec_eqtl}
        '''

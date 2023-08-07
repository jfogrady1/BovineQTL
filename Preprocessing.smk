import pandas as pd
import os
import subprocess

rule all:
    input:
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/multiqc_report.html',
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.out', individual=config["samples"]),
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt',
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt',
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_clean.txt', matrix = config["groups"]),
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/{matrix}_matrix_TPM.txt', matrix = config["groups"]),
        multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data", ".bed", ".bim", ".fam"),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/Admixture_plot.pdf",
        multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".eigenvec", ".eigenval"),
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_liftOver_updated.map",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Master_final.txt",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/PLOT_AlleleFrequencies.pdf"

rule fastqc:
    input:
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/{config["samples"][wildcards.individual]}_{{N}}.fq.gz', N = (1,2))
    output:
        reads = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        12
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/'

rule multiqc:
    input:
        reads = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/{individual}_{N}_fastqc.zip', individual = config['samples'], N = (1,2))
    output:
        report='/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/
        """

rule Alignment:
    input:
        genome = config["star_index"],
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/{config["samples"][wildcards.individual]}_{{N}}.fq.gz', N=(1,2)),
    params:
        prefix = lambda wildcards: f'{config["samples"][wildcards.individual]}'
    output:
        aligned = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.final.out',
        interlog = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.out'
    threads: 20

    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''


rule featureCounts:
    input:
        bam = expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam', sample = config["samples"]),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/GCF_002263795.1_ARS-UCD1.2_genomic.gff"
    output:
        count_matrix = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt'
    threads: 37
    shell:
        '''
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g Dbxref
        '''

rule cleanup:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.sh",
        count_matrix_in = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt',
        final_script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.py",
        temporary = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts_temp.txt"
    output:
        count_matrix = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts2.txt',
        count_matrix_2 = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts_clean.txt',
        cleaned = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt'
    shell:
        '''
        tail -n+2 {input.count_matrix_in} > {output.count_matrix}
        cut -f 1,7-130 {output.count_matrix} > {output.count_matrix_2} # number of samples, get rid of fields we do not want

        python3 {input.final_script} {output.count_matrix_2} {input.temporary} {output.cleaned}
        '''

rule TPM normalisation:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/TPM_normalisation.R",
        annotation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_annotation_MF2.csv",
        count_matrix_in = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt',
    output:
        count_matrix_all = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/ALL_matrix_clean.txt',
        count_matrix_control = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/CONTROL_matrix_clean.txt',
        count_matrix_infected = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/INFECTED_matrix_clean.txt',
        count_TPM_all = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/ALL_matrix_TPM.txt',
        count_TPM_control = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/CONTROL_matrix_TPM.txt',
        count_TPM_infected = '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/INFECTED_matrix_TPM.txt',
        
    shell:
        '''
        Rscript {input.script} {input.count_matrix_in} {input.annotation} {output.count_matrix_all} {output.count_matrix_control} {output.count_matrix_infected} {output.count_TPM_all} {output.count_TPM_control} {output.count_TPM_infected}
        '''

##### Now onto the Genotype data


rule pruning:
    input:
        plink_files = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/SNP_data", ".bed", ".bim", ".fam")
    output:
        filtered = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data", ".bed", ".bim", ".fam"),
        pruned = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".bed", ".bim", ".fam")
    shell:
        '''
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/SNP_data --maf 0.1 --geno 0.95 --mind 0.95 --hwe 0.000001 --allow-extra-chr --double-id --autosome --make-bed --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data --indep-pairwise 1000 5 0.2 --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/pruned_snpset
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data --extract /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/pruned_snpset.prune.in --make-bed --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned
        '''

rule admixture:
    input:
       bed_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.bed",
       admixture_exc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/admixture",
    output:
       admixture_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q"
    shell:
        '''
        {input.admixture_exc} --cv {input.bed_file} 2

        mv /home/workspace/jogrady/eqtl_study/eqtl_nextflow/SNP_Pruned* /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/
        '''

rule admixture_plot:
    input:
       admixture_file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q",
       plot = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Admixture_plot.R"
    output:
       pdf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/Admixture_plot.pdf"
    shell:
        '''
        Rscript {input.plot} {input.admixture_file} /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE Admixture_plot
        '''

rule pca_plot:
    input:
       filtered = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data", ".bed", ".bim", ".fam"),
       pruned = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/pruned_snpset.prune.in"
    output:
       pdf = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned", ".eigenvec", ".eigenval")
    shell:
        '''
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/filtered_SNP_data --double-id --allow-extra-chr --extract {input.pruned} --pca --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned 
        '''

rule liftover:
    input:
       chain = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/bosTau6ToBosTau9.over.chain",
       axiom_master = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/Axiom_GW_Bos_SNP_1.na35.annot.csv",
       script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/LiftOver_remapping.R"

    output: 
       liftover = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_liftOver_updated.map"

    shell:
       """
       Rscript {input.script} {input.chain} {input.axiom_master} {output.liftover}
       """

rule remove_palindromic:
    input:
        axiom_master = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/Axiom_GW_Bos_SNP_1.na35.annot.csv",
        liftover = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_liftOver_updated.map",
        schnabel = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/9913_ARS1.2_648875_BOS1_marker_name_180910", ".map", ".REF_ALLELE"),
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Remove_spurious_SNPs.R"
    output:
        removed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Palindromic_SNPs_removed.txt",
        retained_IDs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Final_raw_IDs_4_analysis.txt",
        liftedover_coord = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_LOver_Rob_final.map",
        final_UMD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Master_final.txt"

    shell:
        '''
        Rscript {input.script} {input.axiom_master} {input.liftover} {input.schnabel} {output.removed} {output.retained_IDs} {output.liftedover_coord} {output.final_UMD}
        '''  
rule remapping:
    input:
        retained_IDs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Final_raw_IDs_4_analysis.txt",
        liftedover_coord = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ARS_LOver_Rob_final.map",
        removed = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Palindromic_SNPs_removed.txt",
        schnabel = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/9913_ARS1.2_648875_BOS1_marker_name_180910", ".REF_ALLELE", ".REF"),
        final_UMD = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Master_final.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/Flipping.R"
    output:
        remapped_plink = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/remapped_SNP_data", ".bed", ".bim", ".fam"),
        remapped_autosomes = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/input_data_final", ".bed", ".bim", ".fam"),
        vcf_4_flipping = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP.Data.For.Flipping.vcf",
        flip_decision = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/ID_Flip_decision.txt",
        final_vcf_4_imputation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
        remapped_alleles = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Axiom_Remapped_Flipped_FINAL_A1.vcf",
        flip_coordinates = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/regions_file_subset.txt"

    shell:
        '''
        # Update coordinates in plink files

        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/SNP_data \
        --keep-allele-order --extract {input.retained_IDs} --update-chr {input.liftedover_coord} 2 1 --update-map {input.liftedover_coord} 3 1 \
        --make-bed --allow-extra-chr --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/remapped_SNP_data
        
        # Restrict to autosomes
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/remapped_SNP_data --keep-allele-order \
        --chr 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 --allow-extra-chr --make-bed \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/input_data_final

        # Remove genotypes and samples which are missing
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/input_data_final --keep-allele-order --mind 0.05 --make-bed \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_1
        
        
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_1 --keep-allele-order --geno 0.05 --make-bed \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_2
        
        
        # For avoidance of any doubt, set the major allele as reference allele and recode as VCF
        plink --cow --bfile /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_data_missing_2 --a1-allele {input.schnabel[0]} --recode vcf \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP.Data.For.Flipping
        
        # Now perform flipping
        Rscript {input.script} {input.final_UMD} {output.vcf_4_flipping} {input.schnabel[1]} {output.flip_decision} {output.final_vcf_4_imputation} {output.remapped_alleles} {output.flip_coordinates}
        
        
        # Redirect files for final 
        cat {output.remapped_alleles} >> {output.final_vcf_4_imputation}
        '''  

rule flip_check:
    input:
        final_vcf_4_imputation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/SNP_Data_Flipped.vcf",
        imputation_master = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/CattleACDr2.variants.phased.filtered.lowMissingnessIDs.unrel.biallelic.GQ25.CR75.annotated.noMW.GQ25.CR75.IMPUTED", ".vcf.gz",".vcf.gz.tbi"),
        European_samples = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/SNP_data/European_samples.txt",
        flip_coordinates = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/regions_file_subset.txt",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/SNP_data/PLOT_strand_flipping.R"
    output:
        AF_gg = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/PLOT_AlleleFrequencies.pdf",
        AF_ss = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Smooth_scatter.png",
        spurious_snps = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Spurious_Flipped_SNPs.txt",
        spurious_snps_pos = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Spurious_SNPS_to_remove_positions.txt"
    shell:
        '''
        # subset genotyped SNPs and European breeds in WGS_reference
        bcftools view -R {input.flip_coordinates} -S {input.European_samples} {input.imputation_master[0]} -Oz \
        -o /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/WGS_array_European_subsetted.vcf.gz
        
        # use vcftools to calculate frequency
        vcftools --vcf {input.final_vcf_4_imputation} --freq \
        --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Reference_Allele_freq
        
        # Calculate the allele frequencies
        vcftools --gzvcf /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/WGS_array_European_subsetted.vcf.gz \
        --freq --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/European_WGS_Allele_freq
        
        # Identify spurious SNPs and plot allele frequencies
        Rscript {input.script} \
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/European_WGS_Allele_freq \
        /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/flip_check/Reference_Allele_freq.frq \
        {output.AF_gg} {output.AF_ss} {output.spurious_snps} {output.spurious_snps_pos}
        
        '''
import pandas as pd
import os
import subprocess

rule all:
    input:
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/multiqc_report.html',
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.out', individual=config["samples"]),
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/gene_counts.txt',
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt'
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
     
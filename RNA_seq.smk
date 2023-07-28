import pandas as pd
import os
import subprocess

rule all:
    input:
        '/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Fastqc/multiqc_report.html',
        expand('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/{individual}_Log.out', individual=config['samples'])


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

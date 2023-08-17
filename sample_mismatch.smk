rule all:
    input:
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/test_mbv.txt"

rule sample_check:
    input:
        vcf = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED", ".vcf.gz", ".vcf.gz.tbi"),
        bam = multiext("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Alignment/C001_Aligned.sortedByCoord.out", ".bam", ".bam.bai")
    output:
        check = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/test_mbv.txt"
    singularity: "docker://jogrady/qtltools:1.3.1"

    shell:
        '''
        QTLtools mbv --vcf {input.vcf[0]} --bam {input.bam[0]} --reg 1 --out {output.check} 
        '''
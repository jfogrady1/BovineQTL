rule all:
    input:
        "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Volcano_plot.pdf"

rule differential_expression:
    input:
        counts = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/count_matrix_clean.txt",
        coldata = "/home/workspace/jogrady/data/eqtl_study/eqtl_nextflow/RNA_seq/covariate_RNA_seq.txt",
        pca_geno = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenvec",
        admixture = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q",
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/Differential_Expression.R"
    output:
        Volcano = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Volcano_plot.pdf",
        Gprofiler = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Gprofiler_enrichment.pdf",
        PCA = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/PCA_Transcriptomic_Admixture.pdf",
        RAW_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/ALL_results.txt",
        Signif_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Significant_results.txt",
        Gprofiler_results = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Gprofiler_results.txt",
        Kirsten_result = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Kirsten_stats.txt"

    shell:
        '''
        Rscript {input.script} {input.counts} {input.coldata} {input.admixture} {output.Volcano} {output.Gprofiler} {output.PCA} \
        {output.RAW_results} {output.Signif_results} {output.Gprofiler_results} {output.Kirsten_result}
        '''
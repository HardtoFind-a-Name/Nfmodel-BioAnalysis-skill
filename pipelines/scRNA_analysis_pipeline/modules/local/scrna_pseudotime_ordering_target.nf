process SCRNA_PSEUDOTIME_ORDERING_TARGET {
    tag "${safe_label}_pseudotime_ordering_target"
    publishDir "${params.outdir}/${params.pseudotime_target_SN}_scrna_pseudotime_ordering_target/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.target_ordering_monocle_cds.rds'), emit: target_monocle_cds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_pseudotime_ordering_target.R \
      --input-rds ${input_rds} \
      --ordering-gene-file ${params.pseudotime_ordering_gene_file} \
      --ordering-gene-col ${params.pseudotime_ordering_gene_col} \
      --outdir . \
      --min-ordering-genes ${params.pseudotime_min_ordering_genes} \
      --sn ${params.pseudotime_target_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

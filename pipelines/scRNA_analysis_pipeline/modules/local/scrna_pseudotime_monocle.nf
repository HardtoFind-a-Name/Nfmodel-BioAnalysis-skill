process SCRNA_PSEUDOTIME_MONOCLE {
    tag "${safe_label}_pseudotime_monocle"
    publishDir "${params.outdir}/${params.pseudotime_SN}_scrna_pseudotime_monocle/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.monocle_cds.rds'), emit: monocle_cds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_pseudotime_monocle.R \
      --input-rds ${input_rds} \
      --outdir . \
      --hvg-nfeatures ${params.pseudotime_hvg_nfeatures} \
      --min-ordering-genes ${params.pseudotime_min_ordering_genes} \
      --sn ${params.pseudotime_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

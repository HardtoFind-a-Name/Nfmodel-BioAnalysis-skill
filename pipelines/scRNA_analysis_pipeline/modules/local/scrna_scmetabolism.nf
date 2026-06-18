process SCRNA_SCMETABOLISM {
    tag "${safe_label}_scmetabolism"
    publishDir "${params.outdir}/${params.scmetabolism_SN}_scrna_scmetabolism/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.scmetabolism_scored.rds'), emit: scmetabolism_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_scmetabolism.R \
      --input-rds ${input_rds} \
      --outdir . \
      --key-celltype "${key_celltype}" \
      --top-n-pathways ${params.scmetabolism_top_n_pathways} \
      --sn ${params.scmetabolism_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

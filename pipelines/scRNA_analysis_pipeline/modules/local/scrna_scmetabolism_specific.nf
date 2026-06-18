process SCRNA_SCMETABOLISM_SPECIFIC {
    tag "${safe_label}_scmetabolism_specific"
    publishDir "${params.outdir}/${params.scmetabolism_specific_SN}_scrna_scmetabolism_specific/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.scmetabolism_specific_scored.rds'), emit: scmetabolism_specific_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_scmetabolism_specific.R \
      --input-rds ${input_rds} \
      --outdir . \
      --key-celltype "${key_celltype}" \
      --top-n-pathways ${params.scmetabolism_top_n_pathways} \
      --min-fdr ${params.scmetabolism_min_fdr} \
      --min-delta ${params.scmetabolism_min_delta} \
      --highlight-groups "${params.scmetabolism_highlight_groups}" \
      --sn ${params.scmetabolism_specific_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

process SCRNA_GSVA {
    tag "${safe_label}_gsva"
    publishDir "${params.outdir}/${params.gsva_SN}_scrna_gsva/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_gsva.R \
      --input-rds ${input_rds} \
      --outdir . \
      --gmt-files "${params.gsva_gmt_files}" \
      --min-geneset-size ${params.gsva_min_geneset_size} \
      --max-geneset-size ${params.gsva_max_geneset_size} \
      --threshold-t ${params.gsva_threshold_t} \
      --threshold-p ${params.gsva_threshold_p} \
      --sn ${params.gsva_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

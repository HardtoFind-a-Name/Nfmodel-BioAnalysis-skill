process GSEA_RISK_NES {
    tag 'gsea_risk_nes'
    publishDir "${params.outdir}/${params.gsea_risk_nes_SN}_gsea_risk_nes", mode: 'copy', overwrite: true
    input:
    path risk_file
    path count_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/gsea_risk_nes.R \
      --risk-file ${risk_file} --count-matrix ${count_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.gsea_risk_nes_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --gmt-dir ${params.gmt_dir} --gmt-gsea_risk ${params.gmt_gsea_risk} \
      --gsea-sort-by ${params.gsea_sort_by}
    """
}

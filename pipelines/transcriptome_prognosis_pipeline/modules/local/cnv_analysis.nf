process CNV_ANALYSIS {
    tag 'cnv_analysis'
    publishDir "${params.outdir}/${params.cnv_SN}_cnv", mode: 'copy', overwrite: true
    input:
    path risk_file
    path expr_file
    path gene_file
    path immune_input
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/cnv_analysis.R \
      --risk-file ${risk_file} --expr ${expr_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.cnv_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --gene-file ${gene_file} \
      --gistic-matrix ${params.gistic_matrix} \
      --download-dir ${params.cnv_download_dir} \
      --immune-input ${immune_input}
    """
}

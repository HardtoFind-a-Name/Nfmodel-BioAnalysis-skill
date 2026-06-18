process GISTIC_INPUT {
    tag 'gistic_input'
    publishDir "${params.outdir}/${params.gistic_SN}_gistic", mode: 'copy', overwrite: true
    input:
    path risk_file
    path expr_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.txt', optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/gistic_input.R \
      --risk-file ${risk_file} --expr ${expr_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.gistic_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --download-dir ${params.cnv_download_dir}
    """
}

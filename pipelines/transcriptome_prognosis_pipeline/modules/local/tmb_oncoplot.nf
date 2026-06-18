process TMB_ONCOPLOT {
    tag 'tmb_oncoplot'
    publishDir "${params.outdir}/${params.tmb_SN}_tmb", mode: 'copy', overwrite: true
    input:
    path risk_file
    path expr_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/tmb_oncoplot.R \
      --risk-file ${risk_file} --maf-file ${params.maf_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.tmb_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

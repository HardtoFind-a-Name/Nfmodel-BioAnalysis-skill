process IPS_ANALYSIS {
    tag 'ips_analysis'
    publishDir "${params.outdir}/${params.ips_SN}_ips", mode: 'copy', overwrite: true
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
    ${params.rscript} ${projectDir}/scripts/ips_analysis.R \
      --risk-file ${risk_file} --expr ${expr_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.ips_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --ips-file ${params.ips_file}
    """
}

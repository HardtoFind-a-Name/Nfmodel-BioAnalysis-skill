process CIRCOS_PLOT {
    tag 'circos_plot'
    publishDir "${params.outdir}/${params.circos_SN}_circos", mode: 'copy', overwrite: true
    input:
    val gwas_file
    path risk_file
    path gene_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/circos_plot.R \
      --gwas-file ${gwas_file} --risk-file ${risk_file} --gene-file ${gene_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.circos_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

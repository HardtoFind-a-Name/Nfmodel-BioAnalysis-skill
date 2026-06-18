process CIBERSORT_IMMUNE {
    tag 'cibersort_immune'
    publishDir "${params.outdir}/${params.cibersort_SN}_cibersort", mode: 'copy', overwrite: true
    input:
    path risk_file
    path expr_file
    path gene_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    path '11_cnv_immune_input.csv', emit: cnv_immune_input, optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/cibersort_immune.R \
      --risk-file ${risk_file} --expr ${expr_file} --gene-file ${gene_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.cibersort_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

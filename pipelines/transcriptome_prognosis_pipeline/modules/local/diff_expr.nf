process DIFF_EXPR {
    tag 'diff_expr'
    publishDir "${params.outdir}/${params.diff_expr_SN}_diff_expr", mode: 'copy', overwrite: true
    input:
    path gene_file
    path expr_file
    path group_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/diff_expr.R \
      --gene-file ${gene_file} --expr ${expr_file} --group ${group_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.diff_expr_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

process LASSO_SELECT {
    tag 'lasso_select'

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/01_lasso_select", mode: 'copy', overwrite: true

    input:
    path expr_matrix
    path survival_table
    path gene_list

    output:
    path '01.lasso_selected_genes_lambda_min.csv', emit: lasso_genes
    path '02.lasso_final_genes_coef.csv', emit: lasso_coef
    path '00.lasso_input_merged.csv', optional: true
    path '03.lasso_cv_curve.pdf', optional: true
    path '03.lasso_cv_curve.png', optional: true
    path '04.lasso_coef_path.pdf', optional: true
    path '04.lasso_coef_path.png', optional: true
    path '05.lasso_summary.txt', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/lasso_select.R \
      --expr ${expr_matrix} \
      --surv ${survival_table} \
      --gene-list ${gene_list} \
      --outdir . \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

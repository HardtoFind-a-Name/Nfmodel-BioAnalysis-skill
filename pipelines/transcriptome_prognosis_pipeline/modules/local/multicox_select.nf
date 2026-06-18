process MULTICOX_SELECT {
    tag "${model_id}_multicox"

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/${model_id}/multicox", mode: 'copy', overwrite: true

    input:
    path expr_matrix
    path survival_table
    path gene_list
    val model_id

    output:
    path '03.multiCox_final_genes_coef.csv', emit: multicox_coef
    path '04.multiCox_gene_list.csv', emit: multicox_genes
    path '01.multiCox_all_results.csv', optional: true
    path '02.multiCox_ph_test.csv', optional: true
    path '05.multiCox_forest.pdf', optional: true
    path '05.multiCox_forest.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/multicox_select.R \
      --expr ${expr_matrix} \
      --surv ${survival_table} \
      --gene-list ${gene_list} \
      --model-id ${model_id} \
      --outdir . \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

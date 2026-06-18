process STEPCOX_SELECT {
    tag "${model_id}_stepcox"

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/${model_id}/stepcox", mode: 'copy', overwrite: true

    input:
    path expr_matrix
    path survival_table
    path gene_list
    val model_id

    output:
    path '01.stepcox_final_genes_coef.csv', emit: stepcox_coef
    path '02.stepcox_all_results.csv', optional: true
    path '03.stepcox_ph_test.csv', optional: true
    path '04.stepcox_forest.pdf', optional: true
    path '04.stepcox_forest.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/stepcox_select.R \
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

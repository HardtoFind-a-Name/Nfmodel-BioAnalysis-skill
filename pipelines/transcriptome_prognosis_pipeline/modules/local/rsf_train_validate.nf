process RSF_TRAIN_VALIDATE {
    tag "${model_id}_rsf"

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/${model_id}", mode: 'copy', overwrite: true

    input:
    path expr_matrix
    path survival_table
    path gene_list
    path validation_sheet
    val model_id

    output:
    path 'train_summary.csv', emit: train_summary
    path 'validation_summary.csv', emit: validation_summary
    path '04.rsf_model.rds', optional: true
    path '01.feature_importance.csv', optional: true
    path '01.feature_importance.pdf', optional: true
    path '01.feature_importance.png', optional: true
    path '02.train_risk_score.csv', optional: true
    path '03.rsf_grid_search.csv', optional: true
    path 'train_135y/*.pdf', optional: true
    path 'train_135y/*.png', optional: true
    path 'train_357y/*.pdf', optional: true
    path 'train_357y/*.png', optional: true
    path '*/01.risk_score.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    def alt_flag = params.run_time_set_alt ? "--time-roc-days-alt ${params.time_roc_points_days_alt}" : "--time-roc-days-alt none"
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/rsf_train_validate.R \
      --expr ${expr_matrix} \
      --surv ${survival_table} \
      --gene-list ${gene_list} \
      --validation-sheet ${validation_sheet} \
      --model-id ${model_id} \
      --outdir . \
      --time-roc-days ${params.time_roc_points_days} \
      ${alt_flag} \
      --ntree-range ${params.rsf_ntree_range} \
      --mtry-range ${params.rsf_mtry_range} \
      --nodesize-range ${params.rsf_nodesize_range} \
      --seed ${params.rsf_seed} \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

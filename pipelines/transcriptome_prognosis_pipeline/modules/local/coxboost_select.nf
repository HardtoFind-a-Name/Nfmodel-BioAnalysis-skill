process COXBOOST_SELECT {
    tag 'coxboost_select'

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/A6_coxboost_coef/coxboost", mode: 'copy', overwrite: true

    input:
    path expr_matrix
    path survival_table
    path gene_list

    output:
    path '02.coxboost_selected_features.csv', emit: coxboost_coef
    path '01.coxboost_model.rds', optional: true
    path '03.coxboost_grid_search.csv', optional: true
    path '04.coxboost_summary.txt', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/coxboost_select.R \
      --expr ${expr_matrix} \
      --surv ${survival_table} \
      --gene-list ${gene_list} \
      --outdir . \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --maxstepno-range ${params.coxboost_maxstepno_range} \
      --penalty-range ${params.coxboost_penalty_range} \
      --auc-threshold ${params.model_auc_threshold} \
      --time-roc-days ${params.time_roc_points_days} \
      --seed ${params.rsf_seed}
    """
}

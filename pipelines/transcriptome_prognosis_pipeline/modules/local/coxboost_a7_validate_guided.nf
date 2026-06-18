process COXBOOST_A7_VALIDATE_GUIDED {
    tag "${model_id}_${cohort_id}"

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/${model_id}/${cohort_id}", mode: 'copy', overwrite: true

    input:
    tuple val(cohort_id), val(time_roc_days), path(v_expr), path(v_surv), path(t_expr), path(t_surv), path(gene_list)
    val model_id

    output:
    path 'train_summary.csv', emit: train_summary
    path 'validation_summary.csv', emit: validation_summary
    path '01.coxboost_model.rds', optional: true
    path '02.coxboost_final_genes_coef.csv', optional: true
    path '02.train_risk_score.csv', optional: true
    path '03.coxboost_grid_search.csv', optional: true
    path 'train_plots/*.pdf', optional: true
    path 'train_plots/*.png', optional: true
    path 'validate/*.pdf', optional: true
    path 'validate/*.png', optional: true
    path 'validate/01.risk_score.csv', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/coxboost_a7_validate_guided.R \
      --train-expr ${t_expr} \
      --train-surv ${t_surv} \
      --gene-list ${gene_list} \
      --cohort-id ${cohort_id} \
      --validation-expr ${v_expr} \
      --validation-surv ${v_surv} \
      --time-roc-days ${time_roc_days} \
      --maxstepno-range ${params.coxboost_maxstepno_range} \
      --penalty-range ${params.coxboost_penalty_range} \
      --auc-threshold ${params.model_auc_threshold} \
      --seed ${params.rsf_seed} \
      --model-id ${model_id} \
      --outdir . \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

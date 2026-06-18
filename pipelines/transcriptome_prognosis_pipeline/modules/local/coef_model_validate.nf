process COEF_MODEL_VALIDATE {
    tag "${model_id}_${cohort_id}"

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/${model_id}/validate/${cohort_id}", mode: 'copy', overwrite: true

    input:
    tuple val(cohort_id), path(expr_matrix), path(survival_table), val(time_roc_days)
    path coef_file
    val model_id

    output:
    path '06.validation_summary.csv', emit: validation_summary
    path '01.risk_score.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/coef_model_validate.R \
      --cohort-id ${cohort_id} \
      --expr ${expr_matrix} \
      --surv ${survival_table} \
      --coef-file ${coef_file} \
      --model-id ${model_id} \
      --outdir . \
      --time-roc-days ${time_roc_days} \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

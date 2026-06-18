process COEF_MODEL_TRAIN {
    tag "${model_id}_train"

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/${model_id}/train", mode: 'copy', overwrite: true

    input:
    path expr_matrix
    path survival_table
    path coef_file
    val model_id

    output:
    path '99.coef_file.csv', emit: coef_file_out
    path '01.train_risk_score.csv', emit: train_risk
    path '06.train_summary.csv', emit: train_summary
    path '00.model_meta.csv', optional: true
    path '135y/*.pdf', optional: true
    path '135y/*.png', optional: true
    path '357y/*.pdf', optional: true
    path '357y/*.png', optional: true

    script:
    def alt_flag = params.run_time_set_alt ? "--time-roc-days-alt ${params.time_roc_points_days_alt}" : "--time-roc-days-alt none"
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/coef_model_train.R \
      --expr ${expr_matrix} \
      --surv ${survival_table} \
      --coef-file ${coef_file} \
      --model-id ${model_id} \
      --outdir . \
      --time-roc-days ${params.time_roc_points_days} \
      ${alt_flag} \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

process CLINICAL_NOMOGRAM {
    tag 'clinical_nomogram'

    publishDir "${params.outdir}/${params.multi_model_SN}_clinical", mode: 'copy', overwrite: true

    input:
    path survival_table
    path clinical_table
    path risk_table

    output:
    path '01.clinical_selected.csv', optional: true
    path '02.survival_selected.csv', optional: true
    path '03.risk_selected.csv', optional: true
    path '04.merged_raw.csv', optional: true
    path '05.merged_cleaned_before_complete_case.csv', optional: true
    path '06.analysis_dataset_complete_case.csv', emit: analysis_dataset
    path '07.univariate_overall_test.csv', optional: true
    path '07.univariate_cox_all.csv', optional: true
    path '08.univariate_cox_display.csv', optional: true
    path '08.univariate_significant_features.csv', optional: true
    path '09.univariate_cox_ph_test.csv', optional: true
    path '10.univariate_cox_forest.pdf', optional: true
    path '10.univariate_cox_forest.png', optional: true
    path '11.multivariate_cox_all.csv', optional: true
    path '11.multivariate_overall_test.csv', optional: true
    path '12.multivariate_cox_display.csv', optional: true
    path '13.multivariate_cox_ph_test.csv', optional: true
    path '13.multivariate_cox_ph_test.pdf', optional: true
    path '13.multivariate_cox_ph_test.png', optional: true
    path '14.multivariate_cox_forest.pdf', optional: true
    path '14.multivariate_cox_forest.png', optional: true
    path '15.nomogram.pdf', optional: true
    path '15.nomogram.png', optional: true
    path '15.nomogram_regplot_observation_prediction.csv', optional: true
    path '16.calibration_curve.pdf', optional: true
    path '16.calibration_curve.png', optional: true
    path '17.timeROC_auc.csv', emit: auc_table
    path '17.timeROC.pdf', optional: true
    path '17.timeROC.png', optional: true
    path '18.dca_curve.pdf', optional: true
    path '18.dca_curve.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/clinical_nomogram.R \
      --surv ${survival_table} \
      --clinical ${clinical_table} \
      --risk ${risk_table} \
      --outdir . \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --time-roc-days ${params.time_roc_points_days} \
      --calibration-bootstrap ${params.nomogram_bootstrap}
    """
}

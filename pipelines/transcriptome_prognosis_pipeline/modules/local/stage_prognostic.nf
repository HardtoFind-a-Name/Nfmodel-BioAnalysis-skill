process PROGNOSTIC {
    tag 'prognostic'
    publishDir "${params.outdir}/${params.stage_SN}_prognostic", mode: 'copy', overwrite: true
    input:
    path risk_file
    path surv_file
    path clinical_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    path '*.txt', optional: true
    script:
    def var_flag = params.stage_variables ? "--variables ${params.stage_variables}" : ""
    """
    ${params.rscript} ${projectDir}/scripts/stage_prognostic.R \
      --risk-file ${risk_file} --surv ${surv_file} --clinical ${clinical_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.stage_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --time-set ${params.stage_time_set} --dca-year ${params.stage_dca_year} \
      --min-level-n ${params.stage_min_level_n} \
      ${var_flag}
    """
}

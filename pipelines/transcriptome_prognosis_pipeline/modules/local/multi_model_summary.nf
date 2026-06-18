process MULTI_MODEL_SUMMARY {
    tag 'multi_model_summary'

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/99_summary", mode: 'copy', overwrite: true

    input:
    val dummy

    output:
    path '00.multi_model_comparison.csv', emit: comparison
    path '01.pass_fail_report.csv', emit: pass_fail

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/multi_model_summary.R \
      --input-dir ${params.outdir}/${params.multi_model_SN}_multi_model \
      --outdir . \
      --auc-threshold ${params.model_auc_threshold} \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

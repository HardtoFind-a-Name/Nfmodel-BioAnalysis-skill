process UNICOX_SCREEN {
    tag 'unicox_screen'

    publishDir "${params.outdir}/${params.multi_model_SN}_multi_model/00_unicox_screen", mode: 'copy', overwrite: true

    input:
    path expr_matrix
    path survival_table
    path candidate_genes

    output:
    path '03.unicox_ph_pass_genes.csv', emit: unicox_genes
    path '01.unicox_all_results.csv', optional: true
    path '02.unicox_ph_test.csv', optional: true
    path '00.unicox_input_merged.csv', optional: true
    path '04.unicox_forest.pdf', optional: true
    path '04.unicox_forest.png', optional: true
    path '06.screening_summary.txt', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/multi_model/unicox_screen.R \
      --expr ${expr_matrix} \
      --surv ${survival_table} \
      --candidate-genes ${candidate_genes} \
      --outdir . \
      --unicox-p ${params.unicox_p_cutoff} \
      --ph-p ${params.ph_p_cutoff} \
      --sn ${params.multi_model_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

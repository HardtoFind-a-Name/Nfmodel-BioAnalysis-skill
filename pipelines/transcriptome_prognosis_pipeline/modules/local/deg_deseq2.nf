process DEG_DESEQ2 {
    tag 'deg_deseq2'

    publishDir "${params.outdir}/${params.deg_SN}_deg", mode: 'copy', overwrite: true

    input:
    path count_matrix
    path group_table
    val min_count
    val min_prop
    val vst_blind
    val logdir
    val rscript_cmd

    output:
    path '02.deg_all.csv', emit: deg_all
    path '01.tcga_luad_vst_expr.csv', emit: vst_expr

    script:
    """
    ${rscript_cmd} ${projectDir}/scripts/deg_deseq2.R \
      --count ${count_matrix} \
      --group ${group_table} \
      --outdir . \
      --min-count ${min_count} \
      --min-prop ${min_prop} \
      --vst-blind ${vst_blind} \
      --sn ${params.deg_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${logdir}
    """
}

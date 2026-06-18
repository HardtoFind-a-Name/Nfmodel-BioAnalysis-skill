process DEG_FILTER_PLOT {
    tag 'deg_filter_plot'

    publishDir "${params.outdir}/${params.deg_SN}_deg", mode: 'copy', overwrite: true

    input:
    path deg_all
    path vst_expr
    path group_table
    val p_cutoff
    val logfc_cutoff
    val use_pvalue
    val logdir
    val rscript_cmd

    output:
    path '03.deg_sig.csv', emit: deg_sig
    path '04.deg_volcano.png', optional: true
    path '04.deg_volcano.pdf', optional: true
    path '05.deg_heatmap.png', optional: true
    path '05.deg_heatmap.pdf', optional: true
    path '05.deg_heatmap.note.txt', optional: true

    script:
    """
    ${rscript_cmd} ${projectDir}/scripts/deg_filter_plot.R \
      --deg-all ${deg_all} \
      --vst ${vst_expr} \
      --group ${group_table} \
      --outdir . \
      --padj ${p_cutoff} \
      --logfc ${logfc_cutoff} \
      --use-fallback ${use_pvalue} \
      --sn ${params.deg_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${logdir}
    """
}

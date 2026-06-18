process SCRNA_KEY_CELL_SELECT {
    tag 'scrna_key_cell_select'
    publishDir "${params.outdir}/${params.keycell_SN}_scrna_key_cell_select", mode: 'copy', overwrite: true

    input:
    path proportions_file
    path enrichment_file
    path keygene_diff_file

    output:
    path 'key_celltypes.csv', emit: key_celltypes
    path '*.csv', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_key_cell_select.R \
      --proportions-file ${proportions_file} \
      --enrichment-file ${enrichment_file} \
      --keygene-diff-file ${keygene_diff_file} \
      --outdir . \
      --key-celltypes "${params.key_celltypes}" \
      --keycell-enrichment-fdr ${params.keycell_enrichment_fdr} \
      --keycell-min-abs-prop-diff ${params.keycell_min_abs_prop_diff} \
      --keycell-min-pct-expr ${params.keycell_min_pct_expr} \
      --keycell-gene-diff-fdr ${params.keycell_gene_diff_fdr} \
      --keycell-min-abs-log2fc ${params.keycell_min_abs_log2fc} \
      --keycell-min-support-genes ${params.keycell_min_support_genes} \
      --sn ${params.keycell_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

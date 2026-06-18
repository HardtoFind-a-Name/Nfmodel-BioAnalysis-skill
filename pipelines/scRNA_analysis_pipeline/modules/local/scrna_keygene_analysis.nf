process SCRNA_KEYGENE_ANALYSIS {
    tag 'scrna_keygene_analysis'
    publishDir "${params.outdir}/${params.keygene_SN}_scrna_keygene_analysis", mode: 'copy', overwrite: true

    input:
    path input_rds

    output:
    path '01_keygene_expr_dotplot_data.csv', emit: dotplot_data
    path '02_keygene_group_diff.csv', emit: group_diff
    path '01.keygene_annotated.rds', emit: keygene_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_keygene_analysis.R \
      --input-rds ${input_rds} \
      --outdir . \
      --target-gene-file ${params.target_gene_file} \
      --target-gene-col ${params.target_gene_col} \
      --sn ${params.keygene_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

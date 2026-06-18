process SCRNA_CELL_ENRICHMENT {
    tag 'scrna_cell_enrichment'
    publishDir "${params.outdir}/${params.enrichment_SN}_scrna_cell_enrichment", mode: 'copy', overwrite: true

    input:
    path input_rds

    output:
    path '01_celltype_counts.csv', emit: counts
    path '02_celltype_proportions.csv', emit: proportions
    path '03_celltype_enrichment_tests.csv', emit: enrichment_tests
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_cell_enrichment.R \
      --input-rds ${input_rds} \
      --outdir . \
      --group-col ${params.group_col} \
      --celltype-col ${params.celltype_col} \
      --sn ${params.enrichment_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

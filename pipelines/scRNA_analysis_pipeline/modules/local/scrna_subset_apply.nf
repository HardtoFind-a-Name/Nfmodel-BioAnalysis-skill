process SCRNA_SUBSET_APPLY {
    tag "${safe_label}_subset_apply"
    publishDir "${params.outdir}/${params.subset_apply_SN}_scrna_subset_apply/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds), path(mapping_file)

    output:
    tuple val(safe_label), val(key_celltype), path('01.subset_annotated.rds'), emit: subset_annotated_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_subset_apply.R \
      --input-rds ${input_rds} \
      --mapping-file ${mapping_file} \
      --outdir . \
      --target-gene-file ${params.target_gene_file} \
      --target-gene-col ${params.target_gene_col} \
      --sn ${params.subset_apply_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

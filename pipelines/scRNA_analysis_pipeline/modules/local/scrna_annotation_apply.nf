process SCRNA_ANNOTATION_APPLY {
    tag 'scrna_annotation_apply'
    publishDir "${params.outdir}/${params.anno_apply_SN}_scrna_annotation_apply", mode: 'copy', overwrite: true

    input:
    path input_rds
    path mapping_file

    output:
    path '01.seurat_annotated.rds', emit: annotated_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_annotation_apply.R \
      --input-rds ${input_rds} \
      --mapping-file ${mapping_file} \
      --outdir . \
      --target-gene-file ${params.target_gene_file} \
      --target-gene-col ${params.target_gene_col} \
      --annotation-marker-reference-file ${params.annotation_marker_reference_file} \
      --annotation-min-support-markers ${params.annotation_min_support_markers} \
      --annotation-allow-low-support ${params.annotation_allow_low_support} \
      --sn ${params.anno_apply_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

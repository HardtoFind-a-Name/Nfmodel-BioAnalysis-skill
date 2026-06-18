process SCRNA_ANNOTATION_PREPARE {
    tag 'scrna_annotation_prepare'
    publishDir "${params.outdir}/${params.anno_prepare_SN}_scrna_annotation_prepare", mode: 'copy', overwrite: true

    input:
    path input_rds

    output:
    path '01.annotation_prepare.rds', emit: prepare_rds
    path '02.annotation_prepare_results.rds', emit: prepare_results
    path '08_cluster_celltype_mapping_template.csv', emit: mapping_template
    path '09_literature_marker_reference_template.csv', emit: literature_marker_template
    path '10_annotation_agent_context.json', emit: agent_context
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_annotation_prepare.R \
      --input-rds ${input_rds} \
      --outdir . \
      --annotation-cluster-col ${params.annotation_cluster_col} \
      --annotation-resolution "${params.annotation_resolution}" \
      --target-gene-file ${params.target_gene_file} \
      --target-gene-col ${params.target_gene_col} \
      --marker-min-pct ${params.marker_min_pct} \
      --marker-logfc-threshold ${params.marker_logfc_threshold} \
      --scrna-cohort-id "${params.scRNA_cohort_id}" \
      --disease-name "${params.disease_name}" \
      --cancer-name "${params.cancer_name}" \
      --annotation-literature-years ${params.annotation_literature_years} \
      --sn ${params.anno_prepare_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

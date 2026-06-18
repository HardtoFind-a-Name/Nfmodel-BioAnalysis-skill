process SCRNA_SUBSET_ANNOTATION_PREPARE {
    tag "${safe_label}_subset_annotation_prepare"
    publishDir "${params.outdir}/${params.subset_anno_prepare_SN}_scrna_subset_annotation_prepare/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.subset_annotation_prepare.rds'), emit: prepare_rds
    path '*.csv', optional: true
    path '*.json', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_subset_annotation_prepare.R \
      --input-rds ${input_rds} \
      --outdir . \
      --key-celltype "${key_celltype}" \
      --safe-label "${safe_label}" \
      --scrna-cohort-id "${params.scRNA_cohort_id}" \
      --disease-name "${params.disease_name}" \
      --cancer-name "${params.cancer_name}" \
      --marker-min-pct ${params.marker_min_pct} \
      --marker-logfc-threshold ${params.marker_logfc_threshold} \
      --sn ${params.subset_anno_prepare_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

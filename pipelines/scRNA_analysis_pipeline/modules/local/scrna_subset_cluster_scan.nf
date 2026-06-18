process SCRNA_SUBSET_CLUSTER_SCAN {
    tag "${safe_label}_subset_cluster_scan"
    publishDir "${params.outdir}/${params.cluster_SN}_subset_cluster_scan/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.seurat_qc_reclustered_final.rds'), emit: clustered_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_cluster_scan.R \
      --input-rds ${input_rds} \
      --outdir . \
      --scan-label subset_${safe_label} \
      --preferred-reduction ${params.preferred_reduction} \
      --fallback-reduction ${params.fallback_reduction} \
      --cluster-scan-use-dims-override "${params.cluster_scan_use_dims_override}" \
      --resolution-grid ${params.subset_resolution_grid} \
      --final-resolution ${params.subset_final_resolution} \
      --cluster-algorithm ${params.cluster_algorithm} \
      --umap-n-neighbors ${params.umap_n_neighbors} \
      --umap-min-dist ${params.umap_min_dist} \
      --sample-col ${params.sample_col} \
      --group-col ${params.group_col} \
      --sn ${params.cluster_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

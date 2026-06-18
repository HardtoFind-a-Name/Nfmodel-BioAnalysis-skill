process SCRNA_CLUSTER_SCAN {
    tag "${scan_label}_cluster_scan"
    publishDir "${params.outdir}/${params.cluster_SN}_${scan_label}_cluster_scan", mode: 'copy', overwrite: true

    input:
    path input_rds
    val scan_label
    val resolution_grid
    val final_resolution

    output:
    path '01.seurat_qc_reclustered_final.rds', emit: clustered_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_cluster_scan.R \
      --input-rds ${input_rds} \
      --outdir . \
      --scan-label ${scan_label} \
      --preferred-reduction ${params.preferred_reduction} \
      --fallback-reduction ${params.fallback_reduction} \
      --cluster-scan-use-dims-override "${params.cluster_scan_use_dims_override}" \
      --resolution-grid ${resolution_grid} \
      --final-resolution ${final_resolution} \
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

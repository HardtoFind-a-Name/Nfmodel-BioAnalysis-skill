process SCRNA_QC_INTEGRATION {
    tag 'scrna_qc_integration'
    publishDir "${params.outdir}/${params.qc_SN}_scrna_qc_integration", mode: 'copy', overwrite: true

    input:
    path input_rds

    output:
    path '01.seurat_precluster.rds', emit: precluster_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_qc_integration.R \
      --input-rds ${input_rds} \
      --outdir . \
      --sample-col ${params.sample_col} \
      --group-col ${params.group_col} \
      --min-features-keep ${params.min_features_keep} \
      --max-features-keep ${params.max_features_keep} \
      --min-counts-keep ${params.min_counts_keep} \
      --max-counts-keep ${params.max_counts_keep} \
      --max-percent-mt ${params.max_percent_mt} \
      --run-doublet-removal ${params.run_doublet_removal} \
      --nfeatures-hvg ${params.nfeatures_hvg} \
      --npcs ${params.npcs} \
      --use-dims ${params.use_dims} \
      --harmony-batch-var ${params.harmony_batch_var} \
      --seed ${params.seed} \
      --sn ${params.qc_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

process SCRNA_PREPARE_INPUT {
    tag 'scrna_prepare_input'
    publishDir "${params.outdir}/${params.input_SN}_scrna_prepare_input", mode: 'copy', overwrite: true

    input:
    path input_rds

    output:
    path '01.input_seurat.rds', emit: seurat_rds
    path '*.csv', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_prepare_input.R \
      --input-rds ${input_rds} \
      --input-format ${params.input_format} \
      --outdir . \
      --project-name ${params.project_id} \
      --sample-col ${params.sample_col} \
      --origin-col ${params.origin_col} \
      --group-col ${params.group_col} \
      --normal-origin ${params.normal_origin_value} \
      --tumor-origin ${params.tumor_origin_value} \
      --min-cells-per-gene ${params.min_cells_per_gene} \
      --min-features-create ${params.min_features_create} \
      --sn ${params.input_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

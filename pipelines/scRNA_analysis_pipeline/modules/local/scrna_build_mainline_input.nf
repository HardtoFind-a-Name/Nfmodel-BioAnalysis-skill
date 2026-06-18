process SCRNA_BUILD_MAINLINE_INPUT {
    tag "${params.scRNA_cohort_id}"
    publishDir "${params.outdir}/${params.raw_build_SN}_scrna_build_mainline_input", mode: 'copy', overwrite: true

    input:
    path raw_source_ready

    output:
    path '01.mainline_input.rds', emit: mainline_input
    path '*.csv', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_build_mainline_input.R \
      --raw-source-manifest ${raw_source_ready} \
      --outdir . \
      --keep-origin "${params.keep_origin}" \
      --exclude-origin "${params.exclude_origin}" \
      --normal-origins "${params.normal_origins}" \
      --tumor-origins "${params.tumor_origins}" \
      --origin-col ${params.origin_col} \
      --sample-col ${params.sample_col} \
      --annotation-cell-id-candidates "${params.annotation_cell_id_candidates}" \
      --strip-10x-barcode-suffix ${params.strip_10x_barcode_suffix} \
      --tenx-assay "${params.tenx_assay}" \
      --sn ${params.raw_build_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

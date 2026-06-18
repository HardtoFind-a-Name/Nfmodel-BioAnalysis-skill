process SCRNA_UNPACK_RAW_UMI {
    tag "${params.scRNA_cohort_id}"
    publishDir "${params.outdir}/${params.raw_unpack_SN}_scrna_unpack_raw_umi", mode: 'copy', overwrite: true

    input:
    path raw_source_manifest

    output:
    path '01.raw_source_after_umi.csv', emit: raw_source_after_umi
    path '*.csv', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_unpack_raw_umi.R \
      --raw-source-manifest ${raw_source_manifest} \
      --outdir . \
      --decompress-max-depth ${params.decompress_max_depth} \
      --sn ${params.raw_unpack_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

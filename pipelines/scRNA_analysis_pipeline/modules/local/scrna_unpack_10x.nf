process SCRNA_UNPACK_10X {
    tag "${params.scRNA_cohort_id}"
    publishDir "${params.outdir}/${params.raw_10x_SN}_scrna_unpack_10x", mode: 'copy', overwrite: true

    input:
    path raw_source_manifest

    output:
    path '01.raw_source_ready.csv', emit: raw_source_ready
    path '*.csv', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_unpack_10x.R \
      --raw-source-manifest ${raw_source_manifest} \
      --outdir . \
      --unpack-root "${params.scRNA_raw_dir}/unpacked_10x" \
      --raw-10x-prefer "${params.raw_10x_prefer}" \
      --sn ${params.raw_10x_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

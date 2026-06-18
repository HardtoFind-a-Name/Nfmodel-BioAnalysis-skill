process SCRNA_DETECT_RAW_SOURCE {
    tag "${params.scRNA_cohort_id}"
    publishDir "${params.outdir}/${params.raw_detect_SN}_scrna_detect_raw_source", mode: 'copy', overwrite: true

    input:
    path resolved_files

    output:
    path '01.raw_source_manifest.csv', emit: raw_source_manifest
    path '*.csv', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_detect_raw_source.R \
      --resolved-file ${resolved_files} \
      --outdir . \
      --raw-source-type "${params.raw_source_type}" \
      --raw-10x-pattern "${params.raw_10x_pattern}" \
      --raw-10x-prefer "${params.raw_10x_prefer}" \
      --sn ${params.raw_detect_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

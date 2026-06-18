process SCRNA_RESOLVE_RAW_FILES {
    tag "${params.scRNA_cohort_id}"
    publishDir "${params.outdir}/${params.raw_resolve_SN}_scrna_resolve_raw_files", mode: 'copy', overwrite: true

    input:
    path manifest

    output:
    path '01.resolved_raw_files.csv', emit: resolved_files

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_resolve_raw_files.R \
      --manifest ${manifest} \
      --raw-dir ${params.scRNA_raw_dir} \
      --raw-umi-file "${params.scRNA_raw_umi_file}" \
      --annotation-file "${params.scRNA_annotation_file}" \
      --raw-umi-pattern "${params.raw_umi_pattern}" \
      --annotation-pattern "${params.annotation_pattern}" \
      --outdir . \
      --sn ${params.raw_resolve_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

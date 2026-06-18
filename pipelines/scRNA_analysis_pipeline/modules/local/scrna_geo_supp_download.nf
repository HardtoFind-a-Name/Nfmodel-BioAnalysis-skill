process SCRNA_GEO_SUPP_DOWNLOAD {
    tag "${params.scRNA_cohort_id}"
    publishDir "${params.outdir}/${params.raw_download_SN}_scrna_geo_supp_download", mode: 'copy', overwrite: true

    output:
    path '01.geo_supp_files.csv', emit: manifest

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_geo_supp_download.R \
      --scrna-cohort-id ${params.scRNA_cohort_id} \
      --database-dir ${params.scRNA_database_dir} \
      --raw-dir ${params.scRNA_raw_dir} \
      --geo-filter-regex "${params.geo_filter_regex}" \
      --outdir . \
      --sn ${params.raw_download_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

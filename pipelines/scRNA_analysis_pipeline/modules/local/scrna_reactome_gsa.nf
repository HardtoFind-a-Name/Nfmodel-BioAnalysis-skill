process SCRNA_REACTOME_GSA {
    tag "${safe_label}_reactome_gsa"
    publishDir "${params.outdir}/${params.reactome_gsa_SN}_scrna_reactome_gsa/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.reactomegsa_results.rds'), emit: reactome_results
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_reactome_gsa.R \
      --input-rds ${input_rds} \
      --outdir . \
      --key-celltype "${key_celltype}" \
      --max-pathways ${params.reactome_max_pathways} \
      --p-cutoff ${params.reactome_p_cutoff} \
      --sn ${params.reactome_gsa_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

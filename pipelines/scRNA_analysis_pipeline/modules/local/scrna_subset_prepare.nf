process SCRNA_SUBSET_PREPARE {
    tag "${safe_label}_subset_prepare"
    publishDir "${params.outdir}/${params.subset_prepare_SN}_scrna_subset_prepare/${safe_label}", mode: 'copy', overwrite: true

    input:
    tuple val(safe_label), val(key_celltype), path(input_rds)

    output:
    tuple val(safe_label), val(key_celltype), path('01.subset_precluster.rds'), emit: subset_precluster_rds
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_subset_prepare.R \
      --input-rds ${input_rds} \
      --outdir . \
      --subset-celltypes "${key_celltype}" \
      --subset-npcs ${params.subset_npcs} \
      --subset-use-dims "${params.subset_use_dims}" \
      --subset-nfeatures-hvg ${params.subset_nfeatures_hvg} \
      --subset-run-jackstraw ${params.subset_run_jackstraw} \
      --seed ${params.seed} \
      --sn ${params.subset_prepare_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

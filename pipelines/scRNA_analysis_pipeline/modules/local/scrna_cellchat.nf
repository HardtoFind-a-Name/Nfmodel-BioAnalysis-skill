process SCRNA_CELLCHAT {
    tag 'scrna_cellchat'
    publishDir "${params.outdir}/${params.cellchat_SN}_scrna_cellchat", mode: 'copy', overwrite: true

    input:
    path input_rds

    output:
    path '01_cellchat_group_status.csv', emit: status
    path '*', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/scrna_cellchat.R \
      --input-rds ${input_rds} \
      --outdir . \
      --group-col ${params.group_col} \
      --celltype-col ${params.celltype_col} \
      --species ${params.cellchat_species} \
      --min-cells ${params.cellchat_min_cells} \
      --sn ${params.cellchat_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

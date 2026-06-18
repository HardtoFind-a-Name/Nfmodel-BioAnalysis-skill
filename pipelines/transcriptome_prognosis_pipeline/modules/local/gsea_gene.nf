process GSEA_GENE {
    tag 'gsea_gene'
    publishDir "${params.outdir}/${params.gsea_gene_SN}_gsea_gene", mode: 'copy', overwrite: true
    input:
    path risk_file
    path expr_file
    path gene_file
    val train_id
    output:
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    script:
    """
    ${params.rscript} ${projectDir}/scripts/gsea_gene.R \
      --risk-file ${risk_file} --expr ${expr_file} \
      --train-id ${train_id} --outdir . \
      --sn ${params.gsea_gene_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --gene-file ${gene_file} \
      --gmt-dir ${params.gmt_dir} --gmt-gsea_gene ${params.gmt_gsea_gene} \
      --gsea-sort-by ${params.gsea_sort_by}
    """
}

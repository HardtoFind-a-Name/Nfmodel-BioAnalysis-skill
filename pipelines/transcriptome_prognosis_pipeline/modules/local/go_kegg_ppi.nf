process GO_KEGG_PPI {
    tag 'go_kegg_ppi'

    publishDir "${params.outdir}/${params.go_kegg_SN}_go_kegg_ppi", mode: 'copy', overwrite: true

    input:
    path candidate_genes
    path string_interactions_input

    output:
    path '01_gene_id_conversion.csv', optional: true
    path '02_GO_all_results.csv', emit: go_all
    path '03_GO_significant_results.csv', optional: true
    path '04.GO_barplot.pdf', optional: true
    path '04.GO_barplot.png', optional: true
    path '06_KEGG_all_results.csv', emit: kegg_all
    path '07_KEGG_significant_results.csv', optional: true
    path '08.KEGG_barplot.pdf', optional: true
    path '08.KEGG_barplot.png', optional: true
    path '10_PPI_node_metrics.csv', emit: ppi_metrics
    path '11_PPI_edges_used_for_plot.csv', optional: true
    path '12.PPI_network.pdf', optional: true
    path '12.PPI_network.png', optional: true
    path '13_PPI_network_summary.csv', optional: true
    path '00.string_interactions_input.tsv', optional: true

    script:
    """
    cp ${string_interactions_input} 00.string_interactions_input.tsv
    ${params.rscript} ${projectDir}/scripts/go_kegg_ppi.R \
      --candidate-genes ${candidate_genes} \
      --string-tsv ${string_interactions_input} \
      --outdir . \
      --sn ${params.go_kegg_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir} \
      --go-top-n-each ${params.go_top_n_each} \
      --kegg-top-n ${params.kegg_top_n} \
      --adj-cutoff ${params.go_kegg_adj_cutoff} \
      --fallback-p ${params.go_kegg_fallback_p} \
      --kegg-gmt ${params.kegg_gmt} \
      --kegg-gmt-id-type ${params.kegg_gmt_id_type} \
      --ppi-score-cutoff ${params.ppi_score_cutoff} \
      --ppi-hub-top-n ${params.ppi_hub_top_n}
    """
}

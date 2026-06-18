process GSEA_CANDIDATE_GENES {
    tag 'gsea_candidate_genes'

    publishDir "${params.outdir}/${params.candidate_SN}_candidate_genes", mode: 'copy', overwrite: true

    input:
    path deg_all
    path deg_sig
    path metascape_input

    output:
    path '10.NRGs_from_GSEA_leading_edge_and_Metascape.csv', emit: candidate_genes
    path '01.gsea_go_bp_all.csv', emit: gsea_all
    path '02.gsea_go_bp_sig.csv', optional: true
    path '03.gsea_go_bp_sig_dotplot.pdf', optional: true
    path '04.gsea_go_bp_sig_dotplot.png', optional: true
    path '05.metascape_neuroimmune_selected_terms.csv', optional: true
    path '06.metascape_neuroimmune_gene_list.csv', optional: true
    path '07.metascape_neuroimmune_barplot.pdf', optional: true
    path '08.metascape_neuroimmune_barplot.png', optional: true
    path '09.gsea_leading_edge_gene_set.csv', optional: true
    path '11.NRGs_venn.pdf', optional: true
    path '12.NRGs_venn.png', optional: true
    path '00.metascape_input.csv', optional: true

    script:
    """
    cp ${metascape_input} 00.metascape_input.csv
    ${params.rscript} ${projectDir}/scripts/gsea_candidate_genes.R \
      --deg-all ${deg_all} \
      --deg-sig ${deg_sig} \
      --metascape ${metascape_input} \
      --outdir . \
      --sn ${params.candidate_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

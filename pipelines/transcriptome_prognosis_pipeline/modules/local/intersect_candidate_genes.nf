process INTERSECT_CANDIDATE_GENES {
    tag 'intersect_candidate_genes'

    publishDir "${params.outdir}/${params.candidate_SN}_candidate_genes", mode: 'copy', overwrite: true

    input:
    path deg_sig
    path gene_sets_sheet

    output:
    path '01.candidate_genes.csv', emit: candidate_genes
    path '02.candidate_venn.pdf', optional: true
    path '02.candidate_venn.png', optional: true

    script:
    """
    ${params.rscript} ${projectDir}/scripts/intersect_candidate_genes.R \
      --deg-sig ${deg_sig} \
      --gene-sets-sheet ${gene_sets_sheet} \
      --outdir . \
      --sn ${params.candidate_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

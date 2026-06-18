process TIDE_ANALYSIS {
    tag 'tide_analysis'
    publishDir "${params.outdir}/${params.tide_SN}_tide", mode: 'copy', overwrite: true
    input:
    path risk_file
    path expr_file
    val train_id
    output:
    path '00.*_tide_scores', optional: true
    path '*.csv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true
    script:
    def tid = train_id.toLowerCase().replace('-', '_')
    """
    # Step 1: tidepy (with cache)
    mkdir -p ${params.tidepy_cache_dir}
    if [ ! -f ${params.tidepy_cache_dir}/03.${tid}_tide_scores.tsv ]; then
        python3 ${projectDir}/scripts/tidepy_run.py \
          --expr ${expr_file} \
          --outdir ${params.tidepy_cache_dir} \
          --train-id ${train_id} \
          --tidepy ${params.tidepy_bin}
    fi
    cp ${params.tidepy_cache_dir}/03.${tid}_tide_scores.tsv ./00.${tid}_tide_scores

    # Step 2: R analysis
    ${params.rscript} ${projectDir}/scripts/tide_analysis.R \
      --risk-file ${risk_file} --tide-scores ./00.${tid}_tide_scores \
      --train-id ${train_id} --outdir . \
      --sn ${params.tide_SN} \
      --summary-dir ${params.summary_dir} \
      --logdir ${params.logdir}
    """
}

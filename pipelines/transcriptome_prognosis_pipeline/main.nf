nextflow.enable.dsl = 2

include { DEG_DESEQ2 }                from './modules/local/deg_deseq2'
include { DEG_FILTER_PLOT }           from './modules/local/deg_filter_plot'
include { GSEA_CANDIDATE_GENES }      from './modules/local/gsea_candidate_genes'
include { INTERSECT_CANDIDATE_GENES } from './modules/local/intersect_candidate_genes'
include { GO_KEGG_PPI }               from './modules/local/go_kegg_ppi'
include { CLINICAL_NOMOGRAM }         from './modules/local/clinical_nomogram'
include { GEN_VALIDATION_SHEET }      from './modules/local/generate_validation_sheet'

include { UNICOX_SCREEN }                          from './modules/local/unicox_screen'
include { LASSO_SELECT }                           from './modules/local/lasso_select'
include { MULTICOX_SELECT as MULTICOX_A2A3 }       from './modules/local/multicox_select'
include { MULTICOX_SELECT as MULTICOX_A4A5 }       from './modules/local/multicox_select'
include { STEPCOX_SELECT  as STEPCOX_A3 }          from './modules/local/stepcox_select'
include { STEPCOX_SELECT  as STEPCOX_A5 }          from './modules/local/stepcox_select'
include { COXBOOST_SELECT }                        from './modules/local/coxboost_select'
include { COEF_MODEL_TRAIN as TRAIN_A1 }           from './modules/local/coef_model_train'
include { COEF_MODEL_TRAIN as TRAIN_A2 }           from './modules/local/coef_model_train'
include { COEF_MODEL_TRAIN as TRAIN_A3 }           from './modules/local/coef_model_train'
include { COEF_MODEL_TRAIN as TRAIN_A4 }           from './modules/local/coef_model_train'
include { COEF_MODEL_TRAIN as TRAIN_A5 }           from './modules/local/coef_model_train'
include { COEF_MODEL_TRAIN as TRAIN_A6 }           from './modules/local/coef_model_train'
include { COEF_MODEL_VALIDATE as VALIDATE_A1 }     from './modules/local/coef_model_validate'
include { COEF_MODEL_VALIDATE as VALIDATE_A2 }     from './modules/local/coef_model_validate'
include { COEF_MODEL_VALIDATE as VALIDATE_A3 }     from './modules/local/coef_model_validate'
include { COEF_MODEL_VALIDATE as VALIDATE_A4 }     from './modules/local/coef_model_validate'
include { COEF_MODEL_VALIDATE as VALIDATE_A5 }     from './modules/local/coef_model_validate'
include { COEF_MODEL_VALIDATE as VALIDATE_A6 }     from './modules/local/coef_model_validate'
include { RSF_TRAIN_VALIDATE as RSF_B1 }           from './modules/local/rsf_train_validate'
include { RSF_TRAIN_VALIDATE as RSF_B2 }           from './modules/local/rsf_train_validate'
include { RSF_B3_VALIDATE_GUIDED }                 from './modules/local/rsf_b3_validate_guided'
include { COXBOOST_A7_VALIDATE_GUIDED }            from './modules/local/coxboost_a7_validate_guided'
include { MULTI_MODEL_SUMMARY }                    from './modules/local/multi_model_summary'
include { PROGNOSTIC }   from './modules/local/stage_prognostic'
include { GSEA_GENE }          from './modules/local/gsea_gene'
include { GSEA_RISK_NES }      from './modules/local/gsea_risk_nes'
include { CIBERSORT_IMMUNE }   from './modules/local/cibersort_immune'
include { CNV_ANALYSIS }       from './modules/local/cnv_analysis'
include { GISTIC_INPUT }       from './modules/local/gistic_input'
include { TIDE_ANALYSIS }      from './modules/local/tide_analysis'
include { IPS_ANALYSIS }       from './modules/local/ips_analysis'
include { DIFF_EXPR }          from './modules/local/diff_expr'
include { CIRCOS_PLOT }        from './modules/local/circos_plot'
include { TMB_ONCOPLOT }       from './modules/local/tmb_oncoplot'

workflow {
    // --- Derived paths (dynamic with CLI overrides) ---
    train_id_lower  = params.train_id.toLowerCase().replace('-', '_')
    count_matrix    = "${params.rawdata_dir}/01.${train_id_lower}_count.csv"
    group_table     = "${params.rawdata_dir}/02.${train_id_lower}_group.csv"
    tpm_matrix      = "${params.rawdata_dir}/03.${train_id_lower}_tpmlog2.csv"
    fpkm_matrix     = "${params.rawdata_dir}/04.${train_id_lower}_fpkmlog2.csv"
    survival_table  = "${params.rawdata_dir}/05.${train_id_lower}_survival.csv"
    clinical_table  = "${params.rawdata_dir}/06.${train_id_lower}_clinical.csv"
    expr_matrix     = params.expr_type == 'tpm' ? tpm_matrix : fpkm_matrix

    clean_root       = "${params.user_dir}/data_cleaned"
    training_dir     = "${clean_root}/Train_set/${params.train_id}"

    validation_ids_list = params.validation_ids ? params.validation_ids.split(',').collect { it.trim() } : []
    validation_dirs_str = validation_ids_list.collect { "${clean_root}/Validation_set/${it}" }.join(' ')
    lasso_model_flags = [
        params.run_model_A1,
        params.run_model_A4,
        params.run_model_A5,
        params.run_model_A6,
        params.run_model_A7,
        params.run_model_B2,
        params.run_model_B3
    ]
    any_lasso_model_enabled = lasso_model_flags.any { it }

    needs_candidate_genes = (
        params.run_go_kegg_ppi ||
        params.run_multi_model
    )

    any_candidate_method = (
        params.run_candidate_genes ||
        params.run_intersect_candidate_genes
    )

    if (!params.project_id?.toString()?.trim()) {
        error "`params.project_id` is required."
    }

    if (!params.rscript?.toString()?.trim()) {
        error "Set `NFMODELS_RSCRIPT_CMD` in ../.env, then launch with ./run_pipeline.sh"
    }

    if (!params.python_cmd?.toString()?.trim()) {
        error "Set `NFMODELS_PYTHON_CMD` in ../.env, then launch with ./run_pipeline.sh"
    }

    if (params.run_tide && !params.tidepy_bin?.toString()?.trim()) {
        error "Set `NFMODELS_TIDEPY_BIN` in ../.env or disable --run_tide"
    }

    if (needs_candidate_genes && !any_candidate_method && !params.candidate_genes) {
        error "Provide `params.candidate_genes` or enable `--run_candidate_genes` / `--run_intersect_candidate_genes`."
    }

    if (params.run_candidate_genes && !params.candidate_genes && !params.run_deg) {
        error "`GSEA_CANDIDATE_GENES` requires `--run_deg true` when `params.candidate_genes` is not provided."
    }

    if (params.run_candidate_genes && !params.candidate_genes && !params.metascape_input_file) {
        error "`GSEA_CANDIDATE_GENES` requires `params.metascape_input_file` when `params.candidate_genes` is not provided."
    }

    if (params.run_intersect_candidate_genes && !params.candidate_genes && !params.run_deg) {
        error "`INTERSECT_CANDIDATE_GENES` requires `--run_deg true` when `params.candidate_genes` is not provided."
    }

    if (params.run_intersect_candidate_genes && !params.candidate_genes && !params.gene_sets_sheet) {
        error "`INTERSECT_CANDIDATE_GENES` requires `params.gene_sets_sheet`."
    }

    if (params.run_go_kegg_ppi && !params.string_interactions_file) {
        error "`GO_KEGG_PPI` requires `params.string_interactions_file` or set `--run_go_kegg_ppi false`."
    }

    if (params.run_multi_model && !params.validation_ids && !params.validation_sheet) {
        error "`run_multi_model` requires `--validation_ids` or `--validation_sheet`."
    }

    if (params.run_multi_model && !params.run_lasso_select && any_lasso_model_enabled) {
        error "`LASSO_SELECT` is required for models A1/A4/A5/A6/A7/B2/B3. Set `--run_lasso_select true` or disable those model toggles."
    }

    if (params.run_clinical_nomogram && !params.clinical_table) {
        error "`CLINICAL_NOMOGRAM` requires `params.clinical_table` or set `--run_clinical_nomogram false`."
    }

    // --- Data copy (runs at parse time, before any process) ---
    ["mkdir", "-p", params.rawdata_dir].execute().waitFor()
    ["mkdir", "-p", params.summary_dir].execute().waitFor()
    def trainSrc = "${clean_root}/Train_set/${params.train_id}"
    ["bash", "-c", "cp ${trainSrc}/*.csv ${params.rawdata_dir}/ 2>/dev/null || true"].execute().waitFor()
    if (params.validation_ids) {
        params.validation_ids.split(',').each { vid ->
            def v = vid.trim()
            if (v) {
                ["bash", "-c", "mkdir -p ${params.rawdata_dir}/${v} && cp ${clean_root}/Validation_set/${v}/*.csv ${params.rawdata_dir}/${v}/ 2>/dev/null || true"].execute().waitFor()
            }
        }
    }

    count_ch = params.run_deg ?
        Channel.fromPath(count_matrix, checkIfExists: true) : null
    group_ch = params.run_deg ?
        Channel.fromPath(group_table, checkIfExists: true) : null
    surv_ch = params.run_multi_model ?
        Channel.fromPath(survival_table, checkIfExists: true) : null
    expr_ch = params.run_multi_model ?
        Channel.fromPath(expr_matrix, checkIfExists: true) : null
    candidate_ch = params.candidate_genes ?
        Channel.fromPath(params.candidate_genes, checkIfExists: true) : null
    clinical_ch = params.run_clinical_nomogram ?
        Channel.fromPath(clinical_table, checkIfExists: true) : null

    // Generate validation_sheet if validation_ids provided
    if (params.validation_ids && validation_ids_list.size() > 0) {
        GEN_VALIDATION_SHEET(params.validation_ids, params.train_id)
        validation_sheet_ch = GEN_VALIDATION_SHEET.out.validation_sheet
    } else {
        validation_sheet_ch = Channel.fromPath(params.validation_sheet, checkIfExists: true)
    }

    // --- DEG ---
    if (params.run_deg) {
        DEG_DESEQ2(count_ch, group_ch, params.min_count, params.min_prop, true, params.logdir, params.rscript)
        DEG_FILTER_PLOT(DEG_DESEQ2.out.deg_all, DEG_DESEQ2.out.vst_expr, group_ch, params.deg_p_cutoff, params.deg_logfc_cutoff, params.deg_use_pvalue, params.logdir, params.rscript)
    }

    // --- Candidate genes ---
    if (params.run_candidate_genes && candidate_ch == null) {
        metascape_ch = Channel.fromPath(params.metascape_input_file, checkIfExists: true)
        GSEA_CANDIDATE_GENES(DEG_DESEQ2.out.deg_all, DEG_FILTER_PLOT.out.deg_sig, metascape_ch)
        candidate_ch = GSEA_CANDIDATE_GENES.out.candidate_genes
    }

    if (params.run_intersect_candidate_genes && candidate_ch == null) {
        gene_sets_sheet_ch = Channel.fromPath(params.gene_sets_sheet, checkIfExists: true)
        INTERSECT_CANDIDATE_GENES(DEG_FILTER_PLOT.out.deg_sig, gene_sets_sheet_ch)
        candidate_ch = INTERSECT_CANDIDATE_GENES.out.candidate_genes
    }

    // --- GO/KEGG/PPI ---
    if (params.run_go_kegg_ppi) {
        string_ch = Channel.fromPath(params.string_interactions_file, checkIfExists: true)
        GO_KEGG_PPI(candidate_ch, string_ch)
    }

    // --- Multi-model risk pipeline ---
    if (params.run_multi_model) {
        UNICOX_SCREEN(expr_ch, surv_ch, candidate_ch)
        if (params.run_lasso_select) {
            LASSO_SELECT(expr_ch, surv_ch, UNICOX_SCREEN.out.unicox_genes)
        }

        validation_sets_ch = validation_sheet_ch
            .splitCsv(header: true)
            .map { row ->
                def td = (row.time_roc_days as String)?.trim()?.replace(';', ',')
                def roc_days = (td != null && td != "") ? td : params.time_roc_points_days
                tuple(row.cohort_id as String,
                      file(row.expr_file as String),
                      file(row.surv_file as String),
                      roc_days)
            }

        // Collect summaries from all models
        train_summaries_ch = Channel.empty()
        valid_summaries_ch = Channel.empty()

        // A1: unicox -> LASSO coef
        if (params.run_model_A1) {
            TRAIN_A1(expr_ch, surv_ch, LASSO_SELECT.out.lasso_coef, 'A1_lasso_coef')
            VALIDATE_A1(validation_sets_ch, TRAIN_A1.out.coef_file_out.first(), 'A1_lasso_coef')
            train_summaries_ch = train_summaries_ch.mix(TRAIN_A1.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(VALIDATE_A1.out.validation_summary)
        }

        // A2: unicox -> multiCox coef
        if (params.run_model_A2 || params.run_model_A3) {
            MULTICOX_A2A3(expr_ch, surv_ch, UNICOX_SCREEN.out.unicox_genes, 'A2A3')
        }
        if (params.run_model_A2) {
            TRAIN_A2(expr_ch, surv_ch, MULTICOX_A2A3.out.multicox_coef, 'A2_multiCox_coef')
            VALIDATE_A2(validation_sets_ch, TRAIN_A2.out.coef_file_out.first(), 'A2_multiCox_coef')
            train_summaries_ch = train_summaries_ch.mix(TRAIN_A2.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(VALIDATE_A2.out.validation_summary)
        }

        // A3: unicox -> multiCox -> stepCox coef
        if (params.run_model_A3) {
            STEPCOX_A3(expr_ch, surv_ch, MULTICOX_A2A3.out.multicox_genes, 'A3_stepCox')
            TRAIN_A3(expr_ch, surv_ch, STEPCOX_A3.out.stepcox_coef, 'A3_stepCox_coef')
            VALIDATE_A3(validation_sets_ch, TRAIN_A3.out.coef_file_out.first(), 'A3_stepCox_coef')
            train_summaries_ch = train_summaries_ch.mix(TRAIN_A3.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(VALIDATE_A3.out.validation_summary)
        }

        // A4: unicox -> LASSO -> multiCox coef
        if (params.run_model_A4 || params.run_model_A5) {
            MULTICOX_A4A5(expr_ch, surv_ch, LASSO_SELECT.out.lasso_genes, 'A4A5')
        }
        if (params.run_model_A4) {
            TRAIN_A4(expr_ch, surv_ch, MULTICOX_A4A5.out.multicox_coef, 'A4_lasso_multiCox_coef')
            VALIDATE_A4(validation_sets_ch, TRAIN_A4.out.coef_file_out.first(), 'A4_lasso_multiCox_coef')
            train_summaries_ch = train_summaries_ch.mix(TRAIN_A4.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(VALIDATE_A4.out.validation_summary)
        }

        // A5: unicox -> LASSO -> multiCox -> stepCox coef
        if (params.run_model_A5) {
            STEPCOX_A5(expr_ch, surv_ch, MULTICOX_A4A5.out.multicox_genes, 'A5_stepCox')
            TRAIN_A5(expr_ch, surv_ch, STEPCOX_A5.out.stepcox_coef, 'A5_lasso_multiCox_stepCox')
            VALIDATE_A5(validation_sets_ch, TRAIN_A5.out.coef_file_out.first(), 'A5_lasso_multiCox_stepCox')
            train_summaries_ch = train_summaries_ch.mix(TRAIN_A5.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(VALIDATE_A5.out.validation_summary)
        }

        // A6: unicox -> LASSO -> CoxBoost coef (experimental)
        if (params.run_model_A6) {
            COXBOOST_SELECT(expr_ch, surv_ch, LASSO_SELECT.out.lasso_genes)
            TRAIN_A6(expr_ch, surv_ch, COXBOOST_SELECT.out.coxboost_coef, 'A6_coxboost_coef')
            VALIDATE_A6(validation_sets_ch, TRAIN_A6.out.coef_file_out.first(), 'A6_coxboost_coef')
            train_summaries_ch = train_summaries_ch.mix(TRAIN_A6.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(VALIDATE_A6.out.validation_summary)
        }

        // B1: unicox -> RSF
        if (params.run_model_B1) {
            RSF_B1(expr_ch, surv_ch, UNICOX_SCREEN.out.unicox_genes, validation_sheet_ch, 'B1_unicox_rsf')
            train_summaries_ch = train_summaries_ch.mix(RSF_B1.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(RSF_B1.out.validation_summary)
        }

        // B2: unicox -> LASSO -> RSF
        if (params.run_model_B2) {
            RSF_B2(expr_ch, surv_ch, LASSO_SELECT.out.lasso_genes, validation_sheet_ch, 'B2_lasso_rsf')
            train_summaries_ch = train_summaries_ch.mix(RSF_B2.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(RSF_B2.out.validation_summary)
        }

        // B3: unicox -> LASSO -> RSF (validation-guided grid search, per-cohort)
        if (params.run_model_B3) {
            b3_input_ch = expr_ch
                .combine(surv_ch)
                .combine(LASSO_SELECT.out.lasso_genes)
                .combine(validation_sets_ch)
                .map { t_expr, t_surv, genes, cid, v_expr, v_surv, td ->
                    tuple(cid, td, v_expr, v_surv, t_expr, t_surv, genes)
                }
            RSF_B3_VALIDATE_GUIDED(b3_input_ch, 'B3_rsf_vg')
            train_summaries_ch = train_summaries_ch.mix(RSF_B3_VALIDATE_GUIDED.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(RSF_B3_VALIDATE_GUIDED.out.validation_summary)
        }

        // A7: unicox -> LASSO -> CoxBoost (validation-guided grid search, per-cohort)
        if (params.run_model_A7) {
            a7_input_ch = expr_ch
                .combine(surv_ch)
                .combine(LASSO_SELECT.out.lasso_genes)
                .combine(validation_sets_ch)
                .map { t_expr, t_surv, genes, cid, v_expr, v_surv, td ->
                    tuple(cid, td, v_expr, v_surv, t_expr, t_surv, genes)
                }
            COXBOOST_A7_VALIDATE_GUIDED(a7_input_ch, 'A7_coxboost_vg')
            train_summaries_ch = train_summaries_ch.mix(COXBOOST_A7_VALIDATE_GUIDED.out.train_summary)
            valid_summaries_ch = valid_summaries_ch.mix(COXBOOST_A7_VALIDATE_GUIDED.out.validation_summary)
        }

        // Summary — wait for all train summaries, then read from publishDir
        MULTI_MODEL_SUMMARY(train_summaries_ch.collect())
    }

    // --- Post-modeling analyses (require --risk_file) ---
    if (params.risk_file) {
        risk_ch = Channel.fromPath(params.risk_file, checkIfExists: true)
        gene_for_post = params.gene_file ? Channel.fromPath(params.gene_file, checkIfExists: true) : null
        expr_for_post = params.run_multi_model ? expr_ch : Channel.fromPath(fpkm_matrix, checkIfExists: true)
        count_for_post = params.run_deg ? count_ch : Channel.fromPath(count_matrix, checkIfExists: true)
        surv_for_post = params.run_multi_model ? surv_ch : Channel.fromPath(survival_table, checkIfExists: true)
        clinical_for_post = Channel.fromPath(clinical_table, checkIfExists: true)
        group_for_post = params.run_deg ? group_ch : Channel.fromPath(group_table, checkIfExists: true)

        if (params.run_stage)       { PROGNOSTIC(risk_ch, surv_for_post, clinical_for_post, params.train_id) }
        if (params.run_diff_expr)   { DIFF_EXPR(gene_for_post, expr_for_post, group_for_post, params.train_id) }
        if (params.run_gsea_gene)   { GSEA_GENE(risk_ch, expr_for_post, gene_for_post, params.train_id) }
        if (params.run_gsea_risk_nes) { GSEA_RISK_NES(risk_ch, count_for_post, params.train_id) }
        if (params.run_cibersort)   { CIBERSORT_IMMUNE(risk_ch, expr_for_post, gene_for_post, params.train_id) }
        immune_for_cnv = params.run_cibersort ?
            CIBERSORT_IMMUNE.out.cnv_immune_input :
            Channel.empty()
        if (params.run_cnv)         { CNV_ANALYSIS(risk_ch, expr_for_post, gene_for_post, immune_for_cnv, params.train_id) }
        if (params.run_gistic)      { GISTIC_INPUT(risk_ch, expr_for_post, params.train_id) }
        if (params.run_tide)        { TIDE_ANALYSIS(risk_ch, expr_for_post, params.train_id) }
        if (params.run_ips)         { IPS_ANALYSIS(risk_ch, expr_for_post, params.train_id) }
        if (params.run_circos)      { CIRCOS_PLOT(params.gwas_file, risk_ch, gene_for_post, params.train_id) }
        if (params.run_tmb)         { TMB_ONCOPLOT(risk_ch, expr_for_post, params.train_id) }
    }

    // --- Clinical nomogram (disabled by default) ---
    if (params.run_clinical_nomogram) {
        error "CLINICAL_NOMOGRAM is disabled until multi-model pipeline stabilizes. Set --run_clinical_nomogram false."
    }
}

workflow.onComplete {
    def infoDir = file("${params.run_dir}/pipeline_info"); if (!infoDir.exists()) infoDir.mkdirs()
    def mmdFile = new File("${params.run_dir}/pipeline_info/dag.mmd")
    if (mmdFile.exists() && mmdFile.text.trim().startsWith("flowchart")) {
        def htmlFile = new File("${params.run_dir}/pipeline_info/dag.html")
        htmlFile.text = """\
<!DOCTYPE html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<script src="https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"></script>
<style>body{margin:0;padding:20px;background:#fff}#c{overflow:auto;width:100%;height:95vh;border:1px solid #ddd}</style>
</head><body><div id="c"><pre class="mermaid">
${mmdFile.text}
</pre></div>
<script>mermaid.initialize({startOnLoad:true,theme:'default',flowchart:{useMaxWidth:false}});</script>
</body></html>"""
        log.info "DAG HTML generated: ${htmlFile}"
    }
}

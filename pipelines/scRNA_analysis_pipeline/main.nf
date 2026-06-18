nextflow.enable.dsl = 2

include { SCRNA_GEO_SUPP_DOWNLOAD }    from './modules/local/scrna_geo_supp_download'
include { SCRNA_RESOLVE_RAW_FILES }    from './modules/local/scrna_resolve_raw_files'
include { SCRNA_DETECT_RAW_SOURCE }    from './modules/local/scrna_detect_raw_source'
include { SCRNA_UNPACK_RAW_UMI }       from './modules/local/scrna_unpack_raw_umi'
include { SCRNA_UNPACK_10X }           from './modules/local/scrna_unpack_10x'
include { SCRNA_BUILD_MAINLINE_INPUT } from './modules/local/scrna_build_mainline_input'
include { SCRNA_PREPARE_INPUT }        from './modules/local/scrna_prepare_input'
include { SCRNA_QC_INTEGRATION }       from './modules/local/scrna_qc_integration'
include { SCRNA_CLUSTER_SCAN as SCRNA_MAIN_CLUSTER_SCAN }   from './modules/local/scrna_cluster_scan'
include { SCRNA_SUBSET_CLUSTER_SCAN } from './modules/local/scrna_subset_cluster_scan'
include { SCRNA_ANNOTATION_PREPARE }   from './modules/local/scrna_annotation_prepare'
include { SCRNA_ANNOTATION_APPLY }     from './modules/local/scrna_annotation_apply'
include { SCRNA_CELL_ENRICHMENT }      from './modules/local/scrna_cell_enrichment'
include { SCRNA_KEYGENE_ANALYSIS }     from './modules/local/scrna_keygene_analysis'
include { SCRNA_KEY_CELL_SELECT }      from './modules/local/scrna_key_cell_select'
include { SCRNA_CELLCHAT }             from './modules/local/scrna_cellchat'
include { SCRNA_SUBSET_PREPARE }       from './modules/local/scrna_subset_prepare'
include { SCRNA_SUBSET_APPLY }         from './modules/local/scrna_subset_apply'
include { SCRNA_SUBSET_ANNOTATION_PREPARE } from './modules/local/scrna_subset_annotation_prepare'
include { SCRNA_PSEUDOTIME_MONOCLE }   from './modules/local/scrna_pseudotime_monocle'
include { SCRNA_PSEUDOTIME_ORDERING_TARGET } from './modules/local/scrna_pseudotime_ordering_target'
include { SCRNA_SCMETABOLISM }         from './modules/local/scrna_scmetabolism'
include { SCRNA_SCMETABOLISM_SPECIFIC } from './modules/local/scrna_scmetabolism_specific'
include { SCRNA_GSVA }                 from './modules/local/scrna_gsva'
include { SCRNA_REACTOME_GSA }         from './modules/local/scrna_reactome_gsa'

workflow {
    if (!params.project_id?.toString()?.trim()) {
        error "`params.project_id` is required."
    }
    if (!params.rscript?.toString()?.trim()) {
        error "Set `NFMODELS_RSCRIPT_CMD` or `SCRNA_RSCRIPT_CMD` in ../.env, then launch with ./run_pipeline.sh"
    }

    has_input_rds = params.input_rds?.toString()?.trim()
    run_raw = params.run_raw_import

    if (!has_input_rds && !run_raw) {
        error "Provide `--input_rds` or set `--run_raw_import true --scRNA_cohort_id <GSE...>`."
    }
    if (has_input_rds && run_raw) {
        log.info "Both input_rds and run_raw_import were supplied; using input_rds and skipping raw import."
    }
    if (!has_input_rds && run_raw && !params.scRNA_cohort_id?.toString()?.trim()) {
        error "`scRNA_cohort_id` is required when `run_raw_import=true`."
    }

    downstream_requires_subset_apply = params.run_pseudotime || params.run_scmetabolism ||
        params.run_scmetabolism_specific || params.run_gsva || params.run_reactome_gsa
    subset_requested = params.run_subset_prepare || params.run_subset_apply || downstream_requires_subset_apply
    need_key_cell_select = params.run_key_cell_select || subset_requested
    if (subset_requested && !params.mapping_file) {
        error "Subset and downstream modules require main `--mapping_file` so celltype_manual exists."
    }
    if ((params.run_subset_apply || downstream_requires_subset_apply) && !params.subset_mapping_file && !params.subset_mapping_manifest) {
        error "Subset apply/downstream modules require `--subset_mapping_file` for a single key cell or `--subset_mapping_manifest` for multiple key cells."
    }
    if (params.run_gsva && !params.gsva_gmt_files?.toString()?.trim()) {
        error "`--gsva_gmt_files` is required when `run_gsva=true`."
    }

    if (has_input_rds) {
        input_ch = Channel.fromPath(params.input_rds, checkIfExists: true)
    } else {
        SCRNA_GEO_SUPP_DOWNLOAD()
        SCRNA_RESOLVE_RAW_FILES(SCRNA_GEO_SUPP_DOWNLOAD.out.manifest)
        SCRNA_DETECT_RAW_SOURCE(SCRNA_RESOLVE_RAW_FILES.out.resolved_files)
        SCRNA_UNPACK_RAW_UMI(SCRNA_DETECT_RAW_SOURCE.out.raw_source_manifest)
        SCRNA_UNPACK_10X(SCRNA_UNPACK_RAW_UMI.out.raw_source_after_umi)
        SCRNA_BUILD_MAINLINE_INPUT(SCRNA_UNPACK_10X.out.raw_source_ready)
        input_ch = SCRNA_BUILD_MAINLINE_INPUT.out.mainline_input
    }

    SCRNA_PREPARE_INPUT(input_ch)
    SCRNA_QC_INTEGRATION(SCRNA_PREPARE_INPUT.out.seurat_rds)
    SCRNA_MAIN_CLUSTER_SCAN(
        SCRNA_QC_INTEGRATION.out.precluster_rds,
        'main',
        params.main_resolution_grid,
        params.main_final_resolution
    )
    SCRNA_ANNOTATION_PREPARE(SCRNA_MAIN_CLUSTER_SCAN.out.clustered_rds)

    if (params.mapping_file) {
        mapping_ch = Channel.fromPath(params.mapping_file, checkIfExists: true)
        SCRNA_ANNOTATION_APPLY(SCRNA_ANNOTATION_PREPARE.out.prepare_rds, mapping_ch)

        if (params.run_cell_enrichment || need_key_cell_select) {
            SCRNA_CELL_ENRICHMENT(SCRNA_ANNOTATION_APPLY.out.annotated_rds)
        }

        if (params.run_keygene_analysis || need_key_cell_select) {
            SCRNA_KEYGENE_ANALYSIS(SCRNA_ANNOTATION_APPLY.out.annotated_rds)
        }

        if (need_key_cell_select) {
            SCRNA_KEY_CELL_SELECT(
                SCRNA_CELL_ENRICHMENT.out.proportions,
                SCRNA_CELL_ENRICHMENT.out.enrichment_tests,
                SCRNA_KEYGENE_ANALYSIS.out.group_diff
            )
        }

        if (params.run_cellchat) {
            SCRNA_CELLCHAT(SCRNA_ANNOTATION_APPLY.out.annotated_rds)
        }

        if (subset_requested) {
            key_cell_ch = SCRNA_KEY_CELL_SELECT.out.key_celltypes
                .splitCsv(header: true)
                .filter { row -> row.key_celltype?.toString()?.trim() }
                .map { row -> tuple(row.safe_label.toString(), row.key_celltype.toString()) }

            subset_input_ch = key_cell_ch
                .combine(SCRNA_ANNOTATION_APPLY.out.annotated_rds)
                .map { safe_label, key_celltype, annotated_rds -> tuple(safe_label, key_celltype, annotated_rds) }

            SCRNA_SUBSET_PREPARE(subset_input_ch)
            SCRNA_SUBSET_CLUSTER_SCAN(SCRNA_SUBSET_PREPARE.out.subset_precluster_rds)
            SCRNA_SUBSET_ANNOTATION_PREPARE(SCRNA_SUBSET_CLUSTER_SCAN.out.clustered_rds)

            if (params.run_subset_apply || downstream_requires_subset_apply) {
                if (params.subset_mapping_manifest) {
                    subset_mapping_ch = Channel.fromPath(params.subset_mapping_manifest, checkIfExists: true)
                        .splitCsv(header: true)
                        .map { row -> tuple(row.safe_label.toString(), file(row.mapping_file.toString())) }
                } else {
                    subset_mapping_ch = key_cell_ch.map { safe_label, key_celltype -> tuple(safe_label, file(params.subset_mapping_file)) }
                }

                subset_apply_input_ch = SCRNA_SUBSET_ANNOTATION_PREPARE.out.prepare_rds
                    .join(subset_mapping_ch)
                    .map { safe_label, key_celltype, prepare_rds, mapping_file -> tuple(safe_label, key_celltype, prepare_rds, mapping_file) }
                SCRNA_SUBSET_APPLY(subset_apply_input_ch)

                if (params.run_pseudotime) {
                    SCRNA_PSEUDOTIME_MONOCLE(SCRNA_SUBSET_APPLY.out.subset_annotated_rds)
                    if (params.pseudotime_ordering_gene_file?.toString()?.trim()) {
                        SCRNA_PSEUDOTIME_ORDERING_TARGET(SCRNA_SUBSET_APPLY.out.subset_annotated_rds)
                    }
                }

                if (params.run_scmetabolism) {
                    SCRNA_SCMETABOLISM(SCRNA_SUBSET_APPLY.out.subset_annotated_rds)
                }
                if (params.run_scmetabolism_specific) {
                    SCRNA_SCMETABOLISM_SPECIFIC(SCRNA_SUBSET_APPLY.out.subset_annotated_rds)
                }
                if (params.run_gsva) {
                    SCRNA_GSVA(SCRNA_SUBSET_APPLY.out.subset_annotated_rds)
                }
                if (params.run_reactome_gsa) {
                    SCRNA_REACTOME_GSA(SCRNA_SUBSET_APPLY.out.subset_annotated_rds)
                }
            }
        }
    } else {
        log.info "No mapping_file supplied; workflow will stop after SCRNA_ANNOTATION_PREPARE."
    }
}

workflow.onComplete {
    def infoDir = file("${params.run_dir}/pipeline_info")
    if (!infoDir.exists()) infoDir.mkdirs()
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

# scRNA Analysis Pipeline Parameters

## Required / Common

| Parameter | Default | Description |
| --- | --- | --- |
| `project_id` | `scrna_test` | Run ID under `${project_root}`. |
| `input_rds` | `null` | Existing input RDS: `list(expr_mat, meta)` or Seurat object. Required unless `run_raw_import=true`. |
| `input_format` | `auto` | `auto`, `mainline_rds`, or `seurat_rds`. |
| `mapping_file` | `null` | Filled manual annotation mapping file. If absent, workflow stops after annotation prepare. |
| `target_gene_file` | `none` | Optional prognosis/candidate gene file. Must be manually specified. |
| `target_gene_col` | `gene` | Gene column in `target_gene_file`. |

## Raw Import

| Parameter | Default | Description |
| --- | --- | --- |
| `run_raw_import` | `false` | Download and build `mainline_input.rds` when `input_rds` is not supplied. |
| `scRNA_cohort_id` | `null` | Single-cell cohort ID, usually a GEO GSE accession. |
| `scRNA_database_dir` | `${user_dir}/database/scRNA` | Root directory for scRNA raw data. |
| `scRNA_raw_dir` | `${scRNA_database_dir}/${scRNA_cohort_id}` | Cohort raw-data directory; CLI-overridable. |
| `geo_filter_regex` | `none` | Optional GEOquery supplementary file filter. |
| `scRNA_raw_umi_file` | `none` | Optional explicit raw source file/path. Despite the legacy name, this can be UMI RDS, 10X matrix directory, 10X archive, or 10X h5. |
| `scRNA_annotation_file` | `none` | Optional explicit annotation file name/path. |
| `raw_source_type` | `auto` | Raw source type: `auto`, `umi_rds`, `10x_dir`, `10x_archive`, or `10x_h5`. |
| `raw_umi_pattern` | `raw_UMI\|raw.*matrix\|UMI_matrix\|filtered_feature_bc_matrix\|raw_feature_bc_matrix\|matrix.mtx\|10x\|h5` | Pattern used to identify the raw source if not explicit. |
| `raw_10x_pattern` | `filtered_feature_bc_matrix\|raw_feature_bc_matrix\|matrix.mtx\|10x` | Pattern for recognizing 10X-like source names. |
| `raw_10x_prefer` | `filtered_feature_bc_matrix` | Preferred 10X directory pattern when multiple matrix directories exist. |
| `annotation_pattern` | `annotation\|cell_annotation\|metadata` | Pattern used to identify annotation file if not explicit. |
| `decompress_max_depth` | `5` | Maximum gzip magic-header decompression iterations for UMI RDS. |
| `strip_10x_barcode_suffix` | `false` | If true, allows matching `AAAC-1` expression barcodes to `AAAC` annotation IDs. |
| `tenx_assay` | `Gene Expression` | Assay selected from multi-assay 10X inputs. |
| `keep_origin` | `nLung,tLung` | Sample origins retained for mainline analysis. |
| `exclude_origin` | `tL/B` | Sample origins excluded after retention filtering. |
| `normal_origins` | `nLung` | Origins mapped to `normal`. |
| `tumor_origins` | `tLung` | Origins mapped to `tumor`. |
| `annotation_cell_id_candidates` | `Index,barcode_sample,Barcode` | Annotation columns tried for cell ID matching. |

## Metadata

| Parameter | Default | Description |
| --- | --- | --- |
| `sample_col` | `Sample` | Sample metadata column. |
| `origin_col` | `Sample_Origin` | Source/origin metadata column. |
| `group_col` | `group` | Group column for Tumor/Normal or other comparisons. |
| `celltype_col` | `celltype_manual` | Manual cell type metadata column after annotation apply. |
| `normal_origin_value` | `nLung` | Origin value mapped to `normal`. |
| `tumor_origin_value` | `tLung` | Origin value mapped to `tumor`. |

## QC / Integration

| Parameter | Default |
| --- | --- |
| `min_cells_per_gene` | `3` |
| `min_features_create` | `0` |
| `min_features_keep` | `200` |
| `max_features_keep` | `4000` |
| `min_counts_keep` | `500` |
| `max_counts_keep` | `20000` |
| `max_percent_mt` | `15` |
| `run_doublet_removal` | `true` |
| `nfeatures_hvg` | `3000` |
| `npcs` | `50` |
| `use_dims` | `1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20` |
| `harmony_batch_var` | `Sample` |
| `seed` | `1234` |

## Cluster Scan

| Parameter | Default |
| --- | --- |
| `preferred_reduction` | `harmony` |
| `fallback_reduction` | `pca` |
| `main_resolution_grid` | `0.1,0.2,0.3,0.4,0.5,0.6` |
| `main_final_resolution` | `0.2` |
| `subset_resolution_grid` | `0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0` |
| `subset_final_resolution` | `0.4` |
| `cluster_scan_use_dims_override` | empty |
| `cluster_algorithm` | `1` |
| `umap_n_neighbors` | `30` |
| `umap_min_dist` | `0.3` |

## Annotation / Downstream

| Parameter | Default | Description |
| --- | --- | --- |
| `annotation_cluster_col` | `seurat_clusters` | Main annotation cluster column. |
| `annotation_resolution` | empty | Optional main resolution override for annotation prepare. |
| `marker_min_pct` | `0.25` | Marker minimum pct. |
| `marker_logfc_threshold` | `0.25` | Marker logFC threshold. |
| `annotation_literature_years` | `10` | Suggested literature search lookback for agent context. |
| `annotation_marker_reference_file` | `null` | Filled `09_literature_marker_reference_filled.csv` used by marker-support audit. |
| `annotation_min_support_markers` | `2` | Minimum expressed, literature-supported markers required per cluster. |
| `annotation_allow_low_support` | `false` | If true, low-support clusters are retained and flagged instead of failing. |
| `run_cell_enrichment` | `true` | Run cell proportion/enrichment after main annotation apply. |
| `run_keygene_analysis` | `true` | Run target/key gene expression and group-difference summaries after main annotation apply. |
| `run_key_cell_select` | `true` | Select key cell types from enrichment and keygene evidence. |
| `key_celltypes` | empty | Manual comma-separated key cell override; otherwise uses automatic selection. |
| `keycell_enrichment_fdr` | `0.05` | Enrichment adjusted p-value cutoff for automatic key cell selection. |
| `keycell_min_abs_prop_diff` | `0.05` | Minimum absolute cell proportion difference. |
| `keycell_min_pct_expr` | `10` | Minimum keygene expression percentage. |
| `keycell_gene_diff_fdr` | `0.05` | Keygene group-difference adjusted p-value cutoff. |
| `keycell_min_abs_log2fc` | `0.25` | Minimum absolute keygene log2FC. |
| `keycell_min_support_genes` | `1` | Minimum number of key genes supporting a key cell. |
| `run_cellchat` | `false` | Run CellChat after main annotation apply. |
| `cellchat_species` | `human` | CellChat DB species. |
| `cellchat_min_cells` | `10` | Minimum cells per celltype in CellChat. |

## Subset / Pseudotime / Pathway Modules

| Parameter | Default | Description |
| --- | --- | --- |
| `run_subset_prepare` | `false` | Build subset precluster object from main annotated object. |
| `run_subset_apply` | `false` | Apply filled subset cluster mapping. |
| `subset_mapping_file` | `null` | Filled subset mapping file shortcut for a single key cell. |
| `subset_mapping_manifest` | `null` | CSV manifest for multiple key cells: `key_celltype,safe_label,mapping_file`. |
| `subset_npcs` | `30` | PCA dimensions computed in subset prepare. |
| `subset_use_dims` | `1,2,...,20` | Subset dims written into object metadata for subset cluster scan. |
| `subset_nfeatures_hvg` | `2000` | HVG count for subset prepare. |
| `subset_run_jackstraw` | `false` | Optional JackStraw for subset PCA. |
| `run_pseudotime` | `false` | Run Monocle pseudotime on subset annotated object. |
| `pseudotime_ordering_gene_file` | `null` | Optional target ordering gene file. |
| `pseudotime_ordering_gene_col` | `gene` | Ordering gene column. |
| `pseudotime_hvg_nfeatures` | `3000` | HVG ordering candidate count. |
| `pseudotime_min_ordering_genes` | `50` | Minimum ordering genes. |
| `run_scmetabolism` | `false` | Run scMetabolism scoring on subset annotated object. |
| `run_scmetabolism_specific` | `false` | Run key-celltype specificity analysis. |
| `scmetabolism_top_n_pathways` | `10` | Top pathway count. |
| `scmetabolism_min_fdr` | `0.05` | FDR cutoff for specificity. |
| `scmetabolism_min_delta` | `0` | Mean-score delta cutoff. |
| `scmetabolism_highlight_groups` | `T/NK cells,Monocytes,Macrophages,B cells,Epithelial cells` | Groups highlighted in specificity plots. |
| `run_gsva` | `false` | Run GSVA on subset annotated object. |
| `gsva_gmt_files` | `null` | Comma-separated GMT files; required when GSVA is enabled. |
| `gsva_min_geneset_size` | `5` | GSVA gene-set lower size filter. |
| `gsva_max_geneset_size` | `500` | GSVA gene-set upper size filter. |
| `gsva_threshold_t` | `2` | GSVA limma t threshold. |
| `gsva_threshold_p` | `0.05` | GSVA p-value threshold. |
| `run_reactome_gsa` | `false` | Run ReactomeGSA on subset annotated object. |
| `reactome_max_pathways` | `10` | Max pathways in ReactomeGSA heatmap. |
| `reactome_p_cutoff` | `0.05` | ReactomeGSA significance cutoff. |

Fixed internal contracts: subset reads `celltype_manual`, subset apply writes `subset_celltype_manual`, pathway/pseudotime modules use each per-key-cell subset annotated object, and group comparisons read `group`. These are not CLI parameters.

## Project Runtime Environment

Runtime binary locations are not hardcoded in the pipeline. Configure them in `../.env`; `run_pipeline.sh` loads that file automatically before launching Nextflow:

```bash
./run_pipeline.sh --input_rds /path/to/input.rds -profile local
```

| Environment variable | Used for | Default in example |
| --- | --- | --- |
| `NFMODELS_NEXTFLOW_BIN` | Shared command used to launch Nextflow | `nextflow` |
| `NFMODELS_RSCRIPT_CMD` | Shared R command used by Nextflow processes | `conda run -n R4.3.3 Rscript` |
| `SCRNA_NEXTFLOW_BIN` | Optional scRNA-specific Nextflow override | inherits `NFMODELS_NEXTFLOW_BIN` |
| `SCRNA_RSCRIPT_CMD` | Optional scRNA-specific R override | inherits `NFMODELS_RSCRIPT_CMD` |

`NFMODELS_RSCRIPT_CMD` and any pipeline-specific override should remain Conda-managed R commands in this workspace.

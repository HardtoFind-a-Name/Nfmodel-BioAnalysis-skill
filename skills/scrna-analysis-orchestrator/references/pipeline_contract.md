# scRNA Pipeline Contract

The scRNA adapter launches `pipelines/scRNA_analysis_pipeline/run_pipeline.sh` and coordinates two manual annotation gates with `scrna-celltype-annotator`.

## Core Inputs

- `--input_rds`: existing mainline RDS or Seurat RDS.
- `--input_format`: `auto`, `mainline_rds`, or `seurat_rds`.
- `--mapping_file`: reviewed main annotation mapping file.
- `--annotation_marker_reference_file`: reviewed literature marker reference file for main annotation support audit.
- `--target_gene_file`: optional, manually supplied only.
- `--target_gene_col`: default `gene`.

## Main Annotation Gate

1. Run without `mapping_file` to stop after `SCRNA_ANNOTATION_PREPARE`.
2. Hand these files to `scrna-celltype-annotator`:
   - `08_cluster_celltype_mapping_template.csv`
   - `09_literature_marker_reference_template.csv`
   - `10_annotation_agent_context.json`
   - marker tables such as `02_all_markers.csv`, `03_top10_markers_by_cluster.csv`, and `04_top20_markers_by_cluster.csv`
3. The annotator produces reviewed candidate files:
   - `08_cluster_celltype_mapping_filled.csv`
   - `09_literature_marker_reference_filled.csv`
   - `annotation_literature_evidence.md`
   - `annotation_decision_log.csv`
4. User review remains required before the files are passed back to NF.

## Main Apply And Subset Prepare

- Run with `--mapping_file` and `--annotation_marker_reference_file` to apply `celltype_manual`.
- Use `--run_subset_prepare` or downstream subset toggles when key-cell/subtype analysis is needed.
- The adapter can emit `scrna_annotation_handoff.json` recording expected prepare outputs and annotator outputs.

## Subset Annotation Gate

- Subset prepare writes per-key-cell directories under `07b_scrna_subset_annotation_prepare/<safe_label>/`.
- Each subset directory is handed to `scrna-celltype-annotator` using subset-specific templates:
  - `08_subset_cluster_celltype_mapping_template.csv`
  - `09_subset_literature_marker_reference_template.csv`
  - `10_subset_annotation_agent_context.json`
- For multiple key cells, create `subset_mapping_manifest` with columns: `key_celltype,safe_label,mapping_file`.
- Continue with `--run_subset_apply true --subset_mapping_manifest <manifest>` after subset mapping review.

## Target Genes

If `target_gene_file` is absent, skip key gene plots and module-score outputs. Do not infer files from transcriptome model directories.

## Boundaries

The adapter generates commands and handoff metadata only. It does not perform literature search, invent cell labels, or auto-approve annotation files.

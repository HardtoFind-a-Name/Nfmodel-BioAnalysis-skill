# scRNA Analysis Pipeline

Independent single-cell RNA-seq analysis pipeline with Nextflow and R.

## Workflow

```text
Existing input_rds
  -> Prepare input (00)
  -> QC + integration / precluster object (01)
  -> Main cluster scan (02)
  -> Annotation prepare (03)
      -> stop here if --mapping_file is not supplied
      -> Annotation apply (04)
          -> Cell enrichment (05, optional)
          -> CellChat (06, optional)

Raw import path
  -> GEO supplementary download (00a)
  -> Resolve raw source + annotation files (00b)
  -> Detect raw source type: umi_rds / 10x_dir / 10x_archive / 10x_h5 (00c)
  -> Unpack UMI gzip layers when needed (00d)
  -> Unpack/locate 10X matrix directory when needed (00e)
  -> Build mainline_input.rds (00f)
  -> Prepare input (00)
```

The pipeline is independent from `transcriptome_prognosis_pipeline`. It can link to prognosis results only through explicit file parameters such as `--target_gene_file`.

## Input

Use one of two entry points:

- Existing `--input_rds`: supports `list(expr_mat, meta)` RDS or Seurat object RDS.
- Raw GEO import: set `--run_raw_import true --scRNA_cohort_id <GSE...>`.

Raw import downloads GEO supplementary files into `${user_dir}/database/scRNA/${scRNA_cohort_id}` by default. Override with `--scRNA_raw_dir` if needed. Raw source and annotation files are not hardcoded; specify them by CLI if automatic pattern matching is ambiguous. 10X/UMI compatibility is handled before the main analysis line by building a standard `mainline_input.rds`.

The pipeline writes standard `.rds` objects and does not require `qs`.

## Agent-assisted Manual Annotation

Run without `--mapping_file` first. The workflow stops after `SCRNA_ANNOTATION_PREPARE` and publishes:

```text
results/03_scrna_annotation_prepare/08_cluster_celltype_mapping_template.csv
results/03_scrna_annotation_prepare/09_literature_marker_reference_template.csv
results/03_scrna_annotation_prepare/10_annotation_agent_context.json
```

Use the draft skill at `/data/nas1/public_JOB062/skill-construct/scrna-celltype-annotator` to perform PubMed-first marker review and fill:

```text
08_cluster_celltype_mapping_filled.csv
09_literature_marker_reference_filled.csv
annotation_literature_evidence.md
annotation_decision_log.csv
```

Rerun apply with the filled mapping and marker-reference table. By default, each cluster must have at least two expressed markers supported by references:

```bash
./run_pipeline.sh \
  --project_id scrna_demo \
  --input_rds /path/to/input.rds \
  --mapping_file /path/to/08_cluster_celltype_mapping_filled.csv \
  --annotation_marker_reference_file /path/to/09_literature_marker_reference_filled.csv \
  -profile local
```

## Raw Import Examples

```bash
./run_pipeline.sh \
  --run_raw_import true \
  --scRNA_cohort_id GSE131907 \
  -profile local
```

Automatically detect UMI RDS or 10X source:

```bash
./run_pipeline.sh \
  --run_raw_import true \
  --scRNA_cohort_id GSE131907 \
  --raw_source_type auto \
  -profile local
```

Explicit 10X directory or archive:

```bash
./run_pipeline.sh \
  --run_raw_import true \
  --scRNA_cohort_id GSE131907 \
  --raw_source_type 10x_dir \
  --scRNA_raw_umi_file /path/to/filtered_feature_bc_matrix \
  --scRNA_annotation_file /path/to/annotation.tsv.gz \
  -profile local
```

Explicit UMI RDS:

```bash
./run_pipeline.sh \
  --run_raw_import true \
  --scRNA_cohort_id GSE131907 \
  --raw_source_type umi_rds \
  --scRNA_raw_umi_file /path/to/raw_UMI_matrix.rds.gz \
  --scRNA_annotation_file /path/to/cell_annotation.txt.gz \
  -profile local
```

## Cluster Scan

`SCRNA_CLUSTER_SCAN` inherits `reduction_used` and `dims_used` from upstream object metadata. Adjust `main_resolution_grid` and `main_final_resolution` for the main branch. Use `cluster_scan_use_dims_override` only when intentionally overriding upstream dims. The same module is aliased for future subset clustering.


## Key Cells, Subsets, And Downstream

After main annotation apply, the pipeline can run keygene analysis and select key cell types from cell enrichment plus keygene expression/group-difference evidence. Override automatic selection with comma-separated `--key_celltypes` when needed.

Each key cell is processed as an independent branch:

```text
key_celltypes.csv
  -> subset prepare per key cell
  -> subset cluster scan per key cell
  -> subset annotation prepare per key cell
  -> agent second-pass annotation
  -> subset apply per key cell
  -> downstream per key cell
```

For one key cell, `--subset_mapping_file` is a shortcut. For multiple key cells, provide `--subset_mapping_manifest` with columns `key_celltype,safe_label,mapping_file`:

```bash
./run_pipeline.sh \
  --input_rds /path/to/input.rds \
  --mapping_file /path/to/main_mapping_filled.csv \
  --annotation_marker_reference_file /path/to/main_marker_refs_filled.csv \
  --key_celltypes "NK cells,T cells" \
  --run_subset_apply true \
  --subset_mapping_manifest /path/to/subset_mapping_manifest.csv \
  --run_pseudotime true \
  --run_scmetabolism_specific true \
  -profile local
```

GSVA requires user-supplied GMT files:

```bash
--run_gsva true --gsva_gmt_files /path/to/hallmark.gmt,/path/to/kegg.gmt
```

Fixed downstream contracts are not CLI parameters: subset reads `celltype_manual`, subset apply writes `subset_celltype_manual`, pathway and pseudotime modules consume the subset annotated object, and group comparisons read the upstream `group` column.

## Target Genes

Prognosis or candidate genes are never hardcoded. Provide them manually:

```bash
--target_gene_file /path/to/01.feature_importance.csv \
--target_gene_col gene
```

If `--target_gene_file` is not supplied, target gene plots are skipped and the main scRNA workflow continues.

## Runtime

Runtime commands are configured in `../.env` and loaded automatically by `run_pipeline.sh`.

```bash
./run_pipeline.sh --input_rds /path/to/input.rds -profile local
```

Nextflow processes use `params.rscript`, which is read from `SCRNA_RSCRIPT_CMD` first, then `NFMODELS_RSCRIPT_CMD`. Keep this as a Conda-managed R command in this workspace; do not use bare `Rscript`.

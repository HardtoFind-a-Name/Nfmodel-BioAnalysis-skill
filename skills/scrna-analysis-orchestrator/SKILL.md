---
name: scrna-analysis-orchestrator
description: Generate, check, and optionally run staged commands for NFmodels single-cell RNA analysis through pipelines/scRNA_analysis_pipeline, including handoff to scrna-celltype-annotator for main and subset manual cell type annotation gates. Use for scRNA input preparation, QC, integration, clustering, annotation_prepare, agent-assisted mapping review, annotation_apply, subset annotation, cell enrichment, CellChat, and optional downstream analyses.
---

# scRNA Analysis Orchestrator

Use this skill as the single-cell adapter for `pipelines/scRNA_analysis_pipeline`.

## Role

- Convert a reviewed scRNA request into a `pipelines/scRNA_analysis_pipeline/run_pipeline.sh` command.
- Enforce the staged annotation flow and hand off mapping templates to `scrna-celltype-annotator`.
- Check core parameters: `input_rds` or raw import settings, `mapping_file`, `annotation_marker_reference_file`, `subset_mapping_file`/`subset_mapping_manifest`, `target_gene_file`, and downstream subset toggles.
- Reuse the project-level `.env` through `run_pipeline.sh`; never ask the user to manually source `.env`.

## Workflow

1. Confirm `nfmodels-environment-check --profile scrna` has passed.
2. Confirm the input mode.
   - Use `input_rds` for an existing `list(expr_mat, meta)` RDS or Seurat RDS.
   - Use raw import parameters only when the user explicitly asks to build input from raw scRNA files.
3. Generate a staged command and optional `scrna_annotation_handoff.json` with `scripts/generate_scrna_command.py`.
4. If `mapping_file` is absent, run through main `annotation_prepare` only; this is expected and not a failure. Hand off `08_cluster_celltype_mapping_template.csv`, `09_literature_marker_reference_template.csv`, and `10_annotation_agent_context.json` to `scrna-celltype-annotator`.
5. After agent annotation and user review, regenerate the command with `mapping_file` and `annotation_marker_reference_file` to continue `annotation_apply`, cell enrichment, key-cell selection, and optional CellChat.
6. If subset analysis is requested, run through subset annotation prepare, then hand off each subset template directory to `scrna-celltype-annotator`; after user review, continue with `subset_mapping_file` or `subset_mapping_manifest`.
7. Include `target_gene_file` only when the user manually provides it. Do not infer prognosis model paths.

## Boundaries

- Do not run the bulk transcriptome pipeline from this skill.
- Do not hardcode R or Nextflow paths in commands.
- Do not invent `target_gene_top_n`, `target_gene_rank_col`, or `target_gene_rank_desc`; those parameters are not part of the scRNA public interface.
- Do not perform literature-backed cell type decisions inside this adapter; delegate that work to `scrna-celltype-annotator` and treat its output as user-review material.

## Resources

- `scripts/generate_scrna_command.py`: generate staged launch commands, preflight blocker reports, and annotation handoff JSON.
- `references/pipeline_contract.md`: scRNA stages, required inputs, and annotation handoff rules.
- `scrna-celltype-annotator`: separate skill used to fill and validate main/subset manual annotation mapping files.

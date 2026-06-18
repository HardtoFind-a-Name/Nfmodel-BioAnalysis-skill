---
name: scrna-celltype-annotator
description: Literature-backed two-pass scRNA cell type annotation from cluster marker tables, annotation templates, and NFmodels scRNA agent context files. Use when Codex needs to fill or review scRNA cluster mapping tables, generate marker-reference evidence, or annotate main/subset scRNA clusters using PubMed-first literature search.
---

# scRNA Celltype Annotator

Use this skill to annotate scRNA clusters from NFmodels scRNA pipeline outputs. Always combine cluster marker evidence with literature evidence; do not assign a cell type from memory alone.

## Inputs

Expected files from the pipeline:

- `08_cluster_celltype_mapping_template.csv` or subset equivalent
- `09_literature_marker_reference_template.csv` or subset equivalent
- `10_annotation_agent_context.json` or subset equivalent
- marker tables such as `02_all_markers.csv`, `03_top10_markers_by_cluster.csv`, `04_top20_markers_by_cluster.csv`

## Workflow

1. Read `references/annotation-workflow.md`.
2. If literature search is needed, use `literature-search-review` with PubMed-first search. Record search date, query, included papers, PMID/DOI/URL, and abstract/full-text basis.
3. Fill the marker reference table with `manual_markers,references`. References must include PMID, DOI, or URL.
4. Fill the mapping table with `cluster,manual_celltype,manual_markers,annotation_flag,annotation_confidence,notes`.
5. Run `scripts/validate_mapping.py` before returning the files. Fix all validation failures unless the user explicitly allows low-support annotations.

## Output Contract

Produce these files in the requested output directory:

- `08_cluster_celltype_mapping_filled.csv`
- `09_literature_marker_reference_filled.csv`
- `annotation_literature_evidence.md`
- `annotation_decision_log.csv`

For subset annotation, keep the same schema and use subset-specific filenames if the pipeline context requests them.

## Annotation Rules

- Use at least two literature-supported markers per cluster by default.
- A marker counts as support only when it is present in the cluster marker/expression evidence and has a non-empty literature reference.
- Use `annotation_flag=manual_review` for uncertain calls and explain the uncertainty in `notes`.
- Use `annotation_flag=low_support` only when low-support annotations are explicitly allowed.
- Never fabricate references, PMIDs, DOIs, marker evidence, or cell type labels.

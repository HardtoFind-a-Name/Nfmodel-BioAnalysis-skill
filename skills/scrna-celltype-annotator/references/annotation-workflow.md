# scRNA Annotation Workflow

## Pass 1: Main Cell Type Annotation

Search concepts:

```text
(scRNA_cohort_id OR disease_name OR cancer_name)
AND (single-cell OR scRNA OR single cell RNA sequencing)
AND (cell type OR marker OR atlas)
```

If the cohort-specific search is sparse, broaden to disease/cancer single-cell atlas papers. Prioritize human disease/cancer scRNA studies, cell atlas papers, and marker validation papers.

## Pass 2: Key Cell Subtype Annotation

Search concepts:

```text
(key_celltype)
AND (disease_name OR cancer_name)
AND (single-cell OR scRNA)
AND (subtype OR state OR marker)
```

Use subset cluster marker genes as additional search terms only after a broad pass.

## Evidence Fields

For each paper keep: title, year, first author, journal, PMID, DOI or URL, dataset/model, main marker/cell type finding, limitation, and whether the basis was abstract-only or full text.

## Mapping Decisions

Fill mapping rows as follows:

- `manual_celltype`: concise biological label.
- `manual_markers`: comma-separated marker genes supporting the label.
- `annotation_flag`: `literature_supported`, `manual_review`, or `low_support`.
- `annotation_confidence`: `high`, `medium`, or `low`.
- `notes`: one-line reason including marker evidence and any ambiguity.

## Validation

A cluster passes by default when at least two markers from `manual_markers` are present in the reference table with non-empty references. If a marker is not in the literature reference table, it does not count as support.

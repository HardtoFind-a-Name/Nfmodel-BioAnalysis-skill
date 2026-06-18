# Bulk Transcriptome Pipeline Capabilities

`pipelines/transcriptome_prognosis_pipeline` is historically named for prognosis modeling, but it supports broader bulk transcriptome analyses.

## Fully Supported

- DEG: `--run_deg`
- candidate genes: `--run_intersect_candidate_genes` or precomputed `--candidate_genes`
- GO/KEGG/PPI: `--run_go_kegg_ppi`
- multi-model prognosis: `--run_multi_model`
- stage/prognostic and nomogram: `--run_stage`
- GSEA: `--run_gsea_gene`, `--run_gsea_risk_nes`
- immune infiltration: `--run_cibersort`
- TIDE: `--run_tide`
- IPS: `--run_ips`
- selected-gene Tumor/Normal differential expression: `--run_diff_expr`

## Supported With Caveats

- CNV depends on TCGA-style resources and optional local matrices.
- GISTIC is currently input generation and handoff, not a full GISTIC2 run.
- Circos and TMB have defaults that may need project-specific review.

## Launch Policy

Generate commands with:

```bash
pipelines/transcriptome_prognosis_pipeline/run_pipeline.sh
```

The launcher loads `.env` (project root); do not add a manual `source .env` step.

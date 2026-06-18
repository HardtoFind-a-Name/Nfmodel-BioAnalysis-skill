---
name: transcriptome-analysis-orchestrator
description: Generate, check, and optionally run commands for NFmodels bulk transcriptome analyses through pipelines/transcriptome_prognosis_pipeline. Use for DEG, candidate genes, GO/KEGG/PPI, multi-model prognosis, nomogram/stage analysis, GSEA, CIBERSORT immune infiltration, TIDE, IPS, CNV, GISTIC input, differential expression, circos, TMB, or other bulk transcriptome pipeline tasks.
---

# Transcriptome Analysis Orchestrator

Use this skill as the bulk transcriptome adapter. It does not parse full project reports or author global analysis plans; the top-level `nfmodels-orchestrator` owns that.

## Role

- Convert a reviewed bulk transcriptome request into a `pipelines/transcriptome_prognosis_pipeline/run_pipeline.sh` command.
- Check required inputs before execution.
- Explain supported modules and current caveats.
- Reuse the project-level `.env` through `run_pipeline.sh`; never ask the user to manually source `.env`.

## Pipeline Capabilities

The pipeline name is historical. Treat it as a general bulk transcriptome pipeline that includes:

- bulk DEG and DEG filtering
- candidate gene generation
- GO/KEGG/PPI
- multi-model prognosis modeling and validation
- stage/prognostic analysis and nomogram
- gene- and risk-based GSEA
- CIBERSORT immune infiltration
- TIDE and IPS
- CNV and GISTIC input generation
- Tumor/Normal differential expression for selected genes
- circos and TMB/oncoplot outputs

Read `references/pipeline_capabilities.md` when module boundaries matter.

## Workflow

1. Confirm `nfmodels-environment-check --profile transcriptome` has passed, using `--require-tidepy` when TIDE is requested.
2. Confirm project and data identifiers: `project_id`, `train_id`, and validation cohorts or a validation sheet when multi-model analysis is requested.
3. Run `scripts/generate_transcriptome_command.py` with the requested analyses and known paths.
4. Review the generated blockers and command.
5. Execute the command only after blockers are resolved and the user wants execution.
6. For post-model analyses, provide `risk_file`; provide `gene_file` when the selected module needs genes.

## Boundaries

- Do not include the old report-to-plan orchestration flow here.
- Do not call `nextflow run main.nf` directly by default.
- Do not hardcode R, Python, Nextflow, or TIDE binary paths in generated commands.
- Keep archived skill files read-only unless the user explicitly asks to inspect history.

## Resources

- `scripts/generate_transcriptome_command.py`: generate a launch command and preflight blocker report.
- `references/pipeline_capabilities.md`: module map and known limitations.

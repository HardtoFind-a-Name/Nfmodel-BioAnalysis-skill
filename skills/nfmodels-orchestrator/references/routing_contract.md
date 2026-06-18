# NFmodels Routing Contract

The top-level orchestrator emits a reviewable analysis plan first, then a route manifest. It must not emit an execution-ready plan before user review.

## Domains

- `transcriptome`: bulk RNA analysis handled by `transcriptome-analysis-orchestrator`.
- `scrna`: single-cell RNA analysis handled by `scrna-analysis-orchestrator`.

## Planning Artifacts

- `01.normalized_request.md`: normalized source request/report.
- `02.analysis_plan.md`: primary user-review artifact and execution contract.
- `03.route_manifest.json`: derived machine-readable route metadata.
- `04.review_status.json`: review gate status.

## Manifest Fields

- `project_id`: requested project identifier or `UNSET_PROJECT_ID`.
- `analysis_plan_path`: path to `02.analysis_plan.md`.
- `source`: request text or normalized file metadata.
- `environment_preflight`: profile and tidepy requirement for `nfmodels-environment-check`.
- `routes`: one entry per detected domain.
- `detected_paths`: absolute-looking paths found in the source text.
- `blockers`: unresolved issues that prevent direct dispatch.
- `review_status`: `pending` unless the user explicitly approves the analysis plan.
- `next_actions`: adapter skills and scripts to run next.

## Boundaries

The total controller should draft the reviewable plan, detect routes, and require environment preflight only after plan approval. It should not choose final model files, install packages, or launch Nextflow directly.

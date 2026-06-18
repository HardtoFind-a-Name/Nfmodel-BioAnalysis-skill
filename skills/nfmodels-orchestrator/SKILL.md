---
name: nfmodels-orchestrator
description: Generate a reviewable NFmodels analysis plan from a user request or .md/.docx report, listing executable analyses, required inputs, outputs, dependencies, pipeline gaps, and manual steps before routing to environment checks and transcriptome/scRNA adapters. Use when NFmodels work needs plan drafting, user review, routing, or dispatch planning before execution.
---

# NFmodels Orchestrator

Use this skill as the top-level controller for NFmodels analysis requests. It first creates a user-reviewable `02.analysis_plan.md`, then routes approved plans to environment checks and domain adapter skills.

## Role

- Accept either a plan/report file (`.md` or `.docx`) or a direct user request.
- Normalize the source into markdown-like text for review.
- Generate `02.analysis_plan.md` listing executable analyses, needed inputs, expected outputs, dependencies, scheduling, execution policy, non-executable steps, and pipeline-vs-plan gaps.
- Detect the requested analysis domain: bulk transcriptome, scRNA, or both.
- Produce a route manifest with adapter skill names, pipeline roots, detected inputs, environment preflight, blockers, and next actions.
- Do not execute Nextflow directly and do not allow environment/adapters to proceed before plan review is approved.

## Workflow

1. Capture the source request.
   - For `.md`, read directly.
   - For `.docx`, convert with `pandoc` if available.
   - For `.doc`, ask the user to convert to `.docx` or `.md`.
2. Run `scripts/route_nfmodels_plan.py` to create planning artifacts.
3. Review `02.analysis_plan.md`; this is the user-facing execution contract.
4. Use `03.route_manifest.json` only after the analysis plan is approved.
5. Run `nfmodels-environment-check` with the manifest profile before execution-oriented adapter work.
6. Dispatch to one or more adapter skills:
   - `transcriptome-analysis-orchestrator` for `pipelines/transcriptome_prognosis_pipeline`
   - `scrna-analysis-orchestrator` for `pipelines/scRNA_analysis_pipeline`
7. Let adapters generate concrete `run_pipeline.sh` commands and domain input blockers.

## Review Policy

- Treat `02.analysis_plan.md` as the primary planning artifact and review contract.
- Treat the route manifest as derived metadata until the user approves the analysis plan.
- Keep runtime binary paths out of skill instructions and generated commands; the pipeline launchers load `.env` (project root).
- Treat `nfmodels-environment-check` as the execution preflight gate after analysis plan approval and before adapter execution.
- Do not route ambiguous data paths silently. Mark them as blockers or assumptions in the manifest.

## Resources

- `scripts/route_nfmodels_plan.py`: normalize a file/request and emit `02.analysis_plan.md`, route manifest, and review status artifacts.
- `references/routing_contract.md`: route manifest fields and domain boundaries.

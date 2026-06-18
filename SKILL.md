---
name: nfmodels-analysis-project
description: Self-contained NFmodels analysis skill project bundling 5 Claude Code compatible skills (orchestrator, environment-check, transcriptome adapter, scRNA adapter, celltype annotator) with 2 Nextflow pipelines (transcriptome_prognosis_pipeline, scRNA_analysis_pipeline). Use when setting up or running NFmodels bioinformatics analysis workflows, generating analysis plans, checking environments, or executing bulk transcriptome and single-cell RNA pipelines.
---

# NFmodels Analysis Project

Self-contained, portable skill project that bundles NFmodels analysis skills and Nextflow pipelines into one directory. Deployable anywhere — no hardcoded paths to the old `NFmodels/` repository.

## Project Structure

```
nfmodels-analysis-project/
├── .env                         ← shared runtime configuration
├── SKILL.md                     ← this file (project meta-skill)
├── README.md                    ← full documentation
├── lib/
│   └── _resolver.py             ← unified path discovery
├── skills/                      ← 5 Claude Code compatible skills
│   ├── nfmodels-orchestrator/
│   ├── nfmodels-environment-check/
│   ├── transcriptome-analysis-orchestrator/
│   ├── scrna-analysis-orchestrator/
│   └── scrna-celltype-annotator/
└── pipelines/                   ← 2 Nextflow pipelines
    ├── transcriptome_prognosis_pipeline/
    └── scRNA_analysis_pipeline/
```

## Skill Inventory

| Skill | Purpose | Entry Point |
|---|---|---|
| `nfmodels-orchestrator` | Parse analysis requests, generate reviewable `02.analysis_plan.md`, route to domain adapters | `skills/nfmodels-orchestrator/scripts/route_nfmodels_plan.py` |
| `nfmodels-environment-check` | Validate `.env`, runtime binaries, R packages | `skills/nfmodels-environment-check/scripts/check_nfmodels_env.py` |
| `transcriptome-analysis-orchestrator` | Generate bulk transcriptome pipeline commands | `skills/transcriptome-analysis-orchestrator/scripts/generate_transcriptome_command.py` |
| `scrna-analysis-orchestrator` | Generate staged scRNA pipeline commands | `skills/scrna-analysis-orchestrator/scripts/generate_scrna_command.py` |
| `scrna-celltype-annotator` | Literature-backed cell type annotation for scRNA clusters | `skills/scrna-celltype-annotator/scripts/validate_mapping.py` |

## Path Resolution

All skill scripts use a unified resolver (`lib/_resolver.py`) to locate the project root:

1. `--nfmodels-dir` CLI argument (explicit)
2. `NFMODELS_ROOT` environment variable
3. Auto-detection from `lib/_resolver.py` location

No `parents[3]` magic numbers. Deploy this project anywhere.

## Quick Start

```bash
# Clone or copy this project
cd nfmodels-analysis-project

# Configure runtime
cp .env .env.local
# Edit .env.local with your Nextflow/R/Python paths
# Then set NFMODELS_ROOT to this directory:
export NFMODELS_ROOT=$(pwd)

# Check environment
python3 skills/nfmodels-environment-check/scripts/check_nfmodels_env.py --profile all

# Generate an analysis plan
python3 skills/nfmodels-orchestrator/scripts/route_nfmodels_plan.py \
  --request "Run DEG and multi-model prognosis on TCGA-LUAD" \
  --project-id demo_project \
  --out-dir /path/to/planning/

# Generate pipeline commands
python3 skills/transcriptome-analysis-orchestrator/scripts/generate_transcriptome_command.py \
  --project-id demo_project --train-id TCGA-LUAD
```

## Boundaries

- Skills reference each other by name (e.g., `nfmodels-orchestrator` dispatches to `nfmodels-environment-check`).
- Pipelines are launched via `pipelines/<name>/run_pipeline.sh` which sources `.env` from the project root.
- No `agents/openai.yaml` — this project targets Claude Code and other agents that use SKILL.md frontmatter.
- Runtime binary paths live in `.env` only, never hardcoded in scripts or skill instructions.

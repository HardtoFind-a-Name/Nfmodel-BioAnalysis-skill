---
name: nfmodels-environment-check
description: Check NFmodels runtime environment before pipeline execution. Use when validating .env (project root), diagnosing missing Nextflow/Python/Rscript/tidepy settings, checking transcriptome or scRNA R package dependencies, generating recommended environment variables, or running a preflight gate before nfmodels-orchestrator dispatches transcriptome/scRNA adapters.
---

# NFmodels Environment Check

Use this skill as the shared preflight gate for NFmodels pipelines. It can run standalone, and `nfmodels-orchestrator` should run it before dispatching execution to domain adapters.

## Role

- Check `.env` (project root) and required runtime variables.
- Recommend values for missing variables from a bounded search of PATH and common Conda env locations.
- Check R package availability with the configured `NFMODELS_RSCRIPT_CMD`; never use bare `Rscript`.
- Treat `NFMODELS_TIDEPY_BIN` as optional globally, but as required when TIDE analysis is planned.
- Report missing R packages without installing anything.

## Workflow

1. Run `scripts/check_nfmodels_env.py` with `--profile all`, `--profile transcriptome`, `--profile scrna`, or `--profile runtime`.
2. Review blockers and warnings.
3. If variables are missing, use `--write-recommended-env` to create `.env.recommended`, then manually review before changing `.env`.
4. If R packages are missing, use the grouped package report to install the required CRAN/Bioconductor/GitHub packages in the configured R environment.
5. Re-run the check before pipeline execution.

## Orchestrator Placement

- Let `nfmodels-orchestrator` call this skill after route review and before adapter execution.
- Keep transcriptome/scRNA adapters thin: they may request this preflight but should not duplicate environment-check logic.
- Run this skill directly when the user asks only to inspect or repair NFmodels runtime configuration.

## Resources

- `scripts/check_nfmodels_env.py`: environment and R package checker.
- `references/r_package_manifest.json`: static dependency manifest used by default.
- `references/environment_contract.md`: variable meanings and preflight policy.

# NFmodels Environment Contract

## Required Variables

- `NFMODELS_NEXTFLOW_BIN`: command or absolute path used by launchers to run Nextflow.
- `NFMODELS_RSCRIPT_CMD`: Conda-managed R command used inside Nextflow processes. Do not use bare `Rscript` in this workspace.
- `NFMODELS_PYTHON_CMD`: Python command used by helper modules.

## Optional Variables

- `NFMODELS_TIDEPY_BIN`: TIDE executable. Optional globally, required only for TIDE analysis.
- `SCRNA_NEXTFLOW_BIN`: optional scRNA-specific override; inherits `NFMODELS_NEXTFLOW_BIN` when absent.
- `SCRNA_RSCRIPT_CMD`: optional scRNA-specific override; inherits `NFMODELS_RSCRIPT_CMD` when absent.

## Policy

The checker reports missing configuration and writes `.env.recommended` only when requested. It must not overwrite `.env` or install packages automatically.

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ENV_FILE="${PROJECT_DIR}/.env"

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "NFmodels environment file not found: ${ENV_FILE}" >&2
  exit 1
fi

source "${ENV_FILE}"

NEXTFLOW_CMD="${SCRNA_NEXTFLOW_BIN:-${NFMODELS_NEXTFLOW_BIN:-}}"
if [[ -z "${NEXTFLOW_CMD}" ]]; then
  echo "Set NFMODELS_NEXTFLOW_BIN or SCRNA_NEXTFLOW_BIN in ${ENV_FILE}" >&2
  exit 1
fi

if [[ -z "${SCRNA_RSCRIPT_CMD:-${NFMODELS_RSCRIPT_CMD:-}}" ]]; then
  echo "Set NFMODELS_RSCRIPT_CMD or SCRNA_RSCRIPT_CMD in ${ENV_FILE}" >&2
  exit 1
fi

cd "${SCRIPT_DIR}"
exec ${NEXTFLOW_CMD} run main.nf "$@"

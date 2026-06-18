#!/usr/bin/env python3
"""
Unified path resolution for NFmodels Analysis Project.

Discovery chain (highest priority first):
    1. --nfmodels-dir CLI argument (explicit)
    2. NFMODELS_ROOT environment variable
    3. Auto-detect from lib/_resolver.py location (project root = lib/../)
    4. Fallback: current working directory
"""

from __future__ import annotations

import os
from pathlib import Path

_ENV_VAR = "NFMODELS_ROOT"


def _detect_from_lib() -> Path | None:
    """Detect project root from this file's location."""
    return Path(__file__).resolve().parent.parent


def resolve_project_root(cli_arg: str | None = None) -> Path:
    """Resolve the NFmodels Analysis Project root directory.

    Args:
        cli_arg: Optional --nfmodels-dir value from argparse.

    Returns:
        Absolute Path to the project root.

    Raises:
        FileNotFoundError: If root cannot be located and .env is missing.
    """
    if cli_arg:
        root = Path(cli_arg).expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"Project root does not exist: {root}")
        return root

    env_val = os.environ.get(_ENV_VAR)
    if env_val:
        root = Path(env_val).expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"NFMODELS_ROOT does not exist: {root}")
        return root

    detected = _detect_from_lib()
    if detected is not None and detected.is_dir():
        return detected

    cwd = Path.cwd().resolve()
    if (cwd / ".env").exists() or (cwd / "pipelines").is_dir():
        return cwd

    raise FileNotFoundError(
        "Cannot locate NFmodels Analysis Project root. "
        f"Set {_ENV_VAR} environment variable or pass --nfmodels-dir explicitly."
    )


def resolve_env(project_root: Path) -> Path:
    """Return the .env file path."""
    env_file = project_root / ".env"
    if not env_file.exists():
        raise FileNotFoundError(f".env not found in {project_root}. "
                                "Create one with NFMODELS_NEXTFLOW_BIN, NFMODELS_RSCRIPT_CMD, etc.")
    return env_file


def resolve_pipeline_root(project_root: Path, domain: str) -> Path:
    """Return the pipeline directory for a given domain.

    Args:
        project_root: Resolved project root directory.
        domain: 'transcriptome' or 'scrna'.

    Returns:
        Absolute Path to the pipeline directory.

    Raises:
        ValueError: If domain is unknown.
        FileNotFoundError: If pipeline directory does not exist.
    """
    pipelines = {
        "transcriptome": project_root / "pipelines" / "transcriptome_prognosis_pipeline",
        "scrna": project_root / "pipelines" / "scRNA_analysis_pipeline",
    }
    if domain not in pipelines:
        raise ValueError(f"Unknown domain: {domain!r}. Use 'transcriptome' or 'scrna'.")
    p = pipelines[domain]
    if not p.is_dir():
        raise FileNotFoundError(f"Pipeline not found: {p}")
    return p


def resolve_launcher(project_root: Path, domain: str) -> Path:
    """Return the run_pipeline.sh path for a given domain."""
    return resolve_pipeline_root(project_root, domain) / "run_pipeline.sh"

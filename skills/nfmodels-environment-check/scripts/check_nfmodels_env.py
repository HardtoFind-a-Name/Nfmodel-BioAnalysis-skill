#!/usr/bin/env python3
"""Check NFmodels runtime variables and R package dependencies."""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

NFMODELS_DIR = Path(__file__).resolve().parents[3]
# Override NFMODELS_DIR via _resolver when available
try:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "lib"))
    from _resolver import resolve_project_root  # noqa: E402
    NFMODELS_DIR = resolve_project_root()
except ImportError:
    pass
MANIFEST_PATH = Path(__file__).resolve().parents[1] / "references" / "r_package_manifest.json"

REQUIRED_VARS = {
    "NFMODELS_NEXTFLOW_BIN": "Nextflow launcher command",
    "NFMODELS_RSCRIPT_CMD": "Conda-managed Rscript command",
    "NFMODELS_PYTHON_CMD": "Python helper command",
}
OPTIONAL_VARS = {
    "NFMODELS_TIDEPY_BIN": "Optional tidepy command for TIDE analysis",
}
ALIAS_VARS = {
    "SCRNA_NEXTFLOW_BIN": "Optional scRNA Nextflow override",
    "SCRNA_RSCRIPT_CMD": "Optional scRNA Rscript override",
}
KNOWN_ENV_ROOT = Path("/data/nas2/software/miniconda3/envs")
KNOWN_CONDA = Path("/data/nas2/software/miniconda3/bin/conda")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def dump_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def source_env_file(env_file: Path) -> tuple[dict[str, str], list[str]]:
    if not env_file.exists():
        return {}, [f"Environment file not found: {env_file}"]
    cmd = f"set -a; source {shlex.quote(str(env_file))}; env -0"
    result = subprocess.run(["bash", "-lc", cmd], capture_output=True)
    if result.returncode != 0:
        stderr = result.stderr.decode(errors="replace").strip()
        return {}, [f"Failed to source {env_file}: {stderr}"]
    env: dict[str, str] = {}
    for item in result.stdout.split(b"\0"):
        if not item or b"=" not in item:
            continue
        key, value = item.split(b"=", 1)
        env[key.decode(errors="replace")] = value.decode(errors="replace")
    return env, []


def first_token(command: str) -> str:
    try:
        parts = shlex.split(command)
    except ValueError:
        return ""
    return parts[0] if parts else ""


def command_available(command: str, env: dict[str, str]) -> tuple[bool, str | None]:
    token = first_token(command)
    if not token:
        return False, None
    if os.path.isabs(token):
        path = Path(token)
        return path.exists() and os.access(path, os.X_OK), str(path) if path.exists() else None
    found = shutil.which(token, path=env.get("PATH") or os.environ.get("PATH"))
    return bool(found), found


def unique(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item and item not in seen:
            seen.add(item)
            out.append(item)
    return out


def bin_candidates(name: str) -> list[str]:
    candidates: list[str] = []
    found = shutil.which(name)
    if found:
        candidates.append(found)
    if KNOWN_ENV_ROOT.exists():
        for path in sorted(KNOWN_ENV_ROOT.glob(f"*/bin/{name}")):
            if path.exists() and os.access(path, os.X_OK):
                candidates.append(str(path))
    return unique(candidates)


def recommend_var(var_name: str) -> str | None:
    if var_name == "NFMODELS_NEXTFLOW_BIN":
        preferred = KNOWN_ENV_ROOT / "nextflow" / "bin" / "nextflow"
        if preferred.exists():
            return str(preferred)
        candidates = bin_candidates("nextflow")
        return candidates[0] if candidates else None
    if var_name == "NFMODELS_RSCRIPT_CMD":
        r433 = KNOWN_ENV_ROOT / "R4.3.3" / "bin" / "Rscript"
        if r433.exists() and (KNOWN_CONDA.exists() or shutil.which("conda")):
            return "conda run -n R4.3.3 Rscript"
        candidates = bin_candidates("Rscript")
        if candidates:
            return candidates[0]
        return None
    if var_name == "NFMODELS_PYTHON_CMD":
        found = shutil.which("python3") or shutil.which("python")
        return "python3" if shutil.which("python3") else found
    if var_name == "NFMODELS_TIDEPY_BIN":
        preferred = KNOWN_ENV_ROOT / "InPETM" / "bin" / "tidepy"
        if preferred.exists():
            return str(preferred)
        candidates = bin_candidates("tidepy")
        return candidates[0] if candidates else None
    return None


def variable_report(env: dict[str, str], require_tidepy: bool) -> tuple[list[dict], list[str], list[str]]:
    rows: list[dict] = []
    blockers: list[str] = []
    warnings: list[str] = []
    all_vars = {**REQUIRED_VARS, **OPTIONAL_VARS, **ALIAS_VARS}
    for name, description in all_vars.items():
        value = env.get(name, "")
        required = name in REQUIRED_VARS or (name == "NFMODELS_TIDEPY_BIN" and require_tidepy)
        available = False
        resolved = None
        if value:
            available, resolved = command_available(value, env)
        recommendation = recommend_var(name)
        status = "ok" if value and available else "missing" if not value else "unavailable"
        if required and status != "ok":
            blockers.append(f"{name} is {status}.")
        elif status != "ok" and name in OPTIONAL_VARS:
            warnings.append(f"{name} is {status}; required only for TIDE analysis.")
        rows.append({
            "name": name,
            "description": description,
            "required": required,
            "configured": value,
            "status": status,
            "resolved_path": resolved,
            "recommendation": recommendation,
        })
    if not env.get("SCRNA_NEXTFLOW_BIN") and env.get("NFMODELS_NEXTFLOW_BIN"):
        warnings.append("SCRNA_NEXTFLOW_BIN is unset and will inherit NFMODELS_NEXTFLOW_BIN.")
    if not env.get("SCRNA_RSCRIPT_CMD") and env.get("NFMODELS_RSCRIPT_CMD"):
        warnings.append("SCRNA_RSCRIPT_CMD is unset and will inherit NFMODELS_RSCRIPT_CMD.")
    return rows, blockers, warnings


def packages_for_profile(manifest: dict, profile: str) -> list[dict]:
    group_names = manifest["profiles"].get(profile, [])
    by_package: dict[str, dict] = {}
    for group in group_names:
        for item in manifest["groups"].get(group, []):
            pkg = item["package"]
            entry = by_package.setdefault(pkg, {"package": pkg, "source": item.get("source", "unknown"), "groups": []})
            entry["groups"].append(group)
            if entry["source"] == "unknown" and item.get("source"):
                entry["source"] = item["source"]
    return [by_package[name] for name in sorted(by_package, key=str.lower)]


def scan_r_scripts(nfmodels_dir: Path) -> list[dict]:
    roots = [nfmodels_dir / "pipelines" / "transcriptome_prognosis_pipeline" / "scripts", nfmodels_dir / "pipelines" / "scRNA_analysis_pipeline" / "scripts"]
    patterns = [
        re.compile(r"\b(?:library|require|requireNamespace)\s*\(\s*[\"']?([A-Za-z][A-Za-z0-9_.]*)"),
        re.compile(r"\b([A-Za-z][A-Za-z0-9_.]*)::"),
    ]
    ignore = {"base", "stats", "utils", "graphics", "grDevices", "methods", "tools", "parallel", "grid", "pkg"}
    found: set[str] = set()
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*.R"):
            text = path.read_text(encoding="utf-8", errors="ignore")
            for pattern in patterns:
                for match in pattern.finditer(text):
                    package = match.group(1)
                    if package not in ignore:
                        found.add(package)
    return [{"package": package, "source": "detected_by_scan", "groups": ["scan_r_scripts"]} for package in sorted(found, key=str.lower)]


def check_r_packages(rscript_cmd: str, packages: list[dict], env: dict[str, str]) -> tuple[list[dict], list[str]]:
    if not packages:
        return [], []
    package_names = [item["package"] for item in packages]
    r_code = "\n".join([
        'pkgs <- strsplit(Sys.getenv("NFMODELS_CHECK_PACKAGES"), ",", fixed = TRUE)[[1]]',
        'pkgs <- pkgs[nzchar(pkgs)]',
        'versions <- installed.packages()[, "Version"]',
        'for (pkg in pkgs) {',
        '  version <- if (pkg %in% names(versions)) versions[[pkg]] else ""',
        '  cat(pkg, "\\t", version, "\\n", sep = "")',
        '}',
    ])
    run_env = os.environ.copy()
    run_env.update(env)
    run_env["NFMODELS_CHECK_PACKAGES"] = ",".join(package_names)
    command = f"{rscript_cmd} -e {shlex.quote(r_code)}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True, env=run_env)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        return [], [f"R package check failed with configured NFMODELS_RSCRIPT_CMD: {stderr}"]
    installed = {}
    for line in result.stdout.splitlines():
        if "\t" not in line:
            continue
        package, version = line.split("\t", 1)
        installed[package] = version.strip()
    rows: list[dict] = []
    meta = {item["package"]: item for item in packages}
    for package in package_names:
        version = installed.get(package, "")
        rows.append({
            "package": package,
            "installed": bool(version),
            "version": version or None,
            "source": meta[package].get("source", "unknown"),
            "groups": meta[package].get("groups", []),
        })
    return rows, []


def write_recommended_env(path: Path, var_rows: list[dict]) -> None:
    values = {row["name"]: row["configured"] or row["recommendation"] for row in var_rows}
    lines = [
        "# Recommended NFmodels runtime environment generated by nfmodels-environment-check.",
        "# Review this file before copying values into .env (project root).",
        "",
    ]
    ordered = ["NFMODELS_NEXTFLOW_BIN", "NFMODELS_RSCRIPT_CMD", "NFMODELS_PYTHON_CMD", "NFMODELS_TIDEPY_BIN"]
    for name in ordered:
        if values.get(name):
            lines.append(f"export {name}={shlex.quote(values[name])}")
    if values.get("NFMODELS_NEXTFLOW_BIN"):
        lines.append('export SCRNA_NEXTFLOW_BIN="${SCRNA_NEXTFLOW_BIN:-${NFMODELS_NEXTFLOW_BIN}}"')
    if values.get("NFMODELS_RSCRIPT_CMD"):
        lines.append('export SCRNA_RSCRIPT_CMD="${SCRNA_RSCRIPT_CMD:-${NFMODELS_RSCRIPT_CMD}}"')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_report(args: argparse.Namespace) -> dict:
    nfmodels_dir = Path(args.nfmodels_dir).resolve()
    env_file = Path(args.env_file).resolve() if args.env_file else nfmodels_dir / ".env"
    env, env_warnings = source_env_file(env_file)
    var_rows, blockers, warnings = variable_report(env, args.require_tidepy)
    warnings = env_warnings + warnings

    manifest = load_json(Path(args.manifest))
    package_requests = packages_for_profile(manifest, args.profile)
    if args.scan_r_scripts:
        scanned = scan_r_scripts(nfmodels_dir)
        existing = {item["package"] for item in package_requests}
        package_requests.extend([item for item in scanned if item["package"] not in existing])
        warnings.append("Experimental --scan-r-scripts was used; prefer the static manifest for production checks.")

    package_rows: list[dict] = []
    package_warnings: list[str] = []
    if args.skip_r_packages or args.profile == "runtime":
        package_warnings.append("R package check skipped.")
    else:
        rscript_cmd = env.get("NFMODELS_RSCRIPT_CMD", "")
        if not rscript_cmd:
            package_warnings.append("R package check skipped because NFMODELS_RSCRIPT_CMD is not configured.")
        else:
            package_rows, package_warnings = check_r_packages(rscript_cmd, package_requests, env)
    warnings.extend(package_warnings)

    missing_packages = [row for row in package_rows if not row["installed"]]
    if missing_packages:
        blockers.append(f"Missing {len(missing_packages)} R package(s) for profile {args.profile}.")

    recommended_env = None
    if args.write_recommended_env:
        recommended_env = str((nfmodels_dir / ".env.recommended").resolve())
        write_recommended_env(Path(recommended_env), var_rows)

    status = "ok" if not blockers else "blocked"
    return {
        "status": status,
        "profile": args.profile,
        "nfmodels_dir": str(nfmodels_dir),
        "env_file": str(env_file),
        "variables": var_rows,
        "blockers": blockers,
        "warnings": warnings,
        "r_packages": {
            "manifest": str(Path(args.manifest).resolve()),
            "checked": bool(package_rows),
            "requested_count": len(package_requests),
            "missing_count": len(missing_packages),
            "missing": missing_packages,
            "installed": [row for row in package_rows if row["installed"]],
        },
        "recommended_env_file": recommended_env,
    }


def print_summary(report: dict) -> None:
    print(f"NFmodels environment check: {report['status']}")
    print(f"Profile: {report['profile']}")
    print(f"Environment file: {report['env_file']}")
    print("")
    print("Runtime variables:")
    for row in report["variables"]:
        marker = "OK" if row["status"] == "ok" else "MISSING" if row["status"] == "missing" else "UNAVAILABLE"
        req = "required" if row["required"] else "optional"
        print(f"- {marker} {row['name']} ({req})")
        if row["configured"]:
            print(f"  configured: {row['configured']}")
        if row["recommendation"] and row["status"] != "ok":
            print(f"  recommendation: {row['recommendation']}")
    if report["blockers"]:
        print("\nBlockers:")
        for item in report["blockers"]:
            print(f"- {item}")
    if report["warnings"]:
        print("\nWarnings:")
        for item in report["warnings"]:
            print(f"- {item}")
    pkg = report["r_packages"]
    if pkg["checked"]:
        print(f"\nR packages: {pkg['requested_count']} checked, {pkg['missing_count']} missing")
        if pkg["missing"]:
            by_source: dict[str, list[str]] = {}
            for row in pkg["missing"]:
                by_source.setdefault(row["source"], []).append(row["package"])
            for source, packages in sorted(by_source.items()):
                print(f"- {source}: {', '.join(sorted(packages, key=str.lower))}")
    if report.get("recommended_env_file"):
        print(f"\nRecommended env written: {report['recommended_env_file']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nfmodels-dir", default=str(NFMODELS_DIR), help="NFmodels project directory")
    parser.add_argument("--env-file", help="Environment file; default: <nfmodels-dir>/.env")
    parser.add_argument("--manifest", default=str(MANIFEST_PATH), help="Static R package manifest JSON")
    parser.add_argument("--profile", choices=["runtime", "transcriptome", "scrna", "all"], default="all")
    parser.add_argument("--require-tidepy", action="store_true", help="Treat NFMODELS_TIDEPY_BIN as required")
    parser.add_argument("--skip-r-packages", action="store_true", help="Skip installed.packages() check")
    parser.add_argument("--scan-r-scripts", action="store_true", help="Experimental: augment package list by scanning R scripts")
    parser.add_argument("--write-recommended-env", action="store_true", help="Write <nfmodels-dir>/.env.recommended")
    parser.add_argument("--json-output", help="Optional JSON report output")
    args = parser.parse_args()

    report = build_report(args)
    if args.json_output:
        dump_json(Path(args.json_output), report)
    print_summary(report)
    sys.exit(0 if report["status"] == "ok" else 1)


if __name__ == "__main__":
    main()

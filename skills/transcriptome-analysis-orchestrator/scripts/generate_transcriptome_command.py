#!/usr/bin/env python3
"""Generate a checked command for NFmodels bulk transcriptome pipeline."""

from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "lib"))
from _resolver import resolve_project_root, resolve_pipeline_root  # noqa: E402

ROOT = resolve_project_root()
PIPELINE_ROOT = resolve_pipeline_root(ROOT, "transcriptome")

MAIN_ANALYSES = {"deg", "candidate_genes", "go_kegg_ppi", "multi_model"}
POST_FLAG_MAP = {
    "stage": "run_stage",
    "prognostic": "run_stage",
    "nomogram": "run_stage",
    "gsea_gene": "run_gsea_gene",
    "gsea_risk_nes": "run_gsea_risk_nes",
    "cibersort": "run_cibersort",
    "immune": "run_cibersort",
    "immune_infiltration": "run_cibersort",
    "cnv": "run_cnv",
    "gistic": "run_gistic",
    "tide": "run_tide",
    "ips": "run_ips",
    "diff_expr": "run_diff_expr",
    "circos": "run_circos",
    "tmb": "run_tmb",
}
ALL_POST_FLAGS = sorted(set(POST_FLAG_MAP.values()))


def shell(parts: list[str]) -> str:
    return shlex.join([str(part) for part in parts])


def split_tokens(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip().lower().replace("-", "_") for item in value.split(",") if item.strip()]


def expand_analyses(tokens: list[str]) -> set[str]:
    if not tokens:
        return set(MAIN_ANALYSES)
    expanded: set[str] = set()
    for token in tokens:
        if token == "main":
            expanded.update(MAIN_ANALYSES)
        elif token == "post":
            expanded.update(POST_FLAG_MAP)
        elif token == "all":
            expanded.update(MAIN_ANALYSES)
            expanded.update(POST_FLAG_MAP)
        else:
            expanded.add(token)
    return expanded


def build(args: argparse.Namespace) -> dict:
    analyses = expand_analyses(split_tokens(args.analyses))
    flags = {flag: "false" for flag in ["run_deg", "run_go_kegg_ppi", "run_multi_model", *ALL_POST_FLAGS]}
    flags["run_deg"] = "true" if "deg" in analyses else "false"
    flags["run_go_kegg_ppi"] = "true" if "go_kegg_ppi" in analyses else "false"
    flags["run_multi_model"] = "true" if "multi_model" in analyses else "false"
    for analysis in analyses:
        flag = POST_FLAG_MAP.get(analysis)
        if flag:
            flags[flag] = "true"

    blockers: list[str] = []
    if not args.project_id:
        blockers.append("Missing --project-id.")
    if not args.train_id:
        blockers.append("Missing --train-id.")
    if flags["run_multi_model"] == "true" and not (args.validation_ids or args.validation_sheet):
        blockers.append("run_multi_model requires --validation-ids or --validation-sheet.")
    if flags["run_go_kegg_ppi"] == "true" and not args.string_interactions_file:
        blockers.append("run_go_kegg_ppi requires --string-interactions-file.")

    post_enabled = any(flags[flag] == "true" for flag in ALL_POST_FLAGS)
    if post_enabled and not args.risk_file:
        blockers.append("Post-model analyses require --risk-file.")
    gene_dependent_flags = {"run_gsea_gene", "run_cibersort", "run_cnv", "run_diff_expr", "run_circos"}
    if any(flags[flag] == "true" for flag in gene_dependent_flags) and not args.gene_file:
        blockers.append("Selected post-model analyses require --gene-file.")

    candidate_enabled = "candidate_genes" in analyses
    if candidate_enabled and not args.candidate_genes and args.candidate_route == "intersect" and not args.gene_sets_sheet:
        blockers.append("Intersect candidate gene route requires --gene-sets-sheet unless --candidate-genes is supplied.")

    parts = [str(PIPELINE_ROOT / "run_pipeline.sh"), "-profile", args.profile]
    parts.extend(["--project_id", args.project_id or "UNSET_PROJECT_ID"])
    parts.extend(["--train_id", args.train_id or "UNSET_TRAIN_ID"])
    parts.extend(["--expr_type", args.expr_type])

    for flag in ["run_deg", "run_go_kegg_ppi", "run_multi_model", *ALL_POST_FLAGS]:
        parts.extend([f"--{flag}", flags[flag]])

    if args.validation_ids:
        parts.extend(["--validation_ids", args.validation_ids])
    if args.validation_sheet:
        parts.extend(["--validation_sheet", args.validation_sheet])
    if args.string_interactions_file:
        parts.extend(["--string_interactions_file", args.string_interactions_file])

    if args.candidate_genes:
        parts.extend(["--candidate_genes", args.candidate_genes, "--run_candidate_genes", "false", "--run_intersect_candidate_genes", "false"])
    elif candidate_enabled:
        if args.candidate_route == "intersect":
            parts.extend(["--run_candidate_genes", "false", "--run_intersect_candidate_genes", "true"])
            if args.gene_sets_sheet:
                parts.extend(["--gene_sets_sheet", args.gene_sets_sheet])
        else:
            parts.extend(["--run_candidate_genes", "true", "--run_intersect_candidate_genes", "false"])
            if args.metascape_input_file:
                parts.extend(["--metascape_input_file", args.metascape_input_file])
    else:
        parts.extend(["--run_candidate_genes", "false", "--run_intersect_candidate_genes", "false"])

    if args.risk_file:
        parts.extend(["--risk_file", args.risk_file])
    if args.gene_file:
        parts.extend(["--gene_file", args.gene_file])

    return {
        "pipeline_root": str(PIPELINE_ROOT),
        "launcher": str(PIPELINE_ROOT / "run_pipeline.sh"),
        "analyses": sorted(analyses),
        "blockers": blockers,
        "command": shell(parts),
        "uses_project_env": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", dest="project_id", required=True)
    parser.add_argument("--train-id", dest="train_id", required=True)
    parser.add_argument("--validation-ids", dest="validation_ids")
    parser.add_argument("--validation-sheet", dest="validation_sheet")
    parser.add_argument("--expr-type", dest="expr_type", default="fpkm", choices=["fpkm", "tpm"])
    parser.add_argument("--profile", default="local")
    parser.add_argument("--analyses", default="main", help="Comma list: main, all, deg, candidate_genes, go_kegg_ppi, multi_model, cibersort, tide, etc.")
    parser.add_argument("--candidate-route", choices=["intersect", "gsea"], default="intersect")
    parser.add_argument("--candidate-genes", dest="candidate_genes")
    parser.add_argument("--gene-sets-sheet", dest="gene_sets_sheet")
    parser.add_argument("--metascape-input-file", dest="metascape_input_file")
    parser.add_argument("--string-interactions-file", dest="string_interactions_file")
    parser.add_argument("--risk-file", dest="risk_file")
    parser.add_argument("--gene-file", dest="gene_file")
    parser.add_argument("--output", help="Optional JSON output path")
    args = parser.parse_args()

    payload = build(args)
    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()

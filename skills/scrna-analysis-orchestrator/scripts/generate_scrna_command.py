#!/usr/bin/env python3
"""Generate checked commands and annotation handoff metadata for NFmodels scRNA pipeline."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "lib"))
from _resolver import resolve_project_root, resolve_pipeline_root  # noqa: E402

ROOT = resolve_project_root()
PIPELINE_ROOT = resolve_pipeline_root(ROOT, "scrna")
ANNOTATOR_SKILL = ROOT / "skills" / "scrna-celltype-annotator"


def shell(parts: list[str]) -> str:
    return shlex.join([str(part) for part in parts])


def is_set(value: str | None) -> bool:
    return bool(value and str(value).strip() and str(value).strip().lower() not in {"none", "null"})


def default_run_dir(project_id: str) -> Path:
    user = os.environ.get("USER") or "default_user"
    return Path("/data/nas1") / user / "project" / project_id


def bool_cli(value: bool) -> str:
    return "true" if value else "false"


def stage_name(args: argparse.Namespace, downstream_requires_subset_apply: bool, subset_requested: bool) -> str:
    if args.run_subset_apply or downstream_requires_subset_apply:
        return "subset_apply_downstream"
    if args.mapping_file and subset_requested:
        return "main_apply_subset_prepare"
    if args.mapping_file:
        return "main_annotation_apply"
    return "main_annotation_prepare"


def expected_paths(run_dir: Path) -> dict:
    results = run_dir / "results"
    main_prepare_dir = results / "03_scrna_annotation_prepare"
    main_agent_dir = run_dir / "annotation_agent" / "main"
    subset_prepare_root = results / "07b_scrna_subset_annotation_prepare"
    subset_agent_root = run_dir / "annotation_agent" / "subset"
    return {
        "run_dir": str(run_dir),
        "results_dir": str(results),
        "main_annotation_prepare_dir": str(main_prepare_dir),
        "main_agent_output_dir": str(main_agent_dir),
        "main_annotation_inputs": {
            "mapping_template": str(main_prepare_dir / "08_cluster_celltype_mapping_template.csv"),
            "literature_marker_template": str(main_prepare_dir / "09_literature_marker_reference_template.csv"),
            "agent_context": str(main_prepare_dir / "10_annotation_agent_context.json"),
            "marker_table": str(main_prepare_dir / "02_all_markers.csv"),
            "top10_marker_table": str(main_prepare_dir / "03_top10_markers_by_cluster.csv"),
            "top20_marker_table": str(main_prepare_dir / "04_top20_markers_by_cluster.csv"),
        },
        "main_annotation_outputs": {
            "mapping_file": str(main_agent_dir / "08_cluster_celltype_mapping_filled.csv"),
            "marker_reference_file": str(main_agent_dir / "09_literature_marker_reference_filled.csv"),
            "evidence_md": str(main_agent_dir / "annotation_literature_evidence.md"),
            "decision_log": str(main_agent_dir / "annotation_decision_log.csv"),
            "validation_report": str(main_agent_dir / "annotation_validation_report.csv"),
        },
        "subset_annotation_prepare_root": str(subset_prepare_root),
        "subset_agent_output_root": str(subset_agent_root),
        "subset_mapping_manifest": str(subset_agent_root / "subset_mapping_manifest.csv"),
    }


def annotator_actions(paths: dict, args: argparse.Namespace) -> list[dict]:
    validate_script = str(ANNOTATOR_SKILL / "scripts" / "validate_mapping.py")
    main_inputs = paths["main_annotation_inputs"]
    main_outputs = paths["main_annotation_outputs"]
    actions = [
        {
            "gate": "main_annotation",
            "skill": "scrna-celltype-annotator",
            "status": "waiting_for_agent_and_user_review" if not args.mapping_file else "supplied_or_completed",
            "input_dir": paths["main_annotation_prepare_dir"],
            "expected_inputs": main_inputs,
            "expected_outputs": main_outputs,
            "validation_command": shell([
                "python", validate_script,
                "--template", main_inputs["mapping_template"],
                "--mapping", main_outputs["mapping_file"],
                "--references", main_outputs["marker_reference_file"],
                "--min-support", str(args.annotation_min_support_markers),
                "--out", main_outputs["validation_report"],
                *( ["--allow-low-support"] if args.annotation_allow_low_support else [] ),
            ]),
        }
    ]
    if args.run_subset_prepare or args.run_subset_apply or any([
        args.run_pseudotime, args.run_scmetabolism, args.run_scmetabolism_specific, args.run_gsva, args.run_reactome_gsa
    ]):
        actions.append({
            "gate": "subset_annotation",
            "skill": "scrna-celltype-annotator",
            "status": "waiting_for_subset_prepare_outputs" if not (args.subset_mapping_file or args.subset_mapping_manifest) else "supplied_or_completed",
            "input_root": paths["subset_annotation_prepare_root"],
            "expected_per_subset_inputs": {
                "mapping_template": "<subset_dir>/08_subset_cluster_celltype_mapping_template.csv",
                "literature_marker_template": "<subset_dir>/09_subset_literature_marker_reference_template.csv",
                "agent_context": "<subset_dir>/10_subset_annotation_agent_context.json",
                "marker_table": "<subset_dir>/02_subset_all_markers.csv",
            },
            "expected_per_subset_outputs": {
                "mapping_file": "<subset_agent_output_root>/<safe_label>/08_subset_cluster_celltype_mapping_filled.csv",
                "marker_reference_file": "<subset_agent_output_root>/<safe_label>/09_subset_literature_marker_reference_filled.csv",
                "evidence_md": "<subset_agent_output_root>/<safe_label>/annotation_literature_evidence.md",
                "decision_log": "<subset_agent_output_root>/<safe_label>/annotation_decision_log.csv",
                "validation_report": "<subset_agent_output_root>/<safe_label>/annotation_validation_report.csv",
            },
            "manifest_contract": {
                "path": paths["subset_mapping_manifest"],
                "columns": ["key_celltype", "safe_label", "mapping_file"],
            },
        })
    return actions


def build(args: argparse.Namespace) -> dict:
    blockers: list[str] = []
    notes: list[str] = []
    if not args.project_id:
        blockers.append("Missing --project-id.")
    if not args.input_rds and not args.run_raw_import:
        blockers.append("Provide --input-rds or enable --run-raw-import with raw import settings.")
    if args.target_gene_file and args.target_gene_file.lower() not in {"none", "null"} and not args.target_gene_col:
        blockers.append("--target-gene-col is required when --target-gene-file is supplied.")
    if args.run_raw_import and not args.scrna_cohort_id and not args.scrna_raw_umi_file:
        blockers.append("Raw import requires --scrna-cohort-id or --scrna-raw-umi-file.")

    downstream_requires_subset_apply = any([
        args.run_pseudotime,
        args.run_scmetabolism,
        args.run_scmetabolism_specific,
        args.run_gsva,
        args.run_reactome_gsa,
    ])
    subset_requested = args.run_subset_prepare or args.run_subset_apply or downstream_requires_subset_apply
    if subset_requested and not args.mapping_file:
        blockers.append("Subset prepare/apply/downstream stages require --mapping-file from main annotation review.")
    if (args.run_subset_apply or downstream_requires_subset_apply) and not (args.subset_mapping_file or args.subset_mapping_manifest):
        blockers.append("Subset apply/downstream stages require --subset-mapping-file or --subset-mapping-manifest from subset annotation review.")
    if args.run_gsva and not args.gsva_gmt_files:
        blockers.append("--run-gsva requires --gsva-gmt-files.")

    stage = stage_name(args, downstream_requires_subset_apply, subset_requested)
    if stage == "main_annotation_prepare":
        notes.append("mapping_file is absent; NF should stop after SCRNA_ANNOTATION_PREPARE and hand off templates to scrna-celltype-annotator.")
    if args.mapping_file and not args.annotation_marker_reference_file:
        notes.append("annotation_marker_reference_file is absent; NF marker-support audit will use expressed mapping markers without literature-reference filtering.")
    if not args.target_gene_file or args.target_gene_file.lower() in {"none", "null"}:
        notes.append("target_gene_file is absent; target-gene overlays and module-score outputs should be skipped.")

    run_dir = Path(args.run_dir).resolve() if args.run_dir else default_run_dir(args.project_id or "UNSET_PROJECT_ID")
    paths = expected_paths(run_dir)

    parts = [str(PIPELINE_ROOT / "run_pipeline.sh"), "-profile", args.profile]
    parts.extend(["--project_id", args.project_id or "UNSET_PROJECT_ID"])
    parts.extend(["--input_format", args.input_format])

    if args.input_rds:
        parts.extend(["--input_rds", args.input_rds])
    if args.run_raw_import:
        parts.extend(["--run_raw_import", "true"])
    if args.scrna_cohort_id:
        parts.extend(["--scRNA_cohort_id", args.scrna_cohort_id])
    if args.scrna_raw_umi_file:
        parts.extend(["--scRNA_raw_umi_file", args.scrna_raw_umi_file])
    if args.disease_name:
        parts.extend(["--disease_name", args.disease_name])
    if args.cancer_name:
        parts.extend(["--cancer_name", args.cancer_name])
    if args.mapping_file:
        parts.extend(["--mapping_file", args.mapping_file])
    if args.annotation_marker_reference_file:
        parts.extend(["--annotation_marker_reference_file", args.annotation_marker_reference_file])
    parts.extend(["--annotation_min_support_markers", str(args.annotation_min_support_markers)])
    if args.annotation_allow_low_support:
        parts.extend(["--annotation_allow_low_support", "true"])
    if args.target_gene_file and args.target_gene_file.lower() not in {"none", "null"}:
        parts.extend(["--target_gene_file", args.target_gene_file, "--target_gene_col", args.target_gene_col])
    if args.run_cellchat:
        parts.extend(["--run_cellchat", "true"])
    if args.key_celltypes:
        parts.extend(["--key_celltypes", args.key_celltypes])
    if args.run_subset_prepare:
        parts.extend(["--run_subset_prepare", "true"])
    if args.run_subset_apply:
        parts.extend(["--run_subset_apply", "true"])
    if args.subset_mapping_file:
        parts.extend(["--subset_mapping_file", args.subset_mapping_file])
    if args.subset_mapping_manifest:
        parts.extend(["--subset_mapping_manifest", args.subset_mapping_manifest])
    downstream_flags = {
        "run_pseudotime": args.run_pseudotime,
        "run_scmetabolism": args.run_scmetabolism,
        "run_scmetabolism_specific": args.run_scmetabolism_specific,
        "run_gsva": args.run_gsva,
        "run_reactome_gsa": args.run_reactome_gsa,
    }
    for flag, enabled in downstream_flags.items():
        if enabled:
            parts.extend([f"--{flag}", "true"])
    if args.gsva_gmt_files:
        parts.extend(["--gsva_gmt_files", args.gsva_gmt_files])

    command = shell(parts)
    review_required = stage in {"main_annotation_prepare", "main_apply_subset_prepare"}
    next_gate = "user_review_main_mapping" if stage == "main_annotation_prepare" else (
        "user_review_subset_mapping" if stage == "main_apply_subset_prepare" else "ready_for_pipeline_if_blockers_empty"
    )
    handoff = {
        "stage": stage,
        "review_required": review_required,
        "pipeline_command": command,
        "pipeline_root": str(PIPELINE_ROOT),
        "paths": paths,
        "agent_annotation": annotator_actions(paths, args),
        "next_gate": next_gate,
    }
    return {
        "pipeline_root": str(PIPELINE_ROOT),
        "launcher": str(PIPELINE_ROOT / "run_pipeline.sh"),
        "stage": stage,
        "blockers": blockers,
        "notes": notes,
        "command": command,
        "uses_project_env": True,
        "handoff": handoff,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-id", dest="project_id", required=True)
    parser.add_argument("--run-dir", dest="run_dir")
    parser.add_argument("--input-rds", dest="input_rds")
    parser.add_argument("--input-format", dest="input_format", default="auto", choices=["auto", "mainline_rds", "seurat_rds"])
    parser.add_argument("--mapping-file", dest="mapping_file")
    parser.add_argument("--annotation-marker-reference-file", dest="annotation_marker_reference_file")
    parser.add_argument("--annotation-min-support-markers", dest="annotation_min_support_markers", type=int, default=2)
    parser.add_argument("--annotation-allow-low-support", dest="annotation_allow_low_support", action="store_true")
    parser.add_argument("--target-gene-file", dest="target_gene_file", default="none")
    parser.add_argument("--target-gene-col", dest="target_gene_col", default="gene")
    parser.add_argument("--profile", default="local")
    parser.add_argument("--run-cellchat", dest="run_cellchat", action="store_true")
    parser.add_argument("--run-raw-import", dest="run_raw_import", action="store_true")
    parser.add_argument("--scrna-cohort-id", dest="scrna_cohort_id")
    parser.add_argument("--scrna-raw-umi-file", dest="scrna_raw_umi_file")
    parser.add_argument("--disease-name", dest="disease_name")
    parser.add_argument("--cancer-name", dest="cancer_name")
    parser.add_argument("--run-subset-prepare", dest="run_subset_prepare", action="store_true")
    parser.add_argument("--run-subset-apply", dest="run_subset_apply", action="store_true")
    parser.add_argument("--subset-mapping-file", dest="subset_mapping_file")
    parser.add_argument("--subset-mapping-manifest", dest="subset_mapping_manifest")
    parser.add_argument("--key-celltypes", dest="key_celltypes")
    parser.add_argument("--run-pseudotime", dest="run_pseudotime", action="store_true")
    parser.add_argument("--run-scmetabolism", dest="run_scmetabolism", action="store_true")
    parser.add_argument("--run-scmetabolism-specific", dest="run_scmetabolism_specific", action="store_true")
    parser.add_argument("--run-gsva", dest="run_gsva", action="store_true")
    parser.add_argument("--run-reactome-gsa", dest="run_reactome_gsa", action="store_true")
    parser.add_argument("--gsva-gmt-files", dest="gsva_gmt_files")
    parser.add_argument("--output", help="Optional JSON output path")
    parser.add_argument("--handoff-output", help="Optional handoff JSON output path")
    args = parser.parse_args()

    payload = build(args)
    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
    if args.handoff_output:
        out = Path(args.handoff_output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload["handoff"], indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()

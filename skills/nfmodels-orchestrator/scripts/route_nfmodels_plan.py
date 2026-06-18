#!/usr/bin/env python3
"""Create an NFmodels reviewable analysis plan and route manifest."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "lib"))
from _resolver import resolve_project_root  # noqa: E402

ROOT = resolve_project_root()

DOMAIN_RULES = {
    "transcriptome": [
        r"\bbulk\b", r"transcriptome", r"RNA[- ]?seq", r"DEG", r"differential expression",
        r"candidate gene", r"prognos", r"survival", r"cox", r"lasso", r"CIBERSORT",
        r"immune infiltration", r"TIDE", r"IPS", r"TMB", r"CNV", r"GSEA", r"GISTIC",
    ],
    "scrna": [
        r"scRNA", r"single[- ]cell", r"Seurat", r"cell type", r"annotation_prepare",
        r"annotation_apply", r"CellChat", r"cluster", r"UMAP", r"mapping_file",
    ],
}

ROUTE_META = {
    "transcriptome": {
        "adapter_skill": "transcriptome-analysis-orchestrator",
        "adapter_script": "skills/transcriptome-analysis-orchestrator/scripts/generate_transcriptome_command.py",
        "pipeline_root": "pipelines/transcriptome_prognosis_pipeline",
    },
    "scrna": {
        "adapter_skill": "scrna-analysis-orchestrator",
        "adapter_script": "skills/scrna-analysis-orchestrator/scripts/generate_scrna_command.py",
        "pipeline_root": "pipelines/scRNA_analysis_pipeline",
    },
}

TRANSCRIPTOME_CATALOG = {
    "rawdata_stage": {
        "name": "Cleaned bulk dataset staging", "phase": "phase0", "toggle": "startup_copy", "policy": "fully_executable",
        "inputs": "data_cleaned/Train_set/{train_id}; data_cleaned/Validation_set/{validation_ids}",
        "outputs": "results/00_rawdata/*.csv", "upstream": "-",
        "requires": ["train_id"],
        "rationale": "Stage cleaned expression, group, survival, and clinical files into the run directory.",
    },
    "deg": {
        "name": "Bulk DEG screening", "phase": "phase1", "toggle": "run_deg", "policy": "fully_executable",
        "inputs": "count matrix; group table", "outputs": "results/01_deg/02.deg_all.csv; 03.deg_sig.csv", "upstream": "rawdata_stage",
        "requires": ["train_id"], "rationale": "DESeq2 DEG and filtered DEG outputs.",
    },
    "candidate_genes": {
        "name": "Candidate gene generation", "phase": "phase1", "toggle": "run_intersect_candidate_genes", "policy": "fully_executable",
        "inputs": "DEG significant genes; gene_sets_sheet or candidate_genes", "outputs": "results/02_candidate_genes/01.candidate_genes.csv", "upstream": "deg",
        "requires_any": ["gene_sets_sheet", "candidate_genes"], "rationale": "Default route is intersect candidate genes unless precomputed candidate_genes is provided.",
    },
    "go_kegg_ppi": {
        "name": "GO/KEGG/PPI interpretation", "phase": "phase1", "toggle": "run_go_kegg_ppi", "policy": "fully_executable",
        "inputs": "candidate genes; STRING interactions", "outputs": "results/03_go_kegg_ppi/", "upstream": "candidate_genes",
        "requires_any": ["string_interactions_file", "candidate_genes"], "rationale": "Functional enrichment and PPI network outputs.",
    },
    "multi_model": {
        "name": "Multi-model prognosis modeling", "phase": "phase1", "toggle": "run_multi_model", "policy": "fully_executable",
        "inputs": "expression; survival; candidate genes; validation_ids or validation_sheet", "outputs": "results/04_multi_model/99_summary/*.csv", "upstream": "candidate_genes",
        "requires_any": ["validation_ids", "validation_sheet"], "rationale": "Train and validate Cox/LASSO/CoxBoost/RSF models.",
    },
    "stage": {
        "name": "Stage/prognostic and nomogram", "phase": "phase2", "toggle": "run_stage", "policy": "fully_executable",
        "inputs": "risk_file; survival; clinical table", "outputs": "results/xx_prognostic/", "upstream": "multi_model or supplied risk_file",
        "requires": ["risk_file"], "rationale": "Independent prognostic analysis, nomogram, calibration, and DCA.",
    },
    "gsea_gene": {
        "name": "Gene-based GSEA", "phase": "phase2", "toggle": "run_gsea_gene", "policy": "fully_executable",
        "inputs": "risk_file; expression; gene_file", "outputs": "results/xx_gsea_gene/", "upstream": "multi_model or supplied risk_file/gene_file",
        "requires": ["risk_file", "gene_file"], "rationale": "GSEA around selected model genes.",
    },
    "gsea_risk_nes": {
        "name": "Risk-based GSEA", "phase": "phase2", "toggle": "run_gsea_risk_nes", "policy": "fully_executable",
        "inputs": "risk_file; count matrix", "outputs": "results/xx_gsea_risk_nes/", "upstream": "multi_model or supplied risk_file",
        "requires": ["risk_file"], "rationale": "Risk-group ranked GSEA/NES outputs.",
    },
    "cibersort": {
        "name": "CIBERSORT immune infiltration", "phase": "phase2", "toggle": "run_cibersort", "policy": "fully_executable",
        "inputs": "risk_file; expression; gene_file", "outputs": "results/xx_cibersort/", "upstream": "multi_model or supplied risk_file/gene_file",
        "requires": ["risk_file", "gene_file"], "rationale": "Immune cell fractions and risk/gene associations.",
    },
    "tide": {
        "name": "TIDE immune response prediction", "phase": "phase2", "toggle": "run_tide", "policy": "fully_executable",
        "inputs": "risk_file; expression; NFMODELS_TIDEPY_BIN", "outputs": "results/xx_tide/", "upstream": "multi_model or supplied risk_file",
        "requires": ["risk_file"], "rationale": "TIDE scores and immune response prediction; tidepy checked by environment preflight.",
    },
    "ips": {
        "name": "IPS immune phenotype score", "phase": "phase2", "toggle": "run_ips", "policy": "fully_executable",
        "inputs": "risk_file; IPS file", "outputs": "results/xx_ips/", "upstream": "multi_model or supplied risk_file",
        "requires": ["risk_file"], "rationale": "IPS comparison by risk group.",
    },
    "cnv": {
        "name": "CNV analysis", "phase": "phase2", "toggle": "run_cnv", "policy": "executable_with_caveats",
        "inputs": "risk_file; expression; gene_file; optional gistic_matrix", "outputs": "results/xx_cnv/", "upstream": "cibersort; multi_model or supplied risk/gene files",
        "requires": ["risk_file", "gene_file"], "rationale": "May depend on TCGA-style CNV resources or supplied matrix.",
    },
    "gistic": {
        "name": "GISTIC input generation", "phase": "phase2", "toggle": "run_gistic", "policy": "manual_handoff_after_pipeline_input",
        "inputs": "CNV-like source data", "outputs": "results/xx_gistic/ input files", "upstream": "rawdata_stage",
        "rationale": "Pipeline prepares GISTIC input; running GISTIC2 remains outside the pipeline.",
    },
    "diff_expr": {
        "name": "Selected gene Tumor/Normal expression", "phase": "phase2", "toggle": "run_diff_expr", "policy": "fully_executable",
        "inputs": "gene_file; expression/group data", "outputs": "results/xx_diff_expr/", "upstream": "supplied gene_file or multi_model",
        "requires": ["gene_file"], "rationale": "Tumor/Normal differential plots for selected genes.",
    },
    "circos": {
        "name": "Circos plot", "phase": "phase2", "toggle": "run_circos", "policy": "executable_with_caveats",
        "inputs": "gene_file; GWAS/locus resources", "outputs": "results/xx_circos/", "upstream": "supplied gene_file or multi_model",
        "requires": ["gene_file"], "rationale": "Current defaults may require project-specific biology review.",
    },
    "tmb": {
        "name": "TMB and oncoplot", "phase": "phase2", "toggle": "run_tmb", "policy": "executable_with_caveats",
        "inputs": "risk_file; MAF or fallback mutation source", "outputs": "results/xx_tmb/", "upstream": "multi_model or supplied risk_file",
        "requires": ["risk_file"], "rationale": "Fallback mutation defaults may need project-specific review.",
    },
}

SCRNA_CATALOG = {
    "scrna_prepare_input": {
        "name": "Prepare scRNA input", "phase": "scrna_stage1", "toggle": "input_rds or run_raw_import", "policy": "fully_executable",
        "inputs": "input_rds or raw import settings", "outputs": "prepared Seurat RDS", "upstream": "-",
        "requires_any": ["input_rds", "run_raw_import"], "rationale": "Accept list(expr_mat, meta), Seurat RDS, or explicit raw import mode.",
    },
    "scrna_qc_integration": {
        "name": "QC, doublet removal, normalization, Harmony/UMAP", "phase": "scrna_stage1", "toggle": "always", "policy": "fully_executable",
        "inputs": "prepared Seurat RDS", "outputs": "precluster RDS; QC plots", "upstream": "scrna_prepare_input",
        "rationale": "Main QC and integration stage.",
    },
    "scrna_cluster_scan": {
        "name": "Resolution scan and final clustering", "phase": "scrna_stage1", "toggle": "always", "policy": "fully_executable",
        "inputs": "precluster RDS", "outputs": "clustered RDS; clustree/resolution outputs", "upstream": "scrna_qc_integration",
        "rationale": "Grid scan and final main clustering.",
    },
    "scrna_annotation_prepare": {
        "name": "Annotation prepare", "phase": "scrna_stage1", "toggle": "always", "policy": "fully_executable",
        "inputs": "clustered RDS", "outputs": "markers; SingleR labels; mapping template", "upstream": "scrna_cluster_scan",
        "rationale": "Generates the manual cell type mapping template for review.",
    },
    "scrna_annotation_apply": {
        "name": "Annotation apply", "phase": "scrna_stage2", "toggle": "mapping_file", "policy": "blocked_until_mapping_file",
        "inputs": "filled mapping_file", "outputs": "annotated RDS with celltype_manual", "upstream": "scrna_annotation_prepare + user-filled mapping",
        "requires": ["mapping_file"], "rationale": "Must wait until the user fills and reviews the mapping template.",
    },
    "scrna_cell_enrichment": {
        "name": "Cell proportion and enrichment", "phase": "scrna_stage2", "toggle": "run_cell_enrichment", "policy": "fully_executable_after_annotation",
        "inputs": "annotated RDS", "outputs": "cell proportion and enrichment tables/plots", "upstream": "scrna_annotation_apply",
        "requires": ["mapping_file"], "rationale": "Runs after manual cell type labels exist.",
    },
    "scrna_cellchat": {
        "name": "CellChat", "phase": "scrna_stage2", "toggle": "run_cellchat", "policy": "optional_executable_after_annotation",
        "inputs": "annotated RDS", "outputs": "CellChat communication tables", "upstream": "scrna_annotation_apply",
        "requires": ["mapping_file"], "rationale": "Optional communication analysis after annotation.",
    },
    "scrna_target_genes": {
        "name": "Target gene overlays/module scores", "phase": "scrna_stage2", "toggle": "target_gene_file", "policy": "optional_executable_when_supplied",
        "inputs": "manual target_gene_file", "outputs": "target gene plots/module score summaries", "upstream": "scrna_annotation_prepare or apply",
        "requires": ["target_gene_file"], "rationale": "Skipped when target_gene_file is absent; no hardcoded prognosis model path.",
    },
}

ANALYSIS_ALIASES = {
    "transcriptome": {
        "deg": [r"\bDEG\b", r"differential expression"],
        "candidate_genes": [r"candidate gene", r"candidate_genes"],
        "go_kegg_ppi": [r"GO/KEGG", r"PPI", r"go_kegg"],
        "multi_model": [r"multi[-_ ]?model", r"prognos", r"cox", r"lasso", r"RSF"],
        "stage": [r"stage", r"nomogram", r"列线图"],
        "gsea_gene": [r"gsea_gene", r"gene[- ]based GSEA"],
        "gsea_risk_nes": [r"gsea_risk", r"risk[- ]based GSEA", r"\bNES\b"],
        "cibersort": [r"CIBERSORT", r"immune infiltration", r"免疫浸润"],
        "tide": [r"\bTIDE\b", r"tidepy"],
        "ips": [r"\bIPS\b"],
        "cnv": [r"\bCNV\b"],
        "gistic": [r"GISTIC"],
        "diff_expr": [r"diff_expr", r"Tumor/Normal"],
        "circos": [r"circos"],
        "tmb": [r"\bTMB\b", r"oncoplot"],
    },
    "scrna": {
        "scrna_prepare_input": [r"input_rds", r"raw import", r"Seurat"],
        "scrna_qc_integration": [r"QC", r"integration", r"Harmony", r"UMAP"],
        "scrna_cluster_scan": [r"cluster", r"resolution", r"clustree"],
        "scrna_annotation_prepare": [r"annotation_prepare", r"annotation", r"cell type"],
        "scrna_annotation_apply": [r"annotation_apply", r"mapping_file"],
        "scrna_cell_enrichment": [r"cell enrichment", r"cell proportion", r"比例"],
        "scrna_cellchat": [r"CellChat"],
        "scrna_target_genes": [r"target_gene_file", r"target gene"],
    },
}

DEFAULT_DOMAIN_TASKS = {
    "transcriptome": ["rawdata_stage", "deg", "candidate_genes", "multi_model"],
    "scrna": ["scrna_prepare_input", "scrna_qc_integration", "scrna_cluster_scan", "scrna_annotation_prepare"],
}


def dump_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def normalize_file(path: Path) -> tuple[str, dict]:
    source = path.resolve()
    suffix = source.suffix.lower()
    if suffix == ".doc":
        raise SystemExit("Direct .doc ingestion is not supported. Convert it to .docx or .md first.")
    if suffix == ".md":
        text = source.read_text(encoding="utf-8")
        source_format = "md"
    elif suffix == ".docx":
        pandoc_bin = shutil.which("pandoc")
        if not pandoc_bin:
            raise SystemExit("pandoc is required to convert .docx reports to markdown.")
        result = subprocess.run([pandoc_bin, str(source), "--to", "gfm", "--wrap=none"], check=True, capture_output=True, text=True)
        text = result.stdout
        source_format = "docx"
    else:
        raise SystemExit(f"Unsupported plan/report format: {source.suffix}")
    return text.replace("\r\n", "\n").replace("\r", "\n"), {"type": "file", "path": str(source), "format": source_format}


def detect_domains(text: str) -> list[str]:
    return [domain for domain, patterns in DOMAIN_RULES.items() if any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns)]


def detect_paths(text: str) -> list[str]:
    candidates = re.findall(r"(?<![\w.-])/(?:[^\s`'\"<>])+", text)
    cleaned = []
    for item in candidates:
        value = re.split(r"[，。；、]", item, maxsplit=1)[0].rstrip(".,);]，。；、")
        if value and value not in cleaned:
            cleaned.append(value)
    return cleaned


def extract_field(text: str, key: str) -> str | None:
    pattern = rf"\b{re.escape(key)}\b\s*(?:=|:|：|为|是)?\s*([^\s，;；。]+)"
    match = re.search(pattern, text, flags=re.IGNORECASE)
    return match.group(1).strip() if match else None


def extract_inputs(text: str, project_id: str) -> dict:
    keys = [
        "train_id", "validation_ids", "validation_sheet", "expr_type", "gene_sets_sheet", "candidate_genes",
        "string_interactions_file", "risk_file", "gene_file", "input_rds", "mapping_file", "target_gene_file",
        "target_gene_col", "run_raw_import",
    ]
    values = {key: extract_field(text, key) for key in keys}
    values["project_id"] = project_id or extract_field(text, "project_id") or "UNSET_PROJECT_ID"
    values["expr_type"] = values.get("expr_type") or "fpkm"
    return {key: value for key, value in values.items() if value}


def detect_requested_tasks(text: str, domains: list[str]) -> dict[str, list[str]]:
    requested: dict[str, list[str]] = {}
    for domain in domains:
        hits = []
        for analysis_id, patterns in ANALYSIS_ALIASES[domain].items():
            if any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns):
                hits.append(analysis_id)
        if domain == "transcriptome" and hits:
            if any(item in hits for item in ["deg", "candidate_genes", "go_kegg_ppi", "multi_model"]):
                hits.insert(0, "rawdata_stage")
        if domain == "scrna" and hits:
            ordered = list(DEFAULT_DOMAIN_TASKS["scrna"])
            ordered.extend(item for item in hits if item not in ordered)
            hits = ordered
        requested[domain] = list(dict.fromkeys(hits or DEFAULT_DOMAIN_TASKS[domain]))
    return requested


def env_profile_for_domains(domains: list[str]) -> str:
    domain_set = set(domains)
    if domain_set == {"transcriptome"}:
        return "transcriptome"
    if domain_set == {"scrna"}:
        return "scrna"
    if domain_set:
        return "all"
    return "runtime"


def missing_requirements(task: dict, inputs: dict, requested: dict[str, list[str]] | None = None) -> list[str]:
    requested = requested or {}
    transcriptome_tasks = set(requested.get("transcriptome", []))
    phase1_handoff = "multi_model" in transcriptome_tasks
    missing = []
    for key in task.get("requires", []):
        if inputs.get(key):
            continue
        if phase1_handoff and key in {"risk_file", "gene_file"}:
            continue
        missing.append(key)
    any_group = task.get("requires_any", [])
    if any_group and not any(inputs.get(key) for key in any_group):
        if not (phase1_handoff and any(key in {"risk_file", "gene_file"} for key in any_group)):
            missing.append(" or ".join(any_group))
    return missing


def task_catalog_for_domain(domain: str) -> dict:
    return TRANSCRIPTOME_CATALOG if domain == "transcriptome" else SCRNA_CATALOG


def build_task_rows(requested: dict[str, list[str]], inputs: dict) -> list[dict]:
    rows = []
    task_num = 0
    for domain, task_ids in requested.items():
        catalog = task_catalog_for_domain(domain)
        for analysis_id in task_ids:
            task = catalog[analysis_id]
            missing = missing_requirements(task, inputs, requested)
            policy = task["policy"] if not missing else f"blocked_missing_input: {', '.join(missing)}"
            if not missing and domain == "transcriptome" and task.get("phase") == "phase2" and "multi_model" in requested.get("transcriptome", []):
                if not inputs.get("risk_file") and any(key in task.get("requires", []) for key in ["risk_file", "gene_file"]):
                    policy = "executable_after_phase1_model_handoff"
            rows.append({
                "task_id": f"{task_num:02d}", "domain": domain, "analysis_id": analysis_id, "task_name": task["name"],
                "phase": task["phase"], "inputs": task["inputs"], "outputs": task["outputs"], "upstream": task["upstream"],
                "pipeline_toggle": task["toggle"], "execution_policy": policy, "rationale": task["rationale"],
            })
            task_num += 1
    return rows


def md_table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        out.append("| " + " | ".join(str(cell).replace("\n", "<br>") for cell in row) + " |")
    return "\n".join(out)


def format_list(items: list[str]) -> str:
    return ", ".join(items) if items else "-"


def build_analysis_plan(text: str, manifest: dict, inputs: dict, requested: dict[str, list[str]], task_rows: list[dict]) -> str:
    unsupported = []
    lower = text.lower()
    for keyword in ["wgcna", "spatial", "空间转录组", "甲基化", "proteomics"]:
        if keyword in lower:
            unsupported.append(keyword)

    missing_input_rows = [row for row in task_rows if row["execution_policy"].startswith("blocked_missing_input")]
    caveat_rows = [row for row in task_rows if "caveat" in row["execution_policy"] or "manual_handoff" in row["execution_policy"]]

    lines = [
        "# NFmodels 分析计划",
        "",
        "## Plan Metadata",
        "",
        md_table(["field", "value"], [
            ["plan_id", f"{manifest['project_id']}_plan"],
            ["project_id", manifest["project_id"]],
            ["source_type", manifest["source"]["type"]],
            ["source_path", manifest["source"].get("path") or "direct_request"],
            ["review_status", manifest["review_status"]],
            ["environment_profile", manifest["environment_preflight"]["profile"]],
            ["require_tidepy", str(manifest["environment_preflight"]["require_tidepy"]).lower()],
        ]),
        "",
        "## Source Summary",
        "",
        "该计划由总控根据输入方案/请求生成，用于用户审查。审查通过前不进入环境检查后的 adapter 执行阶段，也不启动 Nextflow pipeline。",
        "",
        "## Detected Domains And Inputs",
        "",
        md_table(["item", "value"], [
            ["domains", format_list([route["domain"] for route in manifest["routes"]])],
            ["detected_paths", format_list(manifest["detected_paths"])],
            ["project_id", inputs.get("project_id", "UNSET_PROJECT_ID")],
            ["train_id", inputs.get("train_id", "")],
            ["validation_ids", inputs.get("validation_ids", "")],
            ["expr_type", inputs.get("expr_type", "fpkm")],
            ["input_rds", inputs.get("input_rds", "")],
            ["mapping_file", inputs.get("mapping_file", "")],
            ["target_gene_file", inputs.get("target_gene_file", "")],
        ]),
        "",
        "## Requested Analyses",
        "",
        md_table(["domain", "analyses"], [[domain, format_list(tasks)] for domain, tasks in requested.items()]),
        "",
        "## Preconditions",
        "",
        "- Bulk transcriptome pipeline uses cleaned training/validation datasets and project-level runtime variables from `NFmodels/.env`.",
        "- scRNA annotation is two-stage: `annotation_prepare` creates the mapping template; `annotation_apply` requires a filled `mapping_file`.",
        "- `target_gene_file` must be supplied manually for scRNA target gene outputs; no prognosis-model directory is inferred.",
        "- Environment preflight must pass before generating execution-ready adapter commands.",
        "",
        "## Task Inventory",
        "",
        md_table(
            ["task_id", "domain", "analysis_id", "task_name", "phase", "inputs", "outputs", "upstream", "pipeline_toggle", "execution_policy"],
            [[row[h] for h in ["task_id", "domain", "analysis_id", "task_name", "phase", "inputs", "outputs", "upstream", "pipeline_toggle", "execution_policy"]] for row in task_rows],
        ),
        "",
        "## Dependency Map",
        "",
        md_table(["task_id", "depends_on", "enables", "logic"], [[row["task_id"], row["upstream"], row["analysis_id"], row["rationale"]] for row in task_rows]),
        "",
        "## Scheduling",
        "",
        md_table(["stage", "tasks", "start_after", "logic"], schedule_rows(task_rows)),
        "",
        "## Execution Policy",
        "",
        md_table(["task_id", "analysis_id", "execution_policy", "pipeline_toggle", "rationale"], [[row["task_id"], row["analysis_id"], row["execution_policy"], row["pipeline_toggle"], row["rationale"]] for row in task_rows]),
        "",
        "## Required User Inputs Before Execution",
        "",
    ]
    if missing_input_rows:
        lines.append(md_table(["task_id", "analysis_id", "missing_or_needed_input"], [[row["task_id"], row["analysis_id"], row["execution_policy"].replace("blocked_missing_input: ", "")] for row in missing_input_rows]))
    else:
        lines.append("- No missing task-level inputs detected from the current request. Environment/package preflight may still block execution.")

    lines.extend([
        "",
        "## Steps Outside Pipeline",
        "",
        "### Manual-only steps",
        "",
        "- Review and approve this `02.analysis_plan.md` before execution.",
        "- Fill the scRNA annotation mapping template before `annotation_apply` if no `mapping_file` is supplied.",
        "- Interpret biological results and write final report manually.",
        "",
        "### Scriptable-but-not-yet-integrated steps",
        "",
    ])
    scriptable = []
    if any(row["analysis_id"] == "gistic" for row in task_rows):
        scriptable.append("Run GISTIC2 itself after pipeline-generated GISTIC input files; current pipeline only prepares input/handoff files.")
    if unsupported:
        scriptable.append("Unsupported requested analyses detected for future routing: " + ", ".join(unsupported))
    lines.append("\n".join(f"- {item}" for item in scriptable) if scriptable else "- None detected.")

    lines.extend([
        "",
        "## Pipeline-vs-Plan Gaps",
        "",
    ])
    gaps = []
    for row in caveat_rows:
        gaps.append(f"{row['analysis_id']}: {row['execution_policy']} - {row['rationale']}")
    if unsupported:
        gaps.append("Some requested analyses are not mapped to current NFmodels adapters: " + ", ".join(unsupported))
    lines.append("\n".join(f"- {item}" for item in gaps) if gaps else "- No pipeline-vs-plan gaps detected beyond normal manual review gates.")

    lines.extend([
        "",
        "## Environment Preflight",
        "",
        md_table(["field", "value"], [
            ["skill", manifest["environment_preflight"]["skill"]],
            ["profile", manifest["environment_preflight"]["profile"]],
            ["require_tidepy", str(manifest["environment_preflight"]["require_tidepy"]).lower()],
            ["runs_before", "adapter command generation/execution after plan approval"],
        ]),
        "",
        "## Review",
        "",
        md_table(["field", "value"], [
            ["review_status", manifest["review_status"]],
            ["approved_by", ""],
            ["approved_version", "1"],
            ["review_notes", ""],
        ]),
        "",
    ])
    return "\n".join(lines)


def schedule_rows(task_rows: list[dict]) -> list[list[str]]:
    by_phase: dict[str, list[str]] = {}
    for row in task_rows:
        by_phase.setdefault(row["phase"], []).append(f"{row['task_id']}:{row['analysis_id']}")
    starts = {
        "phase0": "plan approval",
        "phase1": "phase0 complete",
        "phase2": "phase1 model outputs or supplied risk/gene files",
        "scrna_stage1": "plan approval",
        "scrna_stage2": "annotation_prepare complete and mapping_file supplied when required",
    }
    return [[phase, ", ".join(tasks), starts.get(phase, "previous dependencies complete"), "Run serial dependencies first; independent tasks can be parallelized by adapter/pipeline."] for phase, tasks in by_phase.items()]


def build_manifest(text: str, source: dict, project_id: str, approved: bool, analysis_plan_path: Path | None = None) -> dict:
    domains = detect_domains(text)
    paths = detect_paths(text)
    blockers = []
    if not domains:
        blockers.append("No NFmodels analysis domain detected. Specify transcriptome, scRNA, or both.")
    if not approved:
        blockers.append("Analysis plan review_status is pending; do not enter environment/adapters execution stage.")

    env_profile = env_profile_for_domains(domains)
    require_tidepy = bool(re.search(r"\bTIDE\b|tidepy", text, flags=re.IGNORECASE))
    routes = [{"domain": domain, **ROUTE_META[domain]} for domain in domains]
    next_actions = [
        {
            "domain": "review",
            "skill": "nfmodels-orchestrator",
            "artifact": str(analysis_plan_path) if analysis_plan_path else None,
            "action": "User reviews and approves 02.analysis_plan.md before environment preflight or adapter execution.",
        },
        {
            "domain": "environment",
            "skill": "nfmodels-environment-check",
            "script": "skills/nfmodels-environment-check/scripts/check_nfmodels_env.py",
            "profile": env_profile,
            "require_tidepy": require_tidepy,
            "action": "Run only after analysis_plan review_status is approved.",
        },
    ]
    for route in routes:
        next_actions.append({
            "domain": route["domain"], "skill": route["adapter_skill"], "script": route["adapter_script"],
            "action": "Generate adapter command after plan approval and environment preflight.",
        })
    return {
        "project_id": project_id or "UNSET_PROJECT_ID",
        "source": source,
        "review_status": "approved" if approved else "pending",
        "analysis_plan_path": str(analysis_plan_path) if analysis_plan_path else None,
        "environment_preflight": {"skill": "nfmodels-environment-check", "profile": env_profile, "require_tidepy": require_tidepy},
        "routes": routes,
        "detected_paths": paths,
        "blockers": blockers,
        "next_actions": next_actions,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--plan-file", help="Path to a .md or .docx plan/report file")
    source.add_argument("--request", help="Natural-language NFmodels request")
    parser.add_argument("--project-id", default="UNSET_PROJECT_ID")
    parser.add_argument("--out-dir", required=True, help="Planning output directory")
    parser.add_argument("--approved", action="store_true", help="Mark generated plan/manifest as already reviewed")
    args = parser.parse_args()

    if args.plan_file:
        text, source_meta = normalize_file(Path(args.plan_file))
    else:
        text = args.request.strip()
        source_meta = {"type": "request", "path": None, "format": "text"}

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    normalized_path = out_dir / "01.normalized_request.md"
    analysis_plan_path = out_dir / "02.analysis_plan.md"
    manifest_path = out_dir / "03.route_manifest.json"
    compat_manifest_path = out_dir / "02.route_manifest.json"
    review_path = out_dir / "04.review_status.json"

    normalized_path.write_text(text + ("" if text.endswith("\n") else "\n"), encoding="utf-8")
    source_meta["normalized_path"] = str(normalized_path.resolve())

    manifest = build_manifest(text, source_meta, args.project_id, args.approved, analysis_plan_path.resolve())
    inputs = extract_inputs(text, manifest["project_id"])
    requested = detect_requested_tasks(text, [route["domain"] for route in manifest["routes"]])
    task_rows = build_task_rows(requested, inputs)
    manifest["detected_inputs"] = inputs
    manifest["requested_analyses"] = requested
    manifest["task_count"] = len(task_rows)

    analysis_plan = build_analysis_plan(text, manifest, inputs, requested, task_rows)
    analysis_plan_path.write_text(analysis_plan, encoding="utf-8")
    dump_json(manifest_path, manifest)
    dump_json(compat_manifest_path, manifest)
    dump_json(review_path, {"review_status": manifest["review_status"], "analysis_plan_path": str(analysis_plan_path.resolve())})
    print(str(analysis_plan_path.resolve()))
    print(str(manifest_path.resolve()))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build subset_mapping_manifest.csv from reviewed subset annotation outputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_context(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def find_rows(subset_agent_root: Path) -> list[dict]:
    rows = []
    for mapping in sorted(subset_agent_root.glob("*/08_subset_cluster_celltype_mapping_filled.csv")):
        out_dir = mapping.parent
        context_candidates = [
            out_dir / "10_subset_annotation_agent_context.json",
            out_dir / "agent_context.json",
        ]
        context = {}
        for candidate in context_candidates:
            if candidate.exists():
                context = read_context(candidate)
                break
        safe_label = str(context.get("safe_label") or out_dir.name)
        key_celltype = str(context.get("key_celltype") or safe_label)
        rows.append({
            "key_celltype": key_celltype,
            "safe_label": safe_label,
            "mapping_file": str(mapping.resolve()),
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subset-agent-root", required=True, help="Root containing one subdirectory per safe_label")
    parser.add_argument("--output", required=True, help="Output subset_mapping_manifest.csv")
    args = parser.parse_args()

    rows = find_rows(Path(args.subset_agent_root))
    if not rows:
        raise SystemExit("No reviewed subset mapping files found under subset agent root.")
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["key_celltype", "safe_label", "mapping_file"])
        writer.writeheader()
        writer.writerows(rows)
    print(str(out.resolve()))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

NULLS = {'', 'na', 'nan', 'none', 'null'}

def read_csv(path):
    with open(path, newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))

def split_markers(value):
    if value is None:
        return []
    return [x.strip() for x in value.replace(';', ',').split(',') if x.strip()]

def bad(value):
    return value is None or value.strip().lower() in NULLS

def main():
    ap = argparse.ArgumentParser(description='Validate scRNA annotation mapping support markers and references.')
    ap.add_argument('--mapping', required=True)
    ap.add_argument('--references', required=True)
    ap.add_argument('--template')
    ap.add_argument('--min-support', type=int, default=2)
    ap.add_argument('--allow-low-support', action='store_true')
    ap.add_argument('--out')
    args = ap.parse_args()

    mapping = read_csv(args.mapping)
    refs = read_csv(args.references)
    required_mapping = ['cluster', 'manual_celltype', 'manual_markers', 'annotation_flag', 'annotation_confidence', 'notes']
    required_refs = ['manual_markers', 'references']
    if not mapping:
        raise SystemExit('Mapping file is empty')
    missing = [c for c in required_mapping if c not in mapping[0]]
    if missing:
        raise SystemExit('Mapping missing columns: ' + ', '.join(missing))
    if not refs:
        raise SystemExit('Reference file is empty')
    missing = [c for c in required_refs if c not in refs[0]]
    if missing:
        raise SystemExit('Reference table missing columns: ' + ', '.join(missing))

    ref_by_marker = {}
    for row in refs:
        for marker in split_markers(row.get('manual_markers')):
            ref = (row.get('references') or '').strip()
            if not bad(ref):
                ref_by_marker.setdefault(marker.upper(), []).append(ref)

    expected_clusters = None
    if args.template:
        template = read_csv(args.template)
        if template and 'cluster' in template[0]:
            expected_clusters = {str(r.get('cluster', '')).strip() for r in template}

    rows = []
    errors = []
    seen = set()
    for row in mapping:
        cluster = str(row.get('cluster', '')).strip()
        seen.add(cluster)
        markers = split_markers(row.get('manual_markers'))
        supported = [m for m in markers if m.upper() in ref_by_marker]
        ok = len(supported) >= args.min_support
        flag = (row.get('annotation_flag') or '').strip()
        if bad(row.get('manual_celltype')):
            errors.append(f'cluster {cluster}: manual_celltype is empty')
        if not ok and not (args.allow_low_support or flag == 'low_support'):
            errors.append(f'cluster {cluster}: only {len(supported)} supported marker(s), need {args.min_support}')
        rows.append({
            'cluster': cluster,
            'manual_celltype': row.get('manual_celltype', ''),
            'manual_markers': row.get('manual_markers', ''),
            'supported_markers': ','.join(supported),
            'support_count': str(len(supported)),
            'passed': str(ok).lower(),
        })

    if expected_clusters is not None:
        missing_clusters = sorted(expected_clusters - seen)
        extra_clusters = sorted(seen - expected_clusters)
        if missing_clusters:
            errors.append('mapping missing clusters: ' + ','.join(missing_clusters))
        if extra_clusters:
            errors.append('mapping has extra clusters: ' + ','.join(extra_clusters))

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    if errors:
        for e in errors:
            print('ERROR:', e, file=sys.stderr)
        return 1
    print(f'OK: {len(mapping)} mapping rows passed with min_support={args.min_support}')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

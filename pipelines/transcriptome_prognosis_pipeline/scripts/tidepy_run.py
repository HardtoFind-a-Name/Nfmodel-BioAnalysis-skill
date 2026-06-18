#!/usr/bin/env python3
"""Merged tidepy pipeline: prepare input -> run tidepy -> extract scores."""
import argparse, os, shutil, subprocess, sys
import pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--expr', required=True, help='Expression matrix CSV (gene rows, sample cols)')
    p.add_argument('--outdir', required=True, help='Output directory for tidepy results')
    p.add_argument('--train-id', required=True)
    p.add_argument('--tidepy', default='tidepy', help='Path to tidepy binary')
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tid = args.train_id.lower().replace('-', '_')

    # Step 1: prepare input
    expr = pd.read_csv(args.expr, index_col=0)
    expr.index.name = 'gene'
    tide_input = os.path.join(args.outdir, f'01.{tid}_tide_input.tsv')
    expr.to_csv(tide_input, sep='\t')
    print(f'[1/3] TIDE input written: {tide_input}')

    # Step 2: run tidepy
    tide_full = os.path.join(args.outdir, f'02.{tid}_tide_full.tsv')
    cmd = [args.tidepy, tide_input, '-o', tide_full, '-c', 'NSCLC', '--ignore_norm']
    print(f'[2/3] Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)

    # Step 3: extract scores
    full = pd.read_csv(tide_full, sep='\t')
    cols = ['Responder', 'TIDE', 'Dysfunction', 'Exclusion', 'CAF', 'MDSC', 'MSI Expr Sig']
    available = [c for c in cols if c in full.columns]
    scores = full[available]
    # tidepy doesn't include sample IDs in output; use expression matrix column names
    sample_ids = list(expr.columns)
    if len(sample_ids) == len(scores):
        scores.insert(0, 'sample', sample_ids)
    else:
        print(f'[WARN] Sample count mismatch: expr={len(sample_ids)} vs scores={len(scores)}')
    tide_r = os.path.join(args.outdir, f'03.{tid}_tide_scores.tsv')
    scores.to_csv(tide_r, sep='\t', index=False)
    print(f'[3/3] R-ready scores: {tide_r} ({len(scores)} samples)')

if __name__ == '__main__':
    main()

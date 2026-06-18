# Transcriptome Prognosis Pipeline

Multi-model transcriptome prognosis pipeline with Nextflow and R. Covers DEG screening → multi-model risk modeling → post-modeling bioinformatics analysis.

## Workflow

```
DEG analysis (01) → Candidate genes (02) → GO/KEGG/PPI (03) → Multi-model risk (04)
                                                                    ↓
                                                         Post-modeling (xx_*)
```

## Directory model

- Pipeline root: `NFmodels/transcriptome_prognosis_pipeline`
- Run root: `${project_root}/${project_id}`
- Results: `${run_dir}/results/`
- R logs: `${run_dir}/logs/`
- NF logs: `${run_dir}/NFlogs/`
- Work cache: `${run_dir}/work/`
- Raw data: `${results_dir}/00_rawdata/`

---

## Multi-model risk pipeline

10 model paths in parallel (A1-A7 coefficient models, B1-B3 RSF models):

| Model | Path | Output |
|-------|------|--------|
| A1 | unicox → LASSO coef | riskScore = Σ gene × coef |
| A2 | unicox → multiCox coef | |
| A3 | unicox → multiCox → stepCox | |
| A4 | unicox → LASSO → multiCox | |
| A5 | unicox → LASSO → multiCox → stepCox | |
| A6 | unicox → LASSO → CoxBoost | grid search on training |
| A7 | unicox → LASSO → CoxBoost | validation-guided grid search |
| B1 | unicox → RSF | |
| B2 | unicox → LASSO → RSF | |
| B3 | unicox → LASSO → RSF | validation-guided grid search |

---

## Post-modeling analysis

All gated behind `--risk_file`. Gene list provided via `--gene_file` (any output with a `gene` column: unicox, lasso, stepCox, multiCox, coxboost).

| Module | Description | Key inputs |
|--------|-------------|------------|
| PROGNOSTIC | Independent prognostic factors + nomogram + DCA | risk_file, survival, clinical |
| GSEA_GENE | Per-gene Spearman → GSEA (multi-GMT) | expr, gene_file, GMT |
| GSEA_RISK_NES | Risk-based DESeq2 → GSEA | count_matrix, risk_file, GMT |
| CIBERSORT_IMMUNE | ESTIMATE + CIBERSORT deconvolution + immune-risk | expr, risk_file, gene_file |
| DIFF_EXPR | Prognostic gene Tumor/Normal Wilcoxon | expr, gene_file, group |
| CNV_ANALYSIS | GISTIC thresholded CNV heatmap + expr/immune by CNV | expr, gene_file, GISTIC matrix |
| GISTIC_INPUT | Generate .seg + marker files for GISTIC2 | risk_file |
| TIDE_ANALYSIS | tidepy TIDE scores + immune response prediction | expr, risk_file |
| IPS_ANALYSIS | IPS immunophenoscore from TCIA | risk_file, ips_file |
| CIRCOS_PLOT | Chromosome circos with GWAS SNPs + prognostic genes | gwas_file, risk_file, gene_file |
| TMB_ONCOPLOT | Tumor mutation burden + oncoplot (TCGAmutations fallback) | risk_file, maf_file |

---

## Database layout

```
{user_dir}/database/
├── MSigdb/
│   ├── c2.cp.kegg.v7.4.symbols.gmt
│   ├── c5.go.bp.v7.4.symbols.gmt
│   └── h.all.v7.5.1.symbols.gmt
├── CNV_downloads/
│   ├── TCGA.LUAD.sampleMap/
│   │   └── Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
│   └── snp6.na35.remap.hg38.subset.txt.gz
├── TCIA/
│   └── {train_id}_ClinicalData.tsv
├── TMB/
│   └── LUAD.maf (optional, TCGAmutations fallback available)
├── GWAS/
│   └── gwas-association.tsv
├── metascape/
│   └── 00.FINAL_GO.csv
└── string/
    └── 00.string_interactions.tsv
```

---

## Example: post-modeling only

```bash
./run_pipeline.sh \
  --project_id nf_r_test_3 \
  --risk_file results/04_multi_model/B3_rsf_vg/GSE31210_357y/02.train_risk_score.csv \
  --gene_file results/04_multi_model/01_lasso_select/02.lasso_final_genes_coef.csv \
  --run_deg false --run_intersect_candidate_genes false --run_go_kegg_ppi false --run_multi_model false \
  --run_stage true --run_gsea_gene true --run_gsea_risk_nes true --run_cibersort true \
  --run_cnv true --run_gistic true --run_tide true --run_ips true \
  --run_tmb true --run_circos true --run_diff_expr true \
  --gmt_gsea_gene h.all.v7.5.1.symbols.gmt --gmt_gsea_risk h.all.v7.5.1.symbols.gmt \
  --stage_variables "riskScore,age,gender,race,stage" --stage_time_set 357 \
  -profile local -log NFlogs/nextflow.log
```

---

## Notes

- Post-modeling modules depend on CIBERSORT → CNV (sequential via channel).
- GSEA sort order configurable via `--gsea_sort_by` (pvalue / p.adjust / NES).
- Stage analysis auto-detects continuous vs categorical variables; drops categorical levels with < 5 samples.
- CNV analysis downloads GISTIC2 matrix from UCSC Xena if local file not found.
- TMB uses TCGAmutations R package as automatic MAF fallback.
- Circos queries biomaRt (Ensembl `www` mirror) for gene chromosome positions.

## Runtime Environment

Runtime commands are configured in `../.env` and loaded automatically by `run_pipeline.sh`.

- `NFMODELS_NEXTFLOW_BIN`: command used to launch Nextflow
- `NFMODELS_RSCRIPT_CMD`: R command used inside Nextflow processes
- `NFMODELS_PYTHON_CMD`: Python command used by helper modules
- `NFMODELS_TIDEPY_BIN`: optional tidepy command for TIDE analysis

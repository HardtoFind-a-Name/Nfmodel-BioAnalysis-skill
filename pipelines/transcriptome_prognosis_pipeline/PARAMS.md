# 转录组预后建模流水线 — 参数文档

## 运行方式

```bash
./run_pipeline.sh -profile local \
  --project_id luad_v2 \
  --train_id TCGA-LUAD \
  --validation_ids GSE11969,GSE13213,GSE31210,GSE72094 \
  --expr_type fpkm \
  --deg_SN 01
```

---

## 1. 用户与项目路径

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--user_dir` | string | `/data/nas1/$USER` | 用户根目录，自动从 `$USER` 环境变量获取 |
| `--project_id` | string | `nf_r_test_2` | 项目 ID |
| `--project_root` | string | `{user_dir}/project` | 项目根目录 |
| `--run_dir` | string | `{project_root}/{project_id}` | 运行目录 |
| `--results_dir` | string | `{run_dir}/results` | 结果输出 |
| `--logdir` | string | `{run_dir}/logs` | R 脚本日志目录 |
| `--nf_log_dir` | string | `{run_dir}/NFlogs` | Nextflow 框架日志目录 |
| `--work_dir` | string | `{run_dir}/work` | Nextflow 工作目录 |
| `--rawdata_dir` | string | `{results_dir}/00_rawdata` | 原始数据（从 clean_root 自动复制） |

---

## 2. 数据源

### 2.1 清洗数据

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--clean_root` | string | `{user_dir}/data_cleaned` | 清洗后数据根目录 |
| `--train_id` | string | `TCGA-LUAD` | 训练集 ID，对应 `{clean_root}/Train_set/{train_id}` |
| `--validation_ids` | string | — | 验证集 ID 列表，逗号分隔，如 `GSE11969,GSE13213` |

流水线启动时自动将训练集和验证集从 `clean_root` 复制到 `rawdata_dir`。

### 2.2 表达矩阵

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--expr_type` | string | `fpkm` | 表达矩阵类型：`tpm` 或 `fpkm`。DEG 固定用 count 矩阵 |

### 2.3 原始数据文件（自动生成，无需手动指定）

| 文件 | 路径格式 |
|------|---------|
| count 矩阵 | `00_rawdata/01.{train_id_lower}_count.csv` |
| 分组表 | `00_rawdata/02.{train_id_lower}_group.csv` |
| TPM log2 | `00_rawdata/03.{train_id_lower}_tpmlog2.csv` |
| FPKM log2 | `00_rawdata/04.{train_id_lower}_fpkmlog2.csv` |
| 生存表 | `00_rawdata/05.{train_id_lower}_survival.csv` |
| 临床表 | `00_rawdata/06.{train_id_lower}_clinical.csv` |

---

## 3. 分析序号 (SN)

### 主流程

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--deg_SN` | `01` | DEG 分析输出 → `01_deg/` |
| `--candidate_SN` | `02` | 候选基因 → `02_candidate_genes/` |
| `--go_kegg_SN` | `03` | GO/KEGG/PPI → `03_go_kegg_ppi/` |
| `--multi_model_SN` | `04` | 多模型 → `04_multi_model/` |

### 后建模分析

| 参数 | 默认值 | 输出目录 |
|------|--------|---------|
| `--stage_SN` | `xx` | `{SN}_prognostic/` |
| `--gsea_gene_SN` | `xx` | `{SN}_gsea_gene/` |
| `--gsea_risk_nes_SN` | `xx` | `{SN}_gsea_risk_nes/` |
| `--cibersort_SN` | `xx` | `{SN}_cibersort/` |
| `--cnv_SN` | `xx` | `{SN}_cnv/` |
| `--gistic_SN` | `xx` | `{SN}_gistic/` |
| `--tide_SN` | `xx` | `{SN}_tide/` |
| `--ips_SN` | `xx` | `{SN}_ips/` |
| `--diff_expr_SN` | `xx` | `{SN}_diff_expr/` |
| `--circos_SN` | `xx` | `{SN}_circos/` |
| `--tmb_SN` | `xx` | `{SN}_tmb/` |

---

## 4. 验证集自动筛选 (validation_sheet.csv)

`validation_sheet.csv` 由 Python 脚本自动生成（`scripts/generate_validation_sheet.py`），无需手动创建。输出路径：`{run_dir}/validation_sheet.csv`

---

## 5. 步骤开关

### 主流程

| 参数 | 默认 | 说明 |
|------|------|------|
| `--run_deg` | `true` | DESeq2 差异表达 |
| `--run_intersect_candidate_genes` | `true` | 交集法候选基因 |
| `--run_go_kegg_ppi` | `true` | GO/KEGG/PPI |
| `--run_multi_model` | `true` | 多模型风险建模 |
| `--run_lasso_select` | `true` | 运行 `LASSO_SELECT` |
| `--run_model_A1` ~ `A7` | `true` | Cox 系数模型 A1-A7 |
| `--run_model_B1` ~ `B3` | `true` | 非 Cox 模型 B1-B3 |
| `--run_time_set_alt` | `true` | 训练集 357 年替代时间点 |

### 后建模分析（全部默认关闭，需 `--risk_file`）

| 参数 | 默认 | 说明 |
|------|------|------|
| `--run_stage` | `false` | 独立预后因素分析 + 列线图 |
| `--run_gsea_gene` | `false` | 靶基因 GSEA（Spearman + KEGG/GO/Hallmark） |
| `--run_gsea_risk_nes` | `false` | 风险评分 GSEA（DESeq2 + log2FC） |
| `--run_cibersort` | `false` | CIBERSORT 免疫浸润 |
| `--run_cnv` | `false` | CNV 拷贝数分析 |
| `--run_gistic` | `false` | GISTIC 输入生成 |
| `--run_tide` | `false` | TIDE 免疫治疗反应预测 |
| `--run_ips` | `false` | IPS 免疫表型评分 |
| `--run_diff_expr` | `false` | 预后基因 Tumor/Normal 差异表达 |
| `--run_circos` | `false` | 染色体定位圈图 |
| `--run_tmb` | `false` | TMB 突变负荷 + oncoplot |

---

## 6. 后建模分析参数

### 通用

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--risk_file` | string | `null` | 风险评分文件（必填，触发后建模） |
| `--gene_file` | string | `null` | 预后基因文件（含 `gene` 列）。兼容所有筛选结果：unicox / lasso / stepCox / multiCox / coxboost |

### Stage / Prognostic

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--stage_time_set` | string | `135` | 生存时间点：`135`（1/3/5 年）或 `357`（3/5/7 年） |
| `--stage_dca_year` | string | `3` | DCA 曲线时间点年数 |
| `--stage_min_level_n` | string | `5` | 分类变量最小水平样本数，低于此值排除 |
| `--stage_variables` | string | `null` | 手动指定变量，逗号分隔。不指定则自动检测所有临床变量 |

### GSEA

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--gmt_dir` | string | `{database_root}/MSigdb` | GMT 文件目录 |
| `--gmt_gsea_gene` | string | `c2.cp.kegg...,c5.go.bp...` | GSEA gene-based 的 GMT 文件名，逗号分隔 |
| `--gmt_gsea_risk` | string | `c2.cp.kegg...` | GSEA risk-based 的 GMT 文件名 |
| `--gsea_sort_by` | string | `p.adjust` | 多通路图排序字段：`pvalue` / `p.adjust` / `NES` |

### CNV / GISTIC

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--cnv_download_dir` | string | `{database_root}/CNV_downloads` | CNV 下载缓存目录 |
| `--gistic_matrix` | string | `null` | GISTIC 阈值化矩阵路径。null 时自动从 Xena 下载 |

### TIDE

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--tidepy_bin` | string | `NFMODELS_TIDEPY_BIN` / `tidepy` | tidepy 可执行命令，由 `../.env` 配置 |
| `--tidepy_cache_dir` | string | `{user_dir}/app/tidepy/{train_id}` | TIDE 结果缓存 |

### IPS / TMB / Circos

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--ips_file` | string | `{database_root}/TCIA/{train_id}_ClinicalData.tsv` | TCIA 临床数据。由 `train_id` 自动推导，无需手动指定 |
| `--maf_file` | string | `{database_root}/TMB/LUAD.maf` | MAF 突变文件。不存在时自动从 `TCGAmutations` 包下载 |
| `--gwas_file` | string | `{database_root}/GWAS/gwas-association.tsv` | GWAS 关联文件 |

---

## 7. 模型参数

### DEG / Cox / timeROC

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--deg_p_cutoff` | `0.05` | DEG adjusted p 阈值 |
| `--deg_logfc_cutoff` | `2` | DEG log2FC 阈值 |
| `--deg_use_pvalue` | `false` | DEG 回退用 raw p-value |
| `--unicox_p_cutoff` | `0.05` | 单因素 Cox p 阈值 |
| `--ph_p_cutoff` | `0.05` | PH 检验 p 阈值 |
| `--time_roc_points_days` | `365,1095,1825` | 主时间点 (135) |
| `--time_roc_points_days_alt` | `1095,1825,2555` | 替代时间点 (357) |
| `--model_auc_threshold` | `0.6` | pass/fail AUC 阈值 |

### RSF

| 参数 | 默认值 |
|------|--------|
| `--rsf_ntree_range` | `100,200,300,500,1000` |
| `--rsf_mtry_range` | `1,2,3` |
| `--rsf_nodesize_range` | `5,10,15,20,25,30` |
| `--rsf_seed` | `10` |

### CoxBoost

| 参数 | 默认值 |
|------|--------|
| `--coxboost_maxstepno_range` | `50,100,200,500` |
| `--coxboost_penalty_range` | `10,50,100,200,500` |

---

## 8. 常用运行示例

### 仅后建模分析

```bash
./run_pipeline.sh \
  --project_id nf_r_test_3 \
  --risk_file project/nf_r_test_3/results/04_multi_model/B3_rsf_vg/GSE31210_357y/02.train_risk_score.csv \
  --gene_file project/nf_r_test_3/results/04_multi_model/01_lasso_select/02.lasso_final_genes_coef.csv \
  --run_deg false --run_intersect_candidate_genes false --run_go_kegg_ppi false --run_multi_model false \
  --run_stage true --run_gsea_gene true --run_gsea_risk_nes true --run_cibersort true \
  --run_cnv true --run_gistic true --run_tide true --run_ips true \
  --run_tmb true --run_circos true --run_diff_expr true \
  --gmt_gsea_gene h.all.v7.5.1.symbols.gmt --gmt_gsea_risk h.all.v7.5.1.symbols.gmt \
  --stage_variables "riskScore,age,gender,race,stage" --stage_time_set 357 \
  -profile local -log NFlogs/nextflow.log
```

### SLURM 集群

```bash
./run_pipeline.sh -profile slurm \
  --project_id cluster_run \
  --validation_ids GSE11969,GSE13213,GSE31210,GSE72094
```

---

## 9. 项目级运行环境

运行时命令不在 pipeline 中硬编码。直接使用启动脚本，脚本会自动读取 `../.env`：

```bash
./run_pipeline.sh --project_id luad_v2 -profile local
```

| 环境变量 | 用途 | 默认示例 |
|------|------|------|
| `NFMODELS_NEXTFLOW_BIN` | 启动 Nextflow 的命令 | `nextflow` |
| `NFMODELS_RSCRIPT_CMD` | Nextflow process 内运行 R 的命令 | `conda run -n R4.3.3 Rscript` |
| `NFMODELS_PYTHON_CMD` | helper Python 脚本命令 | `python3` |
| `NFMODELS_TIDEPY_BIN` | TIDE 分析使用的 tidepy 命令 | `tidepy` |

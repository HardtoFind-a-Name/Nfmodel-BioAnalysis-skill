# NFmodels Analysis Project

版本：v0.2.0

## 介绍

`nfmodels-analysis-project` 是一个独立、可移植的生物信息学分析 skill 项目，
将 NFmodels 的转录组和单细胞分析流程封装为 5 个 Claude Code 兼容的 skill，
并捆绑 2 个 Nextflow pipeline，放入统一目录。

与旧版 `nfmodels-analysis-suite`（Codex 单体 skill）不同，本项目：
- 每个 skill 独立运作，可被 Claude Code 的 Skill 工具直接调用
- 通过 `lib/_resolver.py` 实现路径自动发现，无需硬编码 `parents[3]`
- 去除了 `agents/openai.yaml` 等 Codex 专有文件
- Pipeline 资源内置于 `pipelines/` 目录，无需外部 `NFmodels/` 仓库

该 skill 集成了由 R 脚本实现核心分析逻辑的 Nextflow pipeline：
- `pipelines/transcriptome_prognosis_pipeline`：bulk transcriptome 分析
- `pipelines/scRNA_analysis_pipeline`：single-cell RNA 分析

## 目录结构

```text
nfmodels-analysis-project/
├── .env                         ← 共享运行时配置
├── SKILL.md                     ← 项目级 meta-skill
├── README.md                    ← 本文档
│
├── lib/
│   └── _resolver.py             ← 统一路径发现模块
│
├── skills/                      ← 5 个 Claude Code 兼容 skill
│   ├── nfmodels-orchestrator/           ← 总控：分析方案编排 + 路由
│   │   ├── SKILL.md
│   │   ├── references/routing_contract.md
│   │   └── scripts/route_nfmodels_plan.py
│   ├── nfmodels-environment-check/      ← 环境检查
│   │   ├── SKILL.md
│   │   ├── references/environment_contract.md
│   │   ├── references/r_package_manifest.json
│   │   └── scripts/check_nfmodels_env.py
│   ├── transcriptome-analysis-orchestrator/  ← bulk 转录组适配器
│   │   ├── SKILL.md
│   │   ├── references/pipeline_capabilities.md
│   │   └── scripts/generate_transcriptome_command.py
│   ├── scrna-analysis-orchestrator/     ← 单细胞适配器
│   │   ├── SKILL.md
│   │   ├── references/pipeline_contract.md
│   │   ├── scripts/generate_scrna_command.py
│   │   └── scripts/build_subset_mapping_manifest.py
│   └── scrna-celltype-annotator/        ← 细胞类型注释
│       ├── SKILL.md
│       ├── references/annotation-workflow.md
│       └── scripts/validate_mapping.py
│
└── pipelines/                   ← 2 个 Nextflow pipeline
    ├── transcriptome_prognosis_pipeline/
    │   ├── main.nf, nextflow.config, run_pipeline.sh
    │   ├── conf/{base,local,slurm}.config
    │   ├── modules/local/*.nf
    │   └── scripts/*.R, *.py
    └── scRNA_analysis_pipeline/
        ├── main.nf, nextflow.config, run_pipeline.sh
        ├── conf/{base,local,slurm}.config
        ├── modules/local/*.nf
        └── scripts/*.R
```

## 路径发现机制

所有 skill 脚本通过 `lib/_resolver.py` 定位项目根目录，发现优先级：

1. `--nfmodels-dir` CLI 参数（显式指定）
2. `NFMODELS_ROOT` 环境变量
3. 自动探测（从 `lib/_resolver.py` 位置上溯）
4. 回退：当前工作目录

项目可部署到任意位置，无需依赖外部 `NFmodels/` 仓库。

## 环境配置

核心环境文件：项目根目录下的 `.env`

必须配置：
- `NFMODELS_NEXTFLOW_BIN`：Nextflow 启动命令或绝对路径
- `NFMODELS_RSCRIPT_CMD`：Rscript 命令，推荐使用 Conda 环境形式
- `NFMODELS_PYTHON_CMD`：Python 命令

可选配置：
- `NFMODELS_TIDEPY_BIN`：TIDE 分析需要
- `SCRNA_NEXTFLOW_BIN`：scRNA pipeline 专用 Nextflow 覆盖项
- `SCRNA_RSCRIPT_CMD`：scRNA pipeline 专用 Rscript 覆盖项

Pipeline 的 `run_pipeline.sh` 会自动从项目根加载 `.env`。

## Skill 说明

### 1. nfmodels-orchestrator（总控）

入口 skill。从自然语言需求或 `.md`/`.docx` 方案文件生成可审查的分析计划。

主要输出：
- `01.normalized_request.md`
- `02.analysis_plan.md`
- `03.route_manifest.json`
- `04.review_status.json`

审查通过后，分派到 `nfmodels-environment-check` 和对应 domain adapter。

### 2. nfmodels-environment-check（环境检查）

飞行前检查 gate。验证 `.env` 配置、Nextflow/Python/Rscript 可用性、R 包依赖。

使用方式：
```bash
python3 skills/nfmodels-environment-check/scripts/check_nfmodels_env.py --profile all
```

### 3. transcriptome-analysis-orchestrator（转录组适配器）

为 bulk transcriptome pipeline 生成已检查的 `run_pipeline.sh` 命令。支持的模块包括：
DEG、候选基因、GO/KEGG/PPI、多模型预后、GSEA、CIBERSORT、TIDE、IPS、CNV、GISTIC 等。

### 4. scrna-analysis-orchestrator（单细胞适配器）

为 scRNA pipeline 生成分阶段命令，管理两个人工注释 gate：
1. 主细胞类型注释：`annotation_prepare → annotator → annotation_apply`
2. 关键细胞亚群注释：`subset prepare → annotator → subset_apply`

### 5. scrna-celltype-annotator（注释辅助）

基于文献的 scRNA 细胞类型注释。读取 pipeline 生成的 marker 表和模板，
结合 `literature-search-review` skill 进行 PubMed 文献检索，生成可审查的注释文件。

## 使用方式

### 1. 从方案文件或需求启动

```bash
# 生成分析计划
python3 skills/nfmodels-orchestrator/scripts/route_nfmodels_plan.py \
  --plan-file /path/to/analysis_plan.docx \
  --project-id demo_project \
  --out-dir /path/to/planning/
```

### 2. 转录组分析

```bash
# 生成转录组 pipeline 命令
python3 skills/transcriptome-analysis-orchestrator/scripts/generate_transcriptome_command.py \
  --project-id demo_project \
  --train-id TCGA-LUAD \
  --analyses main
```

### 3. 单细胞分析

```bash
# 第一阶段：annotation_prepare
python3 skills/scrna-analysis-orchestrator/scripts/generate_scrna_command.py \
  --project-id demo_project \
  --input-rds /path/to/seurat.rds

# annotation 完成后，第二阶段：annotation_apply
python3 skills/scrna-analysis-orchestrator/scripts/generate_scrna_command.py \
  --project-id demo_project \
  --input-rds /path/to/seurat.rds \
  --mapping-file /path/to/08_cluster_celltype_mapping_filled.csv
```

## 适配系统

推荐运行环境是 Linux（服务器或 HPC 环境），需要：
- Bash 可用
- Nextflow 可用
- Python 3.8+ 可用
- R/Rscript 可用（推荐 Conda 环境 R4.3.3）
- Pipeline 所需 R 包已安装

## 推荐数据结构

```text
/data/nas1/${USER}/
├── database/           ← 跨项目共享数据库
├── data_cleaned/       ← bulk 清洗数据
│   ├── Train_set/
│   └── Validation_set/
└── project/            ← 分析运行目录
    └── <project_id>/
        ├── planning/
        ├── scripts/
        ├── logs/
        ├── work/
        └── results/
```

这只是推荐结构，实际路径可通过 CLI 参数显式指定。

## 主要生成文件

### 方案编排阶段
- `01.normalized_request.md`
- `02.analysis_plan.md`
- `03.route_manifest.json`
- `04.review_status.json`

### 命令生成阶段
- `transcriptome_command.json`
- `scrna_phase*_command.json`
- `scrna_phase*_handoff.json`

### scRNA 注释阶段
- `08_cluster_celltype_mapping_filled.csv`
- `09_literature_marker_reference_filled.csv`
- `annotation_literature_evidence.md`
- `annotation_decision_log.csv`

## 与旧版差异

| 维度 | nfmodels-analysis-suite (旧) | nfmodels-analysis-project (新) |
|---|---|---|
| 架构 | Codex 单体 skill + 内部 Role Router | 5 个独立 Claude Code skill |
| 路径解析 | `parents[3]` 硬编码 | `_resolver` 三层发现链 |
| Agent 描述 | `agents/openai.yaml` | SKILL.md frontmatter |
| Pipeline 位置 | `assets/nfmodels/` | `pipelines/` |
| 可移植性 | 依赖 `NFmodels/` 仓库布局 | 自包含，可部署到任意位置 |

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

核心环境文件位于项目根目录下的 `.env`：

| 变量 | 必需 | 说明 |
|------|------|------|
| `NFMODELS_NEXTFLOW_BIN` | ✅ | Nextflow 启动命令或绝对路径 |
| `NFMODELS_RSCRIPT_CMD` | ✅ | Rscript 命令，推荐使用 Conda 环境形式 |
| `NFMODELS_PYTHON_CMD` | ✅ | Python 命令 |
| `NFMODELS_TIDEPY_BIN` | ❌ | TIDE 分析需要 |
| `SCRNA_NEXTFLOW_BIN` | ❌ | scRNA pipeline 专用 Nextflow 覆盖项 |
| `SCRNA_RSCRIPT_CMD` | ❌ | scRNA pipeline 专用 Rscript 覆盖项 |

Pipeline 的 `run_pipeline.sh` 会自动从项目根加载 `.env`。

用户通常不需要自己运行环境检查命令。可以直接让 agent 使用本项目的环境检查 skill，例如：

```text
请检查 NFmodels 分析环境是否可用。
```

或者：

```text
请检查这个环境是否能跑转录组和单细胞分析；如果缺少 R 包，请列出需要安装的包。
```

agent 会调用 `nfmodels-environment-check` skill 检查 `.env`、Nextflow、Python、Rscript、R 包和可选 tidepy，并报告阻塞项、警告和建议。

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

## 使用方法

推荐的使用方式是：用户用自然语言告诉 agent 要做什么分析、数据在哪里、哪些结果需要审核；agent 根据本项目的 skill 生成分析计划、环境检查结果、pipeline 命令和注释交接文件。用户不需要直接记住所有 CLI 参数。

### 1. 从方案文件或用户需求开始

如果已有方案文件，可以这样告诉 agent：

```text
请读取 /path/to/analysis_plan.docx，整理可以执行的分析计划。项目 ID 是 demo_project。先不要运行 pipeline，先给我审查 02.analysis_plan.md。
```

如果只有自然语言需求，也可以直接说明：

```text
我想对 TCGA-LUAD 做转录组 DEG、候选基因、预后模型和免疫浸润分析，同时对 GSEXXXX 单细胞数据做注释和 CellChat。请先生成分析计划，不要直接运行。
```

agent 会调用 `nfmodels-orchestrator` skill 生成可审查的分析计划，通常包括：

- `01.normalized_request.md`
- `02.analysis_plan.md`
- `03.route_manifest.json`
- `04.review_status.json`

用户审核计划后，再让 agent 继续执行环境检查和 adapter 命令生成。

### 2. 转录组分析需要提供什么

转录组 pipeline 需要清洗后的表达和表型数据。推荐准备：

- 训练集 count 矩阵：用于 DEG、risk-based GSEA 等需要 count 的步骤
- 训练集 FPKM 矩阵：用于建模、免疫浸润、TIDE、部分表达分析等步骤
- 分组信息：Tumor/Normal 或其他比较分组
- 生存信息：预后建模必须提供
- 临床信息：stage/nomogram 等临床模型需要
- 验证集表达和生存信息：多队列模型验证需要
- 候选基因表或 gene sets：候选基因筛选需要
- 外部参考文件：如 STRING、GMT、IPS、CNV、MAF 等，按模块需要提供

可以按默认结构把清洗数据放在 `/data/nas1/${USER}/data_cleaned/` 下，也可以直接告诉 agent 每个文件的具体路径，例如：

```text
项目 ID 是 demo_project，训练集是 TCGA-LUAD。
清洗后的 count 矩阵在 /path/to/train_count.csv。
清洗后的 FPKM 矩阵在 /path/to/train_fpkm.csv。
生存文件在 /path/to/survival.csv。
临床文件在 /path/to/clinical.csv。
验证集包括 GSE11969 和 GSE31210。
STRING 文件在 /data/nas1/public_JOB062/database/string/00.string_interactions.tsv。
请生成 DEG、候选基因、GO/KEGG/PPI 和多模型预后分析命令，先不要运行。
```

agent 会根据这些信息选择 `transcriptome-analysis-orchestrator` skill，并生成对应 pipeline 命令。

### 3. 单细胞分析怎么启动

如果已有 Seurat RDS 或标准 mainline RDS，可以告诉 agent：

```text
项目 ID 是 demo_project。
单细胞输入 RDS 在 /data/nas1/public_JOB062/project/demo_project/results/00_rawdata/demo_seurat.rds。
请运行到 annotation_prepare，生成细胞类型注释模板后暂停，等待人工注释审核。
```

如果是 GEO 或 10X/UMI 原始数据，可以告诉 agent cohort ID 或原始文件路径：

```text
单细胞数据是 GSEXXXX，请从默认 database/scRNA 目录查找原始文件。如果自动识别不出来，我会再提供 raw UMI 或 annotation 文件路径。
```

第一阶段不需要提供 `mapping_file`。pipeline 会生成 marker 表、mapping template、reference template 和 agent context，然后等待注释审核。

### 4. 主细胞类型注释怎么继续

当 `annotation_prepare` 完成后，可以让 agent 读取输出模板并进行注释辅助：

```text
请根据 results/03_scrna_annotation_prepare 里的 marker 表、mapping template 和 agent context，辅助填写主细胞类型注释，并生成 evidence、decision log 和 validation report。
```

审核通过后，继续告诉 agent：

```text
主细胞类型 mapping 已审核通过，文件是 /path/to/08_cluster_celltype_mapping_filled.csv。
文献 marker reference 文件是 /path/to/09_literature_marker_reference_filled.csv。
请继续 annotation_apply，并为 T cells 和 Myeloid cells 做 subset annotation prepare。
```

agent 会调用 `scrna-analysis-orchestrator` 和 `scrna-celltype-annotator` skill 生成后续命令。

### 5. Subset 注释和下游分析怎么继续

subset annotation prepare 完成后，告诉 agent subset mapping 文件在哪里：

```text
subset mapping manifest 在 /path/to/subset_mapping_manifest.csv。
请继续 subset_apply，并打开 pseudotime 和 scMetabolism 分析。
```

### 6. 路径如何提供

推荐的数据结构只是为了让默认配置更容易工作，不是强制结构。实际运行时，用户可以直接用自然语言告诉 agent 关键路径，例如：

```text
我的 gene file 在 /path/to/genes.csv。
我的 risk file 在 /path/to/risk_score.csv。
我的 target gene file 在 /path/to/target_genes.csv，基因列名是 gene。
我的主细胞注释 mapping 在 /path/to/08_cluster_celltype_mapping_filled.csv。
我的 subset mapping manifest 在 /path/to/subset_mapping_manifest.csv。
请根据这些路径生成命令。
```

agent 会把这些路径转换成底层 CLI 参数。只有在需要复现或调试时，用户才需要查看具体命令。

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

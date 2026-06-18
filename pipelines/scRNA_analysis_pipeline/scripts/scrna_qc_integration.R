rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

input_rds <- get_arg("--input-rds")
outdir <- get_arg("--outdir")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")
sample_col <- get_arg("--sample-col", required = FALSE, default = "Sample")
group_col <- get_arg("--group-col", required = FALSE, default = "group")
min_features_keep <- as.numeric(get_arg("--min-features-keep", required = FALSE, default = "200"))
max_features_keep <- as.numeric(get_arg("--max-features-keep", required = FALSE, default = "4000"))
min_counts_keep <- as.numeric(get_arg("--min-counts-keep", required = FALSE, default = "500"))
max_counts_keep <- as.numeric(get_arg("--max-counts-keep", required = FALSE, default = "20000"))
max_percent_mt <- as.numeric(get_arg("--max-percent-mt", required = FALSE, default = "15"))
run_doublet_removal <- as_bool(get_arg("--run-doublet-removal", required = FALSE, default = "true"))
nfeatures_hvg <- as.integer(get_arg("--nfeatures-hvg", required = FALSE, default = "3000"))
npcs <- as.integer(get_arg("--npcs", required = FALSE, default = "50"))
use_dims <- parse_int_vector(get_arg("--use-dims", required = FALSE, default = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20"))
harmony_batch_var <- get_arg("--harmony-batch-var", required = FALSE, default = "Sample")
seed <- as.integer(get_arg("--seed", required = FALSE, default = "1234"))

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_qc_integration")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_qc_integration")
log_message("Input RDS:", input_rds)
log_message("QC thresholds: features", min_features_keep, "-", max_features_keep,
            "counts", min_counts_keep, "-", max_counts_keep, "percent.mt <", max_percent_mt)
set.seed(seed)

obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
cells_before <- ncol(obj)
genes_before <- count_detected_genes(obj)

mt_pattern <- if (any(grepl("^MT-", rownames(obj)))) "^MT-" else if (any(grepl("^mt-", rownames(obj)))) "^mt-" else NA_character_
if (!is.na(mt_pattern)) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  log_message("Mitochondrial pattern:", mt_pattern)
} else {
  obj[["percent.mt"]] <- 0
  log_message("No mitochondrial gene prefix found; percent.mt set to 0")
}
obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
if (sample_col %in% colnames(obj@meta.data)) {
  p <- VlnPlot(obj, features = qc_features, group.by = sample_col, pt.size = 0, ncol = 1)
  save_plot_dual(p, outdir, "01_vlnplot_before_qc_by_sample", width = 10, height = 14)
}
if (group_col %in% colnames(obj@meta.data)) {
  p <- VlnPlot(obj, features = qc_features, group.by = group_col, pt.size = 0, ncol = 1)
  save_plot_dual(p, outdir, "02_vlnplot_before_qc_by_group", width = 8, height = 12)
}

log_message("Running basic QC filter")
obj <- subset(
  obj,
  subset = nFeature_RNA >= min_features_keep &
    nFeature_RNA <= max_features_keep &
    nCount_RNA >= min_counts_keep &
    nCount_RNA <= max_counts_keep &
    percent.mt < max_percent_mt
)
cells_after_basic_qc <- ncol(obj)
doublet_removed <- 0
doublet_method <- "not_run"

if (run_doublet_removal) {
  if (requireNamespace("scDblFinder", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    log_message("Running scDblFinder")
    ok <- tryCatch({
      sce <- Seurat::as.SingleCellExperiment(obj)
      if (sample_col %in% colnames(SummarizedExperiment::colData(sce))) {
        sce <- scDblFinder::scDblFinder(sce, samples = sce[[sample_col]])
      } else {
        sce <- scDblFinder::scDblFinder(sce)
      }
      obj$doublet_class <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
      obj$doublet_score <- SummarizedExperiment::colData(sce)$scDblFinder.score
      TRUE
    }, error = function(e) {
      log_message("[警告] scDblFinder failed:", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) {
      doublet_removed <- sum(obj$doublet_class == "doublet", na.rm = TRUE)
      obj <- subset(obj, subset = doublet_class != "doublet")
      doublet_method <- "scDblFinder"
    }
  } else {
    log_message("[警告] scDblFinder dependencies not available; doublet removal skipped")
    doublet_method <- "skipped_missing_scDblFinder"
  }
}

cells_after <- ncol(obj)
genes_after <- count_detected_genes(obj)
if (cells_after < 100) save_note_and_stop("QC retained fewer than 100 cells.")

qc_stats <- data.frame(
  metric = c("genes_before_qc", "genes_after_qc", "cells_before_qc", "cells_after_basic_qc",
             "doublets_removed", "cells_after_qc", "cell_retention_rate", "doublet_method"),
  value = c(genes_before, genes_after, cells_before, cells_after_basic_qc, doublet_removed,
            cells_after, round(cells_after / cells_before, 4), doublet_method)
)
write.csv(qc_stats, file.path(outdir, "03.qc_filter_stats.csv"), row.names = FALSE)

if (sample_col %in% colnames(obj@meta.data)) {
  sample_stats <- obj@meta.data %>% dplyr::count(.data[[sample_col]], name = "cells_after_qc")
  write.csv(sample_stats, file.path(outdir, "04.sample_cell_counts_after_qc.csv"), row.names = FALSE)
}

log_message("NormalizeData")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
log_message("FindVariableFeatures")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures_hvg, verbose = FALSE)
p_hvg <- LabelPoints(VariableFeaturePlot(obj), points = head(VariableFeatures(obj), 10), repel = TRUE)
save_plot_dual(p_hvg, outdir, "05_variable_features", width = 8, height = 6)

log_message("ScaleData and PCA")
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
save_plot_dual(ElbowPlot(obj, ndims = min(npcs, 50)), outdir, "06_elbowplot", width = 6, height = 5)

reduction_use <- "pca"
if (nzchar(harmony_batch_var) && harmony_batch_var %in% colnames(obj@meta.data) && requireNamespace("harmony", quietly = TRUE)) {
  log_message("RunHarmony using", harmony_batch_var)
  obj <- harmony::RunHarmony(obj, group.by.vars = harmony_batch_var, reduction = "pca", dims.use = use_dims, verbose = FALSE)
  reduction_use <- "harmony"
} else {
  log_message("[警告] Harmony skipped; using PCA for downstream cluster scan")
}

max_dim <- ncol(Embeddings(obj, reduction = reduction_use))
use_dims <- use_dims[use_dims <= max_dim]
if (length(use_dims) == 0) save_note_and_stop("No valid dimensions available for downstream cluster scan.")
obj <- RunUMAP(obj, reduction = reduction_use, dims = use_dims, verbose = FALSE)

if (sample_col %in% colnames(obj@meta.data)) {
  save_plot_dual(DimPlot(obj, reduction = "umap", group.by = sample_col) + theme_bw(),
                 outdir, "07_umap_by_sample", width = 9, height = 7)
}
if (group_col %in% colnames(obj@meta.data)) {
  save_plot_dual(DimPlot(obj, reduction = "umap", group.by = group_col) + theme_bw(),
                 outdir, "08_umap_by_group", width = 8, height = 6)
}

obj@misc$scrna_pipeline <- obj@misc$scrna_pipeline %||% list()
obj@misc$scrna_pipeline$reduction_used <- reduction_use
obj@misc$scrna_pipeline$dims_used <- use_dims
obj@misc$scrna_pipeline$cluster_source <- "main"

saveRDS(obj, file.path(outdir, "01.seurat_precluster.rds"))
log_message("Script completed: scrna_qc_integration")

write_summary(c(
  "scRNA QC and integration completed.",
  paste0("Final cells: ", cells_after),
  paste0("Final genes detected: ", genes_after),
  paste0("Cell retention rate: ", round(cells_after / cells_before, 4)),
  paste0("Doublet method: ", doublet_method, "; removed: ", doublet_removed),
  paste0("Reduction for downstream clustering: ", reduction_use),
  paste0("Dims for downstream clustering: ", paste(use_dims, collapse = ",")),
  "QC filtering, dimensional reduction, batch correction, and UMAP were completed. Final clustering is performed by cluster_scan.",
  paste0("Object file: 01.seurat_precluster.rds\n", describe_rds(file.path(outdir, "01.seurat_precluster.rds"))),
  "",
  "限制与说明",
  "QC thresholds and dimensional choices are carried forward to cluster_scan through object metadata."
))

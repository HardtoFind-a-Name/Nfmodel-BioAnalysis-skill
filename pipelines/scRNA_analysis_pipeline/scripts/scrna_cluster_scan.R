rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
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
scan_label <- get_arg("--scan-label", required = FALSE, default = "main")
preferred_reduction <- get_arg("--preferred-reduction", required = FALSE, default = "harmony")
fallback_reduction <- get_arg("--fallback-reduction", required = FALSE, default = "pca")
cluster_scan_use_dims_override <- get_arg("--cluster-scan-use-dims-override", required = FALSE, default = "")
resolution_grid <- parse_num_vector(get_arg("--resolution-grid", required = FALSE, default = "0.1,0.2,0.3,0.4,0.5,0.6"))
final_resolution <- as.numeric(get_arg("--final-resolution", required = FALSE, default = "0.2"))
cluster_algorithm <- as.integer(get_arg("--cluster-algorithm", required = FALSE, default = "1"))
umap_n_neighbors <- as.integer(get_arg("--umap-n-neighbors", required = FALSE, default = "30"))
umap_min_dist <- as.numeric(get_arg("--umap-min-dist", required = FALSE, default = "0.3"))
sample_col <- get_arg("--sample-col", required = FALSE, default = "Sample")
group_col <- get_arg("--group-col", required = FALSE, default = "group")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_cluster_scan", subdir = scan_label)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

fmt_res <- function(x) gsub("\\.", "p", sprintf("%.1f", x))

log_message("Script started: scrna_cluster_scan")
log_message("Scan label:", scan_label)
log_message("Input RDS:", input_rds)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))

pipeline_meta <- obj@misc$scrna_pipeline %||% list()
meta_reduction <- pipeline_meta$reduction_used %||% NA_character_
if (!is.na(meta_reduction) && meta_reduction %in% Reductions(obj)) {
  use_reduction <- meta_reduction
} else if (preferred_reduction %in% Reductions(obj)) {
  use_reduction <- preferred_reduction
} else if (fallback_reduction %in% Reductions(obj)) {
  use_reduction <- fallback_reduction
} else {
  save_note_and_stop(paste0("No usable reduction found. Metadata reduction=", meta_reduction,
                            "; preferred/fallback=", preferred_reduction, "/", fallback_reduction))
}

override_dims <- parse_int_vector(cluster_scan_use_dims_override, default = NULL)
if (!is.null(override_dims) && length(override_dims) > 0) {
  use_dims <- override_dims
  dim_source <- "cluster_scan_use_dims_override"
} else if (!is.null(pipeline_meta$dims_used) && length(pipeline_meta$dims_used) > 0) {
  use_dims <- as.integer(pipeline_meta$dims_used)
  dim_source <- "object_metadata"
} else {
  save_note_and_stop("Input object does not contain obj@misc$scrna_pipeline$dims_used. Provide cluster_scan_use_dims_override only if intentional.")
}

use_dims <- use_dims[use_dims <= ncol(Embeddings(obj, reduction = use_reduction))]
if (length(use_dims) == 0) save_note_and_stop("No valid dimensions for selected reduction.")
if (!final_resolution %in% resolution_grid) save_note_and_stop("final_resolution must be included in resolution_grid.")

log_message("Cells:", ncol(obj), "Genes:", nrow(obj))
log_message("Reduction:", use_reduction, "Dims:", paste(use_dims, collapse = ","), "Dim source:", dim_source)

obj <- FindNeighbors(obj, reduction = use_reduction, dims = use_dims, verbose = FALSE)
scan_stat_list <- list()
for (res in resolution_grid) {
  log_message("FindClusters resolution =", res)
  obj <- FindClusters(obj, resolution = res, algorithm = cluster_algorithm, verbose = FALSE)
  res_col <- paste0("RNA_snn_res.", res)
  stat_df <- obj@meta.data %>%
    dplyr::count(.data[[res_col]], name = "cell_count") %>%
    dplyr::rename(cluster = 1) %>%
    dplyr::mutate(resolution = res)
  scan_stat_list[[as.character(res)]] <- stat_df
}

scan_stats <- bind_rows(scan_stat_list)
write.csv(scan_stats, file.path(outdir, "01_resolution_scan_cluster_stats.csv"), row.names = FALSE)

if (requireNamespace("clustree", quietly = TRUE)) {
  tryCatch({
    p_tree <- clustree::clustree(obj@meta.data, prefix = "RNA_snn_res.")
    save_plot_dual(p_tree, outdir, "02_clustree_resolution_scan", width = 12, height = 8)
  }, error = function(e) {
    log_message("[警告] clustree plot failed and was skipped:", conditionMessage(e))
  })
} else {
  log_message("[警告] clustree not available; clustree plot skipped")
}

final_col <- paste0("RNA_snn_res.", final_resolution)
obj$seurat_clusters <- factor(obj[[final_col, drop = TRUE]], levels = sort_cluster_ids(obj[[final_col, drop = TRUE]]))
Idents(obj) <- "seurat_clusters"

if (!"umap" %in% Reductions(obj)) {
  obj <- RunUMAP(obj, reduction = use_reduction, dims = use_dims, n.neighbors = umap_n_neighbors, min.dist = umap_min_dist, verbose = FALSE)
  umap_source <- "generated"
} else {
  umap_source <- "reused"
}

for (res in resolution_grid) {
  res_col <- paste0("RNA_snn_res.", res)
  p <- DimPlot(obj, reduction = "umap", group.by = res_col, label = TRUE, repel = TRUE) +
    ggtitle(paste0("UMAP colored by ", res_col)) + theme_bw()
  save_plot_dual(p, outdir, paste0("03_umap_colored_by_res_", fmt_res(res)), width = 9, height = 7)
}

save_plot_dual(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + theme_bw(),
               outdir, "04_umap_by_final_clusters", width = 9, height = 7)
if (sample_col %in% colnames(obj@meta.data)) {
  save_plot_dual(DimPlot(obj, reduction = "umap", group.by = sample_col) + theme_bw(),
                 outdir, "05_umap_by_sample", width = 10, height = 8)
}
if (group_col %in% colnames(obj@meta.data)) {
  save_plot_dual(DimPlot(obj, reduction = "umap", group.by = group_col) + theme_bw(),
                 outdir, "06_umap_by_group", width = 9, height = 7)
}

cluster_counts <- obj@meta.data %>% dplyr::count(seurat_clusters, name = "cell_count")
write.csv(cluster_counts, file.path(outdir, "07_final_cluster_cell_counts.csv"), row.names = FALSE)

obj@misc$scrna_pipeline <- obj@misc$scrna_pipeline %||% list()
obj@misc$scrna_pipeline$cluster_scan_label <- scan_label
obj@misc$scrna_pipeline$cluster_scan_reduction <- use_reduction
obj@misc$scrna_pipeline$cluster_scan_dims <- use_dims
obj@misc$scrna_pipeline$cluster_scan_resolution_grid <- resolution_grid
obj@misc$scrna_pipeline$cluster_scan_final_resolution <- final_resolution

saveRDS(obj, file.path(outdir, "01.seurat_qc_reclustered_final.rds"))

log_message("Script completed: scrna_cluster_scan")
write_summary(c(
  "scRNA cluster resolution scan completed.",
  paste0("Scan label: ", scan_label),
  paste0("Final cells: ", ncol(obj)),
  paste0("Final genes: ", nrow(obj)),
  paste0("Reduction used: ", use_reduction),
  paste0("Dims used: ", paste(use_dims, collapse = ",")),
  paste0("Dim source: ", dim_source),
  paste0("UMAP source: ", umap_source),
  paste0("Resolution grid: ", paste(resolution_grid, collapse = ",")),
  paste0("Final resolution: ", final_resolution),
  paste0("Final cluster count: ", nrow(cluster_counts)),
  "The selected resolution was applied as seurat_clusters for downstream annotation.",
  paste0("Object file: 01.seurat_qc_reclustered_final.rds\n", describe_rds(file.path(outdir, "01.seurat_qc_reclustered_final.rds"))),
  "",
  "限制与说明",
  "Resolution selection is data dependent; dimensions are inherited from upstream object metadata unless explicitly overridden."
))

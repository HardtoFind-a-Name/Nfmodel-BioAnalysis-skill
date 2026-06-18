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
subset_celltypes <- parse_chr_vector(get_arg("--subset-celltypes", required = FALSE, default = "NK cells"))
subset_npcs <- as.integer(get_arg("--subset-npcs", required = FALSE, default = "30"))
subset_use_dims <- parse_int_vector(get_arg("--subset-use-dims", required = FALSE, default = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20"))
subset_nfeatures_hvg <- as.integer(get_arg("--subset-nfeatures-hvg", required = FALSE, default = "2000"))
subset_run_jackstraw <- as_bool(get_arg("--subset-run-jackstraw", required = FALSE, default = "false"))
seed <- as.integer(get_arg("--seed", required = FALSE, default = "1234"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

subset_annotation_col <- "celltype_manual"

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_subset_prepare")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop
set.seed(seed)

log_message("Script started: scrna_subset_prepare")
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!subset_annotation_col %in% colnames(obj@meta.data)) {
  save_note_and_stop("Required upstream annotation column `celltype_manual` is missing. Run SCRNA_ANNOTATION_APPLY first.")
}
if (length(subset_celltypes) == 0) save_note_and_stop("subset_celltypes is empty.")
subset_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[subset_annotation_col]]) %in% subset_celltypes]
if (length(subset_cells) == 0) {
  save_note_and_stop(paste0("No cells found for subset_celltypes: ", paste(subset_celltypes, collapse = ", ")))
}
subset_obj <- subset(obj, cells = subset_cells)
DefaultAssay(subset_obj) <- "RNA"

composition <- subset_obj@meta.data %>%
  dplyr::count(.data[[subset_annotation_col]], name = "cell_count") %>%
  dplyr::arrange(dplyr::desc(cell_count))
write.csv(composition, file.path(outdir, "01_subset_celltype_composition.csv"), row.names = FALSE)

subset_obj <- NormalizeData(subset_obj, verbose = FALSE)
subset_obj <- FindVariableFeatures(subset_obj, nfeatures = subset_nfeatures_hvg, verbose = FALSE)
subset_obj <- ScaleData(subset_obj, verbose = FALSE)
subset_obj <- RunPCA(subset_obj, npcs = subset_npcs, verbose = FALSE)
valid_dims <- subset_use_dims[subset_use_dims <= ncol(Embeddings(subset_obj, "pca"))]
if (length(valid_dims) == 0) save_note_and_stop("No valid subset_use_dims for PCA reduction.")

p_elbow <- ElbowPlot(subset_obj, ndims = min(subset_npcs, 50)) + theme_bw()
save_plot_dual(p_elbow, outdir, "02_subset_pca_elbowplot", width = 8, height = 6)

if (subset_run_jackstraw) {
  subset_obj <- JackStraw(subset_obj, num.replicate = 50, dims = min(subset_npcs, 30), verbose = FALSE)
  subset_obj <- ScoreJackStraw(subset_obj, dims = seq_len(min(subset_npcs, 30)))
  p_jack <- JackStrawPlot(subset_obj, dims = seq_len(min(subset_npcs, 30))) + theme_bw()
  save_plot_dual(p_jack, outdir, "03_subset_jackstrawplot", width = 12, height = 8)
}

umap_neighbors <- max(2, min(30, ncol(subset_obj) - 1))
subset_obj <- RunUMAP(subset_obj, reduction = "pca", dims = valid_dims, n.neighbors = umap_neighbors, verbose = FALSE)
subset_obj@misc$scrna_pipeline <- subset_obj@misc$scrna_pipeline %||% list()
subset_obj@misc$scrna_pipeline$reduction_used <- "pca"
subset_obj@misc$scrna_pipeline$dims_used <- valid_dims
subset_obj@misc$scrna_pipeline$cluster_source <- "subset"
subset_obj@misc$scrna_pipeline$subset_celltypes <- subset_celltypes

save_plot_dual(DimPlot(subset_obj, reduction = "umap", group.by = subset_annotation_col) + theme_bw(),
               outdir, "04_subset_umap_by_parent_celltype", width = 9, height = 7)
saveRDS(subset_obj, file.path(outdir, "01.subset_precluster.rds"))

log_message("Script completed: scrna_subset_prepare")
write_summary(c(
  "scRNA subset prepare completed.",
  paste0("Subset celltypes: ", paste(subset_celltypes, collapse = ", ")),
  paste0("Subset cells: ", ncol(subset_obj)),
  paste0("Subset genes: ", nrow(subset_obj)),
  paste0("Dims used: ", paste(valid_dims, collapse = ",")),
  "The output object is preclustered and ready for SCRNA_SUBSET_CLUSTER_SCAN.",
  paste0("Object file: 01.subset_precluster.rds\n", describe_rds(file.path(outdir, "01.subset_precluster.rds")))
))

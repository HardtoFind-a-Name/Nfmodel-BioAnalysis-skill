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
mapping_file <- get_arg("--mapping-file")
outdir <- get_arg("--outdir")
target_gene_file <- get_arg("--target-gene-file", required = FALSE, default = "")
target_gene_col <- get_arg("--target-gene-col", required = FALSE, default = "gene")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

subset_cluster_col <- "seurat_clusters"
subset_output_celltype_col <- "subset_celltype_manual"

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_subset_apply")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_subset_apply")
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!file.exists(mapping_file)) save_note_and_stop(paste0("Subset mapping file not found: ", mapping_file))
if (!subset_cluster_col %in% colnames(obj@meta.data)) save_note_and_stop("Subset object lacks seurat_clusters.")
mapping_df <- read.csv(mapping_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("cluster", "manual_celltype", "manual_markers", "annotation_flag")
missing_cols <- setdiff(required_cols, colnames(mapping_df))
if (length(missing_cols) > 0) save_note_and_stop(paste0("Missing required mapping columns: ", paste(missing_cols, collapse = ", ")))

mapping_df$cluster <- as.character(mapping_df$cluster)
mapping_df$manual_celltype <- trimws(as.character(mapping_df$manual_celltype))
mapping_df$manual_markers <- as.character(mapping_df$manual_markers)
mapping_df$annotation_flag <- as.character(mapping_df$annotation_flag)
if (anyDuplicated(mapping_df$cluster) > 0) save_note_and_stop("Duplicated cluster IDs in subset mapping file.")
if (any(!nzchar(mapping_df$manual_celltype))) save_note_and_stop("manual_celltype contains empty values.")
all_clusters <- sort_cluster_ids(obj[[subset_cluster_col, drop = TRUE]])
missing_clusters <- setdiff(all_clusters, mapping_df$cluster)
if (length(missing_clusters) > 0) save_note_and_stop(paste0("Subset clusters missing from mapping file: ", paste(missing_clusters, collapse = ", ")))
mapping_df <- mapping_df[match(all_clusters, mapping_df$cluster), , drop = FALSE]

cluster_map <- setNames(mapping_df$manual_celltype, mapping_df$cluster)
flag_map <- setNames(mapping_df$annotation_flag, mapping_df$cluster)
obj@meta.data[[subset_output_celltype_col]] <- unname(cluster_map[as.character(obj[[subset_cluster_col, drop = TRUE]])])
obj$subset_annotation_flag <- unname(flag_map[as.character(obj[[subset_cluster_col, drop = TRUE]])])
Idents(obj) <- subset_output_celltype_col

write.csv(mapping_df, file.path(outdir, "01_subset_cluster_celltype_mapping_applied.csv"), row.names = FALSE)
celltype_summary <- obj@meta.data %>%
  dplyr::count(.data[[subset_output_celltype_col]], name = "n_cells") %>%
  dplyr::arrange(dplyr::desc(n_cells))
write.csv(celltype_summary, file.path(outdir, "02_subset_celltype_cell_counts.csv"), row.names = FALSE)

save_plot_dual(DimPlot(obj, reduction = "umap", group.by = subset_output_celltype_col, label = TRUE, repel = TRUE) + theme_bw(),
               outdir, "03_subset_umap_by_manual_celltype", width = 10, height = 8)

manual_marker_panel <- unique(unlist(lapply(mapping_df$manual_markers, split_markers)))
manual_marker_panel <- manual_marker_panel[manual_marker_panel %in% rownames(obj)]
if (length(manual_marker_panel) > 0) {
  p <- DotPlot(obj, features = manual_marker_panel, group.by = subset_output_celltype_col) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot_dual(p, outdir, "04_subset_marker_dotplot_after_annotation", width = max(10, length(manual_marker_panel) * 0.35), height = 7)
}

target_genes <- tryCatch(read_target_genes(target_gene_file, target_gene_col), error = function(e) save_note_and_stop(conditionMessage(e)))
target_hits <- intersect(target_genes, rownames(obj))
if (length(target_hits) > 0) {
  p <- DotPlot(obj, features = target_hits, group.by = subset_output_celltype_col) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot_dual(p, outdir, "05_subset_target_gene_dotplot", width = max(8, length(target_hits) * 0.35), height = 6)
  write.csv(data.frame(gene = target_hits), file.path(outdir, "05_subset_target_gene_hits.csv"), row.names = FALSE)
  if (length(target_hits) >= 2) {
    obj <- AddModuleScore(obj, features = list(target_hits), name = "Keygene_avg", assay = DefaultAssay(obj))
    colnames(obj@meta.data)[colnames(obj@meta.data) == "Keygene_avg1"] <- "Keygene_avg"
  }
}

saveRDS(obj, file.path(outdir, "01.subset_annotated.rds"))
log_message("Script completed: scrna_subset_apply")
write_summary(c(
  "scRNA subset manual annotation applied.",
  paste0("Subset cells: ", ncol(obj)),
  paste0("Subset genes: ", nrow(obj)),
  paste0("Subset cell types: ", nrow(celltype_summary)),
  paste0("Target genes matched: ", length(target_hits)),
  "Manual subset labels were added to metadata column subset_celltype_manual.",
  paste0("Object file: 01.subset_annotated.rds\n", describe_rds(file.path(outdir, "01.subset_annotated.rds")))
))

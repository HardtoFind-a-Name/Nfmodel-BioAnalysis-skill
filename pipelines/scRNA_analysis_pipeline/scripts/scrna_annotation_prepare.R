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
annotation_cluster_col <- get_arg("--annotation-cluster-col", required = FALSE, default = "seurat_clusters")
annotation_resolution <- get_arg("--annotation-resolution", required = FALSE, default = "")
target_gene_file <- get_arg("--target-gene-file", required = FALSE, default = "")
target_gene_col <- get_arg("--target-gene-col", required = FALSE, default = "gene")
marker_min_pct <- as.numeric(get_arg("--marker-min-pct", required = FALSE, default = "0.25"))
marker_logfc_threshold <- as.numeric(get_arg("--marker-logfc-threshold", required = FALSE, default = "0.25"))
scRNA_cohort_id <- get_arg("--scrna-cohort-id", required = FALSE, default = "")
disease_name <- get_arg("--disease-name", required = FALSE, default = "")
cancer_name <- get_arg("--cancer-name", required = FALSE, default = "")
annotation_literature_years <- as.integer(get_arg("--annotation-literature-years", required = FALSE, default = "10"))

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_annotation_prepare")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_annotation_prepare")
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!"umap" %in% Reductions(obj)) save_note_and_stop("UMAP reduction not found in input object.")

meta_cols <- colnames(obj@meta.data)
cluster_col <- resolve_cluster_col(meta_cols, annotation_cluster_col, annotation_resolution)
if (is.na(cluster_col) || !cluster_col %in% meta_cols) {
  save_note_and_stop(paste0("Target cluster column not found: ", cluster_col))
}
if ("seurat_clusters" %in% meta_cols && cluster_col != "seurat_clusters") {
  obj$seurat_clusters_original <- obj[["seurat_clusters", drop = TRUE]]
}
obj$seurat_clusters <- factor(obj[[cluster_col, drop = TRUE]], levels = sort_cluster_ids(obj[[cluster_col, drop = TRUE]]))
Idents(obj) <- "seurat_clusters"

cluster_size <- obj@meta.data %>%
  dplyr::count(seurat_clusters, name = "n_cells") %>%
  dplyr::arrange(factor(as.character(seurat_clusters), levels = sort_cluster_ids(seurat_clusters)))
write.csv(cluster_size, file.path(outdir, "01_cluster_size.csv"), row.names = FALSE)
log_message("Running FindAllMarkers")
all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = marker_min_pct, logfc.threshold = marker_logfc_threshold, test.use = "wilcox")
if (nrow(all_markers) == 0) save_note_and_stop("FindAllMarkers returned empty result.")
all_markers <- all_markers %>% dplyr::arrange(cluster, dplyr::desc(avg_log2FC), dplyr::desc(pct.1))
write.csv(all_markers, file.path(outdir, "02_all_markers.csv"), row.names = FALSE)

top10 <- all_markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>% dplyr::ungroup()
top20 <- all_markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(avg_log2FC, n = 20, with_ties = FALSE) %>% dplyr::ungroup()
write.csv(top10, file.path(outdir, "03_top10_markers_by_cluster.csv"), row.names = FALSE)
write.csv(top20, file.path(outdir, "04_top20_markers_by_cluster.csv"), row.names = FALSE)

singleR_labels <- data.frame(cluster = as.character(cluster_size$seurat_clusters), singleR_label = NA_character_)
if (requireNamespace("SingleR", quietly = TRUE) &&
    requireNamespace("celldex", quietly = TRUE) &&
    requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  log_message("Running SingleR cluster-level annotation")
  singleR_labels <- tryCatch({
    ref <- celldex::HumanPrimaryCellAtlasData()
    expr <- as.matrix(get_assay_data_compat(obj, assay = DefaultAssay(obj), layer = "data"))
    pred <- SingleR::SingleR(test = expr, ref = ref, labels = ref$label.main, clusters = obj$seurat_clusters)
    data.frame(cluster = rownames(pred), singleR_label = pred$labels, stringsAsFactors = FALSE)
  }, error = function(e) {
    log_message("[警告] SingleR failed:", conditionMessage(e))
    data.frame(cluster = as.character(cluster_size$seurat_clusters), singleR_label = NA_character_)
  })
}
write.csv(singleR_labels, file.path(outdir, "05_singler_cluster_labels.csv"), row.names = FALSE)

marker_panel <- top10 %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(manual_markers = paste(unique(gene), collapse = ","), .groups = "drop") %>%
  dplyr::mutate(cluster = as.character(cluster))
template_context <- cluster_size %>%
  dplyr::transmute(cluster = as.character(seurat_clusters), n_cells = n_cells) %>%
  dplyr::left_join(singleR_labels, by = "cluster") %>%
  dplyr::left_join(marker_panel, by = "cluster")
template <- template_context %>%
  dplyr::transmute(
    cluster = cluster,
    manual_celltype = "",
    manual_markers = manual_markers %||% "",
    annotation_flag = "manual_pending",
    annotation_confidence = "",
    notes = ""
  )
write.csv(template, file.path(outdir, "08_cluster_celltype_mapping_template.csv"), row.names = FALSE)

literature_template <- data.frame(
  manual_markers = unique(unlist(strsplit(paste(na.omit(template$manual_markers), collapse = ","), ","))),
  references = "",
  stringsAsFactors = FALSE
)
literature_template$manual_markers <- trimws(literature_template$manual_markers)
literature_template <- literature_template[!is.na(literature_template$manual_markers) & nzchar(literature_template$manual_markers), , drop = FALSE]
if (nrow(literature_template) == 0) literature_template <- data.frame(manual_markers = character(), references = character())
write.csv(literature_template, file.path(outdir, "09_literature_marker_reference_template.csv"), row.names = FALSE)

json_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\"", x)
  x <- gsub("\n", "\\n", x)
  x
}
write_json_context <- function(path, context, clusters, top_markers) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  cat("{\n", file = con)
  keys <- names(context)
  for (i in seq_along(keys)) {
    comma <- ","
    cat(sprintf('  "%s": "%s"%s\n', keys[i], json_escape(context[[keys[i]]]), comma), file = con)
  }
  cat('  "clusters": [\n', file = con)
  if (nrow(clusters) > 0) {
    for (i in seq_len(nrow(clusters))) {
      comma <- if (i < nrow(clusters)) "," else ""
      cat(sprintf('    {"cluster": "%s", "n_cells": %s, "singleR_label": "%s", "top_markers": "%s"}%s\n',
                  json_escape(clusters$cluster[i]), clusters$n_cells[i], json_escape(clusters$singleR_label[i]), json_escape(clusters$manual_markers[i]), comma), file = con)
    }
  }
  cat('  ],\n', file = con)
  cat('  "top_markers": [\n', file = con)
  if (nrow(top_markers) > 0) {
    keep_cols <- intersect(c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), colnames(top_markers))
    tm <- top_markers[, keep_cols, drop = FALSE]
    for (i in seq_len(nrow(tm))) {
      fields <- paste(sprintf('"%s": "%s"', names(tm), vapply(tm, function(col) json_escape(col[i]), character(1))), collapse = ", ")
      comma <- if (i < nrow(tm)) "," else ""
      cat(sprintf('    {%s}%s\n', fields, comma), file = con)
    }
  }
  cat('  ]\n', file = con)
  cat("}\n", file = con)
}
agent_context <- list(
  annotation_stage = "main_pass_1",
  scRNA_cohort_id = scRNA_cohort_id,
  disease_name = disease_name,
  cancer_name = cancer_name,
  annotation_literature_years = annotation_literature_years,
  mapping_template = "08_cluster_celltype_mapping_template.csv",
  literature_marker_template = "09_literature_marker_reference_template.csv",
  marker_table = "02_all_markers.csv",
  top_marker_table = "03_top10_markers_by_cluster.csv"
)
write_json_context(file.path(outdir, "10_annotation_agent_context.json"), agent_context, template_context, top10)

target_genes <- tryCatch(read_target_genes(target_gene_file, target_gene_col), error = function(e) save_note_and_stop(conditionMessage(e)))
target_hits <- intersect(target_genes, rownames(obj))
if (length(target_hits) > 0) {
  dot <- DotPlot(obj, features = target_hits, group.by = "seurat_clusters") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot_dual(dot, outdir, "09_target_gene_dotplot_by_cluster", width = max(8, length(target_hits) * 0.35), height = 6)
  write.csv(data.frame(gene = target_hits), file.path(outdir, "09_target_gene_hits.csv"), row.names = FALSE)
} else if (length(target_genes) > 0) {
  log_message("[警告] No target genes matched the Seurat object")
}

prepare_res <- list(
  cluster_col_used = cluster_col,
  cluster_size = cluster_size,
  all_markers = all_markers,
  top10_markers = top10,
  top20_markers = top20,
  singleR_labels = singleR_labels,
  mapping_template = template,
  literature_marker_reference_template = literature_template,
  target_gene_hits = target_hits
)
saveRDS(obj, file.path(outdir, "01.annotation_prepare.rds"))
saveRDS(prepare_res, file.path(outdir, "02.annotation_prepare_results.rds"))

log_message("Script completed: scrna_annotation_prepare")
write_summary(c(
  "scRNA annotation preparation completed.",
  paste0("Final cells: ", ncol(obj)),
  paste0("Final genes: ", nrow(obj)),
  paste0("Cluster column used: ", cluster_col),
  paste0("Cluster count: ", nrow(cluster_size)),
  paste0("Markers detected: ", nrow(all_markers)),
  paste0("Target genes matched: ", length(target_hits)),
  "Manual annotation templates and agent context were generated for literature-backed annotation.",
  paste0("Object file: 01.annotation_prepare.rds\n", describe_rds(file.path(outdir, "01.annotation_prepare.rds"))),
  "",
  "限制与说明",
  "Manual cell type labels must be filled into 08_cluster_celltype_mapping_template.csv and supplied as mapping_file for annotation_apply."
))

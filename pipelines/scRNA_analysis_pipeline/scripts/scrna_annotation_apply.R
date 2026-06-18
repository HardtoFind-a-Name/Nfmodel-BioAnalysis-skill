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
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")
target_gene_file <- get_arg("--target-gene-file", required = FALSE, default = "")
target_gene_col <- get_arg("--target-gene-col", required = FALSE, default = "gene")
annotation_marker_reference_file <- get_arg("--annotation-marker-reference-file", required = FALSE, default = "")
annotation_min_support_markers <- as.integer(get_arg("--annotation-min-support-markers", required = FALSE, default = "2"))
annotation_allow_low_support <- as_bool(get_arg("--annotation-allow-low-support", required = FALSE, default = "false"))

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_annotation_apply")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_annotation_apply")
log_message("Input RDS:", input_rds)
log_message("Mapping file:", mapping_file)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!file.exists(mapping_file)) save_note_and_stop(paste0("Mapping file not found: ", mapping_file))
mapping_df <- read.csv(mapping_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("cluster", "manual_celltype", "manual_markers", "annotation_flag", "annotation_confidence", "notes")
missing_cols <- setdiff(required_cols, colnames(mapping_df))
if (length(missing_cols) > 0) save_note_and_stop(paste0("Missing required mapping columns: ", paste(missing_cols, collapse = ", ")))

mapping_df$cluster <- as.character(mapping_df$cluster)
mapping_df$manual_celltype <- trimws(as.character(mapping_df$manual_celltype))
mapping_df$manual_markers <- as.character(mapping_df$manual_markers)
mapping_df$annotation_flag <- as.character(mapping_df$annotation_flag)
if (anyDuplicated(mapping_df$cluster) > 0) save_note_and_stop("Duplicated cluster IDs in mapping file.")
if (any(!nzchar(mapping_df$manual_celltype))) save_note_and_stop("manual_celltype contains empty values.")

all_clusters <- sort_cluster_ids(obj$seurat_clusters)
missing_clusters <- setdiff(all_clusters, mapping_df$cluster)
if (length(missing_clusters) > 0) save_note_and_stop(paste0("Clusters missing from mapping file: ", paste(missing_clusters, collapse = ", ")))
mapping_df <- mapping_df[match(all_clusters, mapping_df$cluster), , drop = FALSE]

reference_markers <- character()
if (!is_null_path(annotation_marker_reference_file)) {
  if (!file.exists(annotation_marker_reference_file)) save_note_and_stop(paste0("Annotation marker reference file not found: ", annotation_marker_reference_file))
  ref_df <- read.csv(annotation_marker_reference_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("manual_markers", "references") %in% colnames(ref_df))) save_note_and_stop("Reference table must contain manual_markers and references columns.")
  valid_ref <- !is.na(ref_df$references) & nzchar(trimws(as.character(ref_df$references)))
  reference_markers <- unique(toupper(unlist(lapply(ref_df$manual_markers[valid_ref], split_markers))))
}
expr_data <- get_assay_data_compat(obj, assay = DefaultAssay(obj), layer = "data")
audit_rows <- list()
for (i in seq_len(nrow(mapping_df))) {
  cl <- mapping_df$cluster[i]
  markers <- split_markers(mapping_df$manual_markers[i])
  markers_in_obj <- markers[markers %in% rownames(obj)]
  if (length(reference_markers) > 0) {
    markers_supported <- markers_in_obj[toupper(markers_in_obj) %in% reference_markers]
  } else {
    markers_supported <- markers_in_obj
  }
  cells <- colnames(obj)[as.character(obj$seurat_clusters) == cl]
  if (length(cells) > 0 && length(markers_supported) > 0) {
    pct_expr <- Matrix::rowMeans(expr_data[markers_supported, cells, drop = FALSE] > 0) * 100
    expressed_supported <- names(pct_expr)[pct_expr > 0]
  } else {
    pct_expr <- numeric()
    expressed_supported <- character()
  }
  support_count <- length(expressed_supported)
  passed <- support_count >= annotation_min_support_markers
  audit_rows[[cl]] <- data.frame(
    cluster = cl,
    manual_celltype = mapping_df$manual_celltype[i],
    manual_markers = paste(markers, collapse = ","),
    markers_in_object = paste(markers_in_obj, collapse = ","),
    literature_supported_markers = paste(markers_supported, collapse = ","),
    expressed_supported_markers = paste(expressed_supported, collapse = ","),
    support_count = support_count,
    min_support_required = annotation_min_support_markers,
    passed = passed,
    stringsAsFactors = FALSE
  )
}
support_audit <- dplyr::bind_rows(audit_rows)
write.csv(support_audit, file.path(outdir, "annotation_marker_support_audit.csv"), row.names = FALSE)
rejected <- support_audit[!support_audit$passed, , drop = FALSE]
write.csv(rejected, file.path(outdir, "annotation_rejected_or_low_support_clusters.csv"), row.names = FALSE)
if (nrow(rejected) > 0) {
  if (!annotation_allow_low_support) {
    save_note_and_stop(paste0("Annotation marker support audit failed for clusters: ", paste(rejected$cluster, collapse = ", "),
                              ". Provide at least ", annotation_min_support_markers, " supported markers or set annotation_allow_low_support=true."))
  }
  low_clusters <- rejected$cluster
  mapping_df$annotation_flag[mapping_df$cluster %in% low_clusters] <- "low_support"
}

cluster_map <- setNames(mapping_df$manual_celltype, mapping_df$cluster)
flag_map <- setNames(mapping_df$annotation_flag, mapping_df$cluster)
obj$celltype_manual <- unname(cluster_map[as.character(obj$seurat_clusters)])
obj$annotation_flag <- unname(flag_map[as.character(obj$seurat_clusters)])
if ("annotation_confidence" %in% colnames(mapping_df)) {
  conf_map <- setNames(mapping_df$annotation_confidence, mapping_df$cluster)
  obj$annotation_confidence <- unname(conf_map[as.character(obj$seurat_clusters)])
}
Idents(obj) <- "celltype_manual"

write.csv(mapping_df, file.path(outdir, "01_cluster_celltype_mapping_applied.csv"), row.names = FALSE)
celltype_summary <- obj@meta.data %>%
  dplyr::count(celltype_manual, name = "n_cells") %>%
  dplyr::arrange(dplyr::desc(n_cells))
write.csv(celltype_summary, file.path(outdir, "02_celltype_cell_counts.csv"), row.names = FALSE)

save_plot_dual(DimPlot(obj, reduction = "umap", group.by = "celltype_manual", label = TRUE, repel = TRUE) + theme_bw(),
               outdir, "03_umap_by_celltype_manual", width = 10, height = 8)
save_plot_dual(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + theme_bw(),
               outdir, "04_umap_by_cluster_after_annotation", width = 9, height = 7)

manual_marker_panel <- unique(unlist(lapply(mapping_df$manual_markers, split_markers)))
manual_marker_panel <- manual_marker_panel[manual_marker_panel %in% rownames(obj)]
if (length(manual_marker_panel) > 0) {
  p <- DotPlot(obj, features = manual_marker_panel, group.by = "celltype_manual") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot_dual(p, outdir, "05_marker_dotplot_after_annotation", width = max(10, length(manual_marker_panel) * 0.3), height = 7)
}

target_genes <- tryCatch(read_target_genes(target_gene_file, target_gene_col), error = function(e) save_note_and_stop(conditionMessage(e)))
target_hits <- intersect(target_genes, rownames(obj))
if (length(target_hits) > 0) {
  p <- DotPlot(obj, features = target_hits, group.by = "celltype_manual") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot_dual(p, outdir, "06_target_gene_dotplot_by_celltype", width = max(8, length(target_hits) * 0.35), height = 6)
  write.csv(data.frame(gene = target_hits), file.path(outdir, "06_target_gene_hits.csv"), row.names = FALSE)
}

saveRDS(obj, file.path(outdir, "01.seurat_annotated.rds"))

log_message("Script completed: scrna_annotation_apply")
write_summary(c(
  "scRNA manual annotation applied.",
  paste0("Final cells: ", ncol(obj)),
  paste0("Final genes: ", nrow(obj)),
  paste0("Annotated cell types: ", nrow(celltype_summary)),
  paste0("Target genes matched: ", length(target_hits)),
  paste0("Marker support audit passed clusters: ", sum(support_audit$passed), "/", nrow(support_audit)),
  "Manual cell type labels were added to metadata column celltype_manual.",
  paste0("Object file: 01.seurat_annotated.rds\n", describe_rds(file.path(outdir, "01.seurat_annotated.rds"))),
  "",
  "限制与说明",
  "Cell type labels depend on the supplied manual mapping file and should be reviewed before downstream biological interpretation."
))

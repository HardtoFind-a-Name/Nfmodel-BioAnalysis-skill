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
key_celltype <- get_arg("--key-celltype")
safe_label <- get_arg("--safe-label")
scRNA_cohort_id <- get_arg("--scrna-cohort-id", required = FALSE, default = "")
disease_name <- get_arg("--disease-name", required = FALSE, default = "")
cancer_name <- get_arg("--cancer-name", required = FALSE, default = "")
marker_min_pct <- as.numeric(get_arg("--marker-min-pct", required = FALSE, default = "0.25"))
marker_logfc_threshold <- as.numeric(get_arg("--marker-logfc-threshold", required = FALSE, default = "0.25"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_subset_annotation_prepare", subdir = safe_label)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_subset_annotation_prepare")
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!"seurat_clusters" %in% colnames(obj@meta.data)) save_note_and_stop("seurat_clusters missing from subset clustered object.")
Idents(obj) <- "seurat_clusters"
cluster_size <- obj@meta.data %>% dplyr::count(seurat_clusters, name = "n_cells") %>% dplyr::arrange(seurat_clusters)
write.csv(cluster_size, file.path(outdir, "01_subset_cluster_size.csv"), row.names = FALSE)
all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = marker_min_pct, logfc.threshold = marker_logfc_threshold, test.use = "wilcox")
if (nrow(all_markers) == 0) save_note_and_stop("FindAllMarkers returned empty result for subset.")
all_markers <- all_markers %>% dplyr::arrange(cluster, dplyr::desc(avg_log2FC), dplyr::desc(pct.1))
write.csv(all_markers, file.path(outdir, "02_subset_all_markers.csv"), row.names = FALSE)
top10 <- all_markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>% dplyr::ungroup()
top20 <- all_markers %>% dplyr::group_by(cluster) %>% dplyr::slice_max(avg_log2FC, n = 20, with_ties = FALSE) %>% dplyr::ungroup()
write.csv(top10, file.path(outdir, "03_subset_top10_markers_by_cluster.csv"), row.names = FALSE)
write.csv(top20, file.path(outdir, "04_subset_top20_markers_by_cluster.csv"), row.names = FALSE)
marker_panel <- top10 %>% dplyr::group_by(cluster) %>% dplyr::summarise(manual_markers = paste(unique(gene), collapse = ","), .groups = "drop") %>% dplyr::mutate(cluster = as.character(cluster))
template <- cluster_size %>%
  dplyr::transmute(cluster = as.character(seurat_clusters)) %>%
  dplyr::left_join(marker_panel, by = "cluster") %>%
  dplyr::transmute(cluster = cluster, manual_celltype = "", manual_markers = manual_markers %||% "",
                   annotation_flag = "manual_pending", annotation_confidence = "", notes = "")
write.csv(template, file.path(outdir, "08_subset_cluster_celltype_mapping_template.csv"), row.names = FALSE)
lit <- data.frame(manual_markers = unique(unlist(strsplit(paste(na.omit(template$manual_markers), collapse = ","), ","))), references = "", stringsAsFactors = FALSE)
lit$manual_markers <- trimws(lit$manual_markers)
lit <- lit[!is.na(lit$manual_markers) & nzchar(lit$manual_markers), , drop = FALSE]
write.csv(lit, file.path(outdir, "09_subset_literature_marker_reference_template.csv"), row.names = FALSE)

json_escape <- function(x) { x <- as.character(x %||% ""); x <- gsub('\\\\', '\\\\\\\\', x); x <- gsub('"', '\\\\"', x); x }
con <- file(file.path(outdir, "10_subset_annotation_agent_context.json"), "w")
cat("{\n", file = con)
cat(sprintf('  "annotation_stage": "subset_pass_2",\n  "key_celltype": "%s",\n  "safe_label": "%s",\n  "scRNA_cohort_id": "%s",\n  "disease_name": "%s",\n  "cancer_name": "%s",\n  "mapping_template": "08_subset_cluster_celltype_mapping_template.csv",\n  "literature_marker_template": "09_subset_literature_marker_reference_template.csv",\n  "marker_table": "02_subset_all_markers.csv"\n', json_escape(key_celltype), json_escape(safe_label), json_escape(scRNA_cohort_id), json_escape(disease_name), json_escape(cancer_name)), file = con)
cat("}\n", file = con)
close(con)
saveRDS(obj, file.path(outdir, "01.subset_annotation_prepare.rds"))

log_message("Script completed: scrna_subset_annotation_prepare")
write_summary(c("scRNA subset annotation preparation completed.", paste0("Key celltype: ", key_celltype), paste0("Subset clusters: ", nrow(cluster_size)), paste0("Markers detected: ", nrow(all_markers))))

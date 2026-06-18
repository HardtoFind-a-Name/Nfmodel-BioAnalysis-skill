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
key_celltype <- get_arg("--key-celltype", required = FALSE, default = "T/NK cells")
top_n_pathways <- as.integer(get_arg("--top-n-pathways", required = FALSE, default = "10"))
min_fdr <- as.numeric(get_arg("--min-fdr", required = FALSE, default = "0.05"))
min_delta <- as.numeric(get_arg("--min-delta", required = FALSE, default = "0"))
highlight_groups <- parse_chr_vector(get_arg("--highlight-groups", required = FALSE, default = "T/NK cells,Monocytes,Macrophages,B cells,Epithelial cells"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

celltype_col <- "subset_celltype_manual"
ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_scmetabolism_specific")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_scmetabolism_specific")
require_package_or_stop("scMetabolism", save_note_and_stop)
require_package_or_stop("AUCell", save_note_and_stop)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!celltype_col %in% colnames(obj@meta.data)) save_note_and_stop("subset_celltype_manual is missing. Run subset apply first.")
if (!key_celltype %in% unique(as.character(obj@meta.data[[celltype_col]]))) save_note_and_stop(paste0("Key celltype not found: ", key_celltype))

obj2 <- tryCatch(
  scMetabolism::sc.metabolism.Seurat(obj = obj, method = "AUCell", imputation = FALSE, metabolism.type = "KEGG"),
  error = function(e) save_note_and_stop(paste0("scMetabolism scoring failed: ", conditionMessage(e)))
)
score_mat <- as.matrix(obj2@assays$METABOLISM$score)
meta <- obj2@meta.data
common_cells <- intersect(colnames(score_mat), rownames(meta))
score_mat <- score_mat[, common_cells, drop = FALSE]
meta <- meta[common_cells, , drop = FALSE]
key_cells <- rownames(meta)[as.character(meta[[celltype_col]]) == key_celltype]
other_cells <- setdiff(rownames(meta), key_cells)
if (length(key_cells) == 0 || length(other_cells) == 0) save_note_and_stop("Need both key and non-key cells for specificity testing.")

stats <- lapply(rownames(score_mat), function(pw) {
  x <- as.numeric(score_mat[pw, key_cells])
  y <- as.numeric(score_mat[pw, other_cells])
  p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
  data.frame(pathway = pw, mean_key = mean(x, na.rm = TRUE), mean_other = mean(y, na.rm = TRUE),
             delta_mean = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE), p_value = p)
})
stats <- dplyr::bind_rows(stats)
stats$fdr <- p.adjust(stats$p_value, method = "BH")
stats <- stats %>% dplyr::arrange(fdr, dplyr::desc(delta_mean))
stats$rank_group <- ifelse(!is.na(stats$fdr) & stats$fdr < min_fdr & stats$delta_mean > min_delta,
                           "significant_up_in_key", "other")
write.csv(stats, file.path(outdir, "01_pathway_ranking_stats_key_vs_others.csv"), row.names = FALSE)

top_tbl <- stats %>% dplyr::filter(rank_group == "significant_up_in_key") %>% dplyr::slice_head(n = top_n_pathways)
if (nrow(top_tbl) < top_n_pathways) {
  top_tbl <- stats %>% dplyr::arrange(dplyr::desc(delta_mean)) %>% dplyr::slice_head(n = top_n_pathways)
}
write.csv(top_tbl, file.path(outdir, "02_top_specific_pathways.csv"), row.names = FALSE)

groups <- unique(as.character(meta[[celltype_col]]))
preferred_groups <- c(highlight_groups, setdiff(groups, highlight_groups))
preferred_groups <- preferred_groups[preferred_groups %in% groups]
plot_list <- lapply(top_tbl$pathway, function(pw) {
  data.frame(pathway = pw, score = as.numeric(score_mat[pw, ]), group = as.character(meta[[celltype_col]]))
})
plot_df <- dplyr::bind_rows(plot_list)
summary_df <- plot_df %>% dplyr::group_by(group, pathway) %>% dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
summary_df$group <- factor(summary_df$group, levels = preferred_groups)
write.csv(summary_df, file.path(outdir, "03_group_mean_scores_top_pathways.csv"), row.names = FALSE)

p_bar <- ggplot(top_tbl, aes(x = delta_mean, y = reorder(pathway, delta_mean))) + geom_col(fill = "#4C78A8") + theme_bw()
save_plot_dual(p_bar, outdir, "04_barplot_top_specific_pathways", width = 8, height = max(5, nrow(top_tbl) * 0.35))
p_heat <- ggplot(summary_df, aes(x = group, y = pathway, fill = mean_score)) +
  geom_tile() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_dual(p_heat, outdir, "05_heatmap_top_specific_pathways", width = 10, height = max(5, nrow(top_tbl) * 0.35))
saveRDS(obj2, file.path(outdir, "01.scmetabolism_specific_scored.rds"))

log_message("Script completed: scrna_scmetabolism_specific")
write_summary(c(
  "scRNA scMetabolism key-celltype specificity analysis completed.",
  paste0("Key celltype: ", key_celltype),
  paste0("Top pathways: ", nrow(top_tbl)),
  paste0("FDR cutoff: ", min_fdr),
  paste0("Delta cutoff: ", min_delta)
))

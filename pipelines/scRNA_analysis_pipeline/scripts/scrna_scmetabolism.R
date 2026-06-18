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
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

celltype_col <- "subset_celltype_manual"
ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_scmetabolism")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_scmetabolism")
require_package_or_stop("scMetabolism", save_note_and_stop)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!celltype_col %in% colnames(obj@meta.data)) save_note_and_stop("subset_celltype_manual is missing. Run subset apply first.")

obj2 <- tryCatch(
  scMetabolism::sc.metabolism.Seurat(obj = obj, method = "AUCell", imputation = FALSE, metabolism.type = "KEGG"),
  error = function(e) save_note_and_stop(paste0("scMetabolism scoring failed: ", conditionMessage(e)))
)
score <- obj2@assays$METABOLISM$score
score_df <- as.data.frame(t(as.matrix(score)), check.names = FALSE)
score_df$cell_id <- rownames(score_df)
score_df[[celltype_col]] <- obj2@meta.data[rownames(score_df), celltype_col, drop = TRUE]
write.csv(score_df, file.path(outdir, "01_scmetabolism_cell_scores.csv"), row.names = FALSE)

pathway_cols <- setdiff(colnames(score_df), c("cell_id", celltype_col))
long_df <- data.frame(
  celltype = rep(score_df[[celltype_col]], times = length(pathway_cols)),
  pathway = rep(pathway_cols, each = nrow(score_df)),
  score = as.numeric(as.matrix(score_df[, pathway_cols, drop = FALSE])),
  stringsAsFactors = FALSE
)
colnames(long_df)[1] <- celltype_col
summary_df <- long_df %>%
  dplyr::group_by(.data[[celltype_col]], pathway) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
write.csv(summary_df, file.path(outdir, "02_scmetabolism_mean_scores_by_celltype.csv"), row.names = FALSE)

if (key_celltype %in% summary_df[[celltype_col]]) {
  top_paths <- summary_df %>%
    dplyr::filter(.data[[celltype_col]] == key_celltype) %>%
    dplyr::arrange(dplyr::desc(mean_score)) %>%
    dplyr::slice_head(n = top_n_pathways)
  write.csv(top_paths, file.path(outdir, "03_scmetabolism_top_pathways_in_key_celltype.csv"), row.names = FALSE)
  plot_df <- summary_df %>% dplyr::filter(pathway %in% top_paths$pathway)
  p <- ggplot(plot_df, aes(x = .data[[celltype_col]], y = pathway, fill = mean_score)) +
    geom_tile() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot_dual(p, outdir, "04_scmetabolism_top_pathway_heatmap", width = 10, height = 7)
}

saveRDS(obj2, file.path(outdir, "01.scmetabolism_scored.rds"))
log_message("Script completed: scrna_scmetabolism")
write_summary(c(
  "scRNA scMetabolism analysis completed.",
  paste0("Cells: ", ncol(obj2)),
  paste0("Pathways: ", nrow(score)),
  paste0("Key celltype: ", key_celltype),
  paste0("Object file: 01.scmetabolism_scored.rds\n", describe_rds(file.path(outdir, "01.scmetabolism_scored.rds")))
))

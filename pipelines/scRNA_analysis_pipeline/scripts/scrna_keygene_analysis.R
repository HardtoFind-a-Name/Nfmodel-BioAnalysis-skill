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
target_gene_file <- get_arg("--target-gene-file", required = FALSE, default = "")
target_gene_col <- get_arg("--target-gene-col", required = FALSE, default = "gene")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

celltype_col <- "celltype_manual"
group_col <- "group"
ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_keygene_analysis")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

empty_outputs <- function(note) {
  dot_empty <- data.frame(
    features.plot = character(),
    id = character(),
    avg.exp = numeric(),
    pct.exp = numeric(),
    avg.exp.scaled = numeric(),
    stringsAsFactors = FALSE
  )
  diff_empty <- data.frame(
    celltype = character(),
    gene = character(),
    control_group = character(),
    disease_group = character(),
    mean_control = numeric(),
    mean_disease = numeric(),
    pct_expr = numeric(),
    log2FC = numeric(),
    p_value = numeric(),
    p_adjust = numeric(),
    stringsAsFactors = FALSE
  )
  violin_empty <- data.frame(
    cell = character(),
    gene = character(),
    expression = numeric(),
    celltype = character(),
    group = character(),
    stringsAsFactors = FALSE
  )
  write.csv(dot_empty, file.path(outdir, "01_keygene_expr_dotplot_data.csv"), row.names = FALSE)
  write.csv(diff_empty, file.path(outdir, "02_keygene_group_diff.csv"), row.names = FALSE)
  write.csv(violin_empty, file.path(outdir, "03_keygene_violin_data.csv"), row.names = FALSE)
  write_summary(c("scRNA keygene analysis skipped.", note))
  quit(save = "no", status = 0)
}

log_message("Script started: scrna_keygene_analysis")
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
missing_cols <- setdiff(c(celltype_col, group_col), colnames(obj@meta.data))
if (length(missing_cols) > 0) save_note_and_stop(paste0("Missing metadata columns: ", paste(missing_cols, collapse = ", ")))
if (is_null_path(target_gene_file)) empty_outputs("target_gene_file is not supplied.")
target_genes <- tryCatch(read_target_genes(target_gene_file, target_gene_col), error = function(e) save_note_and_stop(conditionMessage(e)))
target_hits <- intersect(target_genes, rownames(obj))
if (length(target_hits) == 0) empty_outputs("No target genes matched the Seurat object.")

Idents(obj) <- celltype_col
p_dot <- DotPlot(obj, features = target_hits, group.by = celltype_col, assay = DefaultAssay(obj)) +
  RotatedAxis() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_dual(p_dot, outdir, "01_keygene_expr_dotplot", width = max(8, length(target_hits) * 0.35), height = 6)
dot_data <- p_dot$data
write.csv(dot_data, file.path(outdir, "01_keygene_expr_dotplot_data.csv"), row.names = FALSE)

expr <- get_assay_data_compat(obj, assay = DefaultAssay(obj), layer = "data")
meta <- obj@meta.data
meta[[celltype_col]] <- as.character(meta[[celltype_col]])
meta[[group_col]] <- as.character(meta[[group_col]])
groups <- sort(unique(meta[[group_col]][!is.na(meta[[group_col]]) & nzchar(meta[[group_col]])]))
control_group <- if ("normal" %in% groups) "normal" else groups[1]
disease_group <- if ("tumor" %in% groups) "tumor" else if (length(groups) >= 2) setdiff(groups, control_group)[1] else NA_character_
if (is.na(disease_group)) save_note_and_stop("Keygene group differential analysis requires at least two group values.")

rows <- list()
for (ct in sort(unique(meta[[celltype_col]]))) {
  ct_cells <- rownames(meta)[meta[[celltype_col]] == ct]
  for (gene in target_hits) {
    control_cells <- intersect(ct_cells, rownames(meta)[meta[[group_col]] == control_group])
    disease_cells <- intersect(ct_cells, rownames(meta)[meta[[group_col]] == disease_group])
    x <- if (length(control_cells) > 0) as.numeric(expr[gene, control_cells]) else numeric()
    y <- if (length(disease_cells) > 0) as.numeric(expr[gene, disease_cells]) else numeric()
    pval <- if (length(x) > 0 && length(y) > 0) tryCatch(wilcox.test(y, x)$p.value, error = function(e) NA_real_) else NA_real_
    mean_control <- mean(x, na.rm = TRUE)
    mean_disease <- mean(y, na.rm = TRUE)
    pct_expr <- mean(as.numeric(expr[gene, ct_cells]) > 0, na.rm = TRUE) * 100
    rows[[paste(ct, gene, sep = "__")]] <- data.frame(
      celltype = ct,
      gene = gene,
      control_group = control_group,
      disease_group = disease_group,
      mean_control = mean_control,
      mean_disease = mean_disease,
      pct_expr = pct_expr,
      log2FC = log2((mean_disease + 1e-6) / (mean_control + 1e-6)),
      p_value = pval,
      stringsAsFactors = FALSE
    )
  }
}
diff_df <- dplyr::bind_rows(rows)
diff_df$p_adjust <- p.adjust(diff_df$p_value, method = "BH")
diff_df <- diff_df %>% dplyr::arrange(p_adjust, dplyr::desc(abs(log2FC)))
write.csv(diff_df, file.path(outdir, "02_keygene_group_diff.csv"), row.names = FALSE)

expr_long <- list()
for (gene in target_hits) {
  expr_long[[gene]] <- data.frame(
    cell = colnames(obj),
    gene = gene,
    expression = as.numeric(expr[gene, ]),
    celltype = meta[colnames(obj), celltype_col, drop = TRUE],
    group = meta[colnames(obj), group_col, drop = TRUE],
    stringsAsFactors = FALSE
  )
}
violin_df <- dplyr::bind_rows(expr_long)
write.csv(violin_df, file.path(outdir, "03_keygene_violin_data.csv"), row.names = FALSE)
p_violin <- ggplot(violin_df, aes(x = celltype, y = expression, fill = group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, position = position_dodge(width = 0.9), outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_dual(p_violin, outdir, "03_keygene_violin", width = 12, height = max(5, 3 * ceiling(length(target_hits) / 2)))

if (length(target_hits) >= 2) {
  obj <- AddModuleScore(obj, features = list(target_hits), name = "Keygene_avg", assay = DefaultAssay(obj))
  colnames(obj@meta.data)[colnames(obj@meta.data) == "Keygene_avg1"] <- "Keygene_avg"
  score_df <- obj@meta.data[, c(celltype_col, group_col, "Keygene_avg"), drop = FALSE]
  write.csv(score_df, file.path(outdir, "04_keygene_modulescore_data.csv"), row.names = TRUE)
}
saveRDS(obj, file.path(outdir, "01.keygene_annotated.rds"))

log_message("Script completed: scrna_keygene_analysis")
write_summary(c(
  "scRNA keygene analysis completed.",
  paste0("Target genes supplied: ", length(target_genes)),
  paste0("Target genes matched: ", length(target_hits)),
  paste0("Cell types tested: ", length(unique(meta[[celltype_col]]))),
  paste0("Groups: ", control_group, ", ", disease_group)
))

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
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
group_col <- get_arg("--group-col", required = FALSE, default = "group")
celltype_col <- get_arg("--celltype-col", required = FALSE, default = "celltype_manual")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_cell_enrichment")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_cell_enrichment")
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
missing_cols <- setdiff(c(group_col, celltype_col), colnames(obj@meta.data))
if (length(missing_cols) > 0) save_note_and_stop(paste0("Missing metadata columns: ", paste(missing_cols, collapse = ", ")))

meta <- obj@meta.data
meta[[group_col]] <- as.character(meta[[group_col]])
meta[[celltype_col]] <- as.character(meta[[celltype_col]])
meta <- meta[!is.na(meta[[group_col]]) & nzchar(meta[[group_col]]) &
               !is.na(meta[[celltype_col]]) & nzchar(meta[[celltype_col]]), , drop = FALSE]
if (nrow(meta) == 0) save_note_and_stop("No valid cells after filtering group/celltype metadata.")

counts <- meta %>%
  dplyr::count(.data[[group_col]], .data[[celltype_col]], name = "cell_count") %>%
  dplyr::rename(group = .data[[group_col]], celltype = .data[[celltype_col]])
write.csv(counts, file.path(outdir, "01_celltype_counts.csv"), row.names = FALSE)

props <- counts %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(group_total = sum(cell_count), proportion = cell_count / group_total) %>%
  dplyr::ungroup()
write.csv(props, file.path(outdir, "02_celltype_proportions.csv"), row.names = FALSE)

groups <- sort(unique(meta[[group_col]]))
celltypes <- sort(unique(meta[[celltype_col]]))
test_rows <- list()
if (length(groups) == 2) {
  for (ct in celltypes) {
    tab <- table(meta[[group_col]], meta[[celltype_col]] == ct)
    ft <- fisher.test(tab)
    test_rows[[ct]] <- data.frame(celltype = ct, method = "fisher", p_value = ft$p.value, odds_ratio = unname(ft$estimate))
  }
} else {
  for (ct in celltypes) {
    tab <- table(meta[[group_col]], meta[[celltype_col]] == ct)
    cs <- suppressWarnings(chisq.test(tab))
    test_rows[[ct]] <- data.frame(celltype = ct, method = "chisq", p_value = cs$p.value, odds_ratio = NA_real_)
  }
}
tests <- dplyr::bind_rows(test_rows)
tests$p_adjust <- p.adjust(tests$p_value, method = "BH")
tests <- tests %>% dplyr::arrange(p_adjust)
write.csv(tests, file.path(outdir, "03_celltype_enrichment_tests.csv"), row.names = FALSE)

p_count <- ggplot(props, aes(x = group, y = cell_count, fill = celltype)) +
  geom_col() + theme_bw() + labs(x = "Group", y = "Cell count", fill = "Cell type")
save_plot_dual(p_count, outdir, "04_celltype_count_barplot", width = 9, height = 6)
p_prop <- ggplot(props, aes(x = group, y = proportion, fill = celltype)) +
  geom_col(position = "fill") + theme_bw() + labs(x = "Group", y = "Proportion", fill = "Cell type")
save_plot_dual(p_prop, outdir, "05_celltype_proportion_barplot", width = 9, height = 6)

top_line <- if (nrow(tests) > 0) {
  paste0("Most different cell type: ", tests$celltype[1], ", adjusted p=", signif(tests$p_adjust[1], 3))
} else {
  "No enrichment tests were available."
}

log_message("Script completed: scrna_cell_enrichment")
write_summary(c(
  "scRNA cell type enrichment completed.",
  paste0("Final cells included: ", nrow(meta)),
  paste0("Groups: ", paste(groups, collapse = ", ")),
  paste0("Cell types: ", length(celltypes)),
  top_line,
  "Cell type composition was compared across groups using Fisher tests for two groups or chi-square tests for more than two groups.",
  "",
  "限制与说明",
  "Cell type enrichment is composition based and does not account for patient-level repeated measures unless sample-level modeling is added later."
))

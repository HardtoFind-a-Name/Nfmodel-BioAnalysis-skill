args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggvenn)
})

deg_sig_file    <- get_arg("--deg-sig")
gene_sets_sheet <- get_arg("--gene-sets-sheet")
outdir          <- get_arg("--outdir")
logdir          <- get_arg("--logdir")
sn              <- get_arg("--sn")
summary_dir     <- get_arg("--summary-dir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

logsetup <- setup_logging(logdir, summary_dir, sn, "intersect_candidate_genes")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

save_pdf_plot <- function(plot_obj, path, width, height) {
  ok <- tryCatch({
    ggsave(path, plot = plot_obj, width = width, height = height, device = cairo_pdf)
    TRUE
  }, error = function(e) {
    log_message("cairo_pdf unavailable, fallback to default pdf device: ", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    ggsave(path, plot = plot_obj, width = width, height = height, device = "pdf")
  }
}

log_message("Script started: intersect_candidate_genes.R")

if (!file.exists(deg_sig_file)) save_note_and_stop(paste0("DEG sig file does not exist: ", deg_sig_file))
if (!file.exists(gene_sets_sheet)) save_note_and_stop(paste0("Gene sets sheet does not exist: ", gene_sets_sheet))

sheet <- read.csv(gene_sets_sheet, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("set_id", "file_path")
missing_cols <- setdiff(required_cols, colnames(sheet))
if (length(missing_cols) > 0) {
  save_note_and_stop(paste0("Gene sets sheet missing columns: ", paste(missing_cols, collapse = ", ")))
}
if (nrow(sheet) == 0) save_note_and_stop("Gene sets sheet is empty.")

deg_df <- read.csv(deg_sig_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!("gene_symbol" %in% colnames(deg_df))) {
  save_note_and_stop("DEG sig file must contain a 'gene_symbol' column.")
}
deg_genes <- unique(trimws(as.character(deg_df$gene_symbol)))
deg_genes <- deg_genes[deg_genes != "" & !is.na(deg_genes)]
log_message("DEGs loaded: ", length(deg_genes))

venn_list <- list(DEGs = deg_genes)

for (i in seq_len(nrow(sheet))) {
  sid  <- trimws(sheet$set_id[i])
  fpath <- trimws(sheet$file_path[i])
  if (!file.exists(fpath)) save_note_and_stop(paste0("Gene set file not found: ", fpath, " (set_id: ", sid, ")"))
  gs_df <- read.csv(fpath, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("gene_symbol" %in% colnames(gs_df))) {
    colnames(gs_df)[1] <- "gene_symbol"
  }
  genes <- unique(trimws(as.character(gs_df$gene_symbol)))
  genes <- genes[genes != "" & !is.na(genes)]
  venn_list[[sid]] <- genes
  log_message("Gene set '", sid, "' loaded: ", length(genes), " genes")
}

if (length(venn_list) < 2) {
  save_note_and_stop("Need at least 2 sets (DEGs + 1 gene set) for intersection.")
}

candidates <- Reduce(intersect, venn_list)
log_message("Intersection result: ", length(candidates), " candidate genes")

if (length(candidates) == 0) {
  log_message("WARNING: intersection is empty, writing empty output.")
}

write.csv(
  data.frame(gene_symbol = candidates, stringsAsFactors = FALSE),
  file.path(outdir, "01.candidate_genes.csv"),
  row.names = FALSE
)

default_colors <- c("#E64B35", "#4DBBD5", "#ff00ff", "#00A087", "#3C5488")
set_names <- names(venn_list)
n_sets <- length(set_names)

if (n_sets >= 2 && n_sets <= 5) {
  fill_colors <- default_colors[seq_len(n_sets)]
  p_venn <- ggvenn(
    venn_list,
    set_names,
    fill_color = fill_colors,
    show_percentage = TRUE,
    stroke_alpha = 0.5,
    stroke_size = 0.3,
    text_size = 4,
    stroke_color = "white",
    stroke_linetype = "solid",
    set_name_color = fill_colors,
    set_name_size = 5,
    text_color = "black"
  )

  save_pdf_plot(p_venn, file.path(outdir, "02.candidate_venn.pdf"), width = 6, height = 6)
  ggsave(file.path(outdir, "02.candidate_venn.png"), plot = p_venn, width = 5, height = 5, dpi = 600)
  log_message("Venn diagram saved.")
} else {
  log_message("Skipping Venn diagram: ggvenn supports 2-5 sets, got ", n_sets)
}

set_lines <- paste0("  ", names(venn_list), ": ", lengths(venn_list), " genes")
write_summary(c(
  "=== Intersect candidate genes ===",
  paste0("Gene sets: ", length(venn_list)),
  set_lines,
  paste0("Intersection: ", length(candidates), " candidate genes"),
  "",
  "=== 限制与说明 ===",
  "Candidates are the strict intersection of all provided gene sets.",
  paste0("Venn diagram shown for 2-5 sets (got ", length(venn_list), ")."),
  "Gene sets are loaded from paths listed in the gene-sets-sheet CSV."
))

log_message("intersect_candidate_genes step finished.")
message("intersect_candidate_genes step finished.")

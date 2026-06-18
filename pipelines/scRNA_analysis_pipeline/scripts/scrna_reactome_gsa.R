rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

input_rds <- get_arg("--input-rds")
outdir <- get_arg("--outdir")
key_celltype <- get_arg("--key-celltype", required = FALSE, default = "T/NK cells")
max_pathways <- as.integer(get_arg("--max-pathways", required = FALSE, default = "10"))
p_cutoff <- as.numeric(get_arg("--p-cutoff", required = FALSE, default = "0.05"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

celltype_col <- "subset_celltype_manual"
group_col <- "group"
ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_reactome_gsa")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

normalize_group_labels <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    tolower(x) %in% c("tumor", "tumour", "luad", "disease", "case") ~ "Tumor",
    tolower(x) %in% c("normal", "control", "ctrl", "healthy") ~ "Normal",
    TRUE ~ x
  )
}

log_message("Script started: scrna_reactome_gsa")
require_package_or_stop("ReactomeGSA", save_note_and_stop)
require_package_or_stop("ReactomeGSA.data", save_note_and_stop)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
missing_cols <- setdiff(c(celltype_col, group_col), colnames(obj@meta.data))
if (length(missing_cols) > 0) save_note_and_stop(paste0("Missing metadata columns: ", paste(missing_cols, collapse = ", ")))
if (!key_celltype %in% unique(as.character(obj@meta.data[[celltype_col]]))) save_note_and_stop(paste0("Key celltype not found: ", key_celltype))

cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[celltype_col]]) == key_celltype]
sub <- subset(obj, cells = cells)
sub$group <- normalize_group_labels(sub@meta.data[[group_col]])
sub <- subset(sub, cells = rownames(sub@meta.data)[sub$group %in% c("Normal", "Tumor")])
if (ncol(sub) == 0) save_note_and_stop("No cells remained after filtering to Normal/Tumor groups.")
write.csv(sub@meta.data, file.path(outdir, "00_reactomegsa_cell_metadata.csv"), row.names = TRUE)

reactome_res <- tryCatch(
  ReactomeGSA::analyse_sc_clusters(sub, verbose = FALSE),
  error = function(e) save_note_and_stop(paste0("ReactomeGSA failed: ", conditionMessage(e)))
)
pathway_expression <- ReactomeGSA::pathways(reactome_res)
write.csv(data.frame(pathway = rownames(pathway_expression), pathway_expression, check.names = FALSE),
          file.path(outdir, "01_pathway_expression_matrix.csv"), row.names = FALSE)
result_tables <- ReactomeGSA::get_result(reactome_res)
if (!"pathways" %in% names(result_tables)) save_note_and_stop("ReactomeGSA result does not contain a `pathways` table.")
pathway_res <- as.data.frame(result_tables$pathways, check.names = FALSE)
write.csv(pathway_res, file.path(outdir, "02_reactomegsa_pathway_results.csv"), row.names = FALSE)
p_col_candidates <- grep("adj|fdr|p[._ ]?value|pvalue", colnames(pathway_res), ignore.case = TRUE, value = TRUE)
if (length(p_col_candidates) > 0) {
  p_col <- p_col_candidates[1]
  sig <- pathway_res[!is.na(pathway_res[[p_col]]) & pathway_res[[p_col]] < p_cutoff, , drop = FALSE]
} else {
  p_col <- NA_character_
  sig <- pathway_res[0, , drop = FALSE]
}
write.csv(sig, file.path(outdir, "03_reactomegsa_significant_pathways.csv"), row.names = FALSE)

pdf(file.path(outdir, "04_reactomegsa_heatmap_top_pathways.pdf"), width = 10, height = 10)
ReactomeGSA::plot_gsva_heatmap(reactome_res, max_pathways = max_pathways, margins = c(10, 22), cexCol = 2)
dev.off()
png(file.path(outdir, "04_reactomegsa_heatmap_top_pathways.png"), width = 3000, height = 3000, res = 300)
ReactomeGSA::plot_gsva_heatmap(reactome_res, max_pathways = max_pathways, margins = c(10, 22), cexCol = 2)
dev.off()
saveRDS(list(reactome_res = reactome_res, pathway_expression = pathway_expression, pathway_res = pathway_res, sig_pathway_res = sig),
        file.path(outdir, "01.reactomegsa_results.rds"))

log_message("Script completed: scrna_reactome_gsa")
write_summary(c(
  "scRNA ReactomeGSA analysis completed.",
  paste0("Key celltype: ", key_celltype),
  paste0("Cells: ", ncol(sub)),
  paste0("Pathway result rows: ", nrow(pathway_res)),
  paste0("Significant rows: ", nrow(sig)),
  paste0("p column: ", p_col),
  paste0("p cutoff: ", p_cutoff)
))

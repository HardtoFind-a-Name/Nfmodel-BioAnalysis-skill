rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

input_rds <- get_arg("--input-rds")
outdir <- get_arg("--outdir")
hvg_nfeatures <- as.integer(get_arg("--hvg-nfeatures", required = FALSE, default = "3000"))
min_ordering_genes <- as.integer(get_arg("--min-ordering-genes", required = FALSE, default = "50"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

pseudotime_subcluster_col <- "Subcluster"
subset_output_celltype_col <- "subset_celltype_manual"
group_col <- "group"

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_pseudotime_monocle")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_pseudotime_monocle")
require_package_or_stop("monocle", save_note_and_stop)
require_package_or_stop("Biobase", save_note_and_stop)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!subset_output_celltype_col %in% colnames(obj@meta.data)) save_note_and_stop("subset_celltype_manual is missing. Run subset apply first.")

expr <- as.matrix(get_assay_data_compat(obj, assay = DefaultAssay(obj), layer = "counts"))
meta <- obj@meta.data
meta[[pseudotime_subcluster_col]] <- as.character(meta[[subset_output_celltype_col]])
feature <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr), stringsAsFactors = FALSE)
pd <- Biobase::new("AnnotatedDataFrame", data = meta)
fd <- Biobase::new("AnnotatedDataFrame", data = feature)
cds <- monocle::newCellDataSet(expr, phenoData = pd, featureData = fd, expressionFamily = monocle::negbinomial.size())
cds <- monocle::estimateSizeFactors(cds)
cds <- monocle::estimateDispersions(cds)

ordering_genes <- head(VariableFeatures(obj), hvg_nfeatures)
ordering_genes <- ordering_genes[ordering_genes %in% rownames(expr)]
if (length(ordering_genes) < min_ordering_genes) save_note_and_stop(paste0("Too few ordering genes: ", length(ordering_genes)))
cds <- monocle::setOrderingFilter(cds, ordering_genes)
cds <- monocle::reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- monocle::orderCells(cds)

pt <- Biobase::pData(cds)
pt$cell_id <- rownames(pt)
write.csv(pt, file.path(outdir, "01_pseudotime_cell_metadata.csv"), row.names = FALSE)
write.csv(data.frame(gene = ordering_genes), file.path(outdir, "02_pseudotime_ordering_genes.csv"), row.names = FALSE)

pdf(file.path(outdir, "03_pseudotime_by_state.pdf"), width = 7, height = 6)
print(monocle::plot_cell_trajectory(cds, color_by = "State"))
dev.off()
png(file.path(outdir, "03_pseudotime_by_state.png"), width = 2100, height = 1800, res = 300)
print(monocle::plot_cell_trajectory(cds, color_by = "State"))
dev.off()
pdf(file.path(outdir, "04_pseudotime_by_subcluster.pdf"), width = 8, height = 6)
print(monocle::plot_cell_trajectory(cds, color_by = pseudotime_subcluster_col))
dev.off()
png(file.path(outdir, "04_pseudotime_by_subcluster.png"), width = 2400, height = 1800, res = 300)
print(monocle::plot_cell_trajectory(cds, color_by = pseudotime_subcluster_col))
dev.off()
if (group_col %in% colnames(pt)) {
  pdf(file.path(outdir, "05_pseudotime_by_group.pdf"), width = 7, height = 6)
  print(monocle::plot_cell_trajectory(cds, color_by = group_col))
  dev.off()
  png(file.path(outdir, "05_pseudotime_by_group.png"), width = 2100, height = 1800, res = 300)
  print(monocle::plot_cell_trajectory(cds, color_by = group_col))
  dev.off()
}
saveRDS(cds, file.path(outdir, "01.monocle_cds.rds"))

log_message("Script completed: scrna_pseudotime_monocle")
write_summary(c(
  "scRNA Monocle pseudotime completed.",
  paste0("Cells: ", ncol(obj)),
  paste0("Ordering genes: ", length(ordering_genes)),
  paste0("States: ", length(unique(pt$State))),
  paste0("Object file: 01.monocle_cds.rds\n", describe_rds(file.path(outdir, "01.monocle_cds.rds")))
))

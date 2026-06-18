rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

input_rds <- get_arg("--input-rds")
ordering_gene_file <- get_arg("--ordering-gene-file")
ordering_gene_col <- get_arg("--ordering-gene-col", required = FALSE, default = "gene")
outdir <- get_arg("--outdir")
min_ordering_genes <- as.integer(get_arg("--min-ordering-genes", required = FALSE, default = "50"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

pseudotime_subcluster_col <- "Subcluster"
subset_output_celltype_col <- "subset_celltype_manual"

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_pseudotime_ordering_target")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_pseudotime_ordering_target")
require_package_or_stop("monocle", save_note_and_stop)
require_package_or_stop("Biobase", save_note_and_stop)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
if (!subset_output_celltype_col %in% colnames(obj@meta.data)) save_note_and_stop("subset_celltype_manual is missing. Run subset apply first.")
ordering_genes <- tryCatch(read_target_genes(ordering_gene_file, ordering_gene_col), error = function(e) save_note_and_stop(conditionMessage(e)))
ordering_genes <- intersect(ordering_genes, rownames(obj))
if (length(ordering_genes) < min_ordering_genes) save_note_and_stop(paste0("Too few matched ordering genes: ", length(ordering_genes)))

expr <- as.matrix(get_assay_data_compat(obj, assay = DefaultAssay(obj), layer = "counts"))
meta <- obj@meta.data
meta[[pseudotime_subcluster_col]] <- as.character(meta[[subset_output_celltype_col]])
feature <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr), stringsAsFactors = FALSE)
cds <- monocle::newCellDataSet(expr,
  phenoData = Biobase::new("AnnotatedDataFrame", data = meta),
  featureData = Biobase::new("AnnotatedDataFrame", data = feature),
  expressionFamily = monocle::negbinomial.size())
cds <- monocle::estimateSizeFactors(cds)
cds <- monocle::estimateDispersions(cds)
cds <- monocle::setOrderingFilter(cds, ordering_genes)
cds <- monocle::reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- monocle::orderCells(cds)

pt <- Biobase::pData(cds)
pt$cell_id <- rownames(pt)
write.csv(pt, file.path(outdir, "01_target_ordering_pseudotime_cell_metadata.csv"), row.names = FALSE)
write.csv(data.frame(gene = ordering_genes), file.path(outdir, "02_target_ordering_genes_used.csv"), row.names = FALSE)
pdf(file.path(outdir, "03_target_ordering_pseudotime.pdf"), width = 8, height = 6)
print(monocle::plot_cell_trajectory(cds, color_by = pseudotime_subcluster_col))
dev.off()
png(file.path(outdir, "03_target_ordering_pseudotime.png"), width = 2400, height = 1800, res = 300)
print(monocle::plot_cell_trajectory(cds, color_by = pseudotime_subcluster_col))
dev.off()
saveRDS(cds, file.path(outdir, "01.target_ordering_monocle_cds.rds"))

log_message("Script completed: scrna_pseudotime_ordering_target")
write_summary(c(
  "scRNA Monocle target-gene ordering completed.",
  paste0("Matched ordering genes: ", length(ordering_genes)),
  paste0("Cells: ", ncol(obj)),
  paste0("Object file: 01.target_ordering_monocle_cds.rds\n", describe_rds(file.path(outdir, "01.target_ordering_monocle_cds.rds")))
))

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

input_rds <- get_arg("--input-rds")
outdir <- get_arg("--outdir")
gmt_files <- parse_chr_vector(get_arg("--gmt-files", required = FALSE, default = ""))
min_geneset_size <- as.integer(get_arg("--min-geneset-size", required = FALSE, default = "5"))
max_geneset_size <- as.integer(get_arg("--max-geneset-size", required = FALSE, default = "500"))
threshold_t <- as.numeric(get_arg("--threshold-t", required = FALSE, default = "2"))
threshold_p <- as.numeric(get_arg("--threshold-p", required = FALSE, default = "0.05"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

sample_col <- "Sample"
group_col <- "group"
ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_gsva")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

read_gmt_as_list <- function(gmt_file) {
  gmt_obj <- GSEABase::getGmt(gmt_file)
  gene_sets <- lapply(gmt_obj, GSEABase::geneIds)
  names(gene_sets) <- names(gmt_obj)
  gene_sets
}

filter_gene_sets <- function(gene_sets, genes, min_size, max_size) {
  out <- lapply(gene_sets, function(gs) intersect(unique(gs), genes))
  sizes <- lengths(out)
  out[sizes >= min_size & sizes <= max_size]
}

run_gsva_matrix <- function(expr_matrix, gene_sets) {
  if (exists("gsvaParam", where = asNamespace("GSVA"), mode = "function")) {
    GSVA::gsva(GSVA::gsvaParam(exprData = expr_matrix, geneSets = gene_sets, kcdf = "Gaussian"), verbose = FALSE)
  } else {
    GSVA::gsva(expr = expr_matrix, gset.idx.list = gene_sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
  }
}

log_message("Script started: scrna_gsva")
require_package_or_stop("GSVA", save_note_and_stop)
require_package_or_stop("limma", save_note_and_stop)
require_package_or_stop("GSEABase", save_note_and_stop)
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
missing_cols <- setdiff(c(sample_col, group_col), colnames(obj@meta.data))
if (length(missing_cols) > 0) save_note_and_stop(paste0("Missing metadata columns: ", paste(missing_cols, collapse = ", ")))
if (length(gmt_files) == 0) save_note_and_stop("gsva_gmt_files must be provided when run_gsva=true.")
missing_gmt <- gmt_files[!file.exists(gmt_files)]
if (length(missing_gmt) > 0) save_note_and_stop(paste0("GMT files not found: ", paste(missing_gmt, collapse = ", ")))

expr <- as.matrix(get_assay_data_compat(obj, assay = DefaultAssay(obj), layer = "data"))
meta <- obj@meta.data
sample_group <- meta %>% dplyr::select(sample = dplyr::all_of(sample_col), group = dplyr::all_of(group_col)) %>% dplyr::distinct()
if (anyDuplicated(sample_group$sample) > 0) save_note_and_stop("Each Sample must map to one group for GSVA pseudobulk.")
samples <- unique(as.character(meta[[sample_col]]))
pseudobulk <- sapply(samples, function(s) Matrix::rowMeans(expr[, rownames(meta)[as.character(meta[[sample_col]]) == s], drop = FALSE]))
if (is.null(dim(pseudobulk))) pseudobulk <- matrix(pseudobulk, ncol = 1, dimnames = list(rownames(expr), samples))
write.csv(data.frame(gene = rownames(pseudobulk), pseudobulk, check.names = FALSE), file.path(outdir, "01_gsva_pseudobulk_expression.csv"), row.names = FALSE)

all_results <- list()
for (gmt in gmt_files) {
  label <- safe_file_label(extract_label(gmt))
  log_message("Run GSVA collection:", label)
  gene_sets <- filter_gene_sets(read_gmt_as_list(gmt), rownames(pseudobulk), min_geneset_size, max_geneset_size)
  if (length(gene_sets) == 0) save_note_and_stop(paste0("No valid gene sets after filtering: ", gmt))
  gsva_res <- run_gsva_matrix(pseudobulk, gene_sets)
  write.csv(data.frame(pathway = rownames(gsva_res), gsva_res, check.names = FALSE), file.path(outdir, paste0("02_", label, "_gsva_scores.csv")), row.names = FALSE)
  design_df <- sample_group[match(colnames(gsva_res), sample_group$sample), , drop = FALSE]
  design_df$group <- factor(design_df$group)
  if (nlevels(design_df$group) != 2) save_note_and_stop("GSVA currently requires exactly two groups in `group`.")
  design <- model.matrix(~ 0 + group, data = design_df)
  colnames(design) <- levels(design_df$group)
  contrast <- limma::makeContrasts(contrasts = paste0(colnames(design)[2], "-", colnames(design)[1]), levels = design)
  fit <- limma::lmFit(gsva_res, design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast))
  tt <- limma::topTable(fit2, number = Inf, sort.by = "P")
  tt$pathway <- rownames(tt)
  tt$collection <- label
  tt$significant <- abs(tt$t) >= threshold_t & tt$P.Value <= threshold_p
  write.csv(tt, file.path(outdir, paste0("03_", label, "_limma_results.csv")), row.names = FALSE)
  all_results[[label]] <- tt
}
combined <- dplyr::bind_rows(all_results)
write.csv(combined, file.path(outdir, "04_gsva_combined_limma_results.csv"), row.names = FALSE)

sig <- combined %>% dplyr::filter(significant) %>% dplyr::arrange(P.Value)
write.csv(sig, file.path(outdir, "05_gsva_significant_pathways.csv"), row.names = FALSE)
if (nrow(sig) > 0) {
  p <- ggplot(head(sig, 20), aes(x = t, y = reorder(pathway, t), fill = collection)) + geom_col() + theme_bw()
  save_plot_dual(p, outdir, "06_gsva_top_significant_pathways", width = 10, height = 8)
}

log_message("Script completed: scrna_gsva")
write_summary(c(
  "scRNA GSVA analysis completed.",
  paste0("GMT collections: ", paste(basename(gmt_files), collapse = ", ")),
  paste0("Samples: ", ncol(pseudobulk)),
  paste0("Significant pathways: ", nrow(sig)),
  paste0("t threshold: ", threshold_t),
  paste0("p threshold: ", threshold_p)
))

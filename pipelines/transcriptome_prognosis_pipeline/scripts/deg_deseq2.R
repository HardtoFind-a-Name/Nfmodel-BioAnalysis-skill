args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

suppressPackageStartupMessages({
  library(DESeq2)
})

count_file <- get_arg("--count")
group_file <- get_arg("--group")
outdir <- get_arg("--outdir")
min_count <- as.numeric(get_arg("--min-count"))
min_prop <- as.numeric(get_arg("--min-prop"))
vst_blind <- tolower(get_arg("--vst-blind", required = FALSE, default = "true")) %in% c("true", "t", "1", "yes")
logdir <- get_arg("--logdir")
sn          <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

logsetup <- setup_logging(logdir, summary_dir, sn, "deg")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: deg_deseq2.R")

if (!file.exists(count_file)) {
  save_note_and_stop(paste0("Count file does not exist: ", count_file))
}
if (!file.exists(group_file)) {
  save_note_and_stop(paste0("Group file does not exist: ", group_file))
}

log_message("Reading count file: ", count_file)
count_df <- tryCatch(
  read.csv(count_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) save_note_and_stop(paste0("Failed to read count file: ", conditionMessage(e)))
)

log_message("Reading group file: ", group_file)
group_df <- tryCatch(
  read.csv(group_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) save_note_and_stop(paste0("Failed to read group file: ", conditionMessage(e)))
)

if (nrow(count_df) == 0 || ncol(count_df) == 0) {
  save_note_and_stop("Count matrix is empty.")
}
if (nrow(group_df) == 0 || ncol(group_df) == 0) {
  save_note_and_stop("Group table is empty.")
}

log_message("Raw count dimensions: ", nrow(count_df), " genes x ", ncol(count_df), " samples")
log_message("Raw group dimensions: ", nrow(group_df), " rows x ", ncol(group_df), " columns")

colnames(count_df) <- trimws(colnames(count_df))
rownames(group_df) <- trimws(rownames(group_df))
colnames(group_df) <- trimws(colnames(group_df))

if ("Type" %in% colnames(group_df) && !("group" %in% colnames(group_df))) {
  group_df$group <- group_df$Type
}
if (!("group" %in% colnames(group_df))) {
  save_note_and_stop("Group table must contain either 'group' or 'Type' column.")
}

common_samples <- intersect(colnames(count_df), rownames(group_df))
if (length(common_samples) == 0) {
  save_note_and_stop("No overlapping samples between count matrix and group table.")
}

count_df <- count_df[, common_samples, drop = FALSE]
group_df <- group_df[common_samples, , drop = FALSE]
group_df <- group_df[colnames(count_df), , drop = FALSE]

if (any(duplicated(colnames(count_df)))) {
  save_note_and_stop("Duplicated sample IDs found in count matrix.")
}
if (any(duplicated(rownames(group_df)))) {
  save_note_and_stop("Duplicated sample IDs found in group table.")
}
if (any(duplicated(rownames(count_df)))) {
  save_note_and_stop("Duplicated gene IDs found in count matrix.")
}
if (!all(colnames(count_df) == rownames(group_df))) {
  save_note_and_stop("Sample order mismatch between count matrix and group table.")
}

group_df$group <- trimws(as.character(group_df$group))
group_df$group[group_df$group %in% c("normal", "NORMAL")] <- "Normal"
group_df$group[group_df$group %in% c("tumor", "TUMOR", "Tumour", "tumour")] <- "Tumor"

valid_groups <- c("Normal", "Tumor")
invalid_groups <- setdiff(unique(group_df$group), valid_groups)
if (length(invalid_groups) > 0) {
  save_note_and_stop(paste0("Unsupported group labels found: ", paste(invalid_groups, collapse = ", ")))
}

group_df$group <- factor(group_df$group, levels = valid_groups)
if (nlevels(droplevels(group_df$group)) < 2) {
  save_note_and_stop("At least two groups are required for differential expression analysis.")
}

group_table <- table(group_df$group)
if (any(group_table < 2)) {
  save_note_and_stop("Each group must contain at least two samples for DESeq2.")
}

log_message("Group counts: ", paste(names(group_table), as.integer(group_table), collapse = "; "))

coldata <- data.frame(
  row.names = rownames(group_df),
  group = group_df$group
)

raw_count_mat <- as.matrix(count_df)
storage.mode(raw_count_mat) <- "numeric"

if (any(is.na(raw_count_mat))) {
  save_note_and_stop("Count matrix contains NA values.")
}
if (any(!is.finite(raw_count_mat))) {
  save_note_and_stop("Count matrix contains non-finite values.")
}
if (any(raw_count_mat < 0, na.rm = TRUE)) {
  save_note_and_stop("Count matrix contains negative values.")
}

non_integer_ratio <- mean(raw_count_mat %% 1 != 0, na.rm = TRUE)
if (non_integer_ratio > 0) {
  save_note_and_stop(paste0(
    "Input matrix contains non-integer values (ratio = ", signif(non_integer_ratio, 4),
    "). DESeq2 requires raw integer counts."
  ))
}

count_range <- range(raw_count_mat, na.rm = TRUE)
count_mat <- round(raw_count_mat)
storage.mode(count_mat) <- "integer"

log_message("Input count range: [", count_range[1], ", ", count_range[2], "]")

dds <- tryCatch(
  DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ group
  ),
  error = function(e) save_note_and_stop(paste0("Failed to build DESeqDataSet: ", conditionMessage(e)))
)

n_genes_before_filter <- nrow(dds)
min_samples <- max(2, ceiling(min_prop * ncol(dds)))

keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]
n_genes_after_filter <- nrow(dds)

if (n_genes_after_filter == 0) {
  save_note_and_stop("No genes retained after low-expression filtering.")
}

log_message("Genes before filtering: ", n_genes_before_filter)
log_message("Genes after filtering: ", n_genes_after_filter)
log_message("Filtering rule: count >= ", min_count, " in at least ", min_samples, " samples")

log_message("Running DESeq2")
dds <- tryCatch(
  DESeq(dds),
  error = function(e) save_note_and_stop(paste0("DESeq2 failed: ", conditionMessage(e)))
)

res <- tryCatch(
  results(dds, contrast = c("group", "Tumor", "Normal")),
  error = function(e) save_note_and_stop(paste0("Failed to extract DESeq2 results: ", conditionMessage(e)))
)

res_df <- as.data.frame(res)
res_df$gene_symbol <- rownames(res_df)
res_df <- res_df[!is.na(res_df$pvalue), , drop = FALSE]

if (nrow(res_df) == 0) {
  save_note_and_stop("No valid DEG results were returned.")
}

res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange), na.last = TRUE), , drop = FALSE]
res_df$padj_plot <- res_df$padj

write.csv(
  res_df,
  file.path(outdir, "02.deg_all.csv"),
  quote = FALSE,
  row.names = FALSE
)

log_message("Running VST transformation")
vsd <- tryCatch(
  vst(dds, blind = vst_blind),
  error = function(e) save_note_and_stop(paste0("VST failed: ", conditionMessage(e)))
)

vsd_mat <- assay(vsd)
write.csv(
  as.data.frame(vsd_mat),
  file.path(outdir, "01.tcga_luad_vst_expr.csv"),
  quote = FALSE
)

write_summary(c(
  "=== DEG base analysis ===",
  paste0("Matched samples: ", length(common_samples)),
  paste0("Group: ", paste(names(group_table), as.integer(group_table), collapse = "; ")),
  paste0("Count range: [", count_range[1], ", ", count_range[2], "]"),
  paste0("Non-integer ratio: ", signif(non_integer_ratio, 4)),
  paste0("Genes before filtering: ", n_genes_before_filter),
  paste0("Filter rule: count >= ", min_count, " in at least ", min_samples, " samples"),
  paste0("Genes after filtering: ", n_genes_after_filter),
  "",
  "=== 限制与说明 ===",
  "DESeq2 requires raw integer counts; non-integer values are flagged.",
  "Filtering removes low-expression genes to improve statistical power.",
  "VST transformation is performed for downstream visualization.",
  "Threshold-dependent DEG files (03/04/05) are generated by deg_filter_plot.R"
))

log_message("DEG base analysis finished.")
message("DEG base analysis finished.")

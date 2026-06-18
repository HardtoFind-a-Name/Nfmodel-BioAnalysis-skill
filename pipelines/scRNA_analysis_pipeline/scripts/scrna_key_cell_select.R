rm(list = ls())
gc()

suppressPackageStartupMessages({ library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

props_file <- get_arg("--proportions-file")
enrichment_file <- get_arg("--enrichment-file")
keygene_diff_file <- get_arg("--keygene-diff-file")
outdir <- get_arg("--outdir")
manual_key_celltypes <- parse_chr_vector(get_arg("--key-celltypes", required = FALSE, default = ""))
keycell_enrichment_fdr <- as.numeric(get_arg("--keycell-enrichment-fdr", required = FALSE, default = "0.05"))
keycell_min_abs_prop_diff <- as.numeric(get_arg("--keycell-min-abs-prop-diff", required = FALSE, default = "0.05"))
keycell_min_pct_expr <- as.numeric(get_arg("--keycell-min-pct-expr", required = FALSE, default = "10"))
keycell_gene_diff_fdr <- as.numeric(get_arg("--keycell-gene-diff-fdr", required = FALSE, default = "0.05"))
keycell_min_abs_log2fc <- as.numeric(get_arg("--keycell-min-abs-log2fc", required = FALSE, default = "0.25"))
keycell_min_support_genes <- as.integer(get_arg("--keycell-min-support-genes", required = FALSE, default = "1"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_key_cell_select")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

make_key_table <- function(celltypes, source, reason) {
  celltypes <- unique(celltypes[!is.na(celltypes) & nzchar(celltypes)])
  data.frame(
    key_celltype = celltypes,
    safe_label = vapply(celltypes, safe_file_label, character(1)),
    source = source,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

log_message("Script started: scrna_key_cell_select")
if (length(manual_key_celltypes) > 0) {
  out <- make_key_table(manual_key_celltypes, "manual", "manual key_celltypes CLI override")
  write.csv(out, file.path(outdir, "key_celltypes.csv"), row.names = FALSE)
  write_summary(c("Key cell selection used manual override.", paste0("Key cells: ", paste(out$key_celltype, collapse = ", "))))
  quit(save = "no", status = 0)
}

if (!file.exists(props_file)) save_note_and_stop(paste0("Proportions file not found: ", props_file))
if (!file.exists(enrichment_file)) save_note_and_stop(paste0("Enrichment file not found: ", enrichment_file))
if (!file.exists(keygene_diff_file)) save_note_and_stop(paste0("Keygene diff file not found: ", keygene_diff_file))
props <- read.csv(props_file, stringsAsFactors = FALSE, check.names = FALSE)
enrich <- read.csv(enrichment_file, stringsAsFactors = FALSE, check.names = FALSE)
keydiff <- read.csv(keygene_diff_file, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(keydiff) == 0) {
  out <- make_key_table(character(), "auto", "no keygene diff rows")
  write.csv(out, file.path(outdir, "key_celltypes.csv"), row.names = FALSE)
  write_summary("No key cells selected because keygene analysis was empty.")
  quit(save = "no", status = 0)
}

prop_summary <- props %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(prop_min = min(proportion, na.rm = TRUE), prop_max = max(proportion, na.rm = TRUE),
                   proportion_diff = prop_max - prop_min, .groups = "drop")
if (all(c("normal", "tumor") %in% unique(props$group))) {
  wide <- reshape(props[, c("group", "celltype", "proportion")], idvar = "celltype", timevar = "group", direction = "wide")
  names(wide) <- sub("^proportion\\.", "", names(wide))
  wide$proportion_diff <- wide$tumor - wide$normal
  prop_summary <- wide[, c("celltype", "proportion_diff"), drop = FALSE]
}

enrich_keep <- enrich %>%
  dplyr::filter(!is.na(p_adjust), p_adjust <= keycell_enrichment_fdr) %>%
  dplyr::select(celltype, enrichment_p_adjust = p_adjust)
key_keep <- keydiff %>%
  dplyr::filter(!is.na(p_adjust), p_adjust <= keycell_gene_diff_fdr,
                !is.na(pct_expr), pct_expr >= keycell_min_pct_expr,
                !is.na(log2FC), abs(log2FC) >= keycell_min_abs_log2fc) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(support_gene_count = dplyr::n_distinct(gene), support_genes = paste(unique(gene), collapse = ","), .groups = "drop") %>%
  dplyr::filter(support_gene_count >= keycell_min_support_genes)

candidate <- prop_summary %>%
  dplyr::mutate(abs_prop_diff = abs(proportion_diff)) %>%
  dplyr::filter(abs_prop_diff >= keycell_min_abs_prop_diff) %>%
  dplyr::inner_join(enrich_keep, by = "celltype") %>%
  dplyr::inner_join(key_keep, by = "celltype") %>%
  dplyr::arrange(enrichment_p_adjust, dplyr::desc(abs_prop_diff), dplyr::desc(support_gene_count))
write.csv(candidate, file.path(outdir, "01_key_cell_selection_candidates.csv"), row.names = FALSE)
reason <- if (nrow(candidate) > 0) "passed enrichment, proportion, keygene expression, and keygene group-diff thresholds" else "no celltype passed all automatic thresholds"
out <- make_key_table(candidate$celltype, "auto", reason)
write.csv(out, file.path(outdir, "key_celltypes.csv"), row.names = FALSE)

log_message("Script completed: scrna_key_cell_select")
write_summary(c(
  "scRNA key cell selection completed.",
  paste0("Key cells selected: ", nrow(out)),
  paste0("Key cells: ", paste(out$key_celltype, collapse = ", ")),
  paste0("Thresholds: enrichment_fdr=", keycell_enrichment_fdr,
         ", abs_prop_diff=", keycell_min_abs_prop_diff,
         ", pct_expr=", keycell_min_pct_expr,
         ", gene_diff_fdr=", keycell_gene_diff_fdr,
         ", abs_log2FC=", keycell_min_abs_log2fc)
))

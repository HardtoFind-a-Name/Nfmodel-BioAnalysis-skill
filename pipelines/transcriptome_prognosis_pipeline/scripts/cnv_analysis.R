suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(tibble)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

risk_file    <- get_arg("--risk-file", required = FALSE, default = NULL)  # for interface compat
gene_file    <- get_arg("--gene-file")
expr_file    <- get_arg("--expr")
train_id     <- get_arg("--train-id")
outdir       <- get_arg("--outdir")
logdir       <- get_arg("--logdir")
sn           <- get_arg("--sn")
summary_dir  <- get_arg("--summary-dir")
gistic_file  <- get_arg("--gistic-matrix")
download_dir <- get_arg("--download-dir", required = FALSE, default = NULL)
immune_file  <- get_arg("--immune-input", required = FALSE, default = NULL)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

# --- Logging ---
logsetup <- setup_logging(logdir, summary_dir, sn, "cnv")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop
log_message("Start CNV analysis")

# --- Helpers ---
strip_gene_suffix <- function(x) sub("\\|.*$", "", as.character(x))
standardize_tcga <- function(x) { x <- gsub("\\.", "-", as.character(x)); substr(x, 1, 16) }
sample_type_code <- function(x) substr(standardize_tcga(x), 14, 15)
patient_barcode <- function(x) substr(standardize_tcga(x), 1, 12)
looks_like_tcga <- function(x) grepl("^TCGA[-.][A-Za-z0-9]{2}[-.][A-Za-z0-9]{4}", as.character(x))

read_table_auto <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  pl <- tolower(path)
  if (grepl("\\.(tsv|txt)(\\.gz)?$", pl)) return(readr::read_tsv(path, show_col_types = FALSE))
  if (grepl("\\.csv(\\.gz)?$", pl)) return(readr::read_csv(path, show_col_types = FALSE))
  if (grepl("\\.gz$", pl)) return(readr::read_tsv(path, show_col_types = FALSE))
  stop("Unsupported format: ", path)
}

deduplicate_rows <- function(mat) {
  rn <- strip_gene_suffix(rownames(mat))
  keep <- !duplicated(rn) & !is.na(rn) & nzchar(rn)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- rn[keep]; mat
}

table_to_matrix <- function(df, genes, label = "matrix") {
  df <- as.data.frame(df, check.names = FALSE)
  if (ncol(df) < 2) stop(label, ": input too small")
  fc <- as.character(df[[1]]); oc <- colnames(df)[-1]
  gh <- sum(strip_gene_suffix(fc) %in% genes, na.rm = TRUE)
  th <- mean(looks_like_tcga(oc), na.rm = TRUE)
  if (gh >= 1 || th > 0.5) {
    mat <- data.matrix(df[, -1, drop = FALSE])
    rownames(mat) <- strip_gene_suffix(fc); colnames(mat) <- oc
  } else {
    mat <- t(data.matrix(df[, -1, drop = FALSE]))
    rownames(mat) <- strip_gene_suffix(oc); colnames(mat) <- fc
  }
  storage.mode(mat) <- "numeric"; mat
}

match_barcodes <- function(emat, gmat) {
  eid <- colnames(emat); gid <- colnames(gmat)
  cs <- intersect(eid, gid)
  if (length(cs) > 0) return(list(level = "sample", ids = cs,
    emat = emat[, cs, drop = FALSE], gmat = gmat[, cs, drop = FALSE]))
  ep <- patient_barcode(eid); gp <- patient_barcode(gid)
  cp <- intersect(ep, gp)
  if (length(cp) == 0) return(list(level = "none", ids = character(0)))
  emat2 <- emat[, match(cp, ep), drop = FALSE]
  gmat2 <- gmat[, match(cp, gp), drop = FALSE]
  colnames(emat2) <- cp
  colnames(gmat2) <- cp
  list(level = "patient", ids = cp, emat = emat2, gmat = gmat2)
}

save_both <- function(plot_obj, fb, w = 8, h = 5, dpi = 300) {
  ggsave(paste0(fb, ".pdf"), plot_obj, width = w, height = h, units = "in")
  ggsave(paste0(fb, ".png"), plot_obj, width = w, height = h, units = "in", dpi = dpi)
}

cnv_labels <- c("-2" = "Deep Del", "-1" = "Shallow Del", "0" = "Diploid", "1" = "Gain", "2" = "Amplification")
cnv_colors <- c("-2" = "#2166AC", "-1" = "#67A9CF", "0" = "#BDBDBD", "1" = "#F4A582", "2" = "#B2182B")
cnv_levels <- unname(cnv_labels[c("-2", "-1", "0", "1", "2")])

run_kw <- function(df, vcol) {
  df <- df[!is.na(df[[vcol]]) & !is.na(df$cnv_state), , drop = FALSE]
  if (nrow(df) == 0 || dplyr::n_distinct(df$cnv_state) < 2) return(NA_real_)
  kruskal.test(reformulate("cnv_state", response = vcol), data = df)$p.value
}

# --- Extract prognostic genes from gene selection file ---
log_message("Extracting prognostic genes from: ", gene_file)
if (!file.exists(gene_file)) stop(paste("Gene file not found:", gene_file))
gene_df <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene" %in% colnames(gene_df)) stop("Gene file must have a 'gene' column")
prog_genes <- unique(trimws(as.character(gene_df$gene)))
prog_genes <- prog_genes[!is.na(prog_genes) & nzchar(prog_genes)]
if (length(prog_genes) == 0) stop("No prognostic genes found in gene file")
log_message("Prognostic genes: ", paste(prog_genes, collapse = ", "))

# --- Validate inputs ---
if (!file.exists(expr_file)) stop("Expression file not found: ", expr_file)

# =========================================================
# 1. Expression matrix
# =========================================================
log_message("Step 1: expression matrix")
expr_input <- read_table_auto(expr_file)
exp_mat <- table_to_matrix(expr_input, prog_genes, "expression")
colnames(exp_mat) <- standardize_tcga(colnames(exp_mat))
exp_mat <- exp_mat[, !duplicated(colnames(exp_mat)), drop = FALSE]
exp_mat <- exp_mat[intersect(prog_genes, rownames(exp_mat)), , drop = FALSE]
if (nrow(exp_mat) == 0) stop("No prognostic genes in expression matrix")

# =========================================================
# 2. GISTIC matrix (local file → Xena cache → Xena download)
# =========================================================
log_message("Step 2: GISTIC matrix")

extract_first_xena_table <- function(x, label = "Xena") {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    idx <- which(vapply(x, is.data.frame, logical(1)))
    if (length(idx) == 0) stop(label, ": no data.frame in XenaPrepare result")
    return(x[[idx[1]]])
  }
  stop(label, ": unsupported XenaPrepare result type")
}

train_id_dots <- gsub("-", ".", train_id)

if (!is.null(gistic_file) && nzchar(gistic_file) && file.exists(gistic_file)) {
  log_message("Using local GISTIC file: ", gistic_file)
  gistic_raw <- read_table_auto(gistic_file)
} else {
  xena_cache <- file.path(download_dir, paste0(train_id_dots, ".sampleMap"),
    "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz")
  if (file.exists(xena_cache)) {
    log_message("Using Xena-cached GISTIC file: ", xena_cache)
    gistic_raw <- read_table_auto(xena_cache)
  } else {
    log_message("Downloading GISTIC2 from UCSC Xena for ", train_id)
    library(UCSCXenaTools)
    library(SummarizedExperiment)
    dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
    xena_query <- XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
      XenaFilter(filterCohorts = gsub("^TCGA-", "", train_id)) %>%
      XenaFilter(filterDatasets = "Gistic2_CopyNumber_Gistic2_all_thresholded.*genes") %>%
      XenaQuery()
    if (nrow(xena_query) == 0) stop("No GISTIC2 dataset found on UCSC Xena for ", train_id)
    xena_dl <- XenaDownload(xena_query, destdir = download_dir, force = FALSE)
    xena_prep <- XenaPrepare(xena_dl)
    gistic_raw <- extract_first_xena_table(xena_prep, "GISTIC")
    log_message("GISTIC2 downloaded and prepared")
  }
}

gistic_mat <- table_to_matrix(gistic_raw, prog_genes, "GISTIC")
gistic_mat <- deduplicate_rows(gistic_mat)
gid <- standardize_tcga(colnames(gistic_mat))
if (any(nchar(gid) >= 15 & grepl("^[0-9]{2}$", sample_type_code(gid)))) {
  colnames(gistic_mat) <- gid
} else {
  colnames(gistic_mat) <- patient_barcode(gid)
}
gistic_mat <- gistic_mat[, !duplicated(colnames(gistic_mat)), drop = FALSE]
gistic_mat <- gistic_mat[intersect(prog_genes, rownames(gistic_mat)), , drop = FALSE]
gistic_mat <- round(gistic_mat)
if (nrow(gistic_mat) == 0) stop("No prognostic genes in GISTIC matrix")

# =========================================================
# 3. CNV Heatmap
# =========================================================
log_message("Step 3: CNV heatmap")
cnv_long <- as.data.frame(as.table(gistic_mat), stringsAsFactors = FALSE) %>%
  setNames(c("gene", "sample", "cnv_value")) %>%
  mutate(cnv_value = as.integer(cnv_value),
         cnv_state = factor(unname(cnv_labels[as.character(cnv_value)]), levels = cnv_levels))
ordered <- cnv_long %>% group_by(sample) %>%
  summarise(burden = sum(abs(cnv_value), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(burden)) %>% pull(sample)

hm_mat <- gistic_mat[, ordered, drop = FALSE]
ht <- Heatmap(hm_mat, name = "CNV",
  col = colorRamp2(c(-2, -1, 0, 1, 2), unname(cnv_colors[c("-2", "-1", "0", "1", "2")])),
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  show_row_names = TRUE, row_names_side = "left", rect_gp = gpar(col = "white", lwd = 0.3),
  row_title = "Genes", column_title = train_id)
pdf(file.path(outdir, "01.CNV_heatmap.pdf"), width = 14, height = 4 + 0.35 * nrow(hm_mat)); draw(ht); dev.off()
png(file.path(outdir, "01.CNV_heatmap.png"), width = 14, height = 4 + 0.35 * nrow(hm_mat),
    units = "in", res = 300); draw(ht); dev.off()

# =========================================================
# 4. Expression by CNV
# =========================================================
log_message("Step 4: expression by CNV")
mr <- match_barcodes(exp_mat, gistic_mat)
if (length(mr$ids) == 0) {
  log_message("WARNING: No matched samples — skipping expression-CNV plots")
} else {
  cg <- intersect(rownames(mr$emat), rownames(mr$gmat))
  if (length(cg) == 0) {
    log_message("WARNING: No common genes between expression and GISTIC — skipping")
  } else {
    e2 <- mr$emat[cg, , drop = FALSE]; g2 <- mr$gmat[cg, , drop = FALSE]
    cs <- colnames(e2)

    plot_genes <- intersect(rownames(e2), rownames(g2))
    if (length(plot_genes) == 0) {
      log_message("WARNING: No common genes in matched matrices — skipping")
    } else {
    expr_long <- purrr::map_dfr(plot_genes, function(g) {
    tibble(gene = g, sample = cs,
      expression = as.numeric(e2[g, cs]),
      cnv_value = as.integer(g2[g, cs]),
      cnv_state = factor(unname(cnv_labels[as.character(as.integer(g2[g, cs]))]), levels = cnv_levels))
  })
  expr_kw <- expr_long %>% group_by(gene) %>%
    group_modify(~ tibble(kruskal_p = run_kw(.x, "expression"))) %>%
    ungroup() %>% mutate(kruskal_fdr = p.adjust(kruskal_p, method = "BH"))

  for (g in unique(expr_long$gene)) {
    pdf <- expr_long %>% filter(gene == g)
    p_kw <- expr_kw %>% filter(gene == g) %>% pull(kruskal_p)
    p_fdr <- expr_kw %>% filter(gene == g) %>% pull(kruskal_fdr)
    p <- ggboxplot(pdf, x = "cnv_state", y = "expression", fill = "cnv_state", color = "cnv_state",
      add = "jitter", add.params = list(size = 0.6, alpha = 0.4)) +
      scale_fill_manual(values = unname(cnv_colors[c("-2", "-1", "0", "1", "2")])) +
      scale_color_manual(values = unname(cnv_colors[c("-2", "-1", "0", "1", "2")])) +
      labs(title = paste0(g, " expression by CNV state"),
           subtitle = paste0("Kruskal-Wallis P = ", signif(p_kw, 3), "; FDR = ", signif(p_fdr, 3)),
           x = "CNV state", y = "Expression") +
      theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
    save_both(p, file.path(outdir, paste0("02.expression_by_CNV_", g)), w = 7.5, h = 5)
    }
  }
  }
}

# =========================================================
# 5. Immune by CNV
# =========================================================
log_message("Step 5: immune by CNV")
`%||%` <- function(x, y) if (is.null(x)) y else x
if (!is.null(immune_file) && file.exists(immune_file) && exists("mr") && length(mr$ids) > 0) {
  idf <- read_table_auto(immune_file)
  sc <- if ("sample" %in% colnames(idf)) "sample" else colnames(idf)[1]
  ift <- setdiff(colnames(idf), sc)

  if (length(ift) > 0) {
    idf2 <- idf %>% rename(sample = all_of(sc)) %>%
      mutate(sample = if (mr$level == "patient") patient_barcode(sample) else standardize_tcga(sample)) %>%
      distinct(sample, .keep_all = TRUE)

    gm_display <- if (exists("g2")) g2 else gistic_mat
    cg <- intersect(rownames(exp_mat), rownames(gm_display))
    gm <- gm_display[cg, , drop = FALSE]
    cs <- colnames(gm)

    ilong <- purrr::map_dfr(rownames(gm), function(g) {
      tibble(gene = g, sample = cs,
        cnv_value = as.integer(gm[g, cs]),
        cnv_state = factor(unname(cnv_labels[as.character(as.integer(gm[g, cs]))]), levels = cnv_levels))
    }) %>% inner_join(idf2, by = "sample") %>%
      pivot_longer(cols = all_of(ift), names_to = "immune_feature", values_to = "value") %>%
      mutate(value = suppressWarnings(as.numeric(value))) %>% filter(!is.na(cnv_state), !is.na(value))

    if (nrow(ilong) > 0) {
      p_imm <- ggplot(ilong, aes(x = immune_feature, y = value, fill = cnv_state, color = cnv_state)) +
        geom_boxplot(position = position_dodge(width = 0.78), width = 0.58,
                     outlier.shape = NA, alpha = 0.75, linewidth = 0.35) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.16, dodge.width = 0.78),
                    size = 0.8, alpha = 0.5) +
        facet_wrap(~ gene, scales = "free_y", ncol = 2) +
        scale_fill_manual(values = unname(cnv_colors[c("-2", "-1", "0", "1", "2")])) +
        scale_color_manual(values = unname(cnv_colors[c("-2", "-1", "0", "1", "2")])) +
        labs(title = "Immune infiltration across CNV states", x = NULL, y = "Infiltration",
             fill = "CNV", color = "CNV") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
              strip.background = element_rect(fill = "grey90", color = "grey40"),
              strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
              legend.position = "right")
      save_both(p_imm, file.path(outdir, "03.immune_by_CNV_combined"), w = 14, h = 8)
    }
  }
}

total_genes_cnv <- nrow(gistic_mat)
total_samples_cnv <- ncol(gistic_mat)
amp_count <- sum(gistic_mat == 2, na.rm = TRUE)
del_count <- sum(gistic_mat == -2, na.rm = TRUE)
gain_count <- sum(gistic_mat == 1, na.rm = TRUE)
shallow_del_count <- sum(gistic_mat == -1, na.rm = TRUE)
total_calls <- total_genes_cnv * total_samples_cnv
immune_cnv_status <- if (!is.null(immune_file) && file.exists(immune_file) && exists("ilong") && nrow(ilong) > 0) {
  paste0("Immune-CNV analysis performed with ", length(unique(ilong$immune_feature)), " immune features")
} else "Immune-CNV analysis not performed"
write_summary(c("CNV analysis summary",
  paste0("Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("SN: ", sn),
  paste0("Expression: ", expr_file),
  paste0("GISTIC: ", gistic_file),
  paste0("Immune input: ", if (is.null(immune_file)) "(none)" else immune_file),
  paste0("Prognostic genes analyzed: ", length(prog_genes)),
  paste0("Genes in GISTIC matrix: ", total_genes_cnv),
  paste0("Samples in GISTIC matrix: ", total_samples_cnv),
  paste0("Amplification (2) frequency: ", amp_count, "/", total_calls, " (", sprintf("%.2f", amp_count/total_calls*100), "%)"),
  paste0("Gain (1) frequency: ", gain_count, "/", total_calls, " (", sprintf("%.2f", gain_count/total_calls*100), "%)"),
  paste0("Diploid (0): ", total_calls - amp_count - gain_count - del_count - shallow_del_count, "/", total_calls),
  paste0("Shallow deletion (-1) frequency: ", shallow_del_count, "/", total_calls, " (", sprintf("%.2f", shallow_del_count/total_calls*100), "%)"),
  paste0("Deep deletion (-2) frequency: ", del_count, "/", total_calls, " (", sprintf("%.2f", del_count/total_calls*100), "%)"),
  immune_cnv_status,
  "",
  "限制与说明：",
  "1. CNV分析基于GISTIC2阈值化数据，-2/-1/0/1/2分别表示Deep Del/Shallow Del/Diploid/Gain/Amplification。",
  "2. 表达-CNV关联使用Kruskal-Wallis检验，结果受样本量和CNV事件数影响。",
  "3. 免疫-CNV分析依赖于免疫浸润输入数据，缺失时跳过。",
  "4. 本结果为生物信息学预测，需进一步实验验证。"
))

log_message("CNV analysis complete")

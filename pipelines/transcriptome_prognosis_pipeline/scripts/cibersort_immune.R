suppressPackageStartupMessages({
  library(ggpubr)
  library(IOBR)
  library(RColorBrewer)
  library(reshape2)
  library(psych)
  library(ggplot2)
  library(corrplot)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(estimate)
})

options(stringsAsFactors = FALSE)

# --- CLI args ---
args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

risk_file   <- get_arg("--risk-file")
gene_file   <- get_arg("--gene-file")
expr_file   <- get_arg("--expr")
train_id    <- get_arg("--train-id")
outdir      <- get_arg("--outdir")
logdir      <- get_arg("--logdir")
sn          <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
perm_n <- 1000

logsetup <- setup_logging(logdir, summary_dir, sn, "cibersort")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop
log_message("Start CIBERSORT immune analysis")

# --- Input validation ---
if (!file.exists(expr_file)) stop(paste("Expression matrix not found:", expr_file))
if (!file.exists(risk_file)) stop(paste("Risk file not found:", risk_file))

# --- Read expression matrix ---
log_message("Reading expression matrix")
dat <- read.csv(expr_file, row.names = 1, check.names = FALSE)
dat <- as.data.frame(dat)
colnames(dat) <- gsub("\\.", "-", colnames(dat))

# --- Read risk file ---
log_message("Reading risk file")
risk_df <- read.csv(risk_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(risk_df)) {
  colnames(risk_df)[1] <- "sample"
}
risk_df$sample <- as.character(risk_df$sample)
risk_df$sample <- gsub("\\.", "-", risk_df$sample)
if (!"riskGroup" %in% colnames(risk_df)) stop("Risk file must have 'riskGroup' column")
risk_df$riskGroup <- as.character(risk_df$riskGroup)
risk_df <- risk_df[!is.na(risk_df$sample) & nzchar(risk_df$sample) &
                   !is.na(risk_df$riskGroup) & nzchar(risk_df$riskGroup), ]
risk_df <- risk_df[risk_df$sample %in% colnames(dat), ]
if (nrow(risk_df) == 0) stop("No overlapping samples between expression and risk data")
risk_df <- risk_df[match(intersect(colnames(dat), risk_df$sample), risk_df$sample), ]
dat <- dat[, risk_df$sample, drop = FALSE]
rownames(risk_df) <- risk_df$sample

# --- Extract key genes from gene selection file ---
log_message("Extracting key genes from: ", gene_file)
if (!file.exists(gene_file)) stop(paste("Gene file not found:", gene_file))
gene_df <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene" %in% colnames(gene_df)) stop("Gene file must have a 'gene' column")
hubgene <- unique(trimws(as.character(gene_df$gene)))
hubgene <- hubgene[!is.na(hubgene) & nzchar(hubgene)]
hubgene <- hubgene[hubgene %in% rownames(dat)]
log_message("Key genes: ", if (length(hubgene) > 0) paste(hubgene, collapse = ", ") else "(none)")

# =========================================================
# ESTIMATE
# =========================================================
estimate_input_file <- file.path(outdir, "00_estimate_input.gct")
estimate_score_file <- file.path(outdir, "00_estimate_score.gct")
estimate_score_csv  <- file.path(outdir, "00_estimate_score.csv")
cnv_immune_input_file <- file.path(outdir, "11_cnv_immune_input.csv")

log_message("Running ESTIMATE")
estimate_expr <- dat
estimate_expr$GeneSymbol <- rownames(estimate_expr)
estimate_expr <- estimate_expr[, c("GeneSymbol", colnames(dat)), drop = FALSE]
estimate_expr <- estimate_expr[!duplicated(estimate_expr$GeneSymbol) & !is.na(estimate_expr$GeneSymbol), ]
estimate_write_df <- estimate_expr
colnames(estimate_write_df)[1] <- "NAME"
estimate_write_df$Description <- estimate_write_df$NAME
estimate_write_df <- estimate_write_df[, c("NAME", "Description", colnames(dat)), drop = FALSE]
write.table(estimate_write_df, file = estimate_input_file, sep = "\t", quote = FALSE, row.names = FALSE)

estimate::filterCommonGenes(input.f = estimate_input_file, output.f = estimate_input_file, id = "GeneSymbol")
estimate::estimateScore(input.ds = estimate_input_file, output.ds = estimate_score_file, platform = "illumina")

estimate_score_raw <- read.table(estimate_score_file, skip = 2, header = TRUE, sep = "\t",
                                  check.names = FALSE, stringsAsFactors = FALSE)
estimate_score_df <- estimate_score_raw
rownames(estimate_score_df) <- estimate_score_df[, 1]
estimate_score_df <- estimate_score_df[, -c(1, 2), drop = FALSE]
estimate_score_df <- as.data.frame(t(estimate_score_df), check.names = FALSE)
estimate_score_df$sample <- rownames(estimate_score_df)

score_cols_needed <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
score_cols_exist <- intersect(score_cols_needed, colnames(estimate_score_df))
estimate_score_df$sample <- gsub("\\.", "-", estimate_score_df$sample)
estimate_score_df <- merge(estimate_score_df[, c("sample", score_cols_exist), drop = FALSE],
                           risk_df[, c("sample", "riskGroup"), drop = FALSE], by = "sample")
estimate_score_df$riskGroup <- factor(estimate_score_df$riskGroup, levels = c("Low risk", "High risk"))
write.csv(estimate_score_df, estimate_score_csv, row.names = FALSE)

plot_estimate_box <- function(data, score_col, ylab_text, title_text, out_prefix) {
  ymax <- max(data[[score_col]], na.rm = TRUE); ymin <- min(data[[score_col]], na.rm = TRUE)
  ypad <- 0.20 * (ymax - ymin)
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.5
  p <- ggplot(data, aes(x = riskGroup, y = .data[[score_col]], fill = riskGroup)) +
    geom_violin(trim = TRUE, color = "black", alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif") +
    scale_y_continuous(limits = c(ymin, ymax + ypad), expand = expansion(mult = c(0.05, 0.12))) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("Low risk" = "#4DBBD5", "High risk" = "#E64B35")) +
    theme_bw() + labs(title = title_text, x = "", y = ylab_text, fill = "") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.x = element_text(face = "bold", colour = "black", size = 11),
          legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = margin(12, 18, 10, 16))
  ggsave(file.path(outdir, paste0(out_prefix, ".pdf")), p, width = 6, height = 5.5)
  ggsave(file.path(outdir, paste0(out_prefix, ".png")), p, width = 6, height = 5.5, dpi = 300)
}

if ("ImmuneScore"  %in% score_cols_exist)
  plot_estimate_box(estimate_score_df, "ImmuneScore",  "Immune Score",  "Immune Score by risk group",  "01_immune_score")
if ("StromalScore" %in% score_cols_exist)
  plot_estimate_box(estimate_score_df, "StromalScore", "Stromal Score", "Stromal Score by risk group", "02_stromal_score")
if ("ESTIMATEScore" %in% score_cols_exist)
  plot_estimate_box(estimate_score_df, "ESTIMATEScore", "ESTIMATE Score", "ESTIMATE Score by risk group", "03_estimate_score")
if ("TumorPurity"   %in% score_cols_exist)
  plot_estimate_box(estimate_score_df, "TumorPurity",   "Tumor Purity",   "Tumor Purity by risk group",   "04_tumor_purity")

# =========================================================
# CIBERSORT
# =========================================================
rdata_file <- file.path(outdir, "00_cibersort.Rdata")
if (!file.exists(rdata_file)) {
  log_message("Running CIBERSORT (perm=", perm_n, ") ...")
  cibersort <- deconvo_tme(eset = dat, method = "cibersort", arrays = FALSE, perm = perm_n)
  save(cibersort, dat, risk_df, file = rdata_file)
  log_message("CIBERSORT done")
} else {
  log_message("Loading cached CIBERSORT")
  load(rdata_file)
}
if (!exists("cibersort")) stop("cibersort object not found")
write.csv(cibersort, file.path(outdir, "00_cibersort_full_result.csv"))

# --- P-value distribution ---
pvalue_df <- data.frame(sample = rownames(cibersort), pvalue = cibersort$`P-value_CIBERSORT`,
  pass = ifelse(cibersort$`P-value_CIBERSORT` < 0.05, "P < 0.05", "P >= 0.05"), stringsAsFactors = FALSE)
write.csv(pvalue_df, file.path(outdir, "01_cibersort_pvalue_distribution.csv"), row.names = FALSE)

p_hist <- ggplot(pvalue_df, aes(x = pvalue, fill = pass)) +
  geom_histogram(binwidth = 0.05, color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("P < 0.05" = "#E64B35", "P >= 0.05" = "#4DBBD5")) +
  theme_bw() + labs(title = "CIBERSORT P-value distribution", x = "P-value", y = "Count", fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(outdir, "01_cibersort_pvalue_distribution.pdf"), p_hist, width = 7, height = 5)
ggsave(file.path(outdir, "01_cibersort_pvalue_distribution.png"), p_hist, width = 7, height = 5, dpi = 300)

# --- Sample filter ---
p_bar_plot <- ggplot(pvalue_df, aes(x = pass, fill = pass)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = c("P < 0.05" = "#E64B35", "P >= 0.05" = "#4DBBD5")) +
  theme_bw() + labs(title = "CIBERSORT sample filtering", x = "", y = "Count", fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(outdir, "02_cibersort_sample_filtering.pdf"), p_bar_plot, width = 6, height = 5)
ggsave(file.path(outdir, "02_cibersort_sample_filtering.png"), p_bar_plot, width = 6, height = 5, dpi = 300)

cibersort2 <- cibersort[cibersort$`P-value_CIBERSORT` < 0.05, , drop = FALSE]
if (nrow(cibersort2) == 0) stop("No samples passed CIBERSORT P < 0.05")
cibersort.res <- cibersort2 %>% as.data.frame()
colnames(cibersort.res) <- gsub("_CIBERSORT", "", colnames(cibersort.res))
colnames(cibersort.res) <- gsub("_", " ", colnames(cibersort.res))
rownames(cibersort.res) <- cibersort.res$ID
drop_cols <- grep("(^ID$|P-value|Correlation|RMSE)", colnames(cibersort.res), value = TRUE)
cibersort.res <- cibersort.res[, setdiff(colnames(cibersort.res), drop_cols), drop = FALSE]
cibersort.res <- cibersort.res[, apply(cibersort.res, 2, function(col) !all(col == 0)), drop = FALSE]
write.csv(cibersort.res, file.path(outdir, "03_CIBERSORT_result_filtered.csv"))
risk_df2 <- risk_df[rownames(cibersort.res), , drop = FALSE]

# --- Rename immune cells ---
cell_name_map <- c(
  "B cells naive" = "Naive B cells",
  "B cells memory" = "Memory B cells",
  "Plasma cells" = "Plasma cells",
  "T cells CD8" = "CD8+ T cells",
  "T cells CD4 naive" = "Naive CD4+ T cells",
  "T cells CD4 memory resting" = "Resting memory CD4+ T cells",
  "T cells CD4 memory activated" = "Activated memory CD4+ T cells",
  "T cells follicular helper" = "Follicular helper T cells",
  "T cells regulatory (Tregs)" = "Regulatory T cells (Tregs)",
  "T cells gamma delta" = "Gamma delta T cells",
  "NK cells resting" = "Resting NK cells",
  "NK cells activated" = "Activated NK cells",
  "Monocytes" = "Monocytes",
  "Macrophages M0" = "M0 macrophages",
  "Macrophages M1" = "M1 macrophages",
  "Macrophages M2" = "M2 macrophages",
  "Dendritic cells resting" = "Resting dendritic cells",
  "Dendritic cells activated" = "Activated dendritic cells",
  "Mast cells resting" = "Resting mast cells",
  "Mast cells activated" = "Activated mast cells",
  "Eosinophils" = "Eosinophils",
  "Neutrophils" = "Neutrophils"
)
cell_name_map_dot <- c(
  "B.cells.naive" = "Naive B cells",
  "B.cells.memory" = "Memory B cells",
  "Plasma.cells" = "Plasma cells",
  "T.cells.CD8" = "CD8+ T cells",
  "T.cells.CD4.naive" = "Naive CD4+ T cells",
  "T.cells.CD4.memory.resting" = "Resting memory CD4+ T cells",
  "T.cells.CD4.memory.activated" = "Activated memory CD4+ T cells",
  "T.cells.follicular.helper" = "Follicular helper T cells",
  "T.cells.regulatory..Tregs." = "Regulatory T cells (Tregs)",
  "T.cells.gamma.delta" = "Gamma delta T cells",
  "NK.cells.resting" = "Resting NK cells",
  "NK.cells.activated" = "Activated NK cells",
  "Monocytes" = "Monocytes",
  "Macrophages.M0" = "M0 macrophages",
  "Macrophages.M1" = "M1 macrophages",
  "Macrophages.M2" = "M2 macrophages",
  "Dendritic.cells.resting" = "Resting dendritic cells",
  "Dendritic.cells.activated" = "Activated dendritic cells",
  "Mast.cells.resting" = "Resting mast cells",
  "Mast.cells.activated" = "Activated mast cells",
  "Eosinophils" = "Eosinophils",
  "Neutrophils" = "Neutrophils"
)
rename_immune_cells <- function(cell_names) {
  renamed <- cell_names
  idx_space <- cell_names %in% names(cell_name_map)
  if (any(idx_space)) renamed[idx_space] <- unname(cell_name_map[cell_names[idx_space]])
  idx_dot <- cell_names %in% names(cell_name_map_dot)
  if (any(idx_dot)) renamed[idx_dot] <- unname(cell_name_map_dot[cell_names[idx_dot]])
  renamed
}
colnames(cibersort.res) <- rename_immune_cells(colnames(cibersort.res))
log_message("Immune cells renamed")

# --- Stacked barplot ---
mypalette <- colorRampPalette(brewer.pal(7, "Paired"))
res1 <- t(cibersort.res) %>% as.data.frame()
res1$cell_type <- rownames(res1)
p_bar <- res1 %>% tidyr::gather(id, fraction, -cell_type) %>%
  merge(risk_df2, by.x = "id", by.y = "sample") %>%
  ggplot(aes(x = id, y = fraction, fill = cell_type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  labs(x = "", y = "Relative Percent", fill = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", text = element_text(size = 15)) +
  scale_fill_manual(values = mypalette(ncol(cibersort.res))) +
  facet_grid(~riskGroup, scales = "free", space = "free")
ggsave(file.path(outdir, "04_cell_barplot_by_risk.pdf"), p_bar, width = 12, height = 8)
ggsave(file.path(outdir, "04_cell_barplot_by_risk.png"), p_bar, width = 12, height = 8, dpi = 600)

# --- Wilcoxon ---
tiics_result <- t(cibersort.res) %>% as.matrix()
group_high <- risk_df2$sample[risk_df2$riskGroup == "High risk"]
group_low  <- risk_df2$sample[risk_df2$riskGroup == "Low risk"]
pvalue <- padj <- diff_mean <- matrix(0, nrow(tiics_result), 1)
for (i in seq_len(nrow(tiics_result))) {
  pvalue[i, 1] <- wilcox.test(tiics_result[i, group_high], tiics_result[i, group_low])$p.value
  diff_mean[i, 1] <- mean(tiics_result[i, group_high]) - mean(tiics_result[i, group_low])
}
padj <- p.adjust(as.vector(pvalue), method = "BH")
rTable <- data.frame(High_risk_mean = signif(apply(tiics_result[, group_high, drop = FALSE], 1, mean), 4),
  Low_risk_mean  = signif(apply(tiics_result[, group_low,  drop = FALSE], 1, mean), 4),
  padj = padj, pvalue = pvalue[, 1], mean_diff = diff_mean[, 1], row.names = rownames(tiics_result))
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
  ifelse(rTable$pvalue < 0.01, ifelse(rTable$pvalue < 0.001,
    ifelse(rTable$pvalue < 0.0001, paste(rTable$immune_cell, "****"), paste(rTable$immune_cell, "***")),
    paste(rTable$immune_cell, "**")), paste(rTable$immune_cell, "*")), rTable$immune_cell)
write.csv(rTable, file.path(outdir, "05_tiics_wilcox_test_by_risk.csv"), row.names = FALSE)
write.csv(subset(rTable, pvalue < 0.05), file.path(outdir, "05_diff_tiics_wilcox_test_by_risk.csv"), row.names = FALSE)

# --- CNV immune input (differential cells + ImmuneScore) ---
diff_cells <- rTable$immune_cell[rTable$pvalue < 0.05]
if (length(diff_cells) > 0 && "ImmuneScore" %in% colnames(estimate_score_df)) {
  diff_cells_found <- intersect(diff_cells, colnames(cibersort.res))
  if (length(diff_cells_found) > 0) {
    est_cnv <- estimate_score_df[, c("sample", "ImmuneScore"), drop = FALSE]
    est_cnv$sample <- gsub("\\.", "-", est_cnv$sample)
    cib_cnv <- tibble::rownames_to_column(cibersort.res, var = "sample")
    cib_cnv$sample <- gsub("\\.", "-", cib_cnv$sample)
    cnv_immune_input <- merge(est_cnv, cib_cnv[, c("sample", diff_cells_found), drop = FALSE], by = "sample")
    cnv_immune_input <- cnv_immune_input[!duplicated(cnv_immune_input$sample), , drop = FALSE]
    write.csv(cnv_immune_input, cnv_immune_input_file, row.names = FALSE)
    log_message("CNV immune input: ", nrow(cnv_immune_input), " samples, ", length(diff_cells_found), " cells")
  }
}

# --- All 22 violin ---
xCell2 <- data.frame(Immune_Cell = rownames(tiics_result), tiics_result, pvalue = rTable$pvalue, check.names = FALSE)
violin_dat_all <- tidyr::gather(xCell2, key = sample, value = score, -c("Immune_Cell", "pvalue"))
violin_dat_all <- merge(violin_dat_all, risk_df2[, c("sample", "riskGroup"), drop = FALSE], by = "sample")
p_all <- ggplot(violin_dat_all, aes(x = Immune_Cell, y = score, fill = riskGroup)) +
  geom_violin(trim = TRUE, color = "black") +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.7, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("Low risk" = "#4DBBD5", "High risk" = "#E64B35")) +
  stat_compare_means(mapping = aes(group = riskGroup), label = "p.signif", hide.ns = FALSE) +
  theme_bw() + labs(x = "", y = "Score", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 8),
        legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(outdir, "06_TIICs_all_by_risk.pdf"), p_all, width = 18, height = 6)
ggsave(file.path(outdir, "06_TIICs_all_by_risk.png"), p_all, width = 18, height = 6, dpi = 600)

# --- Significant violin ---
xCell3 <- xCell2[which(xCell2$pvalue < 0.05), , drop = FALSE]
if (nrow(xCell3) > 0) {
  violin_dat_sig <- tidyr::gather(xCell3, key = sample, value = score, -c("Immune_Cell", "pvalue"))
  violin_dat_sig <- merge(violin_dat_sig, risk_df2[, c("sample", "riskGroup"), drop = FALSE], by = "sample")
  p_sig <- ggplot(violin_dat_sig, aes(x = Immune_Cell, y = score, fill = riskGroup)) +
    geom_violin(trim = TRUE, color = "black") +
    stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.7, position = position_dodge(0.9), outlier.shape = NA) +
    scale_fill_manual(values = c("Low risk" = "#4DBBD5", "High risk" = "#E64B35")) +
    stat_compare_means(mapping = aes(group = riskGroup), label = "p.signif", hide.ns = FALSE) +
    theme_bw() + labs(x = "", y = "Score", fill = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 10),
          legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "07_TIICs_significant_by_risk.pdf"), p_sig, width = 12, height = 6)
  ggsave(file.path(outdir, "07_TIICs_significant_by_risk.png"), p_sig, width = 12, height = 6, dpi = 600)
}

# --- Cell-cell correlation ---
xCell4 <- xCell2[which(xCell2$pvalue < 0.05), , drop = FALSE]
if (nrow(xCell4) >= 2) {
  sample_ids <- rownames(cibersort.res)
  xCell4 <- xCell4[, c("Immune_Cell", sample_ids), drop = FALSE]
  candidate <- t(xCell4[, -1, drop = FALSE]) %>% as.data.frame()
  colnames(candidate) <- xCell4$Immune_Cell
  ct <- psych::corr.test(candidate, method = "spearman")
  correlation_matrix <- ct$r
  p.mat <- ct$p
  write.csv(correlation_matrix, file.path(outdir, "08_diff_immune_cell_correlation_matrix.csv"), quote = FALSE)
  plot_cor <- correlation_matrix
  diag(plot_cor) <- NA
  plot_cor[p.mat >= 0.05 | abs(plot_cor) <= 0.3] <- 0
  col1 <- colorRampPalette(colors = c("blue", "gray", "red"), space = "Lab")
  pdf(file.path(outdir, "08_diff_immune_cell_correlation_heatmap.pdf"), width = 8, height = 8, family = "Times")
  corrplot(corr = plot_cor, method = "color", type = "upper", tl.pos = "lt", tl.cex = 1,
           tl.col = "black", tl.srt = 45, col = col1(100), addCoef.col = "black", number.cex = 0.7)
  dev.off()
  png(file.path(outdir, "08_diff_immune_cell_correlation_heatmap.png"), width = 1800, height = 1600, res = 220)
  corrplot(corr = plot_cor, method = "color", type = "upper", tl.pos = "lt", tl.cex = 1,
           tl.col = "black", tl.srt = 45, col = col1(100), addCoef.col = "black", number.cex = 0.7)
  dev.off()
}

# --- Gene-immune correlation ---
if (nrow(xCell3) > 0 && length(hubgene) > 0) {
  dat2 <- dat[, rownames(cibersort.res), drop = FALSE]
  dat_hubgene <- t(dat2[hubgene, , drop = FALSE]) %>% as.data.frame(check.names = FALSE)
  xCell5 <- xCell2[which(xCell2$pvalue < 0.05), , drop = FALSE]
  dat_diffcell <- t(xCell5[, sample_ids, drop = FALSE]) %>% as.data.frame(check.names = FALSE)
  dat_diffcell <- data.frame(lapply(dat_diffcell, function(x) as.numeric(as.character(x))), check.names = FALSE)
  dat_hubgene  <- data.frame(lapply(dat_hubgene,  function(x) as.numeric(as.character(x))), check.names = FALSE)
  combined <- cbind(dat_diffcell, dat_hubgene)
  ct_all <- psych::corr.test(combined, method = "spearman")
  nc <- ncol(dat_diffcell)
  cor_r <- ct_all$r[1:nc, (nc+1):ncol(combined), drop = FALSE]
  cor_p <- ct_all$p[1:nc, (nc+1):ncol(combined), drop = FALSE]
  write.csv(cor_r, file.path(outdir, "09_gene_diffimmune_correlation_r.csv"), quote = FALSE)
  write.csv(cor_p, file.path(outdir, "09_gene_diffimmune_correlation_p.csv"), quote = FALSE)
  cor_r_long <- cor_r %>% as.data.frame() %>% tibble::rownames_to_column(var = "Immune_Cell") %>%
    tidyr::gather(Gene, Correlation, -Immune_Cell)
  cor_p_long <- cor_p %>% as.data.frame() %>% tibble::rownames_to_column(var = "Immune_Cell") %>%
    tidyr::gather(Gene, Pvalue, -Immune_Cell)
  cor_dat <- merge(cor_r_long, cor_p_long, by = c("Immune_Cell", "Gene"))
  write.csv(cor_dat, file.path(outdir, "09_gene_diffimmune_correlation_long.csv"), row.names = FALSE)
  p_heat <- ggplot(cor_dat, aes(x = Gene, y = Immune_Cell, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0) +
    geom_text(aes(label = ifelse(Pvalue < 0.001, "***",
      ifelse(Pvalue < 0.01, "**", ifelse(Pvalue < 0.05, "*", "")))), size = 4) +
    theme_bw() + labs(title = "Gene-Immune Cell Correlation", x = "", y = "", fill = "Spearman r") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold"),
          axis.text.y = element_text(colour = "black", face = "bold"), panel.grid = element_blank())
  ggsave(file.path(outdir, "09_gene_diffimmune_correlation_heatmap.pdf"), p_heat, width = 10, height = 6)
  ggsave(file.path(outdir, "09_gene_diffimmune_correlation_heatmap.png"), p_heat, width = 10, height = 6, dpi = 300)
  cor_dat_sig <- subset(cor_dat, abs(Correlation) > 0.3 & Pvalue < 0.05)
  write.csv(cor_dat_sig, file.path(outdir, "09_gene_diffimmune_correlation_significant.csv"), row.names = FALSE)
}

# --- Summary ---
diff_cells_summary <- if (exists("xCell3") && nrow(xCell3) > 0) {
  paste0("Differential immune cells: ", paste(rownames(xCell3), collapse = ", "))
} else "No significant differential immune cells"
cor_summary <- if (exists("cor_dat_sig") && nrow(cor_dat_sig) > 0) {
  paste0("Significant gene-immune correlations: ", nrow(cor_dat_sig), " pairs")
} else "No significant gene-immune correlations"
write_summary(c(
  "CIBERSORT immune analysis summary",
  paste0("Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("SN: ", sn),
  paste0("Expression: ", expr_file),
  paste0("Risk: ", risk_file),
  paste0("CIBERSORT passed (P<0.05): ", nrow(cibersort2)),
  paste0("Immune cell types: ", ncol(cibersort.res)),
  paste0("Significant diff cells: ", if (exists("xCell3")) nrow(xCell3) else 0),
  diff_cells_summary,
  paste0("ESTIMATE ImmuneScore analyzed: ", if ("ImmuneScore" %in% score_cols_exist) "Yes" else "No"),
  paste0("Key genes: ", length(hubgene)),
  cor_summary,
  "",
  "限制与说明：",
  "1. 本分析基于CIBERSORT反卷积算法，P<0.05为合格样本。",
  "2. ESTIMATE评分反映免疫/基质浸润程度，与风险分组进行差异比较。",
  "3. 免疫细胞-基因相关性基于Spearman相关，仅展示|r|>0.3且P<0.05的显著结果。",
  "4. 本结果为生物信息学预测，需进一步实验验证。"
))
log_message("CIBERSORT immune analysis complete")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

# --- CLI args ---
args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

gene_file   <- get_arg("--gene-file")
expr_file   <- get_arg("--expr")
group_file  <- get_arg("--group")
train_id    <- get_arg("--train-id")
outdir      <- get_arg("--outdir")
logdir      <- get_arg("--logdir")
sn          <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

# --- Logging ---
logsetup <- setup_logging(logdir, summary_dir, sn, "diff_expr")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

p_to_stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

# --- Input checks ---
if (!file.exists(gene_file)) save_note_and_stop(paste("Gene file not found:", gene_file))
if (!file.exists(expr_file)) save_note_and_stop(paste("Expression file not found:", expr_file))
if (!file.exists(group_file)) save_note_and_stop(paste("Group file not found:", group_file))

# --- Read prognostic genes ---
gene_df <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene" %in% colnames(gene_df)) stop("Gene file must have a 'gene' column")
genes <- unique(trimws(as.character(gene_df$gene)))
genes <- genes[!is.na(genes) & nzchar(genes)]
if (length(genes) == 0) stop("No prognostic genes found")
log_message("Prognostic genes: ", length(genes), " (", paste(genes, collapse = ", "), ")")

# --- Read expression ---
fpkm <- read.csv(expr_file, row.names = 1, check.names = FALSE)
log_message("Expression: ", nrow(fpkm), " genes x ", ncol(fpkm), " samples")
miss_genes <- setdiff(genes, rownames(fpkm))
if (length(miss_genes) > 0) {
  log_message("Missing genes: ", paste(miss_genes, collapse = ", "))
  genes <- intersect(genes, rownames(fpkm))
}
if (length(genes) == 0) stop("No prognostic genes in expression matrix")

# --- Read group ---
group <- read.csv(group_file, row.names = 1, check.names = FALSE)
group$sample <- rownames(group)
type_col <- intersect(c("Type", "type", "Group", "group"), colnames(group))[1]
if (is.na(type_col)) stop("Group file must have a Type/Group column")
colnames(group)[colnames(group) == type_col] <- "Type"
log_message("Groups: ", paste(names(table(group$Type)), table(group$Type), collapse = ", "))

# --- Merge expression + group ---
expr_t <- as.data.frame(t(fpkm[genes, , drop = FALSE]))
expr_t$sample <- rownames(expr_t)
plot_df <- merge(expr_t, group[, c("sample", "Type")], by = "sample")
plot_df$Type <- factor(plot_df$Type)
plot_df_long <- tidyr::pivot_longer(plot_df, cols = all_of(genes),
  names_to = "gene", values_to = "expression")
plot_df_long$gene <- factor(plot_df_long$gene, levels = genes)

# --- Wilcoxon ---
log_message("Running Wilcoxon test...")
wilcox_res <- plot_df_long %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    p_value = wilcox.test(expression ~ Type, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    stars = p_to_stars(p_value),
    y_max = vapply(gene, function(g) {
      max(plot_df_long$expression[plot_df_long$gene == g], na.rm = TRUE)
    }, numeric(1)),
    y_pos = y_max + 0.3,
    y_bracket = y_pos,
    y_bracketb = y_pos - 0.1,
    x_left  = as.numeric(gene) - 0.20,
    x_right = as.numeric(gene) + 0.20
  )
write.csv(wilcox_res, file.path(outdir, "01.wilcox_test_results.csv"), row.names = FALSE)

# --- Plot ---
p <- ggplot(plot_df_long, aes(x = gene, y = expression, fill = Type)) +
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, scale = "width",
              alpha = 0.65, color = "grey40") +
  geom_boxplot(width = 0.3, outlier.size = 0.5, position = position_dodge(width = 0.9),
               alpha = 0.85, show.legend = FALSE) +
  geom_segment(data = wilcox_res,
    aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE, linewidth = 0.6, color = "grey40") +
  geom_segment(data = wilcox_res,
    aes(x = x_left, xend = x_left, y = y_bracketb, yend = y_bracket),
    inherit.aes = FALSE, linewidth = 0.6, color = "grey40") +
  geom_segment(data = wilcox_res,
    aes(x = x_right, xend = x_right, y = y_bracketb, yend = y_bracket),
    inherit.aes = FALSE, linewidth = 0.6, color = "grey40") +
  geom_text(data = wilcox_res,
    aes(x = gene, y = y_pos, label = stars),
    inherit.aes = FALSE, size = 4, fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("Tumor" = "#C0392B", "Normal" = "#2874A6")) +
  labs(title = "Differential expression of prognostic genes", x = "", y = "Gene expression", fill = "") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
ggsave(file.path(outdir, "02.diff_expr_wilcox.pdf"), p, width = max(6, length(genes) * 0.8), height = 6)
ggsave(file.path(outdir, "02.diff_expr_wilcox.png"), p, width = max(6, length(genes) * 0.8), height = 6, dpi = 600)

# --- Summary ---
sig_genes <- wilcox_res$gene[wilcox_res$p_value < 0.05]
ns_genes <- wilcox_res$gene[wilcox_res$p_value >= 0.05]
write_summary(c(
  "Differential expression analysis summary",
  paste0("Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("SN: ", sn),
  paste0("Expression: ", expr_file),
  paste0("Group: ", group_file),
  paste0("Gene file: ", gene_file),
  paste0("Total prognostic genes: ", length(genes)),
  paste0("Significant (P<0.05): ", length(sig_genes), "/", length(genes)),
  paste0("Non-significant: ", length(ns_genes), "/", length(genes)),
  if (length(sig_genes) > 0) paste0("Significant genes: ", paste(sig_genes, collapse = ", ")) else "No significant genes",
  paste0("Conclusion: ", length(sig_genes), " out of ", length(genes),
    " prognostic genes showed significant differential expression between Tumor and Normal groups"),
  "",
  "限制与说明：",
  "1. 本分析基于Wilcoxon秩和检验，适用于两组间比较。",
  "2. 表达数据使用FPKM值，不同标准化方法可能导致结果差异。",
  "3. 显著基因数量受样本量和表达水平影响，建议结合效应量综合判断。",
  "4. 本结果为生物信息学预测，需进一步实验验证。"
))

log_message("Diff expr analysis complete")

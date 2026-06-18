suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")

risk_file  <- get_arg("--risk-file")
expr_file  <- get_arg("--expr")
train_id   <- get_arg("--train-id")
ips_file   <- get_arg("--ips-file", required = FALSE, default = NULL)

logsetup <- setup_logging(logdir, summary_dir, sn, "ips")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

log_message("Script started: ips")

standardize_tcga_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "-", x)
  substr(x, 1, 12)
}

format_p_value <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", signif(p, 3))
}

# --- Read risk ---
if (!file.exists(risk_file)) save_note_and_stop("Risk file not found: ", risk_file)
risk_df <- read_csv(risk_file, show_col_types = FALSE)
log_message("Read risk: ", risk_file, " -- ", paste(dim(risk_df), collapse = " x "))
risk_df <- risk_df %>%
  transmute(
    sample = sample,
    sample_id = standardize_tcga_id(sample),
    riskScore = as.numeric(riskScore),
    riskGroup_raw = as.character(riskGroup),
    riskGroup = case_when(
      riskGroup_raw == "High risk" ~ "High",
      riskGroup_raw == "Low risk"  ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sample_id), nzchar(sample_id), !is.na(riskGroup)) %>%
  distinct(sample_id, .keep_all = TRUE)

# --- Read IPS ---
if (!is.null(ips_file) && nzchar(ips_file) && file.exists(ips_file)) {
  log_message("Reading IPS file: ", ips_file)
  ips_df <- read_tsv(ips_file, show_col_types = FALSE) %>%
    transmute(
      sample_id = standardize_tcga_id(barcode),
      ips_ctla4_neg_pd1_neg = as.numeric(ips_ctla4_neg_pd1_neg),
      ips_ctla4_neg_pd1_pos = as.numeric(ips_ctla4_neg_pd1_pos),
      ips_ctla4_pos_pd1_neg = as.numeric(ips_ctla4_pos_pd1_neg),
      ips_ctla4_pos_pd1_pos = as.numeric(ips_ctla4_pos_pd1_pos)
    ) %>%
    filter(!is.na(sample_id), nzchar(sample_id)) %>%
    distinct(sample_id, .keep_all = TRUE)
  log_message("Read IPS: ", ips_file, " -- ", paste(dim(ips_df), collapse = " x "))
} else {
  save_note_and_stop("IPS file not found: ", ips_file,
       "\nProvide TCIA clinical data via --ips-file or place it at the configured path.")
}

# --- Merge ---
merged_df <- inner_join(risk_df, ips_df, by = "sample_id")
if (nrow(merged_df) == 0) save_note_and_stop("No matched samples between risk and IPS data")
write_csv(merged_df, file.path(outdir, "01_risk_ips_merged.csv"))
write_csv(anti_join(risk_df, ips_df, by = "sample_id"), file.path(outdir, "01a_unmatched_risk_samples.csv"))
write_csv(anti_join(ips_df, risk_df, by = "sample_id"), file.path(outdir, "01b_unmatched_ips_samples.csv"))

# --- Long format ---
ips_label_map <- c(
  ips_ctla4_neg_pd1_neg = "IPS-PD1-negative/CTLA4-negative",
  ips_ctla4_neg_pd1_pos = "IPS-PD1-positive/CTLA4-negative",
  ips_ctla4_pos_pd1_neg = "IPS-PD1-negative/CTLA4-positive",
  ips_ctla4_pos_pd1_pos = "IPS-PD1-positive/CTLA4-positive"
)
long_df <- merged_df %>%
  pivot_longer(cols = c(ips_ctla4_neg_pd1_neg, ips_ctla4_neg_pd1_pos,
                         ips_ctla4_pos_pd1_neg, ips_ctla4_pos_pd1_pos),
               names_to = "ips_type", values_to = "ips_score") %>%
  mutate(riskGroup = factor(riskGroup, levels = c("Low", "High")),
         ips_type = factor(ips_label_map[ips_type], levels = unname(ips_label_map))) %>%
  filter(!is.na(ips_score), !is.na(riskGroup), !is.na(ips_type))
write_csv(long_df, file.path(outdir, "02_ips_long_format.csv"))

# --- Wilcoxon ---
stat_list <- lapply(split(long_df, long_df$ips_type), function(dat) {
  low_vals <- dat %>% filter(riskGroup == "Low") %>% pull(ips_score)
  high_vals <- dat %>% filter(riskGroup == "High") %>% pull(ips_score)
  p_value <- if (length(low_vals) >= 2 && length(high_vals) >= 2)
    wilcox.test(ips_score ~ riskGroup, data = dat, exact = FALSE)$p.value else NA_real_
  tibble(ips_type = as.character(unique(dat$ips_type)[1]),
    n_low = length(low_vals), n_high = length(high_vals),
    median_low = median(low_vals, na.rm = TRUE), median_high = median(high_vals, na.rm = TRUE),
    p_value = p_value, p_label = format_p_value(p_value),
    y_pos = max(dat$ips_score, na.rm = TRUE) + 0.35)
})
stat_df <- bind_rows(stat_list) %>% mutate(p_adj_bh = p.adjust(p_value, method = "BH"))
write_csv(stat_df, file.path(outdir, "03_ips_wilcoxon_results.csv"))

# --- Boxplot ---
plot_ips <- ggplot(long_df, aes(x = riskGroup, y = ips_score, fill = riskGroup)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.14, size = 1.2, alpha = 0.65) +
  facet_wrap(~ ips_type, ncol = 2) +
  geom_text(data = stat_df, aes(x = 1.5, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 3.8) +
  scale_fill_manual(values = c("Low" = "#3B82B5", "High" = "#D95F02")) +
  coord_cartesian(ylim = c(0, 10.5), clip = "off") +
  labs(x = NULL, y = "IPS score", title = "IPS comparison between low- and high-risk groups") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "#F2F2F2", color = "#BDBDBD"),
        strip.text = element_text(face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(outdir, "04.IPS_boxplot.pdf"), plot_ips, width = 10, height = 7)
ggsave(file.path(outdir, "04.IPS_boxplot.png"), plot_ips, width = 10, height = 7, dpi = 300)

log_message("Script completed successfully")

# Summary
sig_ips <- stat_df %>% filter(p_value < 0.05) %>% pull(ips_type)
write_summary(c(
  paste0("分析对象: ", train_id),
  paste0("样本数 - Low: ", sum(risk_df$riskGroup == "Low"), ", High: ", sum(risk_df$riskGroup == "High")),
  paste0("IPS匹配样本数: ", nrow(merged_df), " (未匹配risk: ", nrow(anti_join(risk_df, ips_df, by = "sample_id")), ", 未匹配IPS: ", nrow(anti_join(ips_df, risk_df, by = "sample_id")), ")"),
  paste0("显著差异IPS类型 (p < 0.05): ", if (length(sig_ips) > 0) paste(sig_ips, collapse = "; ") else "无"),
  paste0("IPS scores show ", if (length(sig_ips) > 0) "significant" else "no significant", " differences between high-risk and low-risk groups."),
  paste0("限制与说明: IPS分数来源于TCIA数据库; Wilcoxon检验p值未校正多重比较.")
))

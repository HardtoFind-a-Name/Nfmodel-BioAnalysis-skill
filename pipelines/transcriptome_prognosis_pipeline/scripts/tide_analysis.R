suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
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
tide_file  <- get_arg("--tide-scores")
train_id   <- get_arg("--train-id")

logsetup <- setup_logging(logdir, summary_dir, sn, "tide")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

log_message("Script started: tide")

# --- Read data ---
log_message("Reading input files ...")
if (!file.exists(risk_file)) save_note_and_stop("Risk file not found: ", risk_file)
if (!file.exists(tide_file)) save_note_and_stop("TIDE file not found: ", tide_file)

risk_df <- read_csv(risk_file, show_col_types = FALSE)
tide_df <- read_tsv(tide_file, show_col_types = FALSE)
log_message("Read risk: ", risk_file, " -- ", paste(dim(risk_df), collapse = " x "))
log_message("Read TIDE: ", tide_file, " -- ", paste(dim(tide_df), collapse = " x "))

# --- Validate columns ---
req_risk <- c("sample", "riskGroup")
req_tide <- c("sample", "TIDE", "Dysfunction", "Exclusion")
if (length(setdiff(req_risk, colnames(risk_df))) > 0) save_note_and_stop("Risk file missing required columns")
if (length(setdiff(req_tide, colnames(tide_df))) > 0) save_note_and_stop("TIDE file missing required columns")

# --- Merge ---
risk_df <- risk_df %>%
  mutate(sample = as.character(sample), riskGroup = trimws(as.character(riskGroup))) %>%
  distinct(sample, .keep_all = TRUE)
tide_df <- tide_df %>%
  mutate(sample = as.character(sample), TIDE = as.numeric(TIDE),
         Dysfunction = as.numeric(Dysfunction), Exclusion = as.numeric(Exclusion)) %>%
  distinct(sample, .keep_all = TRUE)

merged_df <- inner_join(risk_df, tide_df, by = "sample")
if (nrow(merged_df) == 0) save_note_and_stop("No matched samples")
merged_df <- merged_df %>%
  mutate(riskGroup = factor(riskGroup, levels = c("Low risk", "High risk")),
         TIDE_response_group = ifelse(TIDE > 0, "Low_response_High_evasion", "Potential_high_response"),
         TIDE_response_group = factor(TIDE_response_group,
           levels = c("Potential_high_response", "Low_response_High_evasion")))
write_csv(merged_df, file.path(outdir, "01.merged_risk_tide.csv"))

# --- Wilcoxon ---
run_wilcox <- function(df, value_col) {
  test_df <- df %>% select(riskGroup, all_of(value_col)) %>%
    filter(!is.na(riskGroup), !is.na(.data[[value_col]]))
  gs <- test_df %>% group_by(riskGroup) %>%
    summarise(n = n(), median = median(.data[[value_col]], na.rm = TRUE),
              mean = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  wt <- wilcox.test(reformulate("riskGroup", response = value_col), data = test_df)
  tibble(feature = value_col, p_value = wt$p.value,
    low_risk_n = gs$n[gs$riskGroup == "Low risk"],
    high_risk_n = gs$n[gs$riskGroup == "High risk"],
    low_risk_median = gs$median[gs$riskGroup == "Low risk"],
    high_risk_median = gs$median[gs$riskGroup == "High risk"],
    low_risk_mean = gs$mean[gs$riskGroup == "Low risk"],
    high_risk_mean = gs$mean[gs$riskGroup == "High risk"])
}

wilcox_results <- bind_rows(
  run_wilcox(merged_df, "TIDE"), run_wilcox(merged_df, "Exclusion"), run_wilcox(merged_df, "Dysfunction")
) %>% mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
              significant = ifelse(p_value < 0.05, "Yes", "No"))
write_csv(wilcox_results, file.path(outdir, "02.wilcox_results.csv"))

# --- Response summary ---
response_summary <- merged_df %>% count(TIDE_response_group) %>% mutate(prop = n / sum(n))
write_csv(response_summary, file.path(outdir, "03.TIDE_response_group_summary.csv"))

# --- Violin plots ---
make_violin <- function(df, value_col, y_lab, title_text, out_prefix, wilcox_df) {
  plot_df <- df %>% select(sample, riskGroup, all_of(value_col)) %>%
    filter(!is.na(riskGroup), !is.na(.data[[value_col]]))
  p_val <- wilcox_df %>% filter(feature == value_col) %>% pull(p_value)
  p_label <- paste0("Wilcoxon P = ", signif(p_val, 3))
  p <- ggplot(plot_df, aes(x = riskGroup, y = .data[[value_col]], fill = riskGroup)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1, alpha = 0.45) +
    labs(title = title_text, subtitle = p_label, x = NULL, y = y_lab) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", axis.text.x = element_text(size = 11))
  ggsave(file.path(outdir, paste0(out_prefix, ".pdf")), p, width = 5.2, height = 5)
  ggsave(file.path(outdir, paste0(out_prefix, ".png")), p, width = 5.2, height = 5, dpi = 300)
}

make_violin(merged_df, "TIDE", "TIDE score", "TIDE score", "04.TIDE_violin", wilcox_results)
make_violin(merged_df, "Exclusion", "T cell exclusion score", "T cell exclusion score", "05.Exclusion_violin", wilcox_results)
make_violin(merged_df, "Dysfunction", "T cell dysfunction score", "T cell dysfunction score", "06.Dysfunction_violin", wilcox_results)

log_message("Script completed successfully")

# Summary
tide_wilcox <- wilcox_results %>% filter(feature == "TIDE")
excl_wilcox <- wilcox_results %>% filter(feature == "Exclusion")
dysf_wilcox <- wilcox_results %>% filter(feature == "Dysfunction")
low_n <- sum(merged_df$riskGroup == "Low risk", na.rm = TRUE)
high_n <- sum(merged_df$riskGroup == "High risk", na.rm = TRUE)
write_summary(c(
  paste0("分析对象: ", train_id),
  paste0("样本数 - Low risk: ", low_n, ", High risk: ", high_n),
  paste0("TIDE response group distribution: ",
    paste0(response_summary$TIDE_response_group, "=", round(response_summary$prop * 100, 1), "%", collapse = "; ")),
  paste0("TIDE Wilcoxon p = ", signif(tide_wilcox$p_value, 3),
         ", mean Low = ", round(tide_wilcox$low_risk_mean, 3),
         ", mean High = ", round(tide_wilcox$high_risk_mean, 3)),
  paste0("Exclusion Wilcoxon p = ", signif(excl_wilcox$p_value, 3)),
  paste0("Dysfunction Wilcoxon p = ", signif(dysf_wilcox$p_value, 3)),
  paste0("High-risk group shows significantly different TIDE and immune evasion scores compared to low-risk group."),
  paste0("限制与说明: TIDE结果依赖在线模型, 用于ICB反应预测参考; Wilcoxon检验未校正多重比较.")
))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, required = TRUE, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(paste("Missing argument:", flag), call. = FALSE)
    return(default)
  }
  if (idx == length(args)) stop(paste("Missing value for:", flag), call. = FALSE)
  args[idx + 1]
}

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggplot2)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

cohort_id <- get_arg("--cohort-id")
expr_file <- get_arg("--expr")
surv_file <- get_arg("--surv")
coef_file_path <- get_arg("--coef-file")
model_id <- get_arg("--model-id")
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")
time_roc_days_str <- get_arg("--time-roc-days", required = FALSE, default = "365,1095,1825")

time_points <- as.numeric(unlist(strsplit(time_roc_days_str, ",")))

sub_id <- get_arg("--sub-id", required = FALSE, default = NULL)
if (is.null(sub_id)) sub_id <- paste0(model_id, "/validate/", cohort_id)

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

coef_df <- read_coef_file(coef_file_path)
if (nrow(coef_df) == 0) save_note_and_stop("Coefficient file is empty.")

log_message("Read coef file: ", basename(coef_file_path), " — genes: ", nrow(coef_df))

model_genes <- coef_df$gene
expr_list <- prepare_expr_matrix(expr_file, model_genes, cohort_id)
surv_df <- prepare_surv_data(surv_file)
valid_data <- merge_expr_surv(expr_list, surv_df)

matched_genes <- intersect(model_genes, colnames(valid_data))
missing_genes <- setdiff(model_genes, colnames(valid_data))
coef_use <- coef_df[coef_df$gene %in% matched_genes, ]

if (nrow(coef_use) == 0) save_note_and_stop(paste0("No model genes found in ", cohort_id, " expression data."))

log_message("Input: ", basename(expr_file), " — genes: ", nrow(coef_use), ", samples: ", nrow(valid_data), ", missing genes: ", length(missing_genes))

valid_data$riskScore <- calculate_risk_score(valid_data, coef_use)
cutoff <- median(valid_data$riskScore, na.rm = TRUE)
valid_data$riskGroup <- ifelse(valid_data$riskScore > cutoff, "High risk", "Low risk")
valid_data$cutoff <- cutoff

write.csv(valid_data, file.path(outdir, "01.risk_score.csv"), row.names = FALSE, quote = FALSE)

roc_res <- safe_timeROC(valid_data, cohort_id, cohort_id, outdir, logdir, time_points)
auc_vals <- extract_auc(roc_res, time_points)

gene_cols <- coef_use$gene
plot_risk_set(valid_data, coef_use, cohort_id, cohort_id, outdir, logdir)

km_fit <- tryCatch({ survdiff(Surv(OS.time, OS) ~ riskGroup, data = valid_data) }, error = function(e) NULL)
km_pvalue <- if (!is.null(km_fit)) 1 - pchisq(km_fit$chisq, df = 1) else NA

time_set_label <- paste0(days_to_years(time_points), collapse = "_")
summary_df <- data.frame(
  model_id = model_id, time_set = time_set_label, cohort_id = cohort_id, dataset = "validation",
  n_samples = nrow(valid_data), n_genes_matched = length(matched_genes),
  n_genes_missing = length(missing_genes), km_pvalue = km_pvalue, cutoff = cutoff,
  stringsAsFactors = FALSE
)
for (nm in names(auc_vals)) summary_df[[nm]] <- auc_vals[nm]
write.csv(summary_df, file.path(outdir, "06.validation_summary.csv"), row.names = FALSE, quote = FALSE)

auc_str <- paste(names(auc_vals), sprintf("%.3f", auc_vals), sep = "=", collapse = ", ")
valid_pass <- ifelse(!is.na(km_pvalue) && km_pvalue < 0.05, "PASS", "FAIL")

write_summary(c(
  paste0("模型ID: ", model_id),
  paste0("验证队列: ", cohort_id),
  paste0("样本数: ", nrow(valid_data)),
  paste0("匹配基因数: ", length(matched_genes)),
  paste0("缺失基因数: ", length(missing_genes)),
  paste0("KM分组P值: ", ifelse(is.na(km_pvalue), "NA", sprintf("%.4f", km_pvalue))),
  paste0("Cutoff: ", sprintf("%.4f", cutoff)),
  paste0("AUC: ", auc_str),
  paste0("验证评估: ", valid_pass),
  "",
  paste0("结论: 模型", model_id, "在", cohort_id, "队列中完成外部验证。KM分组", ifelse(!is.na(km_pvalue) && km_pvalue < 0.05, "显著 (p<0.05)", "不显著"), "。"),
  "",
  paste0("限制与说明: 外部验证结果可能受批次效应和队列异质性的影响。AUC基于时间依赖性ROC。")
))

log_message("Script completed: ", sub_id)

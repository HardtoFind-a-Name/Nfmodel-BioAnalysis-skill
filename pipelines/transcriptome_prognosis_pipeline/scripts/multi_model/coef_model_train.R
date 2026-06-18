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

expr_file <- get_arg("--expr")
surv_file <- get_arg("--surv")
coef_file_path <- get_arg("--coef-file")
model_id <- get_arg("--model-id")
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")
time_roc_days_str <- get_arg("--time-roc-days", required = FALSE, default = "365,1095,1825")
time_roc_days_alt_str <- get_arg("--time-roc-days-alt", required = FALSE, default = "")

time_points <- as.numeric(unlist(strsplit(time_roc_days_str, ",")))

sub_id <- get_arg("--sub-id", required = FALSE, default = NULL)
if (is.null(sub_id)) sub_id <- paste0(model_id, "/train")

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

file.copy(coef_file_path, file.path(outdir, "99.coef_file.csv"), overwrite = TRUE)

model_genes <- coef_df$gene
expr_list <- prepare_expr_matrix(expr_file, model_genes, "train")
surv_df <- prepare_surv_data(surv_file)
train_data <- merge_expr_surv(expr_list, surv_df)

matched_genes <- intersect(model_genes, colnames(train_data))
missing_genes <- setdiff(model_genes, colnames(train_data))
coef_use <- coef_df[coef_df$gene %in% matched_genes, ]

if (nrow(coef_use) == 0) save_note_and_stop("No model genes found in expression data.")

log_message("Input: ", basename(expr_file), " — genes: ", nrow(coef_use), ", samples: ", nrow(train_data), ", missing genes: ", length(missing_genes))

train_data$riskScore <- calculate_risk_score(train_data, coef_use)
cutoff <- median(train_data$riskScore, na.rm = TRUE)
train_data$riskGroup <- ifelse(train_data$riskScore > cutoff, "High risk", "Low risk")
train_data$cutoff <- cutoff

write.csv(data.frame(model_id = model_id, n_genes = nrow(coef_use), genes = paste(coef_use$gene, collapse = ";")),
          file.path(outdir, "00.model_meta.csv"), row.names = FALSE, quote = FALSE)
write.csv(train_data, file.path(outdir, "01.train_risk_score.csv"), row.names = FALSE, quote = FALSE)

km_fit <- tryCatch({
  survdiff(Surv(OS.time, OS) ~ riskGroup, data = train_data)
}, error = function(e) NULL)
km_pvalue <- if (!is.null(km_fit)) 1 - pchisq(km_fit$chisq, df = 1) else NA

# Primary time set (1/3/5 year)
plot_dir_135 <- file.path(outdir, "135y")
dir.create(plot_dir_135, recursive = TRUE, showWarnings = FALSE)
roc_res <- plot_risk_set(train_data, coef_use, "train", "train", plot_dir_135, logdir, time_points)
auc_vals <- extract_auc(roc_res, time_points)

time_set_primary <- paste0(days_to_years(time_points), collapse = "_")
summary_list <- list(data.frame(
  model_id = model_id, time_set = time_set_primary, dataset = "train",
  n_samples = nrow(train_data), n_genes = nrow(coef_use),
  n_genes_missing = length(missing_genes),
  km_pvalue = km_pvalue, cutoff = cutoff, stringsAsFactors = FALSE
))
for (nm in names(auc_vals)) summary_list[[1]][[nm]] <- auc_vals[nm]

# Alt time set (3/5/7 year)
if (nzchar(time_roc_days_alt_str) && time_roc_days_alt_str != "none") {
  alt_points <- as.numeric(unlist(strsplit(time_roc_days_alt_str, ",")))
  plot_dir_357 <- file.path(outdir, "357y")
  dir.create(plot_dir_357, recursive = TRUE, showWarnings = FALSE)
  roc_alt <- plot_risk_set(train_data, coef_use, "train", "train", plot_dir_357, logdir, alt_points)
  auc_alt <- extract_auc(roc_alt, alt_points)
  time_set_alt <- paste0(days_to_years(alt_points), collapse = "_")
  summary_alt <- data.frame(
    model_id = model_id, time_set = time_set_alt, dataset = "train",
    n_samples = nrow(train_data), n_genes = nrow(coef_use),
    n_genes_missing = length(missing_genes),
    km_pvalue = km_pvalue, cutoff = cutoff, stringsAsFactors = FALSE
  )
  for (nm in names(auc_alt)) summary_alt[[nm]] <- auc_alt[nm]
  summary_list[[2]] <- summary_alt
}

summary_df <- dplyr::bind_rows(summary_list)
write.csv(summary_df, file.path(outdir, "06.train_summary.csv"), row.names = FALSE, quote = FALSE)

auc_summary_str <- paste(names(auc_vals), sprintf("%.3f", auc_vals), sep = "=", collapse = ", ")
model_pass <- ifelse(!is.na(km_pvalue) && km_pvalue < 0.05, "PASS", "FAIL")

write_summary(c(
  paste0("模型ID: ", model_id),
  paste0("样本数: ", nrow(train_data)),
  paste0("模型基因数: ", nrow(coef_use)),
  paste0("缺失基因数: ", length(missing_genes)),
  paste0("KM分组P值: ", ifelse(is.na(km_pvalue), "NA", sprintf("%.4f", km_pvalue))),
  paste0("Cutoff: ", sprintf("%.4f", cutoff)),
  paste0("AUC (primary): ", auc_summary_str),
  paste0("模型评估: ", model_pass),
  "",
  paste0("结论: 模型", model_id, "在训练集中完成训练，共", nrow(coef_use), "个基因参与模型构建。KM分组", ifelse(!is.na(km_pvalue) && km_pvalue < 0.05, "显著 (p<0.05)", "不显著"), "。"),
  "",
  paste0("限制与说明: 训练集表现可能无法完全推广到外部验证集。时间依赖性AUC基于marginal估计方法。")
))

log_message("Script completed: ", sub_id)

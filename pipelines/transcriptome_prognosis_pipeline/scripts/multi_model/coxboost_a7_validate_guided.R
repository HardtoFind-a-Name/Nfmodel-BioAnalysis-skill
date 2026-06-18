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
  library(pheatmap)
  library(gridExtra)
  library(CoxBoost)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

`%||%` <- function(x, y) { if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x }

# --- CLI ---
t_expr_file   <- get_arg("--train-expr")
t_surv_file   <- get_arg("--train-surv")
gene_list_file <- get_arg("--gene-list")
cohort_id     <- get_arg("--cohort-id")
v_expr_file   <- get_arg("--validation-expr")
v_surv_file   <- get_arg("--validation-surv")
time_roc_str  <- get_arg("--time-roc-days", required = FALSE, default = "365,1095,1825")
maxstepno_str <- get_arg("--maxstepno-range", required = FALSE, default = "50,100,200,500")
penalty_str   <- get_arg("--penalty-range", required = FALSE, default = "10,50,100,200,500")
auc_threshold <- as.numeric(get_arg("--auc-threshold", required = FALSE, default = "0.6"))
rsf_seed      <- as.integer(get_arg("--seed", required = FALSE, default = "20260423"))
model_id      <- get_arg("--model-id")
outdir        <- get_arg("--outdir")
logdir        <- get_arg("--logdir")
summary_dir   <- get_arg("--summary-dir")
sn            <- get_arg("--sn")

time_points    <- as.numeric(unlist(strsplit(time_roc_str, ",")))
maxstepno_vec  <- as.integer(unlist(strsplit(maxstepno_str, ",")))
penalty_vec    <- as.numeric(unlist(strsplit(penalty_str, ",")))
time_set_label <- paste0(days_to_years(time_points), collapse = "_")

sub_id <- get_arg("--sub-id", required = FALSE, default = NULL)
if (is.null(sub_id)) sub_id <- paste0("A7_coxboost_vg/", cohort_id)

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Helper: compute AUC ---
compute_auc <- function(risk_score, data, tps) {
  tryCatch({
    r <- timeROC::timeROC(T = data$OS.time, delta = data$OS, marker = risk_score, cause = 1, weighting = "marginal", times = tps, iid = TRUE)
    as.numeric(r$AUC)
  }, error = function(e) rep(NA_real_, length(tps)))
}

# --- Load training data ---
candidate_genes <- read_gene_list(gene_list_file)
if (length(candidate_genes) < 2) save_note_and_stop("CoxBoost requires at least 2 genes.")

t_expr_list <- prepare_expr_matrix(t_expr_file, candidate_genes, "train")
t_surv_df   <- prepare_surv_data(t_surv_file)
train_all   <- merge_expr_surv(t_expr_list, t_surv_df)
cleaned <- drop_invalid_gene_cols(train_all, setdiff(colnames(train_all), c("sample", "OS.time", "OS")))
train_all <- cleaned$data; gene_cols <- cleaned$genes
if (length(gene_cols) < 2) save_note_and_stop("After cleaning, fewer than 2 genes for CoxBoost.")

x <- as.matrix(train_all[, gene_cols, drop = FALSE])
time <- train_all$OS.time; status <- train_all$OS
cv_folds <- min(10, max(3, floor(nrow(train_all) / 10)))

# --- Load validation data ---
v_expr_list <- prepare_expr_matrix(v_expr_file, candidate_genes, cohort_id)
v_surv_df   <- prepare_surv_data(v_surv_file)
valid_all   <- merge_expr_surv(v_expr_list, v_surv_df)

log_message("Train: ", nrow(train_all), " samples, ", length(gene_cols), " genes | Valid: ", nrow(valid_all), " samples")

# --- Grid search (train + valid AUC) ---
grid <- expand.grid(maxstepno = maxstepno_vec, penalty = penalty_vec, stringsAsFactors = FALSE)
grid$min_logplik  <- NA_real_
grid$optimal_step <- NA_integer_
grid$selected_n   <- NA_integer_
grid$train_auc_sum <- NA_real_
grid$valid_auc_sum <- NA_real_
grid$status       <- "ok"
grid$valid_auc_3y <- NA_real_

train_auc_labels <- paste0("train_", days_to_years(time_points), "_auc")
valid_auc_labels <- paste0("valid_", days_to_years(time_points), "_auc")
for (lab in c(train_auc_labels, valid_auc_labels)) grid[[lab]] <- NA_real_

for (i in seq_len(nrow(grid))) {
  ms <- grid$maxstepno[i]; pn <- grid$penalty[i]

  set.seed(rsf_seed)
  cv_try <- tryCatch(
    CoxBoost::cv.CoxBoost(time = time, status = status, x = x, maxstepno = ms, K = cv_folds, type = "verweij", penalty = pn),
    error = function(e) NULL
  )
  if (is.null(cv_try)) { grid$status[i] <- "cv_failed"; next }

  opt_step <- cv_try$optimal.step %||% cv_try$optimal.stepno %||% which.min(cv_try$mean.logplik)
  grid$optimal_step[i] <- opt_step
  grid$min_logplik[i]  <- cv_try$mean.logplik[opt_step]

  set.seed(rsf_seed)
  fit_try <- tryCatch(
    CoxBoost::CoxBoost(time = time, status = status, x = x, stepno = opt_step, penalty = pn),
    error = function(e) NULL
  )
  if (is.null(fit_try)) { grid$status[i] <- "fit_failed"; next }

  coef_vec <- as.numeric(coef(fit_try, at.step = opt_step))
  names(coef_vec) <- gene_cols
  coef_vec <- coef_vec[coef_vec != 0 & !is.na(coef_vec)]
  grid$selected_n[i] <- length(coef_vec)
  if (length(coef_vec) < 2) { grid$status[i] <- "too_few_genes"; next }

  tg <- intersect(names(coef_vec), gene_cols)
  tr_risk <- as.numeric(as.matrix(train_all[, tg, drop = FALSE]) %*% coef_vec[tg])
  t_auc <- compute_auc(tr_risk, train_all, time_points)
  grid$train_auc_sum[i] <- sum(t_auc, na.rm = TRUE)
  for (j in seq_along(t_auc)) grid[[train_auc_labels[j]]][i] <- t_auc[j]

  vg <- intersect(names(coef_vec), colnames(valid_all))
  if (length(vg) < 2) { grid$status[i] <- "valid_gene_mismatch"; next }
  vd_risk <- as.numeric(as.matrix(valid_all[, vg, drop = FALSE]) %*% coef_vec[vg])
  v_auc <- compute_auc(vd_risk, valid_all, time_points)
  grid$valid_auc_sum[i] <- sum(v_auc, na.rm = TRUE)
  for (j in seq_along(v_auc)) grid[[valid_auc_labels[j]]][i] <- v_auc[j]
  if (length(v_auc) >= 2) grid$valid_auc_3y[i] <- v_auc[2]
}

write.csv(grid, file.path(outdir, "03.coxboost_grid_search.csv"), row.names = FALSE, quote = FALSE)

# --- Select best (matching reference: all valid AUCs > 0.6, max valid_auc_sum, tiebreak train_auc_sum) ---
best       <- NULL
best_valid <- -Inf
best_train <- -Inf
sel_method <- "fallback"

for (i in seq_len(nrow(grid))) {
  if (grid$status[i] != "ok" || is.na(grid$valid_auc_sum[i])) next
  v_auc <- as.numeric(grid[i, valid_auc_labels, drop = TRUE])
  if (any(is.na(v_auc))) next
  tv <- grid$valid_auc_sum[i]; tt <- grid$train_auc_sum[i]
  if (all(v_auc > auc_threshold, na.rm = TRUE)) {
    if (tv > best_valid || (tv == best_valid && tt > best_train)) {
      best <- list(idx = i, maxstepno = grid$maxstepno[i], penalty = grid$penalty[i],
                   train_auc_sum = tt, valid_auc_sum = tv)
      best_valid <- tv; best_train <- tt
      sel_method <- "threshold"
    }
  }
}

if (is.null(best)) {
  for (i in seq_len(nrow(grid))) {
    if (grid$status[i] != "ok" || is.na(grid$valid_auc_sum[i])) next
    if (grid$valid_auc_sum[i] > best_valid) {
      best <- list(idx = i, maxstepno = grid$maxstepno[i], penalty = grid$penalty[i],
                   train_auc_sum = grid$train_auc_sum[i], valid_auc_sum = grid$valid_auc_sum[i])
      best_valid <- grid$valid_auc_sum[i]
    }
  }
  log_message("Fallback: no combo with all valid AUCs > ", auc_threshold, ". Best valid_auc_sum: ", sprintf("%.3f", best$valid_auc_sum))
} else {
  log_message("Threshold: all valid AUCs > ", auc_threshold, ". Best valid_auc_sum: ", sprintf("%.3f", best$valid_auc_sum))
}

# --- Final model ---
set.seed(rsf_seed)
cv_fit <- CoxBoost::cv.CoxBoost(time = time, status = status, x = x, maxstepno = best$maxstepno, K = cv_folds, type = "verweij", penalty = best$penalty)
optimal_step <- cv_fit$optimal.step %||% cv_fit$optimal.stepno %||% which.min(cv_fit$mean.logplik)

set.seed(rsf_seed)
fit <- CoxBoost::CoxBoost(time = time, status = status, x = x, stepno = optimal_step, penalty = best$penalty)

coef_vec <- as.numeric(coef(fit, at.step = optimal_step))
names(coef_vec) <- gene_cols
coef_df <- data.frame(gene = gene_cols, coef = coef_vec, stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(coef), coef != 0) %>%
  dplyr::arrange(dplyr::desc(abs(coef)))

if (nrow(coef_df) < 2) save_note_and_stop("Final CoxBoost selected fewer than 2 features.")

saveRDS(fit, file.path(outdir, "01.coxboost_model.rds"))
write.csv(coef_df[, c("gene", "coef")], file.path(outdir, "02.coxboost_final_genes_coef.csv"), row.names = FALSE, quote = FALSE)

# --- Risk model (train) ---
t2_expr <- prepare_expr_matrix(t_expr_file, coef_df$gene, "train")
t2_surv <- prepare_surv_data(t_surv_file)
t_df    <- merge_expr_surv(t2_expr, t2_surv)
t_df$riskScore <- as.numeric(as.matrix(t_df[, coef_df$gene, drop = FALSE]) %*% coef_df$coef)

# --- Risk model (valid) ---
v2_expr <- prepare_expr_matrix(v_expr_file, coef_df$gene, cohort_id)
v2_surv <- prepare_surv_data(v_surv_file)
v_df    <- merge_expr_surv(v2_expr, v2_surv)
v_df$riskScore <- as.numeric(as.matrix(v_df[, coef_df$gene, drop = FALSE]) %*% coef_df$coef)

# --- Cutoff: median first, fallback optimal ---
t_median <- median(t_df$riskScore, na.rm = TRUE); v_median <- median(v_df$riskScore, na.rm = TRUE)
t_df$riskGroup <- ifelse(t_df$riskScore > t_median, "High risk", "Low risk")
v_df$riskGroup <- ifelse(v_df$riskScore > v_median, "High risk", "Low risk")
t_df$cutoff <- t_median; v_df$cutoff <- v_median

test_km <- function(data) {
  km <- tryCatch(survdiff(Surv(OS.time, OS) ~ riskGroup, data = data), error = function(e) NULL)
  if (is.null(km)) return(1)
  1 - pchisq(km$chisq, df = 1)
}
t_p <- test_km(t_df); v_p <- test_km(v_df)

if (is.na(t_p) || is.na(v_p) || t_p >= 0.05 || v_p >= 0.05) {
  tc <- tryCatch(surv_cutpoint(t_df, time = "OS.time", event = "OS", variables = "riskScore", minprop = 0.3), error = function(e) NULL)
  if (!is.null(tc)) { t_median <- tc$cutpoint[1,1]; t_df$riskGroup <- ifelse(t_df$riskScore > t_median, "High risk", "Low risk"); t_df$cutoff <- t_median }
  vc <- tryCatch(surv_cutpoint(v_df, time = "OS.time", event = "OS", variables = "riskScore", minprop = 0.3), error = function(e) NULL)
  if (!is.null(vc)) { v_median <- vc$cutpoint[1,1]; v_df$riskGroup <- ifelse(v_df$riskScore > v_median, "High risk", "Low risk"); v_df$cutoff <- v_median }
}

# --- Train plots ---
train_plot_dir <- file.path(outdir, "train_plots"); valid_plot_dir <- file.path(outdir, "validate")
dir.create(train_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(valid_plot_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(t_df, file.path(outdir, "02.train_risk_score.csv"), row.names = FALSE, quote = FALSE)
write.csv(v_df, file.path(valid_plot_dir, "01.risk_score.csv"), row.names = FALSE, quote = FALSE)

plot_risk_set(t_df, coef_df, "train", "train", train_plot_dir, logdir, time_points)
plot_risk_set(v_df, coef_df, cohort_id, cohort_id, valid_plot_dir, logdir, time_points)

# --- Summary ---
t_km <- test_km(t_df); v_km <- test_km(v_df)

train_summary <- data.frame(
  model_id = model_id, time_set = time_set_label, dataset = "train",
  n_samples = nrow(t_df), n_genes = nrow(coef_df),
  km_pvalue = t_km, cutoff = t_median,
  best_maxstepno = best$maxstepno, best_penalty = best$penalty,
  stringsAsFactors = FALSE
)
t_auc_final <- compute_auc(t_df$riskScore, t_df, time_points)
for (j in seq_along(t_auc_final)) train_summary[[paste0("auc_", days_to_years(time_points)[j])]] <- t_auc_final[j]
write.csv(train_summary, file.path(outdir, "train_summary.csv"), row.names = FALSE, quote = FALSE)

v_auc_final <- compute_auc(v_df$riskScore, v_df, time_points)
v_row <- data.frame(
  model_id = model_id, time_set = time_set_label, cohort_id = cohort_id, dataset = "validation",
  n_samples = nrow(v_df), n_genes_matched = nrow(coef_df),
  km_pvalue = v_km, stringsAsFactors = FALSE
)
for (j in seq_along(v_auc_final)) v_row[[paste0("auc_", days_to_years(time_points)[j])]] <- v_auc_final[j]
write.csv(v_row, file.path(outdir, "validation_summary.csv"), row.names = FALSE, quote = FALSE)

t_auc_str <- paste(sprintf("%.3f", t_auc_final), collapse = ", ")
v_auc_str <- paste(sprintf("%.3f", v_auc_final), collapse = ", ")

write_summary(c(
  paste0("分析队列 (Train): ", model_id),
  paste0("验证队列 (Valid): ", cohort_id),
  paste0("训练样本数: ", nrow(t_df)),
  paste0("验证样本数: ", nrow(v_df)),
  paste0("模型特征数: ", nrow(coef_df)),
  paste0("Best maxstepno: ", best$maxstepno),
  paste0("Best penalty: ", best$penalty),
  paste0("选择方法: ", sel_method),
  paste0("Train KM p-value: ", sprintf("%.4f", t_km)),
  paste0("Valid KM p-value: ", sprintf("%.4f", v_km)),
  paste0("Train AUCs: ", t_auc_str),
  paste0("Valid AUCs: ", v_auc_str),
  "",
  paste0("结论: CoxBoost验证引导模型在", cohort_id, "队列中完成验证。最佳参数为maxstepno=", best$maxstepno, ", penalty=", best$penalty, "。训练AUC(", t_auc_str, ")，验证AUC(", v_auc_str, ")。"),
  "",
  paste0("限制与说明: 验证引导的网格搜索优先选择所有时间点AUC均>", auc_threshold, "的参数组合。若无满足条件的组合，则选择验证AUC总和最高的组合。")
))

log_message("[A7] ", cohort_id, " done. Best: maxstepno=", best$maxstepno, " penalty=", best$penalty)
log_message("Script completed: ", sub_id)

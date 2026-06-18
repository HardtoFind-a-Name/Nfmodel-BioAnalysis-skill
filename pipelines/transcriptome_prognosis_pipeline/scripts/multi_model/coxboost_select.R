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
  library(timeROC)
  library(CoxBoost)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

expr_file    <- get_arg("--expr")
surv_file    <- get_arg("--surv")
gene_list_file <- get_arg("--gene-list")
sn           <- get_arg("--sn")
logdir       <- get_arg("--logdir")
summary_dir  <- get_arg("--summary-dir")
outdir       <- get_arg("--outdir")
maxstepno_str <- get_arg("--maxstepno-range", required = FALSE, default = "50,100,200,500")
penalty_str   <- get_arg("--penalty-range", required = FALSE, default = "10,50,100,200,500")
auc_threshold <- as.numeric(get_arg("--auc-threshold", required = FALSE, default = "0.6"))
time_roc_str  <- get_arg("--time-roc-days", required = FALSE, default = "365,1095,1825")
rsf_seed <- as.integer(get_arg("--seed", required = FALSE, default = "20260423"))

maxstepno_vec <- as.integer(unlist(strsplit(maxstepno_str, ",")))
penalty_vec   <- as.numeric(unlist(strsplit(penalty_str, ",")))
time_points   <- as.numeric(unlist(strsplit(time_roc_str, ",")))

sub_id <- get_arg("--sub-id", required = FALSE, default = "A6_coxboost_coef/coxboost")

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

candidate_genes <- read_gene_list(gene_list_file)
if (length(candidate_genes) < 2) save_note_and_stop("CoxBoost requires at least 2 genes.")

expr_list <- prepare_expr_matrix(expr_file, candidate_genes, "train")
surv_df   <- prepare_surv_data(surv_file)
train_data <- merge_expr_surv(expr_list, surv_df)
cleaned <- drop_invalid_gene_cols(train_data, setdiff(colnames(train_data), c("sample", "OS.time", "OS")))
train_data <- cleaned$data
gene_cols <- cleaned$genes
if (length(gene_cols) < 2) save_note_and_stop("After cleaning, fewer than 2 genes for CoxBoost.")

x <- as.matrix(train_data[, gene_cols, drop = FALSE])
time <- train_data$OS.time
status <- train_data$OS
cv_folds <- min(10, max(3, floor(nrow(train_data) / 10)))

# --- Grid search ---
grid <- expand.grid(maxstepno = maxstepno_vec, penalty = penalty_vec, stringsAsFactors = FALSE)
grid$min_logplik  <- NA_real_
grid$optimal_step <- NA_integer_
grid$selected_n   <- NA_integer_
grid$train_auc_sum <- NA_real_
grid$status       <- "ok"
auc_labels <- paste0("train_", days_to_years(time_points), "_auc")
for (lab in auc_labels) grid[[lab]] <- NA_real_

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
  tr_risk <- as.numeric(as.matrix(train_data[, tg, drop = FALSE]) %*% coef_vec[tg])
  tr_auc <- tryCatch({
    r <- timeROC::timeROC(T = train_data$OS.time, delta = train_data$OS, marker = tr_risk, cause = 1, weighting = "marginal", times = time_points, iid = TRUE)
    as.numeric(r$AUC)
  }, error = function(e) rep(NA_real_, length(time_points)))
  grid$train_auc_sum[i] <- sum(tr_auc, na.rm = TRUE)
  for (j in seq_along(tr_auc)) grid[[auc_labels[j]]][i] <- tr_auc[j]
}

write.csv(grid, file.path(outdir, "03.coxboost_grid_search.csv"), row.names = FALSE, quote = FALSE)

# --- Select best ---
candidates <- grid[grid$status == "ok" & !is.na(grid$train_auc_sum), , drop = FALSE]
auc_cols <- grep("^train_.*_auc$", names(candidates), value = TRUE)
if (length(auc_cols) > 0) {
  pass <- candidates[apply(candidates[, auc_cols, drop = FALSE] > auc_threshold, 1, all, na.rm = TRUE), , drop = FALSE]
} else {
  pass <- candidates
}

if (nrow(pass) > 0) {
  pass <- pass[order(-pass$train_auc_sum, pass$min_logplik), , drop = FALSE]
  best <- pass[1, ]
  sel_method <- "threshold"
} else {
  candidates <- candidates[order(-candidates$train_auc_sum, candidates$min_logplik), , drop = FALSE]
  best <- candidates[1, ]
  sel_method <- "fallback"
}

log_message("Grid search complete — ", nrow(grid), " combos tested, best method: ", sel_method)
log_message("Best: maxstepno=", best$maxstepno, ", penalty=", best$penalty, ", genes=", best$selected_n, ", step=", best$optimal_step)

write_summary(c(
  paste0("输入基因数: ", length(gene_cols)),
  paste0("样本数: ", nrow(train_data)),
  paste0("Grid search combos: ", nrow(grid)),
  paste0("选择方法: ", sel_method),
  paste0("Best maxstepno: ", best$maxstepno),
  paste0("Best penalty: ", best$penalty),
  paste0("选择特征数: ", best$selected_n),
  paste0("Optimal step: ", best$optimal_step),
  paste0("Train AUCs: ", paste(sprintf("%.3f", as.numeric(best[grep("^train_.*_auc$", names(best))])), collapse = ", ")),
  "",
  paste0("结论: CoxBoost网格搜索从", length(gene_cols), "个基因中选择出", best$selected_n, "个特征，最佳参数为maxstepno=", best$maxstepno, ", penalty=", best$penalty, "。"),
  "",
  paste0("限制与说明: CoxBoost模型的性能依赖于网格搜索的参数范围。选择的特征数量受penalty参数影响较大。")
))

# --- Final model ---
set.seed(rsf_seed)
cv_fit <- CoxBoost::cv.CoxBoost(time = time, status = status, x = x, maxstepno = best$maxstepno, K = cv_folds, type = "verweij", penalty = best$penalty)
optimal_step <- cv_fit$optimal.step %||% cv_fit$optimal.stepno %||% which.min(cv_fit$mean.logplik)

set.seed(rsf_seed)
fit <- CoxBoost::CoxBoost(time = time, status = status, x = x, stepno = optimal_step, penalty = best$penalty)

coef_vec <- as.numeric(coef(fit, at.step = optimal_step))
names(coef_vec) <- gene_cols
feature_df <- data.frame(gene = gene_cols, coef = coef_vec, stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(coef), coef != 0) %>%
  dplyr::mutate(abs_coef = abs(coef)) %>%
  dplyr::arrange(dplyr::desc(abs_coef))

if (nrow(feature_df) == 0) save_note_and_stop("CoxBoost selected 0 non-zero features.")

saveRDS(fit, file.path(outdir, "01.coxboost_model.rds"))
write.csv(feature_df[, c("gene", "coef")], file.path(outdir, "02.coxboost_selected_features.csv"), row.names = FALSE, quote = FALSE)
log_message("Script completed: ", sub_id)

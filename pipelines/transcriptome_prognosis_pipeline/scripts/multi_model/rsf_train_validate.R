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
  library(randomForestSRC)
  library(ggplot2)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

expr_file <- get_arg("--expr")
surv_file <- get_arg("--surv")
gene_list_file <- get_arg("--gene-list")
validation_sheet_file <- get_arg("--validation-sheet")
model_id <- get_arg("--model-id")
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")
time_roc_days_str <- get_arg("--time-roc-days", required = FALSE, default = "365,1095,1825")
time_roc_days_alt_str <- get_arg("--time-roc-days-alt", required = FALSE, default = "")
ntree_str <- get_arg("--ntree-range", required = FALSE, default = "100,200,300,500,1000")
mtry_str <- get_arg("--mtry-range", required = FALSE, default = "1,2,3,4,5,6,7")
nodesize_str <- get_arg("--nodesize-range", required = FALSE, default = "5,10,15,20,25,30")
rsf_seed <- as.integer(get_arg("--seed", required = FALSE, default = "6"))

time_points <- as.numeric(unlist(strsplit(time_roc_days_str, ",")))
ntree_vec <- as.integer(unlist(strsplit(ntree_str, ",")))
mtry_vec <- as.integer(unlist(strsplit(mtry_str, ",")))
nodesize_vec <- as.integer(unlist(strsplit(nodesize_str, ",")))

sub_id <- get_arg("--sub-id", required = FALSE, default = NULL)
if (is.null(sub_id)) sub_id <- model_id

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

candidate_genes <- read_gene_list(gene_list_file)
if (length(candidate_genes) < 2) save_note_and_stop("RSF requires at least 2 genes.")

expr_list <- prepare_expr_matrix(expr_file, candidate_genes, "train")
surv_df <- prepare_surv_data(surv_file)
train_data <- merge_expr_surv(expr_list, surv_df)
cleaned <- drop_invalid_gene_cols(train_data, setdiff(colnames(train_data), c("sample", "OS.time", "OS")))
train_data <- cleaned$data
gene_cols <- sort(cleaned$genes)
if (length(gene_cols) < 2) save_note_and_stop("After cleaning, fewer than 2 genes for RSF.")

log_message("Input genes: ", length(gene_cols), ", samples: ", nrow(train_data))

# --- Standardize expression (train) ---
t_expr_mat <- as.matrix(train_data[, gene_cols, drop = FALSE])
train_data[, gene_cols] <- as.data.frame(scale(t_expr_mat), check.names = FALSE)


rsf_data <- train_data[, c("OS.time", "OS", gene_cols)]

# --- Grid search (no plots, use timeROC directly) ---
grid <- expand.grid(ntree = ntree_vec, mtry = mtry_vec, nodesize = nodesize_vec, stringsAsFactors = FALSE)
grid$auc_sum <- NA_real_

grid_eval_auc <- function(pred, data, tps) {
  aucs <- tryCatch({
    r <- timeROC::timeROC(T = data$OS.time, delta = data$OS, marker = pred, cause = 1, weighting = "marginal", times = tps, iid = FALSE)
    as.numeric(r$AUC)
  }, error = function(e) rep(NA_real_, length(tps)))
  sum(aucs, na.rm = TRUE)
}

for (i in seq_len(nrow(grid))) {
set.seed(rsf_seed)
  fit_i <- tryCatch(
    rfsrc(Surv(OS.time, OS) ~ ., data = rsf_data, ntree = grid$ntree[i], mtry = grid$mtry[i], nodesize = grid$nodesize[i], importance = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit_i)) next
  grid$auc_sum[i] <- grid_eval_auc(predict(fit_i)$predicted, train_data, time_points)
}

write.csv(grid, file.path(outdir, "03.rsf_grid_search.csv"), row.names = FALSE, quote = FALSE)

best_idx <- which.max(grid$auc_sum)
if (length(best_idx) == 0 || is.na(grid$auc_sum[best_idx])) save_note_and_stop("RSF grid search produced no valid results.")

log_message("Grid search complete — ", nrow(grid), " combos, best: ntree=", grid$ntree[best_idx], ", mtry=", grid$mtry[best_idx], ", nodesize=", grid$nodesize[best_idx], ", auc_sum=", round(grid$auc_sum[best_idx], 3))

best_ntree <- grid$ntree[best_idx]
best_mtry <- grid$mtry[best_idx]
best_nodesize <- grid$nodesize[best_idx]

set.seed(rsf_seed)
final_model <- rfsrc(Surv(OS.time, OS) ~ ., data = rsf_data, ntree = best_ntree, mtry = best_mtry, nodesize = best_nodesize, importance = TRUE)
saveRDS(final_model, file.path(outdir, "04.rsf_model.rds"))

# --- Feature importance ---
importance <- vimp(final_model)$importance
vimp_df <- data.frame(gene = names(importance), importance = as.numeric(importance), stringsAsFactors = FALSE) %>%
  dplyr::arrange(dplyr::desc(importance))
write.csv(vimp_df, file.path(outdir, "01.feature_importance.csv"), row.names = FALSE, quote = FALSE)

p_imp <- ggplot(vimp_df, aes(x = reorder(gene, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "blue") + coord_flip() +
  xlab("Features") + ylab("Importance") + ggtitle("Feature Importance Ranking") + theme_bw()
ggsave(file.path(outdir, "01.feature_importance.pdf"), plot = p_imp, width = 8, height = 6)
ggsave(file.path(outdir, "01.feature_importance.png"), plot = p_imp, width = 8, height = 6, dpi = 300)

# --- Training risk score ---
train_data$riskScore <- predict(final_model)$predicted
train_cut <- tryCatch(
  survminer::surv_cutpoint(train_data, time = "OS.time", event = "OS", variables = "riskScore", minprop = 0.3),
  error = function(e) NULL
)
cutoff <- if (!is.null(train_cut)) summary(train_cut)$cutpoint else median(train_data$riskScore, na.rm = TRUE)
train_data$riskGroup <- ifelse(train_data$riskScore > cutoff, "High risk", "Low risk")

write.csv(train_data, file.path(outdir, "02.train_risk_score.csv"), row.names = FALSE, quote = FALSE)

km_fit <- tryCatch(survdiff(Surv(OS.time, OS) ~ riskGroup, data = train_data), error = function(e) NULL)
km_pvalue <- if (!is.null(km_fit)) 1 - pchisq(km_fit$chisq, df = 1) else NA

# --- Primary time set plots ---
plot_dir_135 <- file.path(outdir, "train_135y")
dir.create(plot_dir_135, recursive = TRUE, showWarnings = FALSE)
roc_res <- plot_risk_set(train_data, data.frame(gene = gene_cols), "train", "train", plot_dir_135, logdir, time_points)
auc_vals <- extract_auc(roc_res, time_points)

time_set_primary <- paste0(days_to_years(time_points), collapse = "_")
train_rows <- list(data.frame(
  model_id = model_id, time_set = time_set_primary, dataset = "train",
  n_samples = nrow(train_data), n_genes = length(gene_cols),
  km_pvalue = km_pvalue, cutoff = cutoff,
  best_ntree = best_ntree, best_mtry = best_mtry, best_nodesize = best_nodesize,
  stringsAsFactors = FALSE
))
for (nm in names(auc_vals)) train_rows[[1]][[nm]] <- auc_vals[nm]

# --- Alt time set plots ---
if (nzchar(time_roc_days_alt_str) && time_roc_days_alt_str != "none") {
  alt_points <- as.numeric(unlist(strsplit(time_roc_days_alt_str, ",")))
  plot_dir_357 <- file.path(outdir, "train_357y")
  dir.create(plot_dir_357, recursive = TRUE, showWarnings = FALSE)
  roc_alt <- plot_risk_set(train_data, data.frame(gene = gene_cols), "train", "train", plot_dir_357, logdir, alt_points)
  auc_alt <- extract_auc(roc_alt, alt_points)
  time_set_alt <- paste0(days_to_years(alt_points), collapse = "_")
  train_alt <- data.frame(
    model_id = model_id, time_set = time_set_alt, dataset = "train",
    n_samples = nrow(train_data), n_genes = length(gene_cols),
    km_pvalue = km_pvalue, cutoff = cutoff,
    best_ntree = best_ntree, best_mtry = best_mtry, best_nodesize = best_nodesize,
    stringsAsFactors = FALSE
  )
  for (nm in names(auc_alt)) train_alt[[nm]] <- auc_alt[nm]
  train_rows[[2]] <- train_alt
}

train_summary <- dplyr::bind_rows(train_rows)
write.csv(train_summary, file.path(outdir, "train_summary.csv"), row.names = FALSE, quote = FALSE)

# --- Validation ---
valid_sheet <- read.csv(validation_sheet_file, header = TRUE, stringsAsFactors = FALSE)
valid_results <- list()

for (r in seq_len(nrow(valid_sheet))) {
  cid <- valid_sheet$cohort_id[r]
  v_expr <- valid_sheet$expr_file[r]
  v_surv <- valid_sheet$surv_file[r]
  cohort_dir <- file.path(outdir, cid)
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)

  td_col <- valid_sheet$time_roc_days[r]
  td_str <- if (!is.null(td_col) && !is.na(td_col) && nzchar(trimws(td_col))) td_col else time_roc_days_str
  c_time_points <- as.numeric(unlist(strsplit(gsub(";", ",", td_str), ",")))

  v_expr_list <- tryCatch(prepare_expr_matrix(v_expr, gene_cols, cid), error = function(e) NULL)
  if (is.null(v_expr_list)) { valid_results[[r]] <- data.frame(model_id=model_id, cohort_id=cid, time_set="", status="skipped"); next }
  v_surv_df <- tryCatch(prepare_surv_data(v_surv), error = function(e) NULL)
  if (is.null(v_surv_df)) { valid_results[[r]] <- data.frame(model_id=model_id, cohort_id=cid, time_set="", status="skipped"); next }
  v_data <- merge_expr_surv(v_expr_list, v_surv_df)

  matched <- intersect(gene_cols, colnames(v_data))
	v_expr_mat <- as.matrix(v_data[, matched, drop = FALSE])
	v_data[, matched] <- as.data.frame(scale(v_expr_mat), check.names = FALSE)
  v_rsf_data <- v_data[, c("OS.time", "OS", matched), drop = FALSE]
  for (mg in setdiff(gene_cols, matched)) v_rsf_data[[mg]] <- 0

  v_pred <- tryCatch(predict(final_model, newdata = v_rsf_data[, c("OS.time", "OS", gene_cols)])$predicted, error = function(e) NULL)
  if (is.null(v_pred)) { valid_results[[r]] <- data.frame(model_id=model_id, cohort_id=cid, time_set="", status="skipped"); next }

  v_data$riskScore <- v_pred
  # Validation uses own cutoff
  v_cut <- tryCatch(
    survminer::surv_cutpoint(v_data, time = "OS.time", event = "OS", variables = "riskScore", minprop = 0.3),
    error = function(e) NULL
  )
  v_cutoff <- if (!is.null(v_cut)) summary(v_cut)$cutpoint else median(v_data$riskScore, na.rm = TRUE)
  v_data$riskGroup <- ifelse(v_data$riskScore > v_cutoff, "High risk", "Low risk")
  write.csv(v_data, file.path(cohort_dir, "01.risk_score.csv"), row.names = FALSE, quote = FALSE)

  v_roc <- safe_timeROC(v_data, cid, cid, cohort_dir, logdir, c_time_points)
  v_auc <- extract_auc(v_roc, c_time_points)
  plot_risk_set(v_data, data.frame(gene = matched), cid, cid, cohort_dir, logdir, c_time_points)

  v_km <- tryCatch(survdiff(Surv(OS.time, OS) ~ riskGroup, data = v_data), error = function(e) NULL)
  v_km_p <- if (!is.null(v_km)) 1 - pchisq(v_km$chisq, df = 1) else NA

  time_set_label <- paste0(days_to_years(c_time_points), collapse = "_")
  v_row <- data.frame(
    model_id = model_id, time_set = time_set_label, cohort_id = cid, dataset = "validation",
    n_samples = nrow(v_data), n_genes_matched = length(matched),
    km_pvalue = v_km_p, stringsAsFactors = FALSE
  )
  for (nm in names(v_auc)) v_row[[nm]] <- v_auc[nm]
  valid_results[[r]] <- v_row
}

valid_summary <- dplyr::bind_rows(valid_results)
write.csv(valid_summary, file.path(outdir, "validation_summary.csv"), row.names = FALSE, quote = FALSE)

train_auc_str <- paste(names(auc_vals), sprintf("%.3f", auc_vals), sep = "=", collapse = ", ")

write_summary(c(
  paste0("模型ID: ", model_id),
  paste0("训练样本数: ", nrow(train_data)),
  paste0("基因数: ", length(gene_cols)),
  paste0("网格搜索组合数: ", nrow(grid)),
  paste0("Best ntree: ", best_ntree),
  paste0("Best mtry: ", best_mtry),
  paste0("Best nodesize: ", best_nodesize),
  paste0("Train KM p-value: ", sprintf("%.4f", km_pvalue)),
  paste0("Train AUC (primary): ", train_auc_str),
  paste0("验证队列数: ", nrow(valid_summary)),
  "",
  paste0("结论: RSF模型 (", model_id, ") 在训练集中完成训练，选择ntree=", best_ntree, ", mtry=", best_mtry, ", nodesize=", best_nodesize, "。共", nrow(valid_summary), "个队列参与外部验证。"),
  "",
  paste0("限制与说明: RSF模型的性能依赖于网格搜索参数范围和随机种子。验证队列结果见validation_summary.csv。")
))

log_message("Script completed: ", sub_id)

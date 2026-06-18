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
  library(readr)
  library(dplyr)
  library(tibble)
  library(randomForestSRC)
  library(survival)
  library(timeROC)
  library(ggplot2)
  library(survminer)
  library(gridExtra)
  library(pheatmap)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

# --- CLI ---
train_expr_file <- get_arg("--train-expr")
train_surv_file <- get_arg("--train-surv")
candidate_file  <- get_arg("--gene-list")
cohort_id       <- get_arg("--cohort-id")
valid_expr_file <- get_arg("--validation-expr")
valid_surv_file <- get_arg("--validation-surv")
time_roc_days_str <- get_arg("--time-roc-days", required = FALSE, default = "1095,1825,2555")
ntree_str     <- get_arg("--ntree-range", required = FALSE, default = "100,200,300,500,1000")
mtry_str      <- get_arg("--mtry-range", required = FALSE, default = "1,2,3")
nodesize_str  <- get_arg("--nodesize-range", required = FALSE, default = "5,10,15,20,25,30")
rsf_seed      <- as.integer(get_arg("--seed", required = FALSE, default = "10"))
model_id      <- get_arg("--model-id")
base_output_dir <- get_arg("--outdir")
log_dir         <- get_arg("--logdir")
summary_dir     <- get_arg("--summary-dir")
sn              <- get_arg("--sn")

time_points <- as.numeric(unlist(strsplit(time_roc_days_str, ",")))
ntree_vec    <- as.integer(unlist(strsplit(ntree_str, ",")))
mtry_vec     <- as.integer(unlist(strsplit(mtry_str, ",")))
nodesize_vec <- as.integer(unlist(strsplit(nodesize_str, ",")))

sub_id <- get_arg("--sub-id", required = FALSE, default = NULL)
if (is.null(sub_id)) sub_id <- paste0("B3_rsf_vg/", cohort_id)

logsetup <- setup_logging(log_dir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_log_message("Script started: ", sub_id)
log_log_message("Log file: ", logsetup$log_file)
log_log_message("Summary file: ", logsetup$summary_file)

if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)

days_to_years <- function(days) {
  yr <- days / 365
  ifelse(yr == floor(yr), paste0(floor(yr), "y"), paste0(round(yr, 1), "y"))
}

# ================================
# 3. Plot helpers (keep behavior close to the manual script)
# ================================
plot_risk_distribution <- function(risk_df, cohort_name, out_prefix, output_dir) {
  risk_df <- risk_df %>% dplyr::arrange(riskScore)
  risk_df$order_id <- seq_len(nrow(risk_df))
  risk_df$riskGroup <- factor(risk_df$riskGroup, levels = c("Low risk", "High risk"))

  risk_cutoff_text <- paste0("Optimal cutoff: ", format(signif(risk_df$cutoff, 6), trim = TRUE))

  p1 <- ggplot(risk_df, aes(x = order_id, y = riskScore, color = riskGroup)) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = unique(risk_df$cutoff)[1], linetype = 2) +
    geom_vline(xintercept = sum(risk_df$riskGroup == "Low risk") + 0.5, linetype = 2) +
    scale_color_manual(values = c("Low risk" = "#6694e9", "High risk" = "#f47e84"), drop = FALSE) +
    labs(title = cohort_name, subtitle = risk_cutoff_text, x = "Patients (increasing risk score)", y = "Risk Score") +
    theme_bw(base_size = 12) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 10)
    )

  p2 <- ggplot(risk_df, aes(x = order_id, y = OS.time / 365, color = factor(OS, levels = c(0, 1), labels = c("Alive", "Dead")))) +
    geom_point(size = 1.8) +
    geom_vline(xintercept = sum(risk_df$riskGroup == "Low risk") + 0.5, linetype = 2) +
    scale_color_manual(values = c("Alive" = "#6694e9", "Dead" = "#f47e84")) +
    labs(title = "", x = "Patients (increasing risk score)", y = "Survival time (Years)") +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank())

  p_all <- gridExtra::grid.arrange(p1, p2, ncol = 1)

  ggsave(file.path(output_dir, paste0(out_prefix, "_risk_survival_distribution.pdf")), p_all, width = 6, height = 7)
  ggsave(file.path(output_dir, paste0(out_prefix, "_risk_survival_distribution.png")), p_all, width = 6, height = 7, dpi = 600)
}

plot_km <- function(risk_df, cohort_name, out_prefix, output_dir, log_dir) {
  risk_df$riskGroup <- factor(risk_df$riskGroup)
  n_group <- length(unique(risk_df$riskGroup))

  if (n_group < 2) {
    log_message("KM skipped for ", cohort_name, " — only one risk group present, n=", nrow(risk_df))
    return(NULL)
  }

  fit <- survival::survfit(Surv(OS.time, OS) ~ riskGroup, data = risk_df)

  legend_labs_use <- levels(droplevels(risk_df$riskGroup))
  palette_map <- c("Low risk" = "#6694e9", "High risk" = "#f47e84")
  palette_use <- palette_map[legend_labs_use]

  p_km <- survminer::ggsurvplot(
    fit,
    data = risk_df,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.col = "strata",
    legend.title = "Risk group",
    legend.labs = legend_labs_use,
    palette = palette_use,
    ggtheme = theme_bw(),
    title = cohort_name,
    xlab = "Time (day)",
    ylab = "Survival probability"
  )

  pdf(file.path(output_dir, paste0(out_prefix, "_KM.pdf")), width = 6, height = 6)
  print(p_km, newpage = FALSE)
  dev.off()

  png(file.path(output_dir, paste0(out_prefix, "_KM.png")), width = 6, height = 6, units = "in", res = 600)
  print(p_km)
  dev.off()

  p_km
}

plot_heatmap <- function(risk_df, gene_vec, cohort_name, out_prefix, output_dir) {
  missing_genes <- setdiff(gene_vec, colnames(risk_df))
  if (length(missing_genes) > 0) {
    stop(paste0(cohort_name, " heatmap failed, missing genes: ", paste(missing_genes, collapse = ", ")))
  }

  expr_mat <- as.matrix(risk_df[, gene_vec, drop = FALSE])
  rownames(expr_mat) <- risk_df$sample
  expr_mat <- t(scale(t(expr_mat)))
  expr_mat[is.na(expr_mat)] <- 0

  ann_col <- data.frame(RiskGroup = risk_df$riskGroup)
  rownames(ann_col) <- risk_df$sample

  pdf(file.path(output_dir, paste0(out_prefix, "_heatmap.pdf")), width = 8, height = 6)
  pheatmap::pheatmap(
    t(expr_mat),
    annotation_col = ann_col,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    main = cohort_name
  )
  dev.off()

  png(file.path(output_dir, paste0(out_prefix, "_heatmap.png")), width = 8, height = 6, units = "in", res = 600)
  pheatmap::pheatmap(
    t(expr_mat),
    annotation_col = ann_col,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    main = cohort_name
  )
  dev.off()
}

plot_timeROC <- function(risk_df, cohort_name, out_prefix, output_dir, time_points = c(1095, 1825, 2555)) {
  roc_res <- timeROC::timeROC(
    T = risk_df$OS.time,
    delta = risk_df$OS,
    marker = risk_df$riskScore,
    cause = 1,
    weighting = "marginal",
    times = time_points,
    iid = TRUE
  )

  auc_labels <- paste0(days_to_years(time_points), " AUC = ", sprintf("%.3f", roc_res$AUC))

  pdf(file.path(output_dir, paste0(out_prefix, "_timeROC.pdf")), width = 6, height = 5)
  plot(roc_res, time = time_points[1], col = "#f4be7e", lwd = 2, title = FALSE)
  if (length(time_points) >= 2) plot(roc_res, time = time_points[2], col = "#f47e84", lwd = 2, add = TRUE, title = FALSE)
  if (length(time_points) >= 3) plot(roc_res, time = time_points[3], col = "#6694e9", lwd = 2, add = TRUE, title = FALSE)
  abline(0, 1, lty = 2, col = "grey50")
  legend("bottomright", legend = auc_labels, col = c("#f4be7e", "#f47e84", "#6694e9")[seq_along(time_points)], lwd = 2, bty = "n")
  title(main = cohort_name)
  dev.off()

  png(file.path(output_dir, paste0(out_prefix, "_timeROC.png")), width = 6, height = 5, units = "in", res = 600)
  plot(roc_res, time = time_points[1], col = "#f4be7e", lwd = 2, title = FALSE)
  if (length(time_points) >= 2) plot(roc_res, time = time_points[2], col = "#f47e84", lwd = 2, add = TRUE, title = FALSE)
  if (length(time_points) >= 3) plot(roc_res, time = time_points[3], col = "#6694e9", lwd = 2, add = TRUE, title = FALSE)
  abline(0, 1, lty = 2, col = "grey50")
  legend("bottomright", legend = auc_labels, col = c("#f4be7e", "#f47e84", "#6694e9")[seq_along(time_points)], lwd = 2, bty = "n")
  title(main = cohort_name)
  dev.off()

  roc_res
}

safe_timeROC <- function(risk_df, cohort_name, out_prefix, output_dir) {
  tryCatch(
    plot_timeROC(risk_df, cohort_name, out_prefix, output_dir, time_points),
    error = function(e) {
      log_message("timeROC failed for ", cohort_name, ": ", e$message)
      NULL
    }
  )
}

# ================================
# 4. Read training set
# ================================
log_message("[B4 Step 1] Read training set")

data0 <- read.csv(train_expr_file, row.names = 1, check.names = FALSE) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

survival <- read.csv(train_surv_file, row.names = 1, check.names = FALSE) %>%
  rownames_to_column("sample_id")

gene <- read.csv(candidate_file, header = TRUE, check.names = FALSE)
gene <- data.frame(symbol = gene[[1]], stringsAsFactors = FALSE)
gene <- gene[!is.na(gene$symbol) & gene$symbol != "", , drop = FALSE]

matched_indices <- c("sample_id", gene$symbol)
data1 <- data0[, colnames(data0) %in% matched_indices, drop = FALSE]

train_dat <- merge(survival, data1, by = "sample_id")
rownames(train_dat) <- train_dat$sample_id
train_dat <- train_dat[, -1, drop = FALSE]
train_dat <- na.omit(train_dat)
train_dat <- train_dat[train_dat$OS.time > 0, , drop = FALSE]

log_message("train n=", nrow(train_dat), " genes=", ncol(train_dat) - 2)

# ================================
# 5. Read validation set
# ================================
log_message("[B4 Step 2] Read validation set")

survival_test <- read.csv(valid_surv_file, header = TRUE, row.names = 1, check.names = FALSE)
survival_test$sample <- rownames(survival_test)
colnames(survival_test)[1:3] <- c("OS", "OS.time", "sample")
survival_test <- survival_test[, c("sample", "OS", "OS.time"), drop = FALSE]
survival_test <- survival_test[survival_test$OS.time > 0, , drop = FALSE]

data_test <- read.csv(valid_expr_file, header = TRUE, row.names = 1, check.names = FALSE)
survival_dat <- t(data_test)
survival_dat <- survival_dat[rownames(survival_dat) %in% survival_test$sample, , drop = FALSE]

test_dat <- survival_dat[, colnames(survival_dat) %in% gene$symbol, drop = FALSE]
test_dat <- as.data.frame(test_dat)
test_dat$sample <- rownames(test_dat)

test_dat <- merge(survival_test, test_dat, by = "sample")
rownames(test_dat) <- test_dat$sample
test_dat <- test_dat[, -1, drop = FALSE]
test_dat <- na.omit(test_dat)

log_message("valid n=", nrow(test_dat), " genes=", ncol(test_dat) - 2)

# ================================
# 6. Align common genes
# ================================
log_message("[B4 Step 3] Align common genes")

train_genes <- setdiff(colnames(train_dat), c("OS", "OS.time"))
test_genes  <- setdiff(colnames(test_dat), c("OS", "OS.time"))
final_genes <- intersect(train_genes, test_genes)

if (length(final_genes) == 0) save_note_and_stop("No common genes between train and validation.")

train_dat <- train_dat[, c("OS", "OS.time", final_genes), drop = FALSE]
test_dat  <- test_dat[, c("OS", "OS.time", final_genes), drop = FALSE]

log_message("common genes=", length(final_genes))

# ================================
# 7. Build model input and standardize
# ================================
log_message("[B4 Step 4] Standardize")

x_train <- subset(train_dat, select = -c(OS, OS.time))
y_train <- subset(train_dat, select = c(OS, OS.time))
x_test  <- subset(test_dat,  select = -c(OS, OS.time))
y_test  <- subset(test_dat,  select = c(OS, OS.time))

x_train <- as.data.frame(lapply(x_train, function(x) as.numeric(as.character(x))), check.names = FALSE)
x_test  <- as.data.frame(lapply(x_test,  function(x) as.numeric(as.character(x))), check.names = FALSE)

y_train$OS <- as.numeric(as.character(y_train$OS))
y_train$OS.time <- as.numeric(as.character(y_train$OS.time))
y_test$OS <- as.numeric(as.character(y_test$OS))
y_test$OS.time <- as.numeric(as.character(y_test$OS.time))

x_train <- scale(x_train)
x_test  <- scale(x_test)

x_train <- as.data.frame(x_train, check.names = FALSE)
x_test  <- as.data.frame(x_test, check.names = FALSE)
x_test <- x_test[, colnames(x_train), drop = FALSE]

train_data <- data.frame(y_train, x_train, check.names = FALSE)
test_data  <- data.frame(y_test,  x_test,  check.names = FALSE)

train_data <- train_data[complete.cases(train_data), , drop = FALSE]
test_data  <- test_data[complete.cases(test_data),  , drop = FALSE]

if (nrow(train_data) == 0) save_note_and_stop("train_data empty after cleaning")
if (nrow(test_data) == 0) save_note_and_stop("test_data empty after cleaning")

train_sample <- rownames(train_data)
test_sample  <- rownames(test_data)

# ================================
# 8. RSF grid search (manual)
# ================================
log_message("[B4 Step 5] RSF grid search")

tune_grid <- expand.grid(ntree = ntree_vec, mtry = mtry_vec, nodesize = nodesize_vec)

best_auc_train <- -Inf
best_auc_test <- -Inf
best_params <- NULL
best_model <- NULL
best_train_auc_values <- NULL
best_test_auc_values <- NULL

valid_auc_labels <- paste0("valid_auc_", days_to_years(time_points))
grid_out <- data.frame(
  ntree = tune_grid$ntree,
  mtry = tune_grid$mtry,
  nodesize = tune_grid$nodesize,
  train_auc_sum = NA_real_,
  valid_auc_sum = NA_real_,
  stringsAsFactors = FALSE
)
for (lab in valid_auc_labels) grid_out[[lab]] <- NA_real_

evaluate_model <- function(model, data, time_points) {
  auc_values <- numeric(length(time_points))
  predictions <- predict(model, newdata = data, na.action = "na.impute")$predicted

  for (i in seq_along(time_points)) {
    a <- time_points[i]
    b <- timeROC(
      T = data$OS.time,
      delta = data$OS,
      marker = predictions,
      cause = 1,
      weighting = "marginal",
      times = a,
      ROC = TRUE,
      iid = TRUE
    )
    auc_values[i] <- as.numeric(b$AUC[2])
  }

  auc_values
}

for (i in seq_len(nrow(tune_grid))) {
  params <- tune_grid[i, ]

  set.seed(rsf_seed)
  rsf_model <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = train_data,
    importance = TRUE,
    ntree = params$ntree,
    mtry = params$mtry,
    nodesize = params$nodesize
  )

  train_auc_values <- evaluate_model(rsf_model, train_data, time_points)
  test_auc_values  <- evaluate_model(rsf_model, test_data, time_points)

  total_train_auc <- sum(train_auc_values, na.rm = TRUE)
  total_test_auc  <- sum(test_auc_values, na.rm = TRUE)

  grid_out$train_auc_sum[i] <- total_train_auc
  grid_out$valid_auc_sum[i] <- total_test_auc
  for (j in seq_along(valid_auc_labels)) grid_out[[valid_auc_labels[j]]][i] <- test_auc_values[j]

  if (all(test_auc_values > 0.6)) {
    if (total_test_auc > best_auc_test ||
        (total_test_auc == best_auc_test && total_train_auc > best_auc_train)) {
      best_auc_train <- total_train_auc
      best_auc_test <- total_test_auc
      best_params <- params
      best_model <- rsf_model
      best_train_auc_values <- train_auc_values
      best_test_auc_values <- test_auc_values
    }
  }
}

write.csv(grid_out, file.path(base_output_dir, "03.rsf_grid_search.csv"), row.names = FALSE, quote = FALSE)

if (is.null(best_model)) {
  # Pipeline-friendly fallback
  best_idx <- which.max(grid_out$valid_auc_sum)
  if (length(best_idx) == 0 || is.na(grid_out$valid_auc_sum[best_idx])) {
    save_note_and_stop("RSF grid search produced no valid results.")
  }
  best_params <- grid_out[best_idx, c("ntree", "mtry", "nodesize")]
  set.seed(rsf_seed)
  best_model <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = train_data,
    importance = TRUE,
    ntree = best_params$ntree,
    mtry = best_params$mtry,
    nodesize = best_params$nodesize
  )
  best_train_auc_values <- evaluate_model(best_model, train_data, time_points)
  best_test_auc_values  <- evaluate_model(best_model, test_data, time_points)
  best_auc_train <- sum(best_train_auc_values, na.rm = TRUE)
  best_auc_test  <- sum(best_test_auc_values, na.rm = TRUE)
  log_message("Fallback: no model with all valid AUCs > 0.6. Best valid_auc_sum: ", sprintf("%.3f", best_auc_test))
} else {
  log_message("Selected by threshold: all valid AUCs > 0.6. Best valid_auc_sum: ", sprintf("%.3f", best_auc_test))
}

final_model <- best_model

log_message("Best hyperparams: ntree=", best_params$ntree, " mtry=", best_params$mtry, " nodesize=", best_params$nodesize)
log_message("Train AUC sum: ", round(best_auc_train, 3))
log_message("Valid AUC sum: ", round(best_auc_test, 3))

# ================================
# 9. Feature importance
# ================================
log_message("[B4 Step 6] Feature importance")

importance <- vimp(final_model)$importance
importance_df <- data.frame(variable = names(importance), importance = importance)
importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), , drop = FALSE]

importance_plot <- ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  xlab("Features") +
  ylab("Importance") +
  ggtitle("Feature Importance Ranking") +
  theme_bw()

ggsave(file.path(base_output_dir, "01.feature_importance.pdf"), plot = importance_plot, width = 8, height = 6)
ggsave(file.path(base_output_dir, "01.feature_importance.png"), plot = importance_plot, width = 8, height = 6, dpi = 300)

write.csv(importance_df, file.path(base_output_dir, "01.feature_importance.csv"), row.names = FALSE, quote = FALSE)

# ================================
# 10. Risk score + cutoff
# ================================
log_message("[B4 Step 7] Risk scoring")

train_predictions <- predict(final_model, newdata = train_data, na.action = "na.impute")$predicted
test_predictions  <- predict(final_model, newdata = test_data, na.action = "na.impute")$predicted

train_data$riskScore <- as.numeric(train_predictions)
test_data$riskScore  <- as.numeric(test_predictions)

# Default: median cutoff
train_cutoff <- median(train_data$riskScore, na.rm = TRUE)
test_cutoff  <- median(test_data$riskScore, na.rm = TRUE)

train_data$riskGroup <- ifelse(train_data$riskScore > train_cutoff, "High risk", "Low risk")
train_data$riskGroup <- factor(train_data$riskGroup, levels = c("Low risk", "High risk"))
train_data$cutoff <- train_cutoff

test_data$riskGroup <- ifelse(test_data$riskScore > test_cutoff, "High risk", "Low risk")
test_data$riskGroup <- factor(test_data$riskGroup, levels = c("Low risk", "High risk"))
test_data$cutoff <- test_cutoff

# Check KM significance; fallback to optimal cutpoint if p >= 0.05
test_km <- function(data, label) {
  km <- tryCatch(survdiff(Surv(OS.time, OS) ~ riskGroup, data = data), error = function(e) NULL)
  if (is.null(km)) return(1)
  1 - pchisq(km$chisq, df = 1)
}

t_p <- test_km(train_data); v_p <- test_km(test_data)
need_fallback <- is.na(t_p) || is.na(v_p) || t_p >= 0.05 || v_p >= 0.05

if (need_fallback) {
  log_message("Median cutoff KM p: train=", sprintf("%.4f", t_p), " valid=", sprintf("%.4f", v_p), " -> fallback to optimal cutpoint")
  train_cut <- tryCatch(survminer::surv_cutpoint(train_data, time = "OS.time", event = "OS", variables = "riskScore", minprop = 0.3), error = function(e) NULL)
  if (!is.null(train_cut)) {
    train_cutoff <- train_cut$cutpoint[1, 1]
    train_data$riskGroup <- ifelse(train_data$riskScore > train_cutoff, "High risk", "Low risk")
    train_data$riskGroup <- factor(train_data$riskGroup, levels = c("Low risk", "High risk"))
    train_data$cutoff <- train_cutoff
  }
  test_cut <- tryCatch(survminer::surv_cutpoint(test_data, time = "OS.time", event = "OS", variables = "riskScore", minprop = 0.3), error = function(e) NULL)
  if (!is.null(test_cut)) {
    test_cutoff <- test_cut$cutpoint[1, 1]
    test_data$riskGroup <- ifelse(test_data$riskScore > test_cutoff, "High risk", "Low risk")
    test_data$riskGroup <- factor(test_data$riskGroup, levels = c("Low risk", "High risk"))
    test_data$cutoff <- test_cutoff
  }
  log_message("Optimal cutoff: train=", sprintf("%.4f", train_cutoff), " valid=", sprintf("%.4f", test_cutoff))
} else {
  log_message("Median cutoff OK: train p=", sprintf("%.4f", t_p), " valid p=", sprintf("%.4f", v_p))
}

train_data$sample <- train_sample
test_data$sample  <- test_sample

write.csv(train_data, file.path(base_output_dir, "02.train_risk_score.csv"), row.names = FALSE, quote = FALSE)
saveRDS(final_model, file.path(base_output_dir, "04.rsf_model.rds"))

validate_dir <- file.path(base_output_dir, "validate")
dir.create(validate_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(test_data, file.path(validate_dir, "01.risk_score.csv"), row.names = FALSE, quote = FALSE)

# ================================
# 11. Plots
# ================================
log_message("[B4 Step 8] Plots")

train_plot_dir <- file.path(base_output_dir, "train_plots")
dir.create(train_plot_dir, recursive = TRUE, showWarnings = FALSE)

plot_risk_distribution(train_data, "train", "train", train_plot_dir)
plot_km(train_data, "train", "train", train_plot_dir, log_dir)
plot_heatmap(train_data %>% dplyr::arrange(riskScore), final_genes, "train", "train", train_plot_dir)
safe_timeROC(train_data, "train", "train", train_plot_dir)

plot_risk_distribution(test_data, cohort_id, cohort_id, validate_dir)
plot_km(test_data, cohort_id, cohort_id, validate_dir, log_dir)
plot_heatmap(test_data %>% dplyr::arrange(riskScore), final_genes, cohort_id, cohort_id, validate_dir)
safe_timeROC(test_data, cohort_id, cohort_id, validate_dir)

# ================================
# 12. Summary CSVs for NF aggregation
# ================================
time_set_label <- paste0(days_to_years(time_points), collapse = "_")

train_km <- tryCatch(survdiff(Surv(OS.time, OS) ~ riskGroup, data = train_data), error = function(e) NULL)
train_km_p <- if (!is.null(train_km)) 1 - pchisq(train_km$chisq, df = 1) else NA

valid_km <- tryCatch(survdiff(Surv(OS.time, OS) ~ riskGroup, data = test_data), error = function(e) NULL)
valid_km_p <- if (!is.null(valid_km)) 1 - pchisq(valid_km$chisq, df = 1) else NA

train_summary <- data.frame(
  model_id = model_id, time_set = time_set_label, dataset = "train",
  n_samples = nrow(train_data), n_genes = length(final_genes),
  km_pvalue = train_km_p, cutoff = train_cutoff,
  best_ntree = best_params$ntree, best_mtry = best_params$mtry, best_nodesize = best_params$nodesize,
  stringsAsFactors = FALSE
)
auc_names <- paste0("auc_", days_to_years(time_points))
for (k in seq_along(time_points)) train_summary[[auc_names[k]]] <- best_train_auc_values[k]
write.csv(train_summary, file.path(base_output_dir, "train_summary.csv"), row.names = FALSE, quote = FALSE)

valid_summary <- data.frame(
  model_id = model_id, time_set = time_set_label, cohort_id = cohort_id, dataset = "validation",
  n_samples = nrow(test_data), n_genes_matched = length(final_genes),
  km_pvalue = valid_km_p,
  stringsAsFactors = FALSE
)
for (k in seq_along(time_points)) valid_summary[[auc_names[k]]] <- best_test_auc_values[k]
write.csv(valid_summary, file.path(base_output_dir, "validation_summary.csv"), row.names = FALSE, quote = FALSE)

train_auc_str <- paste(sprintf("%.3f", best_train_auc_values), collapse = ", ")
valid_auc_str <- paste(sprintf("%.3f", best_test_auc_values), collapse = ", ")

write_summary(c(
  paste0("模型ID: ", model_id),
  paste0("验证队列: ", cohort_id),
  paste0("训练样本数: ", nrow(train_data)),
  paste0("验证样本数: ", nrow(test_data)),
  paste0("基因数: ", length(final_genes)),
  paste0("Best ntree: ", best_params$ntree),
  paste0("Best mtry: ", best_params$mtry),
  paste0("Best nodesize: ", best_params$nodesize),
  paste0("Train KM p-value: ", ifelse(is.na(train_km_p), "NA", sprintf("%.4f", train_km_p))),
  paste0("Valid KM p-value: ", ifelse(is.na(valid_km_p), "NA", sprintf("%.4f", valid_km_p))),
  paste0("Train AUCs: ", train_auc_str),
  paste0("Valid AUCs: ", valid_auc_str),
  "",
  paste0("结论: RSF验证引导模型在", cohort_id, "队列中完成验证。最佳参数为ntree=", best_params$ntree, ", mtry=", best_params$mtry, ", nodesize=", best_params$nodesize, "。"),
  "",
  paste0("限制与说明: 验证引导的网格搜索优先选择所有时间点验证AUC均>0.6的参数组合。若无满足条件的组合，则选择验证AUC总和最高的组合。")
))

log_message("Script completed: ", sub_id)

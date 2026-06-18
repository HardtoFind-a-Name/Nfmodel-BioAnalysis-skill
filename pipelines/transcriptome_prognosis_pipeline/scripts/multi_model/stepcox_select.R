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
  library(MASS)
  library(ggplot2)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

expr_file <- get_arg("--expr")
surv_file <- get_arg("--surv")
gene_list_file <- get_arg("--gene-list")
model_id <- get_arg("--model-id")
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")

sub_id <- get_arg("--sub-id", required = FALSE, default = NULL)
if (is.null(sub_id)) sub_id <- paste0(model_id, "/stepcox")

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

candidate_genes <- read_gene_list(gene_list_file)
if (length(candidate_genes) < 2) save_note_and_stop("stepCox requires at least 2 genes.")

expr_list <- prepare_expr_matrix(expr_file, candidate_genes, "train")
surv_df <- prepare_surv_data(surv_file)
train_data <- merge_expr_surv(expr_list, surv_df)
cleaned <- drop_invalid_gene_cols(train_data, setdiff(colnames(train_data), c("sample", "OS.time", "OS")))
train_data <- cleaned$data
gene_cols <- cleaned$genes
if (length(gene_cols) < 2) save_note_and_stop("After cleaning, fewer than 2 genes for stepCox.")

log_message("Input genes: ", length(gene_cols), ", samples: ", nrow(train_data))

formula_str <- paste0("Surv(OS.time, OS) ~ ", paste0("`", gene_cols, "`", collapse = " + "))
full_fit <- tryCatch(
  survival::coxph(as.formula(formula_str), data = train_data),
  error = function(e) { save_note_and_stop(paste0("stepCox full fit failed: ", e$message)) }
)

step_fit <- tryCatch(
  MASS::stepAIC(full_fit, direction = "both", trace = FALSE),
  error = function(e) { save_note_and_stop(paste0("stepAIC failed: ", e$message)) }
)

fit_sum <- summary(step_fit)
step_genes <- names(coef(step_fit))
step_genes <- gsub("`", "", step_genes)

if (length(step_genes) == 0) save_note_and_stop("stepCox selected 0 genes.")

res_df <- data.frame(
  gene = step_genes,
  coef = as.numeric(coef(step_fit)),
  HR = fit_sum$conf.int[, "exp(coef)"],
  HR.95L = fit_sum$conf.int[, "lower .95"],
  HR.95H = fit_sum$conf.int[, "upper .95"],
  pvalue = fit_sum$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
) %>% dplyr::arrange(pvalue)

coef_df <- data.frame(gene = res_df$gene, coef = res_df$coef, stringsAsFactors = FALSE)
write.csv(coef_df, file.path(outdir, "01.stepcox_final_genes_coef.csv"), row.names = FALSE, quote = FALSE)
write.csv(res_df, file.path(outdir, "02.stepcox_all_results.csv"), row.names = FALSE, quote = FALSE)

zph <- tryCatch(survival::cox.zph(step_fit), error = function(e) NULL)
if (!is.null(zph)) {
  ph_df <- as.data.frame(zph$table)
  ph_df$term <- rownames(ph_df)
  rownames(ph_df) <- NULL
  write.csv(ph_df, file.path(outdir, "03.stepcox_ph_test.csv"), row.names = FALSE, quote = FALSE)
}

aic_before <- AIC(full_fit)
aic_after <- AIC(step_fit)

log_message("Stepwise selection complete — AIC before: ", round(aic_before, 2), ", after: ", round(aic_after, 2), ", selected: ", length(step_genes), " genes")

if (nrow(res_df) > 0) {
  plot_cox_forest(res_df, outdir, "04.stepcox_forest", paste0(model_id, " stepCox"))
}

write_summary(c(
  paste0("模型ID: ", model_id),
  paste0("输入基因数: ", length(gene_cols)),
  paste0("样本数: ", nrow(train_data)),
  paste0("Stepwise AIC (before): ", round(aic_before, 2)),
  paste0("Stepwise AIC (after): ", round(aic_after, 2)),
  paste0("步进法选择基因数: ", length(step_genes)),
  paste0("Selected genes: ", paste(step_genes, collapse = ", ")),
  "",
  paste0("结论: 逐步回归 (stepAIC, both方向) 从", length(gene_cols), "个基因中筛选出", length(step_genes), "个基因，AIC从", round(aic_before, 2), "降至", round(aic_after, 2), "。"),
  "",
  paste0("限制与说明: 逐步回归基于AIC准则，可能产生过拟合。建议结合临床意义和LASSO/uniCox结果综合判断。")
))

log_message("Script completed: ", sub_id)

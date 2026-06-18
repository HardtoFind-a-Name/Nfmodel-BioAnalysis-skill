rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
})

# CLI parsing
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

# Source utils
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

# Parse args
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")
expr_file   <- get_arg("--expr")
surv_file   <- get_arg("--surv")
candidate_file <- get_arg("--candidate-genes")
unicox_p_cutoff <- as.numeric(get_arg("--unicox-p", required = FALSE, default = "0.05"))
ph_p_cutoff <- as.numeric(get_arg("--ph-p", required = FALSE, default = "0.05"))

sub_id <- get_arg("--sub-id", required = FALSE, default = "00_unicox_screen")

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_file)) save_note_and_stop(paste0("Expression file not found: ", expr_file))
if (!file.exists(surv_file)) save_note_and_stop(paste0("Survival file not found: ", surv_file))
if (!file.exists(candidate_file)) save_note_and_stop(paste0("Candidate gene file not found: ", candidate_file))

candidate_genes <- read_gene_list(candidate_file)
if (length(candidate_genes) == 0) save_note_and_stop("Candidate gene file is empty.")

log_message("Input genes from file: ", length(candidate_genes))

expr_list <- prepare_expr_matrix(expr_file, candidate_genes, "train")
surv_df <- prepare_surv_data(surv_file)
train_data <- merge_expr_surv(expr_list, surv_df)
gene_cols <- setdiff(colnames(train_data), c("sample", "OS.time", "OS"))
cleaned <- drop_invalid_gene_cols(train_data, gene_cols)
train_data <- cleaned$data
gene_cols <- cleaned$genes

if (length(gene_cols) == 0) save_note_and_stop("No valid genes for uniCox.")

log_message("Input: ", basename(expr_file), " — genes: ", length(gene_cols), ", samples: ", nrow(train_data))

write.csv(train_data, file.path(outdir, "00.unicox_input_merged.csv"), row.names = FALSE, quote = FALSE)

get_ci_rule <- function(p_cutoff) {
  if (p_cutoff <= 0.05) list(ci_level = 0.95, ci_label = "95%")
  else if (p_cutoff <= 0.1) list(ci_level = 0.90, ci_label = "90%")
  else if (p_cutoff <= 0.2) list(ci_level = 0.80, ci_label = "80%")
  else stop("unicox_p_cutoff should not exceed 0.2.")
}
ci_rule <- get_ci_rule(unicox_p_cutoff)
uni_ci_level <- ci_rule$ci_level
uni_ci_label <- ci_rule$ci_label
alpha <- 1 - uni_ci_level
z_crit <- qnorm(1 - alpha / 2)

uni_res <- list()
ph_res <- list()

for (g in gene_cols) {
  cox_formula <- as.formula(paste0("Surv(OS.time, OS) ~ `", g, "`"))
  fit <- tryCatch(survival::coxph(cox_formula, data = train_data), error = function(e) NULL)
  if (is.null(fit)) next
  fit_sum <- summary(fit)
  if (nrow(fit_sum$coefficients) == 0) next
  coef_val <- fit_sum$coefficients[1, "coef"]
  se_val <- fit_sum$coefficients[1, "se(coef)"]
  uni_res[[g]] <- data.frame(
    gene = g,
    coef = coef_val,
    HR = exp(coef_val),
    HR.L = exp(coef_val - z_crit * se_val),
    HR.H = exp(coef_val + z_crit * se_val),
    z = fit_sum$coefficients[1, "z"],
    pvalue = fit_sum$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  if (!is.null(zph)) {
    ztab <- as.data.frame(zph$table)
    ztab$term <- rownames(ztab)
    ztab$gene <- g
    rownames(ztab) <- NULL
    ph_res[[g]] <- ztab
  }
}

if (length(uni_res) == 0) save_note_and_stop("uniCox produced no results.")

log_message("UniCox complete — genes tested: ", length(uni_res), ", passing p < ", unicox_p_cutoff, ": ", sum(uni_df$pvalue < unicox_p_cutoff))

uni_df <- dplyr::bind_rows(uni_res) %>% dplyr::arrange(pvalue)
ph_df <- if (length(ph_res) > 0) dplyr::bind_rows(ph_res) else data.frame()

write.csv(uni_df, file.path(outdir, "01.unicox_all_results.csv"), row.names = FALSE, quote = FALSE)
write.csv(ph_df, file.path(outdir, "02.unicox_ph_test.csv"), row.names = FALSE, quote = FALSE)

ph_gene_df <- ph_df %>%
  dplyr::filter(term != "GLOBAL") %>%
  dplyr::transmute(gene = gene, PH.chisq = chisq, PH.p = p)

pass_df <- uni_df %>%
  dplyr::filter(pvalue < unicox_p_cutoff) %>%
  dplyr::inner_join(ph_gene_df, by = "gene") %>%
  dplyr::filter(PH.p > ph_p_cutoff) %>%
  dplyr::arrange(pvalue)

write.csv(pass_df, file.path(outdir, "03.unicox_ph_pass_genes.csv"), row.names = FALSE, quote = FALSE)

forest_df <- uni_df %>% dplyr::filter(pvalue < unicox_p_cutoff)
if (nrow(forest_df) > 0) {
  plot_cox_forest(forest_df, outdir, "04.unicox_forest", paste0("Unicox (", uni_ci_label, " CI)"), ci_label = uni_ci_label)
}

log_message("PH test complete — genes passing PH (p > ", ph_p_cutoff, "): ", nrow(pass_df))
log_message("Forest plot: ", ifelse(nrow(forest_df) > 0, paste0(nrow(forest_df), " genes plotted"), "skipped (0 genes passed p cutoff)"))

write_summary(c(
  paste0("分析队列: ", basename(dirname(expr_file))),
  paste0("输入基因数 (candidate): ", length(candidate_genes)),
  paste0("cleaned后有效基因数: ", length(gene_cols)),
  paste0("样本数: ", nrow(train_data)),
  paste0("uniCox P 阈值: ", unicox_p_cutoff),
  paste0("PH 检验 P 阈值: ", ph_p_cutoff),
  paste0("uniCox 筛选后基因数 (p < ", unicox_p_cutoff, "): ", nrow(uni_df)),
  paste0("PH 检验通过基因数: ", nrow(pass_df)),
  paste0("最终 pass 基因数: ", nrow(pass_df)),
  "",
  paste0("结论: 共", length(gene_cols), "个候选基因经单因素Cox回归和PH检验筛选，在P<", unicox_p_cutoff, "及PH>", ph_p_cutoff, "的标准下，最终保留", nrow(pass_df), "个基因进入后续建模。筛选标准较为", ifelse(unicox_p_cutoff <= 0.05 & ph_p_cutoff <= 0.05, "严格", "宽松"), "。"),
  "",
  paste0("限制与说明: 仅保留在训练集中满足uniCox和PH显著性标准的基因。未考虑基因间的相互作用。PH检验基于 Schoenfeld 残差。")
))

log_message("Script completed: ", sub_id)

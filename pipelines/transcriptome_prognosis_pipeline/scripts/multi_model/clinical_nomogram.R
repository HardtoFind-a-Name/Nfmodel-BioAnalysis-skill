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
  library(survival)
  library(survminer)
  library(rms)
  library(regplot)
  library(forestplot)
  library(timeROC)
  library(rmda)
  library(dplyr)
  library(grid)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

surv_file <- get_arg("--surv")
clinical_file <- get_arg("--clinical")
risk_file <- get_arg("--risk")
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")
time_roc_days <- as.numeric(strsplit(get_arg("--time-roc-days"), ",", fixed = TRUE)[[1]])
calibration_bootstrap <- as.integer(get_arg("--calibration-bootstrap", required = FALSE, default = "1000"))

sub_id <- get_arg("--sub-id", required = FALSE, default = "clinical_nomogram")

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
save_note_and_stop <- logsetup$save_note_and_stop

# Keep backward-compatible write_summary for existing single-line calls
write_summary <- function(...) {
  logsetup$write_summary(paste0(..., collapse = ""))
}

dataset_summary_file <- file.path(logdir, paste0(sn, "_multi_model"), sub_id, "06.dataset_summary.txt")
ensure_dir(dirname(dataset_summary_file))

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

plot_to_pdf_png <- function(plot_fun, pdf_file, png_file, width, height, png_res = 600, pdf_family = "Times") {
  pdf(pdf_file, width = width, height = height, family = pdf_family)
  plot_fun()
  dev.off()

  png(png_file, width = width, height = height, units = "in", res = png_res)
  plot_fun()
  dev.off()
}

normalize_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "null", "NULL", "not reported", "Not Reported",
             "unknown", "Unknown", "[Not Available]")] <- NA
  x
}

extract_feature_from_label <- function(label) {
  if (grepl("^riskScore$", label)) return("riskScore")
  if (grepl("^age$", label)) return("age")
  if (grepl("vs female$", label)) return("gender")
  if (grepl("^Race \\(Reference: white\\)$", label)) return("race")
  if (label %in% c("black", "asian")) return("race")
  if (grepl("^Pathologic stage \\(Reference: Stage I\\)$", label)) return("pathologic.stage")
  if (label %in% c("Stage II", "Stage III", "Stage IV")) return("pathologic.stage")
  NA_character_
}

build_display_table <- function(result_df, feature_vec, overall_df, prefix = "univariate") {
  display_df <- data.frame()

  for (feature_name in feature_vec) {
    tmp_df <- result_df[grepl(paste0("^", gsub("\\.", "\\\\.", feature_name)), result_df$variable), , drop = FALSE]
    if (nrow(tmp_df) == 0) next

    if (feature_name %in% c("riskScore", "age")) {
      tmp_df$display_label <- feature_name
      tmp_df$is_reference <- FALSE
      tmp_df$is_summary <- FALSE
      display_df <- rbind(display_df, tmp_df)
    }

    if (feature_name == "gender") {
      tmp_df$display_label <- sub("^gender", "", tmp_df$variable)
      tmp_df$display_label <- paste0(tmp_df$display_label, " vs female")
      tmp_df$is_reference <- FALSE
      tmp_df$is_summary <- FALSE
      display_df <- rbind(display_df, tmp_df)
    }

    if (feature_name == "race") {
      ref_row <- data.frame(
        variable = "race_reference",
        HR = 1,
        HR.95L = NA,
        HR.95H = NA,
        p.value = NA,
        feature = feature_name,
        display_label = "Race (Reference: white)",
        is_reference = TRUE,
        is_summary = TRUE,
        stringsAsFactors = FALSE
      )
      tmp_df$display_label <- sub("^race", "", tmp_df$variable)
      tmp_df$is_reference <- FALSE
      tmp_df$is_summary <- FALSE
      display_df <- dplyr::bind_rows(display_df, ref_row, tmp_df)
    }

    if (feature_name == "pathologic.stage") {
      ref_row <- data.frame(
        variable = "pathologic.stage_reference",
        HR = 1,
        HR.95L = NA,
        HR.95H = NA,
        p.value = NA,
        feature = feature_name,
        display_label = "Pathologic stage (Reference: Stage I)",
        is_reference = TRUE,
        is_summary = TRUE,
        stringsAsFactors = FALSE
      )
      tmp_df$display_label <- sub("^pathologic.stage", "", tmp_df$variable)
      tmp_df$is_reference <- FALSE
      tmp_df$is_summary <- FALSE
      display_df <- dplyr::bind_rows(display_df, ref_row, tmp_df)
    }
  }

  display_df$overall_p <- vapply(display_df$display_label, extract_feature_from_label, character(1))
  display_df$overall_p <- ifelse(
    is.na(display_df$overall_p),
    NA_real_,
    overall_df$overall_p[match(display_df$overall_p, overall_df$feature)]
  )
  display_df$display_p <- ifelse(display_df$is_summary, display_df$overall_p, display_df$p.value)

  for (i in seq_len(nrow(display_df))) {
    if (isTRUE(display_df$is_summary[i]) && !is.na(display_df$overall_p[i])) {
      p_text <- ifelse(display_df$overall_p[i] < 0.001, "<0.001", sprintf("%.4f", display_df$overall_p[i]))
      display_df$display_label[i] <- paste0(display_df$display_label[i], "; overall P = ", p_text)
    }
  }

  write.csv(
    display_df,
    file.path(outdir, paste0(ifelse(prefix == "univariate", "08", "12"), ".", prefix, "_cox_display.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  display_df
}

draw_forestplot <- function(display_df, title_text, prefix_num, prefix_name) {
  if (nrow(display_df) == 0) return(invisible(NULL))

  display_df$HR_CI <- ifelse(
    is.na(display_df$HR.95L) | is.na(display_df$HR.95H),
    "",
    paste0(
      sprintf("%.3f", display_df$HR), " (",
      sprintf("%.3f", display_df$HR.95L), "-",
      sprintf("%.3f", display_df$HR.95H), ")"
    )
  )

  tabletext <- cbind(
    c("Variable", display_df$display_label),
    c("P value", ifelse(is.na(display_df$display_p), "", ifelse(display_df$display_p < 0.001, "<0.001", sprintf("%.4f", display_df$display_p)))),
    c("Hazard Ratio (95% CI)", display_df$HR_CI)
  )

  plot_fun <- function() {
    forestplot(
      labeltext = tabletext,
      graph.pos = 3,
      is.summary = c(TRUE, display_df$is_summary),
      mean = c(NA, ifelse(display_df$is_reference, NA, display_df$HR)),
      lower = c(NA, display_df$HR.95L),
      upper = c(NA, display_df$HR.95H),
      zero = 1,
      boxsize = 0.18,
      lwd.ci = 2,
      ci.vertices = TRUE,
      ci.vertices.height = 0.08,
      colgap = unit(5, "mm"),
      lineheight = unit(0.9, "cm"),
      graphwidth = unit(0.35, "npc"),
      xlab = "Hazard Ratio",
      title = title_text,
      col = fpColors(box = "#d9b02c", lines = "#4981bb", zero = "gray50"),
      txt_gp = fpTxtGp(
        label = gpar(cex = 0.9),
        ticks = gpar(cex = 0.85),
        xlab = gpar(cex = 1),
        title = gpar(cex = 1.1),
        summary = gpar(cex = 0.95, fontface = "bold")
      )
    )
  }

  height_use <- max(6, 0.35 * nrow(display_df) + 2)
  pdf(file.path(outdir, paste0(prefix_num, ".", prefix_name, ".pdf")), height = height_use, width = 14, onefile = FALSE, family = "Times")
  plot_fun()
  dev.off()

  png(file.path(outdir, paste0(prefix_num, ".", prefix_name, ".png")), height = height_use, width = 14, units = "in", res = 600)
  plot_fun()
  dev.off()
}

cat("", file = log_file)
cat("", file = summary_file)

log_message("Script started: clinical_nomogram.R")
write_summary("Script: clinical_nomogram.R")
write_summary("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
write_summary("Survival file: ", surv_file)
write_summary("Clinical file: ", clinical_file)
write_summary("Risk file: ", risk_file)

if (!file.exists(surv_file)) save_note_and_stop(paste0("Survival file not found: ", surv_file))
if (!file.exists(clinical_file)) save_note_and_stop(paste0("Clinical file not found: ", clinical_file))
if (!file.exists(risk_file)) save_note_and_stop(paste0("Risk file not found: ", risk_file))

survival_df <- read.csv(surv_file, check.names = FALSE, stringsAsFactors = FALSE)
clinical_df <- read.csv(clinical_file, check.names = FALSE, stringsAsFactors = FALSE)
risk_df <- read.csv(risk_file, check.names = FALSE, stringsAsFactors = FALSE)

log_message("Input dimensions - survival: ", paste(dim(survival_df), collapse = " x "))
log_message("Input dimensions - clinical: ", paste(dim(clinical_df), collapse = " x "))
log_message("Input dimensions - risk: ", paste(dim(risk_df), collapse = " x "))

if (nrow(survival_df) == 0) save_note_and_stop("Survival table is empty.")
if (nrow(clinical_df) == 0) save_note_and_stop("Clinical table is empty.")
if (nrow(risk_df) == 0) save_note_and_stop("Risk table is empty.")

colnames(survival_df)[1] <- "sample"
if (!("sample" %in% colnames(survival_df))) save_note_and_stop("sample column not found in survival file.")
if (!("sample_id" %in% colnames(clinical_df))) save_note_and_stop("sample_id column not found in clinical file.")
clinical_df$sample <- clinical_df$sample_id
if (!("OS.time" %in% colnames(survival_df))) save_note_and_stop("OS.time column not found in survival file.")
if (!("OS" %in% colnames(survival_df))) save_note_and_stop("OS column not found in survival file.")
if (!("riskScore" %in% colnames(risk_df))) save_note_and_stop("riskScore column not found in risk file.")

risk_id_col <- intersect(c("sample", "id", "Sample", "ID"), colnames(risk_df))
if (length(risk_id_col) == 0) save_note_and_stop("No sample/id column found in risk file.")
colnames(risk_df)[colnames(risk_df) == risk_id_col[1]] <- "sample"

if (!("age" %in% colnames(clinical_df))) {
  age_alt <- intersect(c("Age", "age_at_index", "age_at_diagnosis", "age_at_index.demographic"), colnames(clinical_df))
  if (length(age_alt) > 0) clinical_df$age <- clinical_df[[age_alt[1]]]
}
if (!("gender" %in% colnames(clinical_df))) {
  gender_alt <- intersect(c("sex", "Sex", "gender.demographic"), colnames(clinical_df))
  if (length(gender_alt) > 0) clinical_df$gender <- clinical_df[[gender_alt[1]]]
}
if (!("race" %in% colnames(clinical_df))) {
  race_alt <- intersect(c("race", "race.demographic", "Race"), colnames(clinical_df))
  if (length(race_alt) > 0) clinical_df$race <- clinical_df[[race_alt[1]]]
}
if (!("pathologic.stage" %in% colnames(clinical_df))) {
  stage_alt <- intersect(c("ajcc_pathologic_stage.diagnoses", "stage", "pathologic_stage"), colnames(clinical_df))
  if (length(stage_alt) > 0) clinical_df$pathologic.stage <- clinical_df[[stage_alt[1]]]
}

need_cols <- c("sample", "age", "gender", "race", "pathologic.stage")
missing_need_cols <- setdiff(need_cols, colnames(clinical_df))
if (length(missing_need_cols) > 0) {
  save_note_and_stop(paste0("Required clinical columns missing: ", paste(missing_need_cols, collapse = ", ")))
}

clinical_use <- clinical_df[, need_cols, drop = FALSE]
survival_use <- survival_df[, c("sample", "OS.time", "OS"), drop = FALSE]
risk_use <- risk_df[, c("sample", "riskScore"), drop = FALSE]

write.csv(clinical_use, file.path(outdir, "01.clinical_selected.csv"), row.names = FALSE, quote = FALSE)
write.csv(survival_use, file.path(outdir, "02.survival_selected.csv"), row.names = FALSE, quote = FALSE)
write.csv(risk_use, file.path(outdir, "03.risk_selected.csv"), row.names = FALSE, quote = FALSE)

merge_sc <- merge(survival_use, clinical_use, by = "sample")
merge_all <- merge(merge_sc, risk_use, by = "sample")
write.csv(merge_all, file.path(outdir, "04.merged_raw.csv"), row.names = FALSE, quote = FALSE)

prog_df <- merge_all
prog_df$OS.time <- suppressWarnings(as.numeric(prog_df$OS.time))
prog_df$OS <- suppressWarnings(as.numeric(prog_df$OS))
prog_df$age <- suppressWarnings(as.numeric(trimws(as.character(prog_df$age))))
prog_df$age[prog_df$age < 0] <- NA
prog_df$riskScore <- suppressWarnings(as.numeric(prog_df$riskScore))

prog_df$gender <- tolower(normalize_missing(prog_df$gender))
prog_df$gender[prog_df$gender %in% c("female", "f")] <- "female"
prog_df$gender[prog_df$gender %in% c("male", "m")] <- "male"
prog_df$gender[!(prog_df$gender %in% c("female", "male"))] <- NA

prog_df$race <- tolower(normalize_missing(prog_df$race))
prog_df$race[prog_df$race %in% c("white")] <- "white"
prog_df$race[prog_df$race %in% c("black or african american", "black", "african american")] <- "black"
prog_df$race[prog_df$race %in% c("asian")] <- "asian"
prog_df$race[prog_df$race %in% c("american indian or alaska native")] <- NA
prog_df$race[!(prog_df$race %in% c("white", "black", "asian"))] <- NA

prog_df$pathologic.stage <- toupper(normalize_missing(prog_df$pathologic.stage))
prog_df$pathologic.stage <- gsub("\\s+", " ", prog_df$pathologic.stage)
prog_df$pathologic.stage[prog_df$pathologic.stage %in% c("STAGE X", "STAGE 0", "STAGE IS", "0")] <- NA
prog_df$pathologic.stage <- gsub("^STAGE IV[ABC]$", "STAGE IV", prog_df$pathologic.stage)
prog_df$pathologic.stage <- gsub("^STAGE III[ABC]$", "STAGE III", prog_df$pathologic.stage)
prog_df$pathologic.stage <- gsub("^STAGE II[ABC]$", "STAGE II", prog_df$pathologic.stage)
prog_df$pathologic.stage <- gsub("^STAGE I[ABC]$", "STAGE I", prog_df$pathologic.stage)
prog_df$pathologic.stage[prog_df$pathologic.stage %in% c("I", "STAGE I")] <- "Stage I"
prog_df$pathologic.stage[prog_df$pathologic.stage %in% c("II", "STAGE II")] <- "Stage II"
prog_df$pathologic.stage[prog_df$pathologic.stage %in% c("III", "STAGE III")] <- "Stage III"
prog_df$pathologic.stage[prog_df$pathologic.stage %in% c("IV", "STAGE IV")] <- "Stage IV"
prog_df$pathologic.stage[!(prog_df$pathologic.stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"))] <- NA

prog_df$gender <- factor(prog_df$gender, levels = c("female", "male"))
prog_df$race <- factor(prog_df$race, levels = c("white", "black", "asian"))
prog_df$pathologic.stage <- factor(prog_df$pathologic.stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

write.csv(prog_df, file.path(outdir, "05.merged_cleaned_before_complete_case.csv"), row.names = FALSE, quote = FALSE)

model_vars <- c("sample", "OS.time", "OS", "age", "gender", "race", "pathologic.stage", "riskScore")
prog_model_df <- prog_df[, model_vars, drop = FALSE]
prog_model_df <- prog_model_df[complete.cases(prog_model_df), , drop = FALSE]
prog_model_df <- prog_model_df[prog_model_df$OS.time > 0, , drop = FALSE]
prog_model_df <- prog_model_df[prog_model_df$OS %in% c(0, 1), , drop = FALSE]
prog_model_df$gender <- droplevels(prog_model_df$gender)
prog_model_df$race <- droplevels(prog_model_df$race)
prog_model_df$pathologic.stage <- droplevels(prog_model_df$pathologic.stage)
rownames(prog_model_df) <- prog_model_df$sample

if (nrow(prog_model_df) == 0) save_note_and_stop("No complete-case samples remain for clinical nomogram modeling.")

write.csv(prog_model_df, file.path(outdir, "06.analysis_dataset_complete_case.csv"), row.names = FALSE, quote = FALSE)
capture.output(summary(prog_model_df), file = dataset_summary_file)

log_message("Complete-case samples: ", nrow(prog_model_df))
log_message("Events: ", sum(prog_model_df$OS == 1), "; Censored: ", sum(prog_model_df$OS == 0))

uni_var_list <- c("riskScore", "age", "gender", "race", "pathologic.stage")
uni_overall_result <- data.frame()
uni_level_result <- data.frame()

for (feature_name in uni_var_list) {
  uni_formula <- as.formula(paste0("Surv(OS.time, OS) ~ ", feature_name))
  uni_fit <- coxph(uni_formula, data = prog_model_df)
  uni_sum <- summary(uni_fit)
  uni_coef <- as.data.frame(uni_sum$coefficients)
  uni_ci <- as.data.frame(uni_sum$conf.int)

  if (feature_name %in% c("riskScore", "age")) {
    overall_p <- uni_coef[1, "Pr(>|z|)"]
  } else {
    uni_drop1 <- drop1(uni_fit, test = "Chisq")
    overall_p <- uni_drop1[feature_name, "Pr(>Chi)"]
  }

  uni_overall_result <- rbind(
    uni_overall_result,
    data.frame(
      feature = feature_name,
      variable_type = ifelse(feature_name %in% c("riskScore", "age"), "continuous", "categorical"),
      overall_p = overall_p,
      stringsAsFactors = FALSE
    )
  )

  uni_level_result <- rbind(
    uni_level_result,
    data.frame(
      variable = rownames(uni_ci),
      HR = uni_ci[, "exp(coef)"],
      HR.95L = uni_ci[, "lower .95"],
      HR.95H = uni_ci[, "upper .95"],
      p.value = uni_coef[, "Pr(>|z|)"],
      feature = feature_name,
      row.names = NULL,
      check.names = FALSE
    )
  )
}

write.csv(uni_overall_result, file.path(outdir, "07.univariate_overall_test.csv"), row.names = FALSE, quote = FALSE)
write.csv(uni_level_result, file.path(outdir, "07.univariate_cox_all.csv"), row.names = FALSE, quote = FALSE)

uni_result_display <- build_display_table(uni_level_result, uni_var_list, uni_overall_result, prefix = "univariate")
sig_feature <- unique(uni_overall_result$feature[!is.na(uni_overall_result$overall_p) & uni_overall_result$overall_p < 0.05])
write.csv(data.frame(feature = sig_feature), file.path(outdir, "08.univariate_significant_features.csv"), row.names = FALSE, quote = FALSE)

uni_ph_result <- data.frame()
for (feature_name in uni_var_list) {
  uni_formula <- as.formula(paste0("Surv(OS.time, OS) ~ ", feature_name))
  uni_fit <- coxph(uni_formula, data = prog_model_df)
  uni_zph <- tryCatch(cox.zph(uni_fit), error = function(e) NULL)

  if (is.null(uni_zph)) {
    uni_ph_result <- rbind(
      uni_ph_result,
      data.frame(chisq = NA, df = NA, p = NA, term = NA, model = feature_name, note = "cox.zph_failed")
    )
    next
  }

  uni_zph_df <- as.data.frame(uni_zph$table)
  uni_zph_df$term <- rownames(uni_zph_df)
  uni_zph_df$model <- feature_name
  uni_zph_df$note <- "ok"
  uni_ph_result <- rbind(uni_ph_result, uni_zph_df)

  plot_to_pdf_png(
    plot_fun = function() print(ggcoxzph(uni_zph)),
    pdf_file = file.path(outdir, paste0("09.PH_univariate_", feature_name, ".pdf")),
    png_file = file.path(outdir, paste0("09.PH_univariate_", feature_name, ".png")),
    width = 8,
    height = 6
  )
}
write.csv(uni_ph_result, file.path(outdir, "09.univariate_cox_ph_test.csv"), row.names = FALSE, quote = FALSE)
draw_forestplot(uni_result_display, "Univariate Cox", "10", "univariate_cox_forest")

if (length(sig_feature) == 0) {
  save_note_and_stop("No variables with univariate overall p < 0.05 for multivariate Cox.")
}

multi_formula <- as.formula(paste0("Surv(OS.time, OS) ~ ", paste(sig_feature, collapse = " + ")))
log_message("Multivariate formula: ", deparse(multi_formula))

multi_fit <- coxph(multi_formula, data = prog_model_df)
multi_sum <- summary(multi_fit)
multi_ci <- as.data.frame(multi_sum$conf.int)
multi_coef <- as.data.frame(multi_sum$coefficients)

multi_result_df <- data.frame(
  variable = rownames(multi_ci),
  HR = multi_ci[, "exp(coef)"],
  HR.95L = multi_ci[, "lower .95"],
  HR.95H = multi_ci[, "upper .95"],
  p.value = multi_coef[, "Pr(>|z|)"],
  row.names = NULL,
  check.names = FALSE
)
write.csv(multi_result_df, file.path(outdir, "11.multivariate_cox_all.csv"), row.names = FALSE, quote = FALSE)

multi_overall_test <- drop1(multi_fit, test = "Chisq")
multi_overall_df <- data.frame(
  feature = rownames(multi_overall_test),
  Df = multi_overall_test[, "Df"],
  AIC = multi_overall_test[, "AIC"],
  LRT = multi_overall_test[, "LRT"],
  overall_p = multi_overall_test[, "Pr(>Chi)"],
  row.names = NULL,
  check.names = FALSE
)
multi_overall_df <- multi_overall_df[multi_overall_df$feature != "<none>", , drop = FALSE]
write.csv(multi_overall_df, file.path(outdir, "11.multivariate_overall_test.csv"), row.names = FALSE, quote = FALSE)

multi_result_display <- build_display_table(multi_result_df, sig_feature, multi_overall_df, prefix = "multivariate")

multi_zph <- cox.zph(multi_fit)
multi_zph_df <- as.data.frame(multi_zph$table)
multi_zph_df$term <- rownames(multi_zph_df)
write.csv(multi_zph_df, file.path(outdir, "13.multivariate_cox_ph_test.csv"), row.names = FALSE, quote = FALSE)
plot_to_pdf_png(
  plot_fun = function() print(ggcoxzph(multi_zph)),
  pdf_file = file.path(outdir, "13.multivariate_cox_ph_test.pdf"),
  png_file = file.path(outdir, "13.multivariate_cox_ph_test.png"),
  width = 10,
  height = 7
)

draw_forestplot(multi_result_display, "Multivariate Cox", "14", "multivariate_cox_forest")

ddist <- datadist(prog_model_df)
options(datadist = "ddist")

nom_fit <- cph(
  multi_formula,
  data = prog_model_df,
  surv = TRUE,
  x = TRUE,
  y = TRUE,
  time.inc = time_roc_days[1]
)
surv_fun <- Survival(nom_fit)

obs_idx <- if (nrow(prog_model_df) >= 2) 2 else 1
obs_data <- prog_model_df[obs_idx, , drop = FALSE]

plot_to_pdf_png(
  plot_fun = function() {
    regplot(
      nom_fit,
      plots = c("density", "boxes"),
      observation = obs_data,
      center = TRUE,
      subticks = TRUE,
      droplines = TRUE,
      title = "",
      points = TRUE,
      failtime = time_roc_days,
      prfail = FALSE,
      rank = "sd",
      clickable = FALSE
    )
  },
  pdf_file = file.path(outdir, "15.nomogram.pdf"),
  png_file = file.path(outdir, "15.nomogram.png"),
  width = 10,
  height = 6,
  png_res = 150
)

obs_lp <- predict(nom_fit, newdata = obs_data, type = "lp")
obs_surv <- vapply(time_roc_days, function(tt) as.numeric(surv_fun(tt, obs_lp)), numeric(1))
obs_result <- data.frame(
  sample = rownames(obs_data),
  riskScore = obs_data$riskScore,
  age = obs_data$age,
  gender = as.character(obs_data$gender),
  race = as.character(obs_data$race),
  pathologic.stage = as.character(obs_data$pathologic.stage),
  linear_predictor = as.numeric(obs_lp),
  surv_1y = if (length(obs_surv) >= 1) obs_surv[1] else NA_real_,
  surv_3y = if (length(obs_surv) >= 2) obs_surv[2] else NA_real_,
  surv_5y = if (length(obs_surv) >= 3) obs_surv[3] else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(obs_result, file.path(outdir, "15.nomogram_regplot_observation_prediction.csv"), row.names = FALSE, quote = FALSE)

calibration_ok <- TRUE
calibration_list <- list()
for (tt in time_roc_days) {
  cal_fit <- tryCatch(
    cph(multi_formula, data = prog_model_df, surv = TRUE, x = TRUE, y = TRUE, time.inc = tt),
    error = function(e) NULL
  )
  if (is.null(cal_fit)) {
    calibration_ok <- FALSE
    break
  }
  cal_obj <- tryCatch(
    calibrate(cal_fit, u = tt, cmethod = "KM", method = "boot", B = calibration_bootstrap),
    error = function(e) NULL
  )
  if (is.null(cal_obj)) {
    calibration_ok <- FALSE
    break
  }
  calibration_list[[as.character(tt)]] <- cal_obj
}

if (calibration_ok && length(calibration_list) > 0) {
  cal_cols <- c("#4DBBD5FF", "#E64B35FF", "#00A087FF")
  plot_to_pdf_png(
    plot_fun = function() {
      par(mar = c(6, 6, 3, 2))
      idx <- 1
      for (nm in names(calibration_list)) {
        cal_obj <- calibration_list[[nm]]
        plot(
          cal_obj,
          add = idx > 1,
          subtitles = FALSE,
          lwd = 2,
          lty = 1,
          errbar.col = cal_cols[idx],
          col = cal_cols[idx],
          xlab = "Nomogram-predicted survival probability",
          ylab = "Observed survival probability",
          xlim = c(0, 1),
          ylim = c(0, 1)
        )
        idx <- idx + 1
      }
      abline(0, 1, lty = 2, lwd = 2, col = "gray50")
      legend(
        "bottomright",
        legend = c(paste0(round(time_roc_days / 365), "-year"), "Ideal")[seq_len(length(calibration_list) + 1)],
        col = c(cal_cols[seq_len(length(calibration_list))], "gray50"),
        lwd = c(rep(2, length(calibration_list)), 2),
        lty = c(rep(1, length(calibration_list)), 2),
        bty = "n"
      )
    },
    pdf_file = file.path(outdir, "16.calibration_curve.pdf"),
    png_file = file.path(outdir, "16.calibration_curve.png"),
    width = 8,
    height = 8
  )
} else {
  log_message("Calibration skipped because cph/calibrate failed for at least one time point.")
}

prog_model_df$nomogram_lp <- predict(multi_fit, newdata = prog_model_df, type = "lp")
prog_model_df$OS.time.year <- prog_model_df$OS.time / 365
roc_times_year <- time_roc_days / 365
roc_res <- timeROC(
  T = prog_model_df$OS.time.year,
  delta = prog_model_df$OS,
  marker = prog_model_df$nomogram_lp,
  cause = 1,
  weighting = "marginal",
  times = roc_times_year,
  iid = TRUE
)

write.csv(
  data.frame(time_year = roc_times_year, AUC = roc_res$AUC),
  file.path(outdir, "17.timeROC_auc.csv"),
  row.names = FALSE,
  quote = FALSE
)

plot_to_pdf_png(
  plot_fun = function() {
    roc_cols <- c("#4DBBD5FF", "#E64B35FF", "#00A087FF")
    plot(roc_res, time = roc_times_year[1], col = roc_cols[1], lwd = 2, title = "", xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity")
    if (length(roc_times_year) >= 2) lines(roc_res$FP[, 2], roc_res$TP[, 2], col = roc_cols[2], lwd = 2)
    if (length(roc_times_year) >= 3) lines(roc_res$FP[, 3], roc_res$TP[, 3], col = roc_cols[3], lwd = 2)
    abline(0, 1, lty = 2, col = "gray50")
    legend(
      "bottomright",
      legend = paste0(paste0(round(roc_times_year), "-year"), ": AUC = ", sprintf("%.3f", roc_res$AUC)),
      col = roc_cols[seq_along(roc_res$AUC)],
      lty = 1,
      lwd = 2,
      bty = "n"
    )
  },
  pdf_file = file.path(outdir, "17.timeROC.pdf"),
  png_file = file.path(outdir, "17.timeROC.png"),
  width = 6,
  height = 6
)

dca_ok <- TRUE
tryCatch({
  prog_model_df$pred_1y_event <- 1 - surv_fun(time_roc_days[1], predict(nom_fit, type = "lp"))
  prog_model_df$event_1y <- ifelse(prog_model_df$OS == 1 & prog_model_df$OS.time <= time_roc_days[1], 1, 0)
  prog_model_df$pathologic.stage_num_dca <- as.numeric(prog_model_df$pathologic.stage)

  nomogram_dca <- decision_curve(
    event_1y ~ pred_1y_event,
    data = prog_model_df,
    family = binomial(link = "logit"),
    thresholds = seq(0.01, 0.99, by = 0.01),
    confidence.intervals = 0.95
  )
  risk_dca <- decision_curve(
    event_1y ~ riskScore,
    data = prog_model_df,
    family = binomial(link = "logit"),
    thresholds = seq(0.01, 0.99, by = 0.01),
    confidence.intervals = 0.95
  )
  stage_dca <- decision_curve(
    event_1y ~ pathologic.stage_num_dca,
    data = prog_model_df,
    family = binomial(link = "logit"),
    thresholds = seq(0.01, 0.99, by = 0.01),
    confidence.intervals = 0.95
  )

  dca_list <- list(nomogram_dca, risk_dca, stage_dca)
  plot_to_pdf_png(
    plot_fun = function() {
      plot_decision_curve(
        dca_list,
        curve.names = c("Nomogram_1y", "riskScore", "pathologic.stage"),
        cost.benefit.axis = FALSE,
        confidence.intervals = FALSE,
        standardize = FALSE,
        legend.position = "bottomright",
        xlab = "Threshold probability"
      )
    },
    pdf_file = file.path(outdir, "18.dca_curve.pdf"),
    png_file = file.path(outdir, "18.dca_curve.png"),
    width = 7,
    height = 6
  )
}, error = function(e) {
  dca_ok <<- FALSE
  log_message("DCA skipped: ", conditionMessage(e))
})

write_summary("sample_n_survival_before_merge: ", nrow(survival_use))
write_summary("sample_n_clinical_before_merge: ", nrow(clinical_use))
write_summary("sample_n_risk_before_merge: ", nrow(risk_use))
write_summary("sample_n_after_survival_clinical_merge: ", nrow(merge_sc))
write_summary("sample_n_after_all_merge: ", nrow(merge_all))
write_summary("sample_n_complete_case: ", nrow(prog_model_df))
write_summary("event_n_complete_case: ", sum(prog_model_df$OS == 1))
write_summary("censor_n_complete_case: ", sum(prog_model_df$OS == 0))
write_summary("variables_in_univariate_model: ", paste(uni_var_list, collapse = " + "))
write_summary("variables_in_multivariate_model: ", paste(sig_feature, collapse = " + "))
write_summary("main_route: A version, pathologic.stage main route")
write_summary("nomogram_time_points_days: ", paste(time_roc_days, collapse = ", "))
write_summary("warning: race keeps only white/black/asian; other values are set to NA according to current rule.")
write_summary("warning: DCA here is approximate binary 1-year event analysis for exploratory use, not strict censored-survival DCA.")
write_summary("multivariate_PH_global_p: ", sprintf("%.6f", multi_zph_df$p[multi_zph_df$term == "GLOBAL"][1]))
write_summary("timeROC_AUC: ", paste(sprintf("%.4f", roc_res$AUC), collapse = ", "))
write_summary("calibration_bootstrap: ", calibration_bootstrap)
write_summary("dca_generated: ", dca_ok)

log_message("Script completed: ", sub_id)

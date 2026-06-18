suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(rms)
  library(regplot)
  library(forestplot)
  library(timeROC)
  library(rmda)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(grid)
  library(pdftools)
  library(qpdf)
})

# =========================================================
# CLI args
# =========================================================
args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

risk_file     <- get_arg("--risk-file")
surv_file     <- get_arg("--surv")
clinical_file <- get_arg("--clinical")
train_id      <- get_arg("--train-id")
outdir        <- get_arg("--outdir")
logdir        <- get_arg("--logdir")
sn            <- get_arg("--sn")
summary_dir   <- get_arg("--summary-dir")
var_str       <- get_arg("--variables", required = FALSE, default = NULL)
time_set      <- get_arg("--time-set", required = FALSE, default = "135")
dca_year      <- as.numeric(get_arg("--dca-year", required = FALSE, default = "3"))
min_level_n   <- as.numeric(get_arg("--min-level-n", required = FALSE, default = "5"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

# Time points
if (time_set == "135") {
  t_days  <- c(365, 1095, 1825)
  t_years <- c(1, 3, 5)
} else if (time_set == "357") {
  t_days  <- c(1095, 1825, 2555)
  t_years <- c(3, 5, 7)
} else {
  stop("--time-set must be 135 or 357, got: ", time_set)
}
dca_days <- dca_year * 365

# Logging
# Logging
logsetup <- setup_logging(logdir, summary_dir, sn, "prognostic")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("[Start] Stage prognostic analysis")
log_message("Time: ", as.character(Sys.time()))
log_message("Risk: ", risk_file)
log_message("Survival: ", surv_file)
log_message("Clinical: ", clinical_file)
log_message("Time set: ", time_set, " (", paste(t_years, collapse = "/"), " years)")

# =========================================================
# 0. Input check
# =========================================================
for (f in c(surv_file, clinical_file, risk_file))
  if (!file.exists(f)) stop("File not found: ", f)

survival_df <- read.csv(surv_file, check.names = FALSE, stringsAsFactors = FALSE)
clinical_df <- read.csv(clinical_file, check.names = FALSE, stringsAsFactors = FALSE)
risk_df     <- read.csv(risk_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("[survival dim] ", paste(dim(survival_df), collapse = " x "), "\n")
cat("[clinical dim] ", paste(dim(clinical_df), collapse = " x "), "\n")
cat("[risk dim] ", paste(dim(risk_df), collapse = " x "), "\n")

# =========================================================
# 1. Column normalization
# =========================================================
colnames(survival_df)[1] <- "sample"
if (!all(c("sample", "OS.time", "OS") %in% colnames(survival_df)))
  stop("Survival file must have sample, OS.time, OS columns")

# Risk: match sample/id column
risk_id_col <- intersect(c("sample", "id", "Sample", "ID"), colnames(risk_df))
if (length(risk_id_col) == 0) stop("No sample/id column in risk file")
colnames(risk_df)[colnames(risk_df) == risk_id_col[1]] <- "sample"
if (!"riskScore" %in% colnames(risk_df)) stop("riskScore not found in risk file")

# Clinical: match sample column
if (!"sample" %in% colnames(clinical_df)) {
  candidates <- c("sample_id", "id", "Sample", "ID", "barcode", "Barcode")
  hit <- intersect(candidates, colnames(clinical_df))[1]
  if (is.na(hit)) stop("No sample column in clinical file")
  clinical_df$sample <- clinical_df[[hit]]
}

# =========================================================
# 2. Auto-detect + auto-clean helpers
# =========================================================
na_vals <- c("", "NA", "N/A", "na", "null", "NULL",
  "not reported", "Not Reported", "unknown", "Unknown",
  "[Not Available]", "[Not Applicable]", "not applicable", "N/A", "none", "None")

auto_detect_type <- function(x) {
  if (is.factor(x)) return("categorical")
  if (is.character(x)) {
    nx <- suppressWarnings(as.numeric(x))
    if (sum(!is.na(nx)) / length(nx) > 0.8) return("continuous")
    return("categorical")
  }
  if (is.numeric(x)) {
    ux <- unique(na.omit(x))
    if (length(ux) <= 5 && all(ux == floor(ux))) return("categorical")
    return("continuous")
  }
  return("categorical")
}

merge_substages <- function(x, vname) {
  # Stage IA/IB/IIA/IIB/IIIA/IIIB/IVA/IVB → Stage I/II/III/IV
  # Also handles T1a/T1b → T1, etc.
  vlower <- tolower(vname)
  xc <- as.character(x)
  if (grepl("stage", vlower)) {
    xc <- gsub("(?i)\\s*[abc]$", "", xc)  # drop trailing A/B/C
    xc <- gsub("(?i)^stage\\s+", "Stage ", xc)
    # Normalize: Stage Iii → Stage III, Stage Ii → Stage II, Stage Iv → Stage IV
    xc <- gsub("(?i)^Stage Iii", "Stage III", xc)
    xc <- gsub("(?i)^Stage Ii", "Stage II", xc)
    xc <- gsub("(?i)^Stage Iv", "Stage IV", xc)
    xc <- gsub("(?i)^Stage I", "Stage I", xc)
  }
  if (grepl("^[tnm]_?stage", vlower)) {
    # T1a/T1b → T1; T2a/T2b → T2; also remove trailing a/b
    xc <- gsub("(?i)([0-9])[a-c]$", "\\1", xc)
  }
  xc
}

drop_x_levels <- function(xf) {
  # Remove "unknown" X levels: TX, NX, MX, Stage X, X
  xc <- as.character(xf)
  is_x <- grepl("^(T|N|M)?X$", xc, ignore.case = TRUE) |
          grepl("^Stage X$", xc, ignore.case = TRUE) |
          grepl("^X$", xc, ignore.case = TRUE)
  if (any(is_x)) {
    cat("  dropping X levels: ", paste(unique(xc[is_x]), collapse = ", "), "\n", sep = "")
    xc[is_x] <- NA
  }
  xc
}

auto_clean <- function(x, vtype, vname = "", min_n = 5) {
  xc <- trimws(as.character(x))
  xc[xc %in% na_vals] <- NA
  if (vtype == "continuous") {
    return(suppressWarnings(as.numeric(xc)))
  } else {
    # Merge substages + drop X levels before factoring
    xc <- merge_substages(xc, vname)
    xf <- droplevels(as.factor(xc))
    xf <- droplevels(as.factor(drop_x_levels(xf)))
    if (nlevels(xf) < 2) return(rep(NA, length(xf)))
    tb <- table(xf)
    small <- names(tb[tb < min_n])
    if (length(small) > 0)
      cat("  dropping levels (n<", min_n, "): ", paste(small, collapse = ", "), "\n", sep = "")
    xf[xf %in% small] <- NA
    xf <- droplevels(xf)
    if (nlevels(xf) < 2) return(rep(NA, length(xf)))
    return(xf)
  }
}

# Now update the call site to pass vname
clean_one <- function(df, v, min_n) {
  vt <- auto_detect_type(df[[v]])
  var_types[[v]] <<- vt
  cat("[", v, "] type=", vt, "\n", sep = "")
  cleaned <- auto_clean(df[[v]], vt, v, min_n)
  if (all(is.na(cleaned))) {
    cat("  -> all NA, dropping variable\n")
    var_types[[v]] <<- "dropped"
    NULL
  } else {
    cleaned
  }
}

# =========================================================
# 3. Determine variables
# =========================================================
riskScore_col <- "riskScore"
clinical_candidates <- setdiff(colnames(clinical_df), c("sample", "sample_id",
  intersect(c("sample_id", "id", "Sample", "ID", "barcode", "Barcode"), colnames(clinical_df))))

if (!is.null(var_str) && nzchar(trimws(var_str))) {
  user_vars <- trimws(unlist(strsplit(var_str, ",")))
  user_vars <- user_vars[nzchar(user_vars)]
  # Ensure riskScore is included
  if (!riskScore_col %in% user_vars) user_vars <- c(riskScore_col, user_vars)
  cat("[User variables] ", paste(user_vars, collapse = ", "), "\n")
  analysis_vars <- user_vars
} else {
  analysis_vars <- c(riskScore_col, clinical_candidates)
  cat("[Auto variables] ", paste(analysis_vars, collapse = ", "), "\n")
}

# =========================================================
# 4. Merge
# =========================================================
survival_use <- survival_df[, c("sample", "OS.time", "OS"), drop = FALSE]
risk_use     <- risk_df[, c("sample", riskScore_col), drop = FALSE]
clin_cols    <- c("sample", intersect(analysis_vars, colnames(clinical_df)))
clinical_use <- clinical_df[, clin_cols, drop = FALSE]

cat("[N pre-merge] surv=", nrow(survival_use), " clin=", nrow(clinical_use),
    " risk=", nrow(risk_use), "\n")

merge_sc  <- merge(survival_use, clinical_use, by = "sample")
merge_all <- merge(merge_sc, risk_use, by = "sample")
cat("[N after merge] ", nrow(merge_all), "\n")
write.csv(merge_all, file.path(outdir, "01.merged_raw.csv"), row.names = FALSE)

# =========================================================
# 5. Data cleaning
# =========================================================
prog_df <- merge_all
prog_df$OS.time <- as.numeric(prog_df$OS.time)
prog_df$OS <- as.numeric(prog_df$OS)
prog_df$OS.time[prog_df$OS.time <= 0] <- NA

var_types <- list()
for (v in analysis_vars) {
  if (v == riskScore_col) {
    var_types[[v]] <- "continuous"
    prog_df[[v]] <- as.numeric(prog_df[[v]])
    next
  }
  if (!v %in% colnames(prog_df)) next
  prog_df[[v]] <- clean_one(prog_df, v, min_level_n)
  if (is.null(prog_df[[v]])) prog_df[[v]] <- NULL
}
analysis_vars <- intersect(analysis_vars, colnames(prog_df))
analysis_vars <- analysis_vars[!sapply(analysis_vars, function(v) isTRUE(var_types[[v]] == "dropped"))]
cat("[Active variables] ", paste(analysis_vars, collapse = ", "), "\n")

write.csv(prog_df, file.path(outdir, "02.merged_cleaned.csv"), row.names = FALSE)

# =========================================================
# 6. Complete case
# =========================================================
mv <- unique(c("sample", "OS.time", "OS", analysis_vars))
prog_model <- prog_df[, intersect(mv, colnames(prog_df)), drop = FALSE]
prog_model <- prog_model[complete.cases(prog_model), ]
prog_model <- prog_model[prog_model$OS.time > 0, ]
prog_model <- prog_model[prog_model$OS %in% c(0, 1), ]

for (v in analysis_vars) {
  if (v %in% colnames(prog_model) && isTRUE(var_types[[v]] == "categorical")) {
    if (is.factor(prog_model[[v]])) prog_model[[v]] <- droplevels(prog_model[[v]])
  }
}
rownames(prog_model) <- prog_model$sample
cat("[Complete case N] ", nrow(prog_model), " events=", sum(prog_model$OS == 1),
    " censor=", sum(prog_model$OS == 0), "\n")
write.csv(prog_model, file.path(outdir, "03.analysis_dataset.csv"), row.names = FALSE)

# =========================================================
# 7. Distribution
# =========================================================
for (v in analysis_vars) {
  if (v %in% colnames(prog_model)) {
    cat("[", v, "]\n", sep = "")
    print(table(prog_model[[v]], useNA = "ifany"))
    cat("[", v, " x OS]\n", sep = "")
    print(table(prog_model[[v]], prog_model$OS, useNA = "ifany"))
  }
}
capture.output(summary(prog_model), file = file.path(outdir, "04.dataset_summary.txt"))

# =========================================================
# 8. Univariate Cox
# =========================================================
uni_overall <- data.frame()
uni_level   <- data.frame()

for (v in analysis_vars) {
  if (!v %in% colnames(prog_model)) next
  uf <- as.formula(paste0("Surv(OS.time, OS) ~ `", v, "`"))
  fit <- tryCatch(coxph(uf, data = prog_model), error = function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  coef <- as.data.frame(s$coefficients)
  ci   <- as.data.frame(s$conf.int)
  vt   <- var_types[[v]]

  if (vt == "continuous") {
    op <- coef[1, "Pr(>|z|)"]
  } else {
    dd <- tryCatch(drop1(fit, test = "Chisq"), error = function(e) NULL)
    op <- if (!is.null(dd) && nrow(dd) > 0) dd[nrow(dd), "Pr(>Chi)"] else NA_real_
  }

  uni_overall <- rbind(uni_overall, data.frame(
    feature = v, variable_type = vt, overall_p = op, stringsAsFactors = FALSE))

  ul <- data.frame(variable = rownames(ci), HR = ci[, "exp(coef)"],
    HR.95L = ci[, "lower .95"], HR.95H = ci[, "upper .95"],
    p.value = coef[, "Pr(>|z|)"], feature = v, row.names = NULL, check.names = FALSE)
  uni_level <- rbind(uni_level, ul)
}

write.csv(uni_overall, file.path(outdir, "05.univariate_overall.csv"), row.names = FALSE)
write.csv(uni_level, file.path(outdir, "05.univariate_level.csv"), row.names = FALSE)

# =========================================================
# 9. Univariate display table (generic)
# =========================================================
build_display_table <- function(level_df, overall_df, model_data, vars) {
  display <- data.frame()
  for (v in vars) {
    tmp <- level_df[level_df$feature == v, , drop = FALSE]
    if (nrow(tmp) == 0) next
    vt <- var_types[[v]]
    nlv <- if (v %in% colnames(model_data) && is.factor(model_data[[v]]))
      nlevels(model_data[[v]]) else 0

    if (vt == "continuous") {
      tmp$display_label <- v
      tmp$is_reference <- FALSE; tmp$is_summary <- FALSE
      display <- rbind(display, tmp)
    } else if (vt == "categorical" && nlv == 2) {
      # 2 levels: show "level2 vs level1", NO summary row, NO overall P
      ref <- levels(model_data[[v]])[1]
      lv <- sub(paste0("^`?", v, "`?"), "", tmp$variable)
      tmp$display_label <- paste0(lv, " vs ", ref)
      tmp$is_reference <- FALSE; tmp$is_summary <- FALSE
      display <- rbind(display, tmp)
    } else if (vt == "categorical" && nlv > 2) {
      # >2 levels: summary row with overall P, level rows show just level name
      ref <- levels(model_data[[v]])[1]
      ref_row <- data.frame(variable = paste0(v, "_ref"), HR = 1, HR.95L = NA, HR.95H = NA,
        p.value = NA, feature = v,
        display_label = paste0(v, " (Ref: ", ref, ")"),
        is_reference = TRUE, is_summary = TRUE, stringsAsFactors = FALSE)
      lv <- sub(paste0("^`?", v, "`?"), "", tmp$variable)
      tmp$display_label <- lv   # just the level name, no "vs" suffix
      tmp$is_reference <- FALSE; tmp$is_summary <- FALSE
      display <- rbind(display, ref_row, tmp)
    }
  }

  if (nrow(display) == 0) return(display)
  display$overall_p <- NA_real_
  for (i in seq_len(nrow(display))) {
    if (display$is_summary[i]) {
      feat <- display$feature[i]
      if (feat %in% overall_df$feature)
        display$overall_p[i] <- overall_df$overall_p[match(feat, overall_df$feature)]
    }
  }
  display$display_p <- ifelse(display$is_summary, display$overall_p, display$p.value)
  for (i in seq_len(nrow(display))) {
    if (display$is_summary[i] && !is.na(display$overall_p[i])) {
      pt <- ifelse(display$overall_p[i] < 0.001, "<0.001", sprintf("%.4f", display$overall_p[i]))
      display$display_label[i] <- paste0(display$display_label[i], "; overall P = ", pt)
    }
  }
  display
}

uni_display <- build_display_table(uni_level, uni_overall, prog_model, analysis_vars)
write.csv(uni_display, file.path(outdir, "06.univariate_display.csv"), row.names = FALSE)

sig_feat <- unique(uni_overall$feature[uni_overall$overall_p < 0.05 & !is.na(uni_overall$overall_p)])
cat("[Uni significant] ", paste(sig_feat, collapse = ", "), "\n")

# =========================================================
# 10. Univariate PH test
# =========================================================
uni_ph <- data.frame()
for (v in analysis_vars) {
  if (!v %in% colnames(prog_model)) next
  uf <- as.formula(paste0("Surv(OS.time, OS) ~ `", v, "`"))
  fit <- coxph(uf, data = prog_model)
  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  if (is.null(zph)) {
    uni_ph <- rbind(uni_ph, data.frame(chisq = NA, df = NA, p = NA,
      term = NA, model = v, note = "cox.zph_failed", stringsAsFactors = FALSE))
  } else {
    zd <- as.data.frame(zph$table)
    zd$term <- rownames(zd); zd$model <- v; zd$note <- "ok"
    uni_ph <- rbind(uni_ph, zd)
    pdf(file.path(outdir, paste0("07.PH_", make.names(v), ".pdf")), width = 8, height = 6, family = "Times")
    print(ggcoxzph(zph))
    dev.off()
    png(file.path(outdir, paste0("07.PH_", make.names(v), ".png")), width = 8, height = 6, units = "in", res = 600)
    print(ggcoxzph(zph))
    dev.off()
  }
}
write.csv(uni_ph, file.path(outdir, "07.PH_test.csv"), row.names = FALSE)

# =========================================================
# 11. Univariate forest
# =========================================================
ufd <- uni_display
ufd$HR_CI <- ifelse(is.na(ufd$HR.95L) | is.na(ufd$HR.95H), "",
  paste0(sprintf("%.3f", ufd$HR), " (", sprintf("%.3f", ufd$HR.95L), "-",
         sprintf("%.3f", ufd$HR.95H), ")"))
utxt <- cbind(
  c("Variable", ufd$display_label),
  c("P value", ifelse(is.na(ufd$display_p), "",
    ifelse(ufd$display_p < 0.001, "<0.001", sprintf("%.4f", ufd$display_p)))),
  c("Hazard Ratio (95% CI)", ufd$HR_CI))

fp_args <- list(labeltext = utxt, graph.pos = 3,
  is.summary = c(TRUE, ufd$is_summary),
  mean = c(NA, ifelse(ufd$is_reference, NA, ufd$HR)),
  lower = c(NA, ufd$HR.95L), upper = c(NA, ufd$HR.95H),
  zero = 1, boxsize = 0.18, lwd.ci = 2, ci.vertices = TRUE, ci.vertices.height = 0.08,
  colgap = unit(5, "mm"), lineheight = unit(0.9, "cm"), graphwidth = unit(0.35, "npc"),
  xlab = "Hazard Ratio", title = "Univariate Cox",
  col = fpColors(box = "#d9b02c", lines = "#4981bb", zero = "gray50"),
  txt_gp = fpTxtGp(label = gpar(cex = 0.9), ticks = gpar(cex = 0.85),
    xlab = gpar(cex = 1), title = gpar(cex = 1.1), summary = gpar(cex = 0.95, fontface = "bold")))

pdf(file.path(outdir, "08.univariate_forest.pdf"),
    height = max(6, 0.35 * nrow(ufd) + 2), width = 14, onefile = FALSE, family = "Times")
do.call(forestplot, fp_args)
dev.off()
png(file.path(outdir, "08.univariate_forest.png"),
    height = max(6, 0.35 * nrow(ufd) + 2), width = 14, units = "in", res = 600)
do.call(forestplot, fp_args)
dev.off()

# =========================================================
# 12. Multivariate Cox
# =========================================================
if (length(sig_feat) == 0) {
  cat("ERROR: No univariate significant variables\n")
  stop("No univariate significant variables")
}

multi_formula <- as.formula(paste0("Surv(OS.time, OS) ~ `",
  paste(sig_feat, collapse = "` + `"), "`"))
cat("[Multi formula] ", deparse(multi_formula), "\n")

multi_fit  <- coxph(multi_formula, data = prog_model)
multi_s    <- summary(multi_fit)
multi_ci   <- as.data.frame(multi_s$conf.int)
multi_coef <- as.data.frame(multi_s$coefficients)
multi_res  <- data.frame(variable = rownames(multi_ci), HR = multi_ci[, "exp(coef)"],
  HR.95L = multi_ci[, "lower .95"], HR.95H = multi_ci[, "upper .95"],
  p.value = multi_coef[, "Pr(>|z|)"], row.names = NULL, check.names = FALSE)
# Map each coefficient term back to its parent variable
multi_res$feature <- NA_character_
for (v in sig_feat) {
  idx <- grepl(paste0("^`?", v, "`?"), multi_res$variable)
  multi_res$feature[idx] <- v
}
write.csv(multi_res, file.path(outdir, "09.multivariate_cox.csv"), row.names = FALSE)

# =========================================================
# 13. Multi-Cox p<0.05 filter -> nomogram variables
# =========================================================
multi_sig <- c()
for (v in sig_feat) {
  vt <- var_types[[v]]
  if (vt == "continuous") {
    m <- grep(paste0("^`?", v, "`?$"), multi_res$variable)
    if (length(m) > 0 && !is.na(multi_res$p.value[m[1]]) && multi_res$p.value[m[1]] < 0.05)
      multi_sig <- c(multi_sig, v)
  } else {
    ps <- multi_res$p.value[grepl(paste0("^`?", v, "`?"), multi_res$variable)]
    if (length(ps) > 0 && any(!is.na(ps) & ps < 0.05))
      multi_sig <- c(multi_sig, v)
  }
}
cat("[Multi-Cox p<0.05] ", paste(multi_sig, collapse = ", "), "\n")
write.csv(data.frame(feature = multi_sig), file.path(outdir, "10.nomogram_variables.csv"), row.names = FALSE)

# Stop if multi-Cox filter leaves <=1 variable or riskScore fails
if (length(multi_sig) <= 1) {
  log_message("ERROR: Only ", length(multi_sig), " variable(s) passed multi-Cox p<0.05. ",
      "Insufficient for nomogram construction.")
  stop("Insufficient multi-Cox significant variables for nomogram")
}
if (!("riskScore" %in% multi_sig)) {
  log_message("ERROR: riskScore did not pass multi-Cox p<0.05 filter. ",
      "The risk model is not an independent prognostic factor.")
  stop("riskScore not significant in multi-Cox model")
}

# Multi-Cox overall test (drop1)
multi_drop1 <- drop1(multi_fit, test = "Chisq")
multi_overall <- data.frame(
  feature = rownames(multi_drop1),
  overall_p = multi_drop1[, "Pr(>Chi)"],
  row.names = NULL, check.names = FALSE)
multi_overall <- multi_overall[multi_overall$feature != "<none>", , drop = FALSE]

# =========================================================
# 14. Multivariate display (multi-Cox p<0.05 variables only)
# =========================================================
multi_display <- build_display_table(multi_res, multi_overall, prog_model, multi_sig)
write.csv(multi_display, file.path(outdir, "11.multivariate_display.csv"), row.names = FALSE)

# =========================================================
# 15. Multivariate forest + PH
# =========================================================
mfd <- multi_display
mfd$HR_CI <- ifelse(is.na(mfd$HR.95L) | is.na(mfd$HR.95H), "",
  paste0(sprintf("%.3f", mfd$HR), " (", sprintf("%.3f", mfd$HR.95L), "-",
         sprintf("%.3f", mfd$HR.95H), ")"))
mtxt <- cbind(
  c("Variable", mfd$display_label),
  c("P value", ifelse(is.na(mfd$display_p), "",
    ifelse(mfd$display_p < 0.001, "<0.001", sprintf("%.4f", mfd$display_p)))),
  c("Hazard Ratio (95% CI)", mfd$HR_CI))

fp_args2 <- list(labeltext = mtxt, graph.pos = 3,
  is.summary = c(TRUE, mfd$is_summary),
  mean = c(NA, ifelse(mfd$is_reference, NA, mfd$HR)),
  lower = c(NA, mfd$HR.95L), upper = c(NA, mfd$HR.95H),
  zero = 1, boxsize = 0.18, lwd.ci = 2, ci.vertices = TRUE, ci.vertices.height = 0.08,
  colgap = unit(5, "mm"), lineheight = unit(0.9, "cm"), graphwidth = unit(0.35, "npc"),
  xlab = "Hazard Ratio", title = "Multivariate Cox",
  col = fpColors(box = "#d9b02c", lines = "#4981bb", zero = "gray50"),
  txt_gp = fpTxtGp(label = gpar(cex = 0.9), ticks = gpar(cex = 0.85),
    xlab = gpar(cex = 1), title = gpar(cex = 1.1), summary = gpar(cex = 0.95, fontface = "bold")))

pdf(file.path(outdir, "12.multivariate_forest.pdf"),
    height = max(6, 0.35 * nrow(mfd) + 2), width = 14, onefile = FALSE, family = "Times")
do.call(forestplot, fp_args2)
dev.off()
png(file.path(outdir, "12.multivariate_forest.png"),
    height = max(6, 0.35 * nrow(mfd) + 2), width = 14, units = "in", res = 600)
do.call(forestplot, fp_args2)
dev.off()

# Multi PH
mzph <- cox.zph(multi_fit)
mzph_df <- as.data.frame(mzph$table)
mzph_df$term <- rownames(mzph_df)
write.csv(mzph_df, file.path(outdir, "13.multivariate_PH.csv"), row.names = FALSE)
pdf(file.path(outdir, "13.multivariate_PH.pdf"), width = 10, height = 7, family = "Times")
print(ggcoxzph(mzph))
dev.off()
png(file.path(outdir, "13.multivariate_PH.png"), width = 10, height = 7, units = "in", res = 600)
print(ggcoxzph(mzph))
dev.off()

# =========================================================
# 16. Nomogram (multi-Cox p<0.05 only)
# =========================================================
# Nomogram formula without backticks (cph compatibility)
nom_vars_clean <- gsub("`", "", multi_sig)
nom_formula <- as.formula(paste0("Surv(OS.time, OS) ~ ",
  paste(nom_vars_clean, collapse = " + ")))
cat("[Nomogram formula] ", deparse(nom_formula), "\n")

ddist <- datadist(prog_model)
options(datadist = "ddist")

nom_fit <- tryCatch(
  cph(nom_formula, data = prog_model, surv = TRUE, x = TRUE, y = TRUE, time.inc = t_days[1]),
  error = function(e) stop("cph failed: ", e$message))
surv_fun <- Survival(nom_fit)

obs_idx <- min(2, nrow(prog_model))
obs_data <- prog_model[obs_idx, , drop = FALSE]

# Nomogram — regplot closes current device and opens its own, so redirect
# the default device instead of pre-opening png()/pdf(). PDF gets a blank
# first page which we strip with pdftools/qpdf. Falls back to rms::nomogram
# when regplot fails (intermittent subscript error in headless mode).

regplot_ok <- tryCatch({
  save_regplot <- function(device_fun, outfile) {
    old_device <- getOption("device")
    options(device = device_fun)
    on.exit(options(device = old_device), add = TRUE)
    regplot(nom_fit, plots = c("density", "boxes"), observation = obs_data,
      center = TRUE, subticks = TRUE, droplines = TRUE, title = "",
      points = TRUE, failtime = t_days, prfail = FALSE, rank = "sd", clickable = FALSE)
    dev.off()
  }
  save_regplot(
    device_fun = function(...) png(
      filename = file.path(outdir, "14.nomogram.png"),
      width = 4000, height = 2400, res = 300, bg = "white"),
    outfile = file.path(outdir, "14.nomogram.png")
  )
  tmp_pdf <- file.path(outdir, "14.nomogram.tmp.pdf")
  save_regplot(
    device_fun = function(...) pdf(file = tmp_pdf, width = 10, height = 6, bg = "white"),
    outfile = tmp_pdf
  )
  info <- pdftools::pdf_info(tmp_pdf)
  pages <- if (info$pages > 1) 2:info$pages else 1
  qpdf::pdf_subset(input = tmp_pdf, pages = pages,
    output = file.path(outdir, "14.nomogram.pdf"))
  unlink(tmp_pdf)
  TRUE
}, error = function(e) {
  cat("[regplot failed, falling back to rms::nomogram] ", e$message, "\n")
  FALSE
})

if (!regplot_ok) {
  fun_list <- lapply(seq_along(t_days), function(i) {
    force(i); local({ ti <- t_days[i]; function(x) surv_fun(ti, x) })
  })
  nomo <- nomogram(nom_fit, fun = fun_list,
    funlabel = paste0(t_years, "-Year Survival"),
    fun.at = seq(0.1, 0.9, by = 0.1))
  pdf(file.path(outdir, "14.nomogram.pdf"), width = 10, height = 6)
  plot(nomo)
  dev.off()
  png(file.path(outdir, "14.nomogram.png"), width = 4000, height = 2400, res = 300)
  plot(nomo)
  dev.off()
}

obs_lp <- predict(nom_fit, newdata = obs_data, type = "lp")
obs_pred <- data.frame(sample = rownames(obs_data), linear_predictor = as.numeric(obs_lp),
  stringsAsFactors = FALSE)
for (i in seq_along(t_days))
  obs_pred[[paste0("surv_", t_years[i], "y")]] <- as.numeric(surv_fun(t_days[i], obs_lp))
for (v in multi_sig)
  if (v %in% colnames(obs_data)) obs_pred[[v]] <- as.character(obs_data[[v]])
write.csv(obs_pred, file.path(outdir, "14.nomogram_prediction.csv"), row.names = FALSE)

# =========================================================
# 17. Calibration
# =========================================================
cal_cols <- c("#4DBBD5FF", "#E64B35FF", "#00A087FF")
pdf(file.path(outdir, "15.calibration.pdf"), width = 8, height = 8, family = "Times")
par(mar = c(6, 6, 3, 2))
for (i in seq_along(t_days)) {
  cf <- cph(nom_formula, data = prog_model, surv = TRUE, x = TRUE, y = TRUE, time.inc = t_days[i])
  co <- calibrate(cf, u = t_days[i], cmethod = "KM", method = "boot", B = 1000)
  if (i == 1) {
    plot(co, subtitles = FALSE, lwd = 2, lty = 1, errbar.col = cal_cols[i], col = cal_cols[i],
      xlab = "Nomogram-predicted survival", ylab = "Observed survival", xlim = c(0, 1), ylim = c(0, 1))
  } else {
    plot(co, add = TRUE, subtitles = FALSE, lwd = 2, lty = 1, errbar.col = cal_cols[i], col = cal_cols[i],
      xlim = c(0, 1), ylim = c(0, 1))
  }
}
abline(0, 1, lty = 2, lwd = 2, col = "gray50")
legend("bottomright",
  legend = c(paste0(t_years, "-year"), "Ideal"),
  col = c(cal_cols[seq_along(t_years)], "gray50"), lwd = 2, lty = c(rep(1, length(t_years)), 2), bty = "n")
dev.off()

png(file.path(outdir, "15.calibration.png"), width = 8, height = 8, units = "in", res = 600)
par(mar = c(6, 6, 3, 2))
for (i in seq_along(t_days)) {
  cf <- cph(nom_formula, data = prog_model, surv = TRUE, x = TRUE, y = TRUE, time.inc = t_days[i])
  co <- calibrate(cf, u = t_days[i], cmethod = "KM", method = "boot", B = 1000)
  if (i == 1) {
    plot(co, subtitles = FALSE, lwd = 2, lty = 1, errbar.col = cal_cols[i], col = cal_cols[i],
      xlab = "Nomogram-predicted survival", ylab = "Observed survival", xlim = c(0, 1), ylim = c(0, 1))
  } else {
    plot(co, add = TRUE, subtitles = FALSE, lwd = 2, lty = 1, errbar.col = cal_cols[i], col = cal_cols[i],
      xlim = c(0, 1), ylim = c(0, 1))
  }
}
abline(0, 1, lty = 2, lwd = 2, col = "gray50")
legend("bottomright",
  legend = c(paste0(t_years, "-year"), "Ideal"),
  col = c(cal_cols[seq_along(t_years)], "gray50"), lwd = 2, lty = c(rep(1, length(t_years)), 2), bty = "n")
dev.off()

# =========================================================
# 18. timeROC
# =========================================================
roc_fit <- coxph(nom_formula, data = prog_model)
prog_model$nomogram_lp <- predict(roc_fit, newdata = prog_model, type = "lp")
prog_model$OS.time.year <- prog_model$OS.time / 365

roc_res <- timeROC(T = prog_model$OS.time.year, delta = prog_model$OS,
  marker = prog_model$nomogram_lp, cause = 1, weighting = "marginal",
  times = t_years, iid = TRUE)

write.csv(data.frame(time_year = t_years, AUC = roc_res$AUC),
  file.path(outdir, "16.timeROC_auc.csv"), row.names = FALSE)

pdf(file.path(outdir, "16.timeROC.pdf"), width = 6, height = 6, family = "Times")
plot(roc_res, time = t_years[1], col = cal_cols[1], lwd = 2, title = "",
  xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity")
if (length(t_years) >= 2) lines(roc_res$FP[, 2], roc_res$TP[, 2], col = cal_cols[2], lwd = 2)
if (length(t_years) >= 3) lines(roc_res$FP[, 3], roc_res$TP[, 3], col = cal_cols[3], lwd = 2)
abline(0, 1, lty = 2, col = "gray50")
legend("bottomright",
  legend = paste0(t_years, "-year: AUC = ", sprintf("%.3f", roc_res$AUC)),
  col = cal_cols[seq_along(t_years)], lty = 1, lwd = 2, bty = "n")
dev.off()

png(file.path(outdir, "16.timeROC.png"), width = 6, height = 6, units = "in", res = 600)
plot(roc_res, time = t_years[1], col = cal_cols[1], lwd = 2, title = "",
  xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity")
if (length(t_years) >= 2) lines(roc_res$FP[, 2], roc_res$TP[, 2], col = cal_cols[2], lwd = 2)
if (length(t_years) >= 3) lines(roc_res$FP[, 3], roc_res$TP[, 3], col = cal_cols[3], lwd = 2)
abline(0, 1, lty = 2, col = "gray50")
legend("bottomright",
  legend = paste0(t_years, "-year: AUC = ", sprintf("%.3f", roc_res$AUC)),
  col = cal_cols[seq_along(t_years)], lty = 1, lwd = 2, bty = "n")
dev.off()

# =========================================================
# 19. DCA curve
# =========================================================
evt_col <- paste0("event_", dca_year, "y")
prog_model[[evt_col]] <- ifelse(prog_model$OS == 1 & prog_model$OS.time <= dca_days, 1, 0)
prog_model$pred_event <- 1 - surv_fun(dca_days, predict(nom_fit, type = "lp"))

dca_list <- list()
dca_names <- character(0)

# Nomogram
dca_list[[1]] <- decision_curve(
  as.formula(paste0(evt_col, " ~ pred_event")),
  data = prog_model, family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.99, by = 0.01), confidence.intervals = 0.95)
dca_names <- c(dca_names, "Nomogram")

# riskScore
dca_list[[length(dca_list) + 1]] <- decision_curve(
  as.formula(paste0(evt_col, " ~ riskScore")),
  data = prog_model, family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.99, by = 0.01), confidence.intervals = 0.95)
dca_names <- c(dca_names, "riskScore")

# Best categorical variable
cat_vars <- analysis_vars[sapply(analysis_vars, function(v)
  isTRUE(var_types[[v]] == "categorical") && v != riskScore_col && v %in% colnames(prog_model))]
if (length(cat_vars) > 0) {
  nlvls <- sapply(cat_vars, function(v) nlevels(prog_model[[v]]))
  best_cat <- cat_vars[which.max(nlvls)]
  prog_model$dca_cat_val <- as.numeric(prog_model[[best_cat]])
  dca_list[[length(dca_list) + 1]] <- decision_curve(
    as.formula(paste0(evt_col, " ~ dca_cat_val")),
    data = prog_model, family = binomial(link = "logit"),
    thresholds = seq(0.01, 0.99, by = 0.01), confidence.intervals = 0.95)
  dca_names <- c(dca_names, best_cat)
}

dca_cols <- c("black", "red", "green3")[seq_along(dca_list)]
dca_ltys <- rep(1, length(dca_list))

pdf(file.path(outdir, "17.dca_curve.pdf"), width = 7, height = 6, family = "Times")
plot_decision_curve(dca_list, curve.names = dca_names, cost.benefit.axis = FALSE,
  confidence.intervals = FALSE, standardize = FALSE, legend.position = "none",
  col = dca_cols, lty = dca_ltys, lwd = rep(2, length(dca_list)),
  xlab = "Threshold probability")
legend("bottomright", legend = c(dca_names, "All", "None"),
  col = c(dca_cols, "gray50", "gray50"), lty = c(dca_ltys, 2, 3),
  lwd = c(rep(2, length(dca_list)), 1, 1), cex = 0.8, bty = "n")
dev.off()

png(file.path(outdir, "17.dca_curve.png"), width = 7, height = 6, units = "in", res = 600)
plot_decision_curve(dca_list, curve.names = dca_names, cost.benefit.axis = FALSE,
  confidence.intervals = FALSE, standardize = FALSE, legend.position = "none",
  col = dca_cols, lty = dca_ltys, lwd = rep(2, length(dca_list)),
  xlab = "Threshold probability")
legend("bottomright", legend = c(dca_names, "All", "None"),
  col = c(dca_cols, "gray50", "gray50"), lty = c(dca_ltys, 2, 3),
  lwd = c(rep(2, length(dca_list)), 1, 1), cex = 0.8, bty = "n")
dev.off()

# =========================================================
# 20. Summary
# =========================================================
write_summary(c(
  "Stage prognostic analysis summary",
  paste0("Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("SN: ", sn),
  paste0("Risk: ", risk_file),
  paste0("Survival: ", surv_file),
  paste0("Clinical: ", clinical_file),
  paste0("Time set: ", time_set, " (", paste(t_years, collapse = "/"), " years)"),
  paste0("DCA year: ", dca_year),
  paste0("Min level N: ", min_level_n),
  paste0("Variables: ", paste(analysis_vars, collapse = ", ")),
  paste0("Complete case N: ", nrow(prog_model)),
  paste0("Events: ", sum(prog_model$OS == 1)),
  paste0("Uni-Cox significant: ", paste(sig_feat, collapse = ", ")),
  paste0("Multi-Cox significant: ", paste(multi_sig, collapse = ", ")),
  paste0("Nomogram variables: ", paste(multi_sig, collapse = ", ")),
  paste0("timeROC AUC: ", paste(paste0(t_years, "y = ", sprintf("%.4f", roc_res$AUC)), collapse = "; ")),
  paste0("DCA models: ", paste(dca_names, collapse = ", ")),
  paste0("Conclusion: The nomogram incorporating ", length(multi_sig),
    " variables showed timeROC AUC of ", paste(sprintf("%.3f", roc_res$AUC), collapse = "/"),
    " at ", paste(t_years, collapse = "/"), "-year, demonstrating prognostic performance."),
  "",
  "限制与说明：",
  "1. 本分析基于单因素与多因素Cox回归筛选变量，纳入变量受样本量和事件数限制。",
  "2. 列线图校准曲线基于Bootstrap重抽样(B=1000)，AUC和DCA需外部验证。",
  "3. DCA分析基于特定时间点的阈值概率，不同阈值下净获益可能不同。",
  "4. 本结果为回顾性分析，需前瞻性队列验证后方可用于临床决策。"
))

cat("[Nomogram formula] ", deparse(nom_formula), "\n")
cat("[timeROC AUC] ", paste(paste0(t_years, "y=", sprintf("%.4f", roc_res$AUC)), collapse = "; "), "\n")
log_message("[Finish] Analysis completed")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")

risk_file     <- get_arg("--risk-file")
expr_file     <- get_arg("--expr")  # unused, kept for interface compatibility
train_id      <- get_arg("--train-id")
download_dir  <- get_arg("--download-dir", required = FALSE, default = file.path(dirname(outdir), "CNV_downloads"))

logsetup <- setup_logging(logdir, summary_dir, sn, "gistic")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

log_message("Script started: gistic")

cohort_id  <- train_id
force_rebuild <- TRUE
marker_url <- "https://api.gdc.cancer.gov/data/9bd7cbce-80f9-449e-8007-ddc9b1e89dfb"
marker_gz  <- file.path(download_dir, "snp6.na35.remap.hg38.subset.txt.gz")


# --- Helpers ---
pick_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
standardize_tcga_barcode <- function(x) { x <- trimws(gsub("\\.", "-", as.character(x))); x }
sample16 <- function(x) substr(standardize_tcga_barcode(x), 1, 16)
normalize_chr <- function(x) { x <- as.character(x); x <- gsub("^chr", "", x, ignore.case = TRUE); x }
is_primary_tumor_tcga <- function(x16) {
  out <- rep(TRUE, length(x16))
  idx <- grepl("^TCGA-", x16) & nchar(x16) >= 15
  out[idx] <- substr(x16[idx], 14, 15) == "01"
  out
}
write_tsv <- function(df, file, col_names = TRUE) {
  write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = col_names, na = "")
}
write_sample_list <- function(x, file) {
  x <- unique(as.character(x)); x <- x[!is.na(x) & x != ""]
  write.table(data.frame(Sample = x), file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# --- Output files ---
seg_all_file   <- file.path(outdir, paste0(cohort_id, "_all_matched.seg"))
seg_high_file  <- file.path(outdir, paste0(cohort_id, "_high_risk.seg"))
seg_low_file   <- file.path(outdir, paste0(cohort_id, "_low_risk.seg"))
seg_all_nh     <- file.path(outdir, paste0(cohort_id, "_all_matched.noheader.seg"))
seg_high_nh    <- file.path(outdir, paste0(cohort_id, "_high_risk.noheader.seg"))
seg_low_nh     <- file.path(outdir, paste0(cohort_id, "_low_risk.noheader.seg"))
marker_out     <- file.path(outdir, "markers_hg38_freqcnvFALSE.txt")
marker_nh      <- file.path(outdir, "markers_hg38_freqcnvFALSE.noheader.txt")
high_list      <- file.path(outdir, "high_risk_samples.txt")
low_list       <- file.path(outdir, "low_risk_samples.txt")
profile_map_f  <- file.path(outdir, "selected_profile_map.txt")

# --- Step 1: Read risk table ---
log_message("Step 1: reading risk table ...")
risk_raw <- read.csv(risk_file, check.names = FALSE, stringsAsFactors = FALSE)
log_message("Read input: ", risk_file, " -- ", nrow(risk_raw), " x ", ncol(risk_raw))
risk_col <- pick_first_existing_col(risk_raw, c("riskGroup", "risk_group", "RiskGroup", "group", "Group"))
if (is.na(risk_col)) save_note_and_stop("Cannot find risk group column")
sample_source_col <- pick_first_existing_col(risk_raw, c("sample", "Sample", "barcode", "Barcode"))
if (is.na(sample_source_col)) sample_source_col <- colnames(risk_raw)[1]

risk_df <- risk_raw %>%
  mutate(sample_raw = .data[[sample_source_col]],
         Sample_full = standardize_tcga_barcode(sample_raw),
         Sample16 = sample16(sample_raw),
         riskGroup_raw = as.character(.data[[risk_col]]),
         riskGroup = case_when(
           grepl("^high$", riskGroup_raw, ignore.case = TRUE) ~ "High",
           grepl("^low$",  riskGroup_raw, ignore.case = TRUE) ~ "Low",
           grepl("^high",  riskGroup_raw, ignore.case = TRUE) ~ "High",
           grepl("^low",   riskGroup_raw, ignore.case = TRUE) ~ "Low",
           TRUE ~ riskGroup_raw)) %>%
  filter(!is.na(Sample16), Sample16 != "", riskGroup %in% c("High", "Low")) %>%
  distinct(Sample16, .keep_all = TRUE)

if (nrow(risk_df) == 0) save_note_and_stop("No valid High/Low samples found")
log_message("Risk samples: ", nrow(risk_df), " (High: ", sum(risk_df$riskGroup == "High"),
     ", Low: ", sum(risk_df$riskGroup == "Low"), ")")

# --- Step 2: Download CNV segments ---
seg_outputs_needed <- c(seg_all_file, seg_high_file, seg_low_file, seg_all_nh, seg_high_nh, seg_low_nh, high_list, low_list, profile_map_f)
if (!force_rebuild && all(file.exists(seg_outputs_needed))) {
  log_message("Step 2 skipped: seg outputs exist")
} else {
  log_message("Step 2: downloading CNV segments ...")
  query_seg <- GDCquery(project = cohort_id, data.category = "Copy Number Variation",
    data.type = "Masked Copy Number Segment", sample.type = "Primary Tumor")
  GDCdownload(query = query_seg, directory = download_dir, files.per.chunk = 20)
  seg_obj <- tryCatch(GDCprepare(query_seg, directory = download_dir, summarizedExperiment = FALSE),
                       error = function(e) GDCprepare(query_seg, directory = download_dir))
  seg_raw <- as.data.frame(seg_obj, check.names = FALSE, stringsAsFactors = FALSE)

  sbc <- pick_first_existing_col(seg_raw, c("Sample", "sample", "Tumor_Sample_Barcode"))
  pic <- pick_first_existing_col(seg_raw, c("GDC_Aliquot", "Aliquot"))
  chr_c <- pick_first_existing_col(seg_raw, c("Chromosome", "chromosome", "chr"))
  sc <- pick_first_existing_col(seg_raw, c("Start", "start"))
  ec <- pick_first_existing_col(seg_raw, c("End", "end", "Stop"))
  nc <- pick_first_existing_col(seg_raw, c("Num_Probes", "NumMarkers", "num_probes"))
  smc <- pick_first_existing_col(seg_raw, c("Segment_Mean", "Seg.CN", "seg.mean"))
  if (any(is.na(c(sbc, pic, chr_c, sc, ec, nc, smc)))) save_note_and_stop("Cannot identify segment columns")

  seg_df <- seg_raw %>%
    transmute(ProfileID = as.character(.data[[pic]]),
              Sample_barcode = standardize_tcga_barcode(.data[[sbc]]),
              Sample16 = sample16(.data[[sbc]]),
              Chromosome = normalize_chr(.data[[chr_c]]),
              Start = suppressWarnings(as.integer(.data[[sc]])),
              End = suppressWarnings(as.integer(.data[[ec]])),
              NumMarkers = suppressWarnings(as.integer(.data[[nc]])),
              Seg.CN = suppressWarnings(as.numeric(.data[[smc]]))) %>%
    filter(!is.na(ProfileID), ProfileID != "", !is.na(Sample_barcode), Sample_barcode != "",
           is_primary_tumor_tcga(Sample16), Chromosome %in% as.character(1:22),
           !is.na(Start), !is.na(End), !is.na(Seg.CN), Start < End) %>%
    mutate(NumMarkers = ifelse(is.na(NumMarkers) | NumMarkers < 1, 1L, NumMarkers)) %>%
    arrange(ProfileID, as.integer(Chromosome), Start, End, desc(NumMarkers)) %>%
    group_by(ProfileID, Chromosome, Start, End) %>% slice(1) %>% ungroup()
  if (nrow(seg_df) == 0) save_note_and_stop("No valid segment rows after filtering")

  profile_stats <- seg_df %>% count(Sample16, Sample_barcode, ProfileID, name = "n_segments") %>% arrange(Sample16, desc(n_segments))
  rep_profiles <- profile_stats %>% group_by(Sample16) %>% arrange(desc(n_segments), ProfileID, .by_group = TRUE) %>% slice(1) %>% ungroup()

  risk_key <- risk_df %>% distinct(Sample16, riskGroup)
  selected_profiles <- risk_key %>% inner_join(rep_profiles, by = "Sample16") %>% distinct()
  if (nrow(selected_profiles) == 0) save_note_and_stop("No matched CNV samples")

  profile_map <- selected_profiles %>%
    mutate(GISTIC_Sample = paste0(Sample16, "__", substr(ProfileID, 1, 12))) %>%
    select(GISTIC_Sample, Sample16, Sample_barcode, ProfileID, riskGroup, n_segments)
  write_tsv(profile_map, profile_map_f)

  seg_matched <- seg_df %>% inner_join(profile_map %>% select(GISTIC_Sample, Sample16, Sample_barcode, ProfileID, riskGroup),
    by = c("ProfileID", "Sample16", "Sample_barcode"))

  seg_out <- seg_matched %>%
    transmute(Sample = GISTIC_Sample, Chromosome = Chromosome, Start = Start,
              End = End, NumMarkers = NumMarkers, Seg.CN = Seg.CN, riskGroup = riskGroup,
              ProfileID = ProfileID, Sample16 = Sample16) %>%
    arrange(Sample, as.integer(Chromosome), Start, End)

  # Overlap detection
  overlap_df <- seg_out %>% group_by(Sample, Chromosome) %>% arrange(Start, End, .by_group = TRUE) %>%
    mutate(prev_end = lag(End)) %>% filter(!is.na(prev_end) & Start <= prev_end) %>% ungroup()
  if (nrow(overlap_df) > 0) {
    bad_profiles <- overlap_df %>% distinct(Sample, ProfileID, Sample16)
    seg_out <- seg_out %>% anti_join(bad_profiles %>% select(Sample), by = "Sample")
    profile_map <- profile_map %>% anti_join(bad_profiles %>% select(Sample) %>% rename(GISTIC_Sample = Sample), by = "GISTIC_Sample")
    overlap_df2 <- seg_out %>% group_by(Sample, Chromosome) %>% arrange(Start, End, .by_group = TRUE) %>%
      mutate(prev_end = lag(End)) %>% filter(!is.na(prev_end) & Start <= prev_end) %>% ungroup()
    if (nrow(overlap_df2) > 0) save_note_and_stop("Overlapping segments persist after dropping problematic profiles")
    log_message("Dropped ", nrow(bad_profiles), " profiles with overlaps")
  }
  write_tsv(profile_map, profile_map_f)

  seg_all  <- seg_out %>% select(Sample, Chromosome, Start, End, NumMarkers, Seg.CN)
  seg_high <- seg_out %>% filter(riskGroup == "High") %>% select(Sample, Chromosome, Start, End, NumMarkers, Seg.CN)
  seg_low  <- seg_out %>% filter(riskGroup == "Low")  %>% select(Sample, Chromosome, Start, End, NumMarkers, Seg.CN)

  write_tsv(seg_all, seg_all_file); write_tsv(seg_all, seg_all_nh, FALSE)
  write_tsv(seg_high, seg_high_file); write_tsv(seg_high, seg_high_nh, FALSE)
  write_tsv(seg_low, seg_low_file); write_tsv(seg_low, seg_low_nh, FALSE)
  write_sample_list(unique(seg_high$Sample), high_list)
  write_sample_list(unique(seg_low$Sample), low_list)
  log_message("Segments -- All: ", nrow(seg_all), " High: ", nrow(seg_high), " Low: ", nrow(seg_low))
}

# --- Step 3: Marker file ---
log_message("Step 3: preparing marker file ...")
if (!file.exists(marker_gz)) download.file(marker_url, destfile = marker_gz, mode = "wb")
marker_raw <- read.delim(gzfile(marker_gz), header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")
mid <- pick_first_existing_col(marker_raw, c("probeid", "Probe.Name", "probe_set_id", "ID", "Name"))
mchr <- pick_first_existing_col(marker_raw, c("chr", "Chromosome", "chromosome"))
mpos <- pick_first_existing_col(marker_raw, c("pos", "Position", "position", "Start"))
mfreq <- pick_first_existing_col(marker_raw, c("freqcnv", "FreqCNV", "freq_cnv"))
if (any(is.na(c(mid, mchr, mpos, mfreq)))) save_note_and_stop("Cannot identify marker columns")

marker_df <- marker_raw %>%
  mutate(MarkerName = as.character(.data[[mid]]),
         Chromosome = normalize_chr(.data[[mchr]]),
         Position = suppressWarnings(as.integer(.data[[mpos]])),
         freq_keep = tolower(as.character(.data[[mfreq]])) %in% c("false", "f", "0")) %>%
  filter(freq_keep, Chromosome %in% as.character(1:22),
         !is.na(MarkerName), MarkerName != "", !is.na(Position), Position > 0) %>%
  select(MarkerName, Chromosome, Position) %>% distinct()
write_tsv(marker_df, marker_out); write_tsv(marker_df, marker_nh, FALSE)
log_message("Markers: ", nrow(marker_df))

log_message("GISTIC input files generated")
log_message("  High risk seg: ", seg_high_nh)
log_message("  Low risk seg:  ", seg_low_nh)
log_message("  Marker file:   ", marker_out)

log_message("Script completed successfully")

n_high <- if (exists("seg_high")) nrow(seg_high) else 0
n_low  <- if (exists("seg_low")) nrow(seg_low) else 0
write_summary(c(
  paste0("分析对象: ", cohort_id),
  paste0("风险分组样本数 - High: ", sum(risk_df$riskGroup == "High"), ", Low: ", sum(risk_df$riskGroup == "Low")),
  paste0("GISTIC2输入文件已生成: SEG文件 (All: ", if (exists("seg_all")) nrow(seg_all) else 0, " segments, High: ", n_high, ", Low: ", n_low, "), Marker文件 (", nrow(marker_df), " probes)"),
  paste0("CNV segments are ready for GISTIC2 analysis, covering high-risk and low-risk groups separately."),
  paste0("限制与说明: GISTIC2需另行运行; segments基于TCGA Masked Copy Number Segment; marker probes过滤了常见CNV区域.")
))

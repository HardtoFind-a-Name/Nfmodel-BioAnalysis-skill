rm(list = ls())
gc()

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_raw_import_utils.R"))

raw_source_manifest <- get_arg("--raw-source-manifest")
outdir <- get_arg("--outdir")
decompress_max_depth <- as.integer(get_arg("--decompress-max-depth", required = FALSE, default = "5"))
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_unpack_raw_umi")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_unpack_raw_umi")
if (!file.exists(raw_source_manifest)) save_note_and_stop(paste0("Raw source manifest not found: ", raw_source_manifest))
manifest <- read.csv(raw_source_manifest, stringsAsFactors = FALSE, check.names = FALSE)
raw_type <- manifest$raw_source_type[manifest$role == "raw_source"][1]
source_path <- manifest$source_path[manifest$role == "raw_source"][1]
annotation_path <- manifest$source_path[manifest$role == "annotation"][1]
if (is.na(raw_type) || is.na(source_path) || !file.exists(source_path)) save_note_and_stop("Raw source manifest is invalid.")
if (is.na(annotation_path) || !file.exists(annotation_path)) save_note_and_stop("Annotation file is missing in raw source manifest.")

if (raw_type != "umi_rds") {
  log_message("Raw source type is", raw_type, "; UMI unpack step passes manifest through unchanged")
  write.csv(manifest, file.path(outdir, "01.raw_source_after_umi.csv"), row.names = FALSE)
  write_summary(c(
    "scRNA raw UMI unpack skipped.",
    paste0("Raw source type: ", raw_type),
    "The manifest was passed through unchanged for downstream 10X handling or mainline construction.",
    "",
    "限制与说明",
    "Only umi_rds sources are recursively gzip-unpacked in this step."
  ))
  quit(save = "no", status = 0)
}

current <- normalizePath(source_path, winslash = "/", mustWork = TRUE)
steps <- data.frame(depth = integer(), input_file = character(), output_file = character(), stringsAsFactors = FALSE)

for (depth in seq_len(decompress_max_depth)) {
  if (!has_gzip_magic(current)) break
  next_path <- file.path(outdir, basename(strip_gzip_suffix(current, depth)))
  log_message("Decompress gzip depth", depth, ":", current, "->", next_path)
  con_in <- gzfile(current, "rb")
  con_out <- file(next_path, "wb")
  ok <- tryCatch({
    repeat {
      buf <- readBin(con_in, what = "raw", n = 1024 * 1024)
      if (length(buf) == 0) break
      writeBin(buf, con_out)
    }
    TRUE
  }, error = function(e) {
    log_message("[致命错误] gzip decompression failed:", conditionMessage(e))
    FALSE
  })
  close(con_in)
  close(con_out)
  if (!ok) save_note_and_stop("gzip decompression failed.")
  steps <- rbind(steps, data.frame(depth = depth, input_file = current, output_file = next_path))
  current <- normalizePath(next_path, winslash = "/", mustWork = TRUE)
}

if (has_gzip_magic(current)) save_note_and_stop(paste0("Raw UMI file is still gzip after max depth: ", decompress_max_depth))

read_ok <- tryCatch({
  obj <- readRDS(current)
  inherits(obj, c("matrix", "dgCMatrix", "data.frame"))
}, error = function(e) {
  log_message("[致命错误] Final raw UMI file is not readable as RDS:", conditionMessage(e))
  FALSE
})
if (!read_ok) save_note_and_stop("Final raw UMI file is not a readable matrix/data.frame RDS.")

if (nrow(steps) == 0) steps <- data.frame(depth = 0, input_file = source_path, output_file = current)
write.csv(steps, file.path(outdir, "02.raw_umi_unpack_steps.csv"), row.names = FALSE)
manifest$source_path[manifest$role == "raw_source"] <- normalizePath(current, winslash = "/", mustWork = TRUE)
manifest$file_name[manifest$role == "raw_source"] <- basename(current)
manifest$note[manifest$role == "raw_source"] <- "UMI RDS recursively unpacked and validated"
write.csv(manifest, file.path(outdir, "01.raw_source_after_umi.csv"), row.names = FALSE)

log_message("Final raw UMI RDS:", current)
log_message("Script completed: scrna_unpack_raw_umi")
write_summary(c(
  "scRNA raw UMI unpack completed.",
  paste0("Final raw UMI RDS: ", basename(current)),
  paste0("Decompression steps: ", max(steps$depth)),
  "The final raw UMI file was validated with readRDS().",
  "",
  "限制与说明",
  "This step only processes umi_rds raw sources; other raw source types are passed through unchanged."
))

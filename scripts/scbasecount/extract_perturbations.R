#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
tcells_dir <- if (length(args) >= 1) args[1] else "../scbasecount - data/tcells"
out_csv <- if (length(args) >= 2) args[2] else file.path(tcells_dir, "perturbation_summary.csv")

files <- list.files(tcells_dir, pattern = "_tcell\\.rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("No *_tcell.rds files found in ", tcells_dir)
}

csv_escape <- function(x) {
  if (is.na(x) || length(x) == 0) return("")
  x <- as.character(x)
  x <- gsub("\"", "\"\"", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

get_first_value <- function(vec) {
  if (is.null(vec)) return(NA_character_)
  vec <- unique(as.character(vec))
  vec <- vec[!is.na(vec) & vec != ""]
  if (length(vec) == 0) return(NA_character_)
  vec[[1]]
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

processed_files <- character(0)
if (file.exists(out_csv)) {
  existing <- tryCatch(read.csv(out_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(existing) && "file" %in% colnames(existing)) {
    processed_files <- unique(existing$file)
  }
}

write_header <- !file.exists(out_csv) || length(processed_files) == 0
con <- file(out_csv, open = "a", encoding = "UTF-8")
if (write_header) {
  writeLines("file,perturbation_raw,disease,tissue,organism", con)
}

total <- length(files)
for (i in seq_along(files)) {
  f <- files[[i]]
  if (basename(f) %in% processed_files) {
    next
  }
  pert <- NA_character_
  disease <- NA_character_
  tissue <- NA_character_
  organism <- NA_character_

  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.null(obj)) {
    md <- tryCatch(obj[[]], error = function(e) NULL)
    if (!is.null(md)) {
      if ("perturbation" %in% colnames(md)) {
        pert <- get_first_value(md$perturbation)
      }
      if ("disease" %in% colnames(md)) {
        disease <- get_first_value(md$disease)
      }
      if ("tissue" %in% colnames(md)) {
        tissue <- get_first_value(md$tissue)
      }
      if ("organism" %in% colnames(md)) {
        organism <- get_first_value(md$organism)
      }
    }
  }

  line <- paste(
    basename(f),
    csv_escape(pert),
    csv_escape(disease),
    csv_escape(tissue),
    csv_escape(organism),
    sep = ","
  )
  writeLines(line, con)

  rm(obj)
  gc()

  if (i %% 50 == 0) {
    message("Processed ", i, " / ", total, " files")
  }
}

close(con)
message("Wrote perturbation summary to: ", out_csv)

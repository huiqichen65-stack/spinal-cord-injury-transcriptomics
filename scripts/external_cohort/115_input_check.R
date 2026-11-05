source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)

project_root <- PROJECT_ROOT
input_dir <- file.path(project_root, "data", "GSE281205")
output_dir <- file.path(project_root, "results", "115_GSE281205_input_QC")

sample_plan <- data.frame(gsm = paste0("GSM", 8614662:8614667), intended_group = c(rep("Young-14dpi",
  3L), rep("Young-28dpi", 3L)), intended_dpi = c(rep(14L, 3L), rep(28L, 3L)), intended_replicate = rep(1:3,
  2L), stringsAsFactors = FALSE)

marker_symbols <- c("Clec7a", "P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Csf1r", "Aif1", "Tyrobp", "Trem2",
  "Apoe", "Lpl", "Cst7", "Gpnmb", "Itgam", "Lyz2", "Ccr2")

WRITE_PER_BARCODE_QC <- FALSE

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir, "\nSave the 21 GEO files there, without changing their names.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(output_dir, "00_115_run_log_v01.log")
cat("", file = log_path)
log_msg <- function(...) {
  x <- paste0(..., collapse = "")
  cat(x, "\n", sep = "")
  cat(x, "\n", file = log_path, append = TRUE, sep = "")
}

save_csv <- function(df, filename) {
  p <- file.path(output_dir, filename)
  utils::write.csv(df, p, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(p)
}

one_value <- function(lines, prefix) {
  hit <- lines[startsWith(lines, prefix)]
  if (!length(hit))
    return(NA_character_)
  trimws(substring(hit[1L], nchar(prefix) + 1L))
}

read_soft_samples <- function(path) {
  if (!file.exists(path))
    stop("Missing SOFT file: ", path)
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  z <- readLines(con, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE[[:space:]]*=[[:space:]]*GSM[0-9]+", z)
  if (!length(starts))
    stop("SOFT file has no ^SAMPLE sections.")
  ends <- c(starts[-1L] - 1L, length(z))
  rows <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    block <- z[starts[i]:ends[i]]
    chr_lines <- block[startsWith(block, "!Sample_characteristics_ch1 = ")]
    chr_values <- substring(chr_lines, nchar("!Sample_characteristics_ch1 = ") + 1L)
    treatment_lines <- grep("^treatment[[:space:]]*:", chr_values, ignore.case = TRUE, value = TRUE)
    strain_lines <- grep("^strain[[:space:]]*:", chr_values, ignore.case = TRUE, value = TRUE)
    tissue_lines <- grep("^tissue[[:space:]]*:", chr_values, ignore.case = TRUE, value = TRUE)
    rel <- block[startsWith(block, "!Sample_relation = BioSample:")]
    rows[[i]] <- data.frame(gsm = trimws(sub("^\\^SAMPLE[[:space:]]*=[[:space:]]*", "", block[1])),
      geo_title = one_value(block, "!Sample_title = "), geo_treatment = if (length(treatment_lines)) {
        trimws(sub("^[^:]+:[[:space:]]*", "", treatment_lines[1]))
      } else NA_character_, geo_strain = if (length(strain_lines)) {
        trimws(sub("^[^:]+:[[:space:]]*", "", strain_lines[1]))
      } else NA_character_, geo_tissue = if (length(tissue_lines)) {
        trimws(sub("^[^:]+:[[:space:]]*", "", tissue_lines[1]))
      } else NA_character_, geo_platform = one_value(block, "!Sample_platform_id = "), biosample_url = if (length(rel)) {
        trimws(sub("^!Sample_relation = BioSample:[[:space:]]*", "", rel[1]))
      } else NA_character_, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

read_series_matrix_metadata <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  z <- readLines(con, warn = FALSE, encoding = "UTF-8")
  fields <- list()
  for (key in c("!Sample_geo_accession", "!Sample_title")) {
    hit <- z[startsWith(z, paste0(key, "\t"))]
    if (length(hit) == 1L) {
      fields[[key]] <- gsub("^\"|\"$", "", strsplit(hit, "\t", fixed = TRUE)[[1]][-1])
    }
  }
  if (is.null(fields[["!Sample_geo_accession"]])) {
    stop("Series matrix has no !Sample_geo_accession: ", basename(path))
  }
  gsms <- fields[["!Sample_geo_accession"]]
  titles <- fields[["!Sample_title"]]
  if (is.null(titles) || length(titles) != length(gsms)) {
    titles <- rep(NA_character_, length(gsms))
  }
  data.frame(gsm = gsms, series_title = titles, series_file = basename(path), stringsAsFactors = FALSE)
}

read_tsv_gz <- function(path, min_cols = 1L) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  x <- utils::read.delim(con, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE,
    check.names = FALSE, fileEncoding = "UTF-8", fill = FALSE)
  if (ncol(x) < min_cols)
    stop("Too few TSV columns: ", basename(path))
  x
}

read_mtx_header <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  first <- readLines(con, n = 1L, warn = FALSE)
  if (length(first) != 1L || !grepl("^%%MatrixMarket[[:space:]]+matrix[[:space:]]+coordinate", first,
    ignore.case = TRUE)) {
    stop("Not a coordinate MatrixMarket file: ", basename(path))
  }
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L)
      stop("MatrixMarket dimensions not found: ", basename(path))
    if (startsWith(line, "%") || !nzchar(trimws(line)))
      next
    vals <- scan(text = line, what = numeric(), quiet = TRUE)
    if (length(vals) != 3L || any(!is.finite(vals)) || any(vals < 0)) {
      stop("Invalid MatrixMarket dimensions: ", basename(path))
    }
    return(setNames(as.numeric(vals), c("n_genes", "n_barcodes", "n_triplets")))
  }
}

safe_median <- function(x) {
  if (!length(x) || all(is.na(x)))
    return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

files_all <- list.files(input_dir, full.names = TRUE, recursive = FALSE)
meta_expected <- c("GSE281205_family.soft.gz", "GSE281205-GPL21103_series_matrix.txt.gz", "GSE281205-GPL24247_series_matrix.txt.gz")

manifest_rows <- list()
expected_file_map <- list()
add_file <- function(gsm, file_kind, expected_name, resolved) {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(gsm = gsm, file_kind = file_kind, expected_name_or_pattern = expected_name,
    found_count = length(resolved), file_name = if (length(resolved) == 1L)
      basename(resolved) else NA_character_, size_bytes = if (length(resolved) == 1L)
      as.numeric(file.info(resolved)$size) else NA_real_, status = if (length(resolved) == 1L)
      "FOUND" else if (!length(resolved))
      "MISSING" else "AMBIGUOUS_MULTIPLE", stringsAsFactors = FALSE)
}

for (nm in meta_expected) {
  fp <- file.path(input_dir, nm)
  add_file("GSE281205", "GEO_METADATA", nm, if (file.exists(fp))
    fp else character(0))
}

for (i in seq_len(nrow(sample_plan))) {
  gsm <- sample_plan$gsm[i]
  group <- sample_plan$intended_group[i]
  repn <- sample_plan$intended_replicate[i]
  prefix <- paste0(gsm, "_", group, "_", repn, "_")
  required <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  one <- setNames(vector("list", length(required)), c("barcodes", "features", "matrix"))
  for (kind in names(one)) {
    suffix <- switch(kind, barcodes = "barcodes.tsv.gz", features = "features.tsv.gz", matrix = "matrix.mtx.gz")
    strict_basename <- paste0(prefix, suffix)

    strict <- file.path(input_dir, strict_basename)
    if (file.exists(strict)) {
      candidates <- strict
    } else {
      patt <- paste0("^", gsm, "_.*_", gsub("\\.", "\\\\.", suffix), "$")
      candidates <- files_all[grepl(patt, basename(files_all), ignore.case = TRUE)]
    }
    one[[kind]] <- if (length(candidates) == 1L)
      candidates else NA_character_
    add_file(gsm, kind, strict_basename, candidates)
  }
  expected_file_map[[gsm]] <- one
}
manifest <- do.call(rbind, manifest_rows)
save_csv(manifest, "01_115_file_manifest_v01.csv")

if (any(manifest$status != "FOUND" & manifest$gsm == "GSE281205")) {
  stop("One or more GEO metadata files are missing. See 01_115_file_manifest_v01.csv")
}

soft <- read_soft_samples(file.path(input_dir, meta_expected[1]))
series <- do.call(rbind, lapply(meta_expected[-1], function(nm) {
  read_series_matrix_metadata(file.path(input_dir, nm))
}))
if (anyDuplicated(soft$gsm)) stop("Duplicated GSM identifiers in SOFT.")
if (anyDuplicated(series$gsm)) stop("Duplicated GSM identifiers in series matrices.")
geo <- merge(soft, series, by = "gsm", all = TRUE, sort = FALSE)
geo <- merge(sample_plan, geo, by = "gsm", all.x = TRUE, sort = FALSE)
geo <- geo[match(sample_plan$gsm, geo$gsm), , drop = FALSE]

match_text <- function(z, phrase) !is.na(z) & grepl(phrase, z, ignore.case = TRUE, fixed = TRUE)
geo$title_matches_intended <- mapply(match_text, geo$geo_title, geo$intended_group)
geo$series_title_matches_intended <- mapply(match_text, geo$series_title, geo$intended_group)
geo$treatment_matches_intended <- mapply(match_text, geo$geo_treatment, geo$intended_group)
geo$treatment_title_conflict <- geo$title_matches_intended & !geo$treatment_matches_intended
geo$metadata_status <- ifelse(is.na(geo$geo_title) | is.na(geo$series_file), "MISSING_GEO_SAMPLE_METADATA",
  ifelse(geo$treatment_title_conflict, "REVIEW_TITLE_VS_TREATMENT_CONFLICT", ifelse(geo$title_matches_intended &
    geo$series_title_matches_intended & geo$treatment_matches_intended, "CONSISTENT", "REVIEW_OTHER_MISMATCH")))
save_csv(geo, "02_115_GEO_sample_metadata_audit_v01.csv")

if (any(geo$treatment_title_conflict)) {
  log_msg("WARNING: GEO title/treatment conflicts detected:")
  for (i in which(geo$treatment_title_conflict)) {
  }
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package Matrix is required. Install with install.packages('Matrix').")
}

qc_rows <- list()
marker_rows <- list()
barcode_dup_rows <- list()
all_barcodes <- character(0)
per_barcode_paths <- character(0)

for (i in seq_len(nrow(sample_plan))) {
  gsm <- sample_plan$gsm[i]
  plan <- sample_plan[i, ]
  ff <- expected_file_map[[gsm]]

  qc <- data.frame(gsm = gsm, intended_group = plan$intended_group, intended_dpi = plan$intended_dpi,
    replicate = plan$intended_replicate, matrix_status = "NOT_RUN", matrix_n_features = NA_integer_,
    matrix_n_barcodes = NA_integer_, matrix_n_triplets_header = NA_real_, n_features_tsv = NA_integer_,
    n_barcodes_tsv = NA_integer_, n_duplicate_gene_ids = NA_integer_, n_duplicate_gene_symbols = NA_integer_,
    n_duplicate_barcodes_within_sample = NA_integer_, n_non_Gene_Expression_features = NA_integer_,
    matrix_sums_match = NA, nonnegative_integer_counts = NA, total_counts = NA_real_, n_barcodes_with_nonzero_counts = NA_integer_,
    median_counts_per_barcode = NA_real_, median_features_per_barcode = NA_real_, median_percent_mito = NA_real_,
    Clec7a_feature_present = NA, n_Clec7a_detected_barcodes = NA_integer_, pct_Clec7a_detected_barcodes = NA_real_,
    error_message = NA_character_, stringsAsFactors = FALSE)

  tryCatch({
    if (any(is.na(unlist(ff))))
      stop("A required 10x file is missing or ambiguous.")
    barcodes <- read_tsv_gz(ff$barcodes, min_cols = 1L)
    features <- read_tsv_gz(ff$features, min_cols = 2L)
    hdr <- read_mtx_header(ff$matrix)

    bcodes <- trimws(as.character(barcodes[[1]]))
    gene_ids <- trimws(as.character(features[[1]]))
    gene_symbols <- trimws(as.character(features[[2]]))
    feature_type <- if (ncol(features) >= 3L) {
      trimws(as.character(features[[3]]))
    } else rep("Gene Expression", nrow(features))

    qc$matrix_n_features <- as.integer(hdr["n_genes"])
    qc$matrix_n_barcodes <- as.integer(hdr["n_barcodes"])
    qc$matrix_n_triplets_header <- as.numeric(hdr["n_triplets"])
    qc$n_features_tsv <- nrow(features)
    qc$n_barcodes_tsv <- nrow(barcodes)
    qc$n_duplicate_gene_ids <- sum(duplicated(gene_ids))
    qc$n_duplicate_gene_symbols <- sum(duplicated(gene_symbols))
    qc$n_duplicate_barcodes_within_sample <- sum(duplicated(bcodes))
    qc$n_non_Gene_Expression_features <- sum(feature_type != "Gene Expression")

    if (anyNA(bcodes) || any(!nzchar(bcodes)))
      stop("Empty barcode(s) in TSV.")
    if (anyNA(gene_ids) || any(!nzchar(gene_ids)))
      stop("Empty gene ID(s) in TSV.")
    if (anyNA(gene_symbols) || any(!nzchar(gene_symbols))) {
      stop("Empty gene symbol(s) in TSV.")
    }
    if (anyDuplicated(bcodes))
      stop("Duplicated barcodes within sample.")
    if (as.integer(hdr["n_genes"]) != nrow(features) || as.integer(hdr["n_barcodes"]) != nrow(barcodes)) {
      stop("MatrixMarket header dimensions do not match features/barcodes TSV.")
    }

    con <- gzfile(ff$matrix, open = "rt")
    mat <- tryCatch(Matrix::readMM(con), finally = close(con))
    if (!identical(as.integer(dim(mat)), as.integer(c(nrow(features), nrow(barcodes))))) {
      stop("Loaded matrix dimensions are inconsistent with TSV files.")
    }
    values <- mat@x
    if (length(values) != as.numeric(hdr["n_triplets"])) {
      stop("Count of stored matrix entries differs from MatrixMarket header.")
    }
    if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-08)) {
      stop("The count matrix has negative, non-finite, or noninteger values.")
    }
    qc$nonnegative_integer_counts <- TRUE
    qc$matrix_sums_match <- TRUE

    ncount <- as.numeric(Matrix::colSums(mat))
    nfeature <- as.numeric(Matrix::colSums(mat > 0))
    mt_rows <- grepl("^mt-", gene_symbols, ignore.case = TRUE) & feature_type == "Gene Expression"
    mito <- if (any(mt_rows)) {
      as.numeric(Matrix::colSums(mat[mt_rows, , drop = FALSE]))
    } else rep(NA_real_, ncol(mat))
    pctmt <- ifelse(ncount > 0 & !is.na(mito), mito/ncount * 100, NA_real_)

    qc$total_counts <- sum(ncount)
    qc$n_barcodes_with_nonzero_counts <- sum(ncount > 0)
    qc$median_counts_per_barcode <- safe_median(ncount)
    qc$median_features_per_barcode <- safe_median(nfeature)
    qc$median_percent_mito <- safe_median(pctmt)

    for (marker in marker_symbols) {
      ix <- which(tolower(gene_symbols) == tolower(marker) & feature_type == "Gene Expression")
      x <- if (length(ix)) {
        as.numeric(Matrix::colSums(mat[ix, , drop = FALSE]))
      } else numeric(0)
      marker_rows[[length(marker_rows) + 1L]] <- data.frame(gsm = gsm, intended_dpi = plan$intended_dpi,
        marker = marker, present_in_features = length(ix) > 0L, feature_rows = length(ix), n_barcodes_with_detected_marker = if (length(ix))
          sum(x > 0) else NA_integer_, percent_barcodes_with_detected_marker = if (length(ix))
          round(100 * mean(x > 0), 4) else NA_real_, total_marker_counts = if (length(ix))
          sum(x) else NA_real_, interpretation = "Raw detection among all barcodes; NOT a cell-state assignment",
        stringsAsFactors = FALSE)
      if (marker == "Clec7a") {
        qc$Clec7a_feature_present <- length(ix) > 0L
        qc$n_Clec7a_detected_barcodes <- if (length(ix))
          sum(x > 0) else NA_integer_
        qc$pct_Clec7a_detected_barcodes <- if (length(ix))
          round(100 * mean(x > 0), 4) else NA_real_
      }
    }

    if (WRITE_PER_BARCODE_QC) {
      fn <- paste0("06_115_per_barcode_QC_", gsm, "_v01.csv")
      save_csv(data.frame(gsm = gsm, barcode = bcodes, nCount_RNA = ncount, nFeature_RNA = nfeature,
        percent_mito = pctmt, stringsAsFactors = FALSE), fn)
      per_barcode_paths <- c(per_barcode_paths, fn)
    }

    barcode_dup_rows[[length(barcode_dup_rows) + 1L]] <- data.frame(gsm = gsm, n_barcodes = length(bcodes),
      n_duplicate_barcodes_within_library = sum(duplicated(bcodes)), stringsAsFactors = FALSE)
    all_barcodes <- c(all_barcodes, paste(gsm, bcodes, sep = "_"))
    qc$matrix_status <- "PASS_TECHNICAL_QC"
    rm(mat, features, barcodes, values, ncount, nfeature, mito, pctmt)
    invisible(gc(verbose = FALSE))

  }, error = function(e) {
    qc$matrix_status <- "FAIL_TECHNICAL_QC"
    qc$error_message <- conditionMessage(e)
    log_msg("  FAIL: ", conditionMessage(e))
    invisible(gc(verbose = FALSE))
  })
  qc_rows[[length(qc_rows) + 1L]] <- qc
}

qc_summary <- do.call(rbind, qc_rows)
save_csv(qc_summary, "03_115_sample_matrix_QC_v01.csv")
if (length(marker_rows)) {
  save_csv(do.call(rbind, marker_rows), "04_115_raw_marker_detection_v01.csv")
}
if (length(barcode_dup_rows)) {
  save_csv(do.call(rbind, barcode_dup_rows), "05_115_barcode_integrity_v01.csv")
}

status <- if (any(qc_summary$matrix_status != "PASS_TECHNICAL_QC")) {
  "NOT_READY_TECHNICAL_MATRIX_FAILURE"
} else if (any(geo$metadata_status != "CONSISTENT")) {
  "TECHNICALLY_READY_BUT_METADATA_REVIEW_REQUIRED"
} else {
  "TECHNICALLY_READY_FOR_NEXT_STAGE_ANNOTATION"
}

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(output_dir, "07_115_REPORT_v01.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "08_115_sessionInfo_v01.txt"), useBytes = TRUE)

outputs <- sort(list.files(output_dir, full.names = TRUE))
index <- data.frame(output_file = basename(outputs), size_bytes = as.numeric(file.info(outputs)$size),
  stringsAsFactors = FALSE)
save_csv(index, "09_115_OUTPUT_INDEX_v01.csv")

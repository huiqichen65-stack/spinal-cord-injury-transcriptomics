rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 1200)
set.seed(114)

ROOT <- PROJECT_ROOT

STEP60_DIR <- file.path(ROOT, "results", "60_GSE195783_Visium_import_coordinate_reconstruction")

EXTRACT_DIR <- file.path(STEP60_DIR, "00_extracted_processed_data_v05")

SAMPLE_MANIFEST_FILE <- file.path(STEP60_DIR, "00_60_GSE195783_final_sample_manifest_v05.csv")

STEP68_DIR <- file.path(ROOT, "results", "68_GSE195783_spatial_section_map_finalize")

SECTION_MANIFEST_FILE <- file.path(STEP68_DIR, "01_68_GSE195783_spatial_section_map_section_manifest_v01.csv")

SPOT_SECTION_MAP_FILE <- file.path(STEP68_DIR, "02_68_GSE195783_spatial_section_map_spot_to_section_v01.csv.gz")

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

RCTD_FILE <- file.path(STEP65_DIR, "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv.gz")

OUT_DIR <- file.path(ROOT, "results", "114_GSE195783_uninjured_section_validation")

VERSION <- "v01"
STEP <- 114L

EXPECTED_SAMPLES <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", 
    "28dpi_rep1", "28dpi_rep2")

EXPECTED_SECTIONS_PER_SAMPLE <- 4L
EXPECTED_TOTAL_SECTIONS <- 32L

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

CANDIDATE_SECTIONS <- c("7dpi_rep2__lower_left", "14dpi_rep1__upper_right", "14dpi_rep2__lower_right", 
    "28dpi_rep1__upper_right", "28dpi_rep2__lower_right")

GENE_SETS <- list(Normal_neuronal = c("Snap25", "Syt1", "Rbfox3", "Tubb3", "Map2", "Slc17a6", "Slc17a7", 
    "Gad1", "Gad2"), Normal_myelin_oligodendrocyte = c("Plp1", "Mbp", "Mog", "Mag", "Cnp", "Mobp", "Mal", 
    "Cldn11", "Sox10"), Injury_myeloid = c("Spp1", "Gpnmb", "Lgals3", "Ctsd", "Lpl", "Apoe", "Trem2", 
    "Cst7", "C1qa", "C1qb", "C1qc", "Tyrobp", "Ctss"), Injury_fibrotic_ECM = c("Col1a1", "Col1a2", "Col3a1", 
    "Fn1", "Dcn", "Lum", "Igfbp6", "Mgp", "Ccn2"), Reactive_glial_support = c("Gfap", "Lcn2", "C3", "Serpina3n", 
    "Vim", "Timp1"))

MIN_GENE_COVERAGE <- 0.50
HIGH_SCORE_QUANTILES <- c(0.75, 0.80)
K_NEIGHBORS_MORAN <- 6L

required_pkgs <- c("Matrix", "ggplot2")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them before running Step 114.")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_csv_utf8 <- write_csv_utf8

read_csv_gz <- function(filename) {
    con <- gzfile(filename, open = "rt")
    on.exit(close(con), add = TRUE)
    read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv_gz <- function(x, filename) {
    con <- gzfile(filename, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(x, file = con, row.names = FALSE, na = "")
}

read_text_or_gz <- function(filename) {
  if (grepl("\\.gz$", filename, ignore.case = TRUE)) {
    con <- gzfile(filename, open = "rt")
  } else {
    con <- file(filename, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

read_tsv_or_gz <- function(filename) {
    con <- if (grepl("\\.gz$", filename, ignore.case = TRUE)) {
        gzfile(filename, open = "rt")
    }
    else {
        file(filename, open = "rt")
    }
    on.exit(close(con), add = TRUE)
    read.delim(con, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, 
        check.names = FALSE)
}

pick_existing <- function(folder, candidates) {
  p <- file.path(folder, candidates)
  p <- p[file.exists(p)]
  if (length(p) == 0L) return(NA_character_)
  p[1]
}

find_matrix_dir <- function(gsm) {
    gsm_dir <- file.path(EXTRACT_DIR, gsm)
    if (!dir.exists(gsm_dir)) {
        stop("Extracted sample directory is missing for ", gsm, ":\n", gsm_dir, "\nRun the existing Step 60 extraction first. Step 114 does not ", 
            "extract ZIP/TAR archives.")
    }
    dirs <- list.dirs(gsm_dir, recursive = TRUE, full.names = TRUE)
    dirs <- dirs[basename(dirs) == "filtered_feature_bc_matrix"]
    if (length(dirs) != 1L) {
        stop(gsm, ": expected exactly one filtered_feature_bc_matrix directory ", "under ", gsm_dir, 
            "; found ", length(dirs), ".")
    }
    dirs[1]
}

read_10x_matrix_direct <- function(matrix_dir) {
  matrix_file <- pick_existing(
    matrix_dir,
    c("matrix.mtx.gz", "matrix.mtx")
  )
  feature_file <- pick_existing(
    matrix_dir,
    c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
  )
  barcode_file <- pick_existing(
    matrix_dir,
    c("barcodes.tsv.gz", "barcodes.tsv")
  )

  if (anyNA(c(matrix_file, feature_file, barcode_file))) {
    stop(
      "Missing one or more 10x matrix files in:\n", matrix_dir,
      "\nRequired: matrix.mtx[.gz], features.tsv[.gz] or genes.tsv[.gz], ",
      "barcodes.tsv[.gz]."
    )
  }

  matrix_con <- if (grepl("\\.gz$", matrix_file, ignore.case = TRUE)) {
    gzfile(matrix_file, open = "rt")
  } else {
    file(matrix_file, open = "rt")
  }
  counts <- tryCatch(
    Matrix::readMM(matrix_con),
    finally = close(matrix_con)
  )

  features <- read_tsv_or_gz(feature_file)
  barcodes <- read_text_or_gz(barcode_file)

  if (nrow(counts) != nrow(features)) {
    stop(
      "Matrix/feature dimension mismatch in ", matrix_dir,
      ": matrix rows=", nrow(counts),
      "; feature rows=", nrow(features), "."
    )
  }

  if (ncol(counts) != length(barcodes)) {
    stop(
      "Matrix/barcode dimension mismatch in ", matrix_dir,
      ": matrix columns=", ncol(counts),
      "; barcodes=", length(barcodes), "."
    )
  }

  if (ncol(features) >= 2L) {
    gene_symbol <- as.character(features[[2]])
    gene_id <- as.character(features[[1]])
  } else {
    gene_symbol <- as.character(features[[1]])
    gene_id <- gene_symbol
  }

  if (ncol(features) >= 3L) {
    feature_type <- as.character(features[[3]])
    gene_expression_rows <- which(
      is.na(feature_type) |
        feature_type == "" |
        feature_type == "Gene Expression"
    )
    if (length(gene_expression_rows) > 0L) {
      counts <- counts[gene_expression_rows, , drop = FALSE]
      gene_symbol <- gene_symbol[gene_expression_rows]
      gene_id <- gene_id[gene_expression_rows]
    }
  }

  gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <-
    gene_id[is.na(gene_symbol) | gene_symbol == ""]

  list(
    counts = methods::as(counts, "dgCMatrix"),
    gene_symbol = gene_symbol,
    barcodes = as.character(barcodes),
    matrix_file = matrix_file,
    feature_file = feature_file,
    barcode_file = barcode_file
  )
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x)
}

safe_quantile <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) 
        return(NA_real_)
    as.numeric(stats::quantile(x, probs = probability, na.rm = TRUE, names = FALSE, type = 7))
}

safe_row_mean_two <- function(a, b) {
  out <- rep(NA_real_, length(a))
  ok_a <- is.finite(a)
  ok_b <- is.finite(b)
  both <- ok_a & ok_b
  out[both] <- (a[both] + b[both]) / 2
  only_a <- ok_a & !ok_b
  only_b <- !ok_a & ok_b
  out[only_a] <- a[only_a]
  out[only_b] <- b[only_b]
  out
}

moran_knn <- function(score, x, y, k = 6L) {
    ok <- is.finite(score) & is.finite(x) & is.finite(y)
    score <- score[ok]
    coords <- cbind(x[ok], y[ok])
    n <- length(score)
    if (n < max(10L, k + 2L)) 
        return(NA_real_)
    centered <- score - mean(score)
    denominator <- sum(centered^2)
    if (!is.finite(denominator) || denominator <= 0) 
        return(NA_real_)
    dmat <- as.matrix(stats::dist(coords))
    diag(dmat) <- Inf
    k_use <- min(k, n - 1L)
    nn <- t(vapply(seq_len(n), function(i) order(dmat[i, ])[seq_len(k_use)], integer(k_use)))
    ii <- rep(seq_len(n), each = k_use)
    jj <- as.vector(t(nn))
    numerator <- sum(centered[ii] * centered[jj])
    weight_sum <- length(ii)
    (n/weight_sum) * numerator/denominator
}

required_inputs <- c(SAMPLE_MANIFEST_FILE, SECTION_MANIFEST_FILE, SPOT_SECTION_MAP_FILE)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0L) {
    stop("Required Step 114 input(s) missing:\n", paste(missing_inputs, collapse = "\n"), "\nPlease complete the existing Step 60 and Step 68 outputs first.")
}

if (!dir.exists(EXTRACT_DIR)) {
    stop("Already-extracted 10x matrices were not found:\n", EXTRACT_DIR, "\nStep 114 does not use ZIP/TAR archives and will not extract them.")
}

sample_manifest <- read.csv(SAMPLE_MANIFEST_FILE, stringsAsFactors = FALSE, check.names = FALSE)

section_manifest <- read.csv(SECTION_MANIFEST_FILE, stringsAsFactors = FALSE, check.names = FALSE)

spot_section_map <- read_csv_gz(SPOT_SECTION_MAP_FILE)

required_manifest_cols <- c("GSM", "sample_label")
required_section_cols <- c("sample_label", "GSM", "timepoint", "section_position_internal", "section_uid_internal")
required_map_cols <- c("barcode", "GSM", "sample_label", "timepoint", "array_x", "array_y", "is_major_component", 
    "section_position_internal", "section_uid_internal")

if (!all(required_manifest_cols %in% colnames(sample_manifest))) {
    stop("Step 60 sample manifest is missing required columns: ", paste(setdiff(required_manifest_cols, 
        colnames(sample_manifest)), collapse = ", "))
}
if (!all(required_section_cols %in% colnames(section_manifest))) {
    stop("Step 68 section manifest is missing required columns: ", paste(setdiff(required_section_cols, 
        colnames(section_manifest)), collapse = ", "))
}
if (!all(required_map_cols %in% colnames(spot_section_map))) {
    stop("Step 68 spot-to-section map is missing required columns: ", paste(setdiff(required_map_cols, 
        colnames(spot_section_map)), collapse = ", "))
}

sample_manifest$GSM <- as.character(sample_manifest$GSM)
sample_manifest$sample_label <- as.character(sample_manifest$sample_label)

if (nrow(sample_manifest) != 8L || !setequal(sample_manifest$sample_label, EXPECTED_SAMPLES) || anyDuplicated(sample_manifest$GSM) > 
    0L || anyDuplicated(sample_manifest$sample_label) > 0L) {
    stop("Step 60 sample manifest does not contain the expected eight unique samples.")
}

if (nrow(section_manifest) != EXPECTED_TOTAL_SECTIONS || anyDuplicated(section_manifest$section_uid_internal) > 
    0L) {
    stop("Step 68 section manifest must contain 32 unique coordinate-defined sections. ", "Observed rows: ", 
        nrow(section_manifest), ".")
}

sections_per_sample <- table(section_manifest$sample_label)
if (!all(sections_per_sample[EXPECTED_SAMPLES] == EXPECTED_SECTIONS_PER_SAMPLE)) {
  stop("Step 68 manifest does not contain exactly four sections per capture area.")
}

if (!all(CANDIDATE_SECTIONS %in% section_manifest$section_uid_internal)) {
    stop("One or more user-proposed candidate section IDs are absent from Step 68:\n", paste(setdiff(CANDIDATE_SECTIONS, 
        section_manifest$section_uid_internal), collapse = "\n"))
}

manifest_map <- sample_manifest[, c("GSM", "sample_label")]
section_map <- unique(section_manifest[, c("GSM", "sample_label")])
check_map <- merge(manifest_map, section_map, by = c("GSM", "sample_label"))
if (nrow(check_map) != 8L) {
  stop("Step 60 and Step 68 GSM/sample_label mappings do not agree.")
}

all_signature_genes <- unique(unlist(GENE_SETS, use.names = FALSE))
all_signature_genes_lower <- tolower(all_signature_genes)

spot_score_list <- list()
coverage_rows <- list()
matrix_audit_rows <- list()

for (sample_i in seq_len(nrow(sample_manifest))) {

  gsm <- sample_manifest$GSM[sample_i]
  sample_label <- sample_manifest$sample_label[sample_i]

  matrix_dir <- find_matrix_dir(gsm)
  tenx <- read_10x_matrix_direct(matrix_dir)

  counts <- tenx$counts
  barcodes <- tenx$barcodes
  gene_symbol <- tenx$gene_symbol

  if (anyDuplicated(barcodes) > 0L) {
    stop("Duplicated barcodes in matrix for ", gsm, ".")
  }

  map_i <- spot_section_map[
    spot_section_map$GSM == gsm,
    ,
    drop = FALSE
  ]

  if (nrow(map_i) == 0L) {
    stop("No Step 68 spot map rows found for ", gsm, ".")
  }

  if (anyDuplicated(map_i$barcode) > 0L) {
    stop("Duplicated barcode rows in Step 68 map for ", gsm, ".")
  }

  map_idx <- match(barcodes, map_i$barcode)
  matched <- !is.na(map_idx)

  match_rate <- mean(matched)

  if (match_rate < 0.95) {
    stop(
      gsm, ": only ", round(100 * match_rate, 2),
      "% of matrix barcodes match the Step 68 spot map; expected >=95%."
    )
  }

  mapped_rows <- map_i[map_idx[matched], , drop = FALSE]
  keep_major <- !is.na(mapped_rows$section_uid_internal) &
    mapped_rows$is_major_component %in% TRUE

  if (sum(keep_major) == 0L) {
    stop("No major-section spots matched for ", gsm, ".")
  }

  keep_matrix_cols <- which(matched)[keep_major]
  mapped_rows <- mapped_rows[keep_major, , drop = FALSE]

  if (anyDuplicated(mapped_rows$barcode) > 0L) {
    stop("Duplicate mapped barcodes after section filtering for ", gsm, ".")
  }

  total_umi_all <- as.numeric(Matrix::colSums(counts))
  detected_features_all <- as.numeric(Matrix::colSums(counts > 0))

  total_umi <- total_umi_all[keep_matrix_cols]
  detected_features <- detected_features_all[keep_matrix_cols]

  target_counts <- matrix(
    0,
    nrow = length(all_signature_genes),
    ncol = length(keep_matrix_cols),
    dimnames = list(all_signature_genes, barcodes[keep_matrix_cols])
  )

  gene_symbol_lower <- tolower(gene_symbol)

  for (g in all_signature_genes) {
    feature_rows <- which(gene_symbol_lower == tolower(g))
    if (length(feature_rows) > 0L) {
      target_counts[g, ] <- as.numeric(
        Matrix::colSums(
          counts[feature_rows, keep_matrix_cols, drop = FALSE]
        )
      )
    }
  }

  feature_genes_lower <- unique(gene_symbol_lower)

  denom <- pmax(total_umi, 1) / 10000
  log_cp10k <- log1p(sweep(target_counts, 2, denom, "/"))

  gene_z <- log_cp10k
  gene_means <- rowMeans(log_cp10k)
  gene_sds <- apply(log_cp10k, 1, stats::sd)

  for (r in seq_len(nrow(gene_z))) {
    if (is.finite(gene_sds[r]) && gene_sds[r] > 0) {
      gene_z[r, ] <- (log_cp10k[r, ] - gene_means[r]) / gene_sds[r]
    } else {
      gene_z[r, ] <- 0
    }
  }

  score_z <- list()
  score_log <- list()
  coverage_for_sample <- list()

  for (set_name in names(GENE_SETS)) {

    genes <- GENE_SETS[[set_name]]
    present <- genes[
      tolower(genes) %in% feature_genes_lower
    ]

    min_required <- max(
      2L,
      as.integer(ceiling(length(genes) * MIN_GENE_COVERAGE))
    )

    enough <- length(present) >= min_required

    if (enough) {

      present <- intersect(present, all_signature_genes)

      score_z[[set_name]] <- colMeans(
        gene_z[present, , drop = FALSE],
        na.rm = TRUE
      )
      score_log[[set_name]] <- colMeans(
        log_cp10k[present, , drop = FALSE],
        na.rm = TRUE
      )
    } else {
      score_z[[set_name]] <- rep(NA_real_, length(keep_matrix_cols))
      score_log[[set_name]] <- rep(NA_real_, length(keep_matrix_cols))
    }

    coverage_for_sample[[set_name]] <- data.frame(
      GSM = gsm,
      sample_label = sample_label,
      gene_set = set_name,
      n_genes_requested = length(genes),
      n_genes_in_feature_table = length(present),
      minimum_required = min_required,
      coverage_pass = enough,
      genes_in_feature_table = paste(present, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  spot_scores <- data.frame(
    barcode = barcodes[keep_matrix_cols],
    GSM = gsm,
    sample_label = sample_label,
    timepoint = as.character(mapped_rows$timepoint),
    array_x = as.numeric(mapped_rows$array_x),
    array_y = as.numeric(mapped_rows$array_y),
    component_id_internal = mapped_rows$component_id_internal,
    section_position_internal = as.character(mapped_rows$section_position_internal),
    section_uid_internal = as.character(mapped_rows$section_uid_internal),
    nCount_raw = total_umi,
    nFeature_raw = detected_features,
    stringsAsFactors = FALSE
  )

  for (set_name in names(GENE_SETS)) {
    z_col <- paste0(set_name, "_score_z")
    log_col <- paste0(set_name, "_score_logCP10K")
    spot_scores[[z_col]] <- as.numeric(score_z[[set_name]])
    spot_scores[[log_col]] <- as.numeric(score_log[[set_name]])
  }

  spot_scores$normal_combined_score_z <- safe_row_mean_two(
    spot_scores$Normal_neuronal_score_z,
    spot_scores$Normal_myelin_oligodendrocyte_score_z
  )

  spot_scores$injury_combined_score_z <- safe_row_mean_two(
    spot_scores$Injury_myeloid_score_z,
    spot_scores$Injury_fibrotic_ECM_score_z
  )

  spot_scores$normal_combined_score_logCP10K <- safe_row_mean_two(
    spot_scores$Normal_neuronal_score_logCP10K,
    spot_scores$Normal_myelin_oligodendrocyte_score_logCP10K
  )

  spot_scores$injury_combined_score_logCP10K <- safe_row_mean_two(
    spot_scores$Injury_myeloid_score_logCP10K,
    spot_scores$Injury_fibrotic_ECM_score_logCP10K
  )

  spot_scores$normal_minus_injury_score_z <-
    spot_scores$normal_combined_score_z -
    spot_scores$injury_combined_score_z

  injury_q75 <- safe_quantile(spot_scores$injury_combined_score_z, 0.75)
  injury_q80 <- safe_quantile(spot_scores$injury_combined_score_z, 0.80)
  normal_q75 <- safe_quantile(spot_scores$normal_combined_score_z, 0.75)
  normal_q80 <- safe_quantile(spot_scores$normal_combined_score_z, 0.80)

  spot_scores$injury_high_q75 <- if (is.finite(injury_q75)) {
    spot_scores$injury_combined_score_z >= injury_q75
  } else NA
  spot_scores$injury_high_q80 <- if (is.finite(injury_q80)) {
    spot_scores$injury_combined_score_z >= injury_q80
  } else NA
  spot_scores$normal_high_q75 <- if (is.finite(normal_q75)) {
    spot_scores$normal_combined_score_z >= normal_q75
  } else NA
  spot_scores$normal_high_q80 <- if (is.finite(normal_q80)) {
    spot_scores$normal_combined_score_z >= normal_q80
  } else NA

  spot_score_list[[sample_label]] <- spot_scores
  coverage_rows[[sample_label]] <- do.call(rbind, coverage_for_sample)

  matrix_audit_rows[[sample_label]] <- data.frame(
    GSM = gsm,
    sample_label = sample_label,
    matrix_dir = matrix_dir,
    n_matrix_barcodes = length(barcodes),
    n_matrix_features = nrow(counts),
    n_matrix_barcodes_matched_to_step68 = sum(matched),
    barcode_match_rate = match_rate,
    n_major_section_spots_scored = nrow(spot_scores),
    n_major_sections_scored = length(unique(spot_scores$section_uid_internal)),
    matrix_mtx_file = tenx$matrix_file,
    features_file = tenx$feature_file,
    barcodes_file = tenx$barcode_file,
    stringsAsFactors = FALSE
  )

  rm(counts, target_counts, log_cp10k, gene_z, tenx)
  invisible(gc())
}

spot_scores_all <- do.call(rbind, spot_score_list)
rownames(spot_scores_all) <- NULL

gene_coverage <- do.call(rbind, coverage_rows)
rownames(gene_coverage) <- NULL

primary_signature_sets <- c("Normal_neuronal", "Normal_myelin_oligodendrocyte", "Injury_myeloid", "Injury_fibrotic_ECM")

bad_primary_coverage <- gene_coverage[gene_coverage$gene_set %in% primary_signature_sets & !gene_coverage$coverage_pass, 
    , drop = FALSE]

if (nrow(bad_primary_coverage) > 0L) {
    write_csv_utf8(bad_primary_coverage, file.path(OUT_DIR, paste0("01b_114_primary_gene_set_coverage_failures_", 
        VERSION, ".csv")))
    stop("At least one primary signature failed the prespecified gene-coverage requirement. ", "See 01b_114_primary_gene_set_coverage_failures_", 
        VERSION, ".csv. ", "Do not interpret candidate-section scores until gene coverage is resolved.")
}

matrix_audit <- do.call(rbind, matrix_audit_rows)
rownames(matrix_audit) <- NULL

observed_section_ids <- unique(spot_scores_all$section_uid_internal)

missing_sections <- setdiff(section_manifest$section_uid_internal, observed_section_ids)

if (length(missing_sections) > 0L) {
    stop("The following Step 68 sections have no matched expression spots:\n", paste(missing_sections, 
        collapse = "\n"))
}

section_spot_counts <- table(spot_scores_all$section_uid_internal)
if (!all(section_spot_counts[section_manifest$section_uid_internal] > 0L)) {
  stop("At least one coordinate-defined section has zero matched spots.")
}

section_summary_rows <- vector("list", nrow(section_manifest))

for (i in seq_len(nrow(section_manifest))) {
    meta <- section_manifest[i, , drop = FALSE]
    uid <- as.character(meta$section_uid_internal)
    samp <- as.character(meta$sample_label)
    d <- spot_scores_all[spot_scores_all$section_uid_internal == uid, , drop = FALSE]
    section_summary_rows[[i]] <- data.frame(sample_label = samp, GSM = as.character(meta$GSM), timepoint = as.character(meta$timepoint), 
        section_position_internal = as.character(meta$section_position_internal), component_id_internal = meta$component_id_internal, 
        section_uid_internal = uid, putative_uninjured_candidate = uid %in% CANDIDATE_SECTIONS, n_spots = nrow(d), 
        mean_UMI_per_spot = safe_mean(d$nCount_raw), median_UMI_per_spot = safe_median(d$nCount_raw), 
        mean_detected_features_per_spot = safe_mean(d$nFeature_raw), median_detected_features_per_spot = safe_median(d$nFeature_raw), 
        normal_neuronal_median_z = safe_median(d$Normal_neuronal_score_z), normal_myelin_median_z = safe_median(d$Normal_myelin_oligodendrocyte_score_z), 
        injury_myeloid_median_z = safe_median(d$Injury_myeloid_score_z), injury_ECM_median_z = safe_median(d$Injury_fibrotic_ECM_score_z), 
        reactive_glial_median_z = safe_median(d$Reactive_glial_support_score_z), normal_combined_median_z = safe_median(d$normal_combined_score_z), 
        injury_combined_median_z = safe_median(d$injury_combined_score_z), normal_minus_injury_median_z = safe_median(d$normal_minus_injury_score_z), 
        normal_combined_mean_logCP10K = safe_mean(d$normal_combined_score_logCP10K), injury_combined_mean_logCP10K = safe_mean(d$injury_combined_score_logCP10K), 
        normal_high_spot_fraction_q75 = safe_mean(as.numeric(d$normal_high_q75)), normal_high_spot_fraction_q80 = safe_mean(as.numeric(d$normal_high_q80)), 
        injury_high_spot_fraction_q75 = safe_mean(as.numeric(d$injury_high_q75)), injury_high_spot_fraction_q80 = safe_mean(as.numeric(d$injury_high_q80)), 
        moran_I_normal_score_knn6 = moran_knn(d$normal_combined_score_z, d$array_x, d$array_y, k = K_NEIGHBORS_MORAN), 
        moran_I_injury_score_knn6 = moran_knn(d$injury_combined_score_z, d$array_x, d$array_y, k = K_NEIGHBORS_MORAN), 
        stringsAsFactors = FALSE)
}

section_summary <- do.call(rbind, section_summary_rows)
rownames(section_summary) <- NULL

section_summary$rank_normal_minus_injury_within_capture <- NA_integer_
section_summary$rank_injury_low_within_capture <- NA_integer_
section_summary$delta_normal_vs_other_sections <- NA_real_
section_summary$delta_injury_vs_other_sections <- NA_real_
section_summary$delta_injury_high_fraction_q75_vs_others <- NA_real_
section_summary$delta_normal_high_fraction_q75_vs_others <- NA_real_
section_summary$candidate_support_components <- NA_integer_
section_summary$candidate_support_pattern <- "not a proposed candidate"

for (samp in EXPECTED_SAMPLES) {

  idx <- which(section_summary$sample_label == samp)
  ss <- section_summary[idx, , drop = FALSE]

  section_summary$rank_normal_minus_injury_within_capture[idx] <-
    rank(ss$normal_minus_injury_median_z, ties.method = "min",
         na.last = "keep")

  section_summary$rank_injury_low_within_capture[idx] <-
    rank(ss$injury_combined_median_z, ties.method = "max",
         na.last = "keep")

  for (j in idx) {
    if (!section_summary$putative_uninjured_candidate[j]) next

    other <- section_summary[
      section_summary$sample_label == samp &
        section_summary$section_uid_internal !=
          section_summary$section_uid_internal[j],
      ,
      drop = FALSE
    ]

    delta_normal <- section_summary$normal_combined_median_z[j] -
      safe_median(other$normal_combined_median_z)

    delta_injury <- section_summary$injury_combined_median_z[j] -
      safe_median(other$injury_combined_median_z)

    delta_injury_fraction <- section_summary$injury_high_spot_fraction_q75[j] -
      safe_median(other$injury_high_spot_fraction_q75)

    delta_normal_fraction <- section_summary$normal_high_spot_fraction_q75[j] -
      safe_median(other$normal_high_spot_fraction_q75)

    support_n <- sum(c(
      is.finite(delta_normal) && delta_normal > 0,
      is.finite(delta_injury) && delta_injury < 0,
      is.finite(delta_injury_fraction) && delta_injury_fraction < 0
    ))

    pattern <- if (support_n == 3L) {
      "concordant_normal_like_direction"
    } else if (support_n == 2L) {
      "partially_supportive_normal_like_direction"
    } else if (
      is.finite(delta_normal) && is.finite(delta_injury) &&
        delta_normal < 0 && delta_injury > 0
    ) {
      "contradictory_injury_like_direction"
    } else if (support_n == 1L) {
      "weak_or_mixed_direction"
    } else {
      "indeterminate"
    }

    section_summary$delta_normal_vs_other_sections[j] <- delta_normal
    section_summary$delta_injury_vs_other_sections[j] <- delta_injury
    section_summary$delta_injury_high_fraction_q75_vs_others[j] <-
      delta_injury_fraction
    section_summary$delta_normal_high_fraction_q75_vs_others[j] <-
      delta_normal_fraction
    section_summary$candidate_support_components[j] <- support_n
    section_summary$candidate_support_pattern[j] <- pattern
  }
}

candidate_validation <- section_summary[section_summary$putative_uninjured_candidate, , drop = FALSE]

candidate_validation$rank_is_top_normal_minus_injury <-
  candidate_validation$rank_normal_minus_injury_within_capture == 4L

candidate_validation$rank_is_lowest_injury <-
  candidate_validation$rank_injury_low_within_capture == 1L

candidate_validation$rank_concordant <- candidate_validation$rank_is_top_normal_minus_injury & candidate_validation$rank_is_lowest_injury

best_section_idx <- unlist(lapply(EXPECTED_SAMPLES, function(samp) {
    idx <- which(section_summary$sample_label == samp)
    if (length(idx) == 0L) 
        return(integer(0))
    idx[which.max(section_summary$normal_minus_injury_median_z[idx])]
}))

section_summary$highest_normal_minus_injury_section <- section_summary$section_uid_internal %in% section_summary$section_uid_internal[best_section_idx]

rctd_section_summary <- data.frame()
rctd_region_proportions <- data.frame()
rctd_available <- file.exists(RCTD_FILE)

if (rctd_available) {
    rctd <- read_csv_gz(RCTD_FILE)
    required_rctd_cols <- c("section_uid_internal", "GSM", "sample_label", "array_x", "array_y", "spatial_region_v03")
    if (all(required_rctd_cols %in% colnames(rctd))) {
        rctd <- rctd[!is.na(rctd$section_uid_internal) & rctd$section_uid_internal != "", , drop = FALSE]
        normal_cols <- intersect(c("Neuron", "Oligodendrocyte"), colnames(rctd))
        injury_cols <- intersect(c("Activated inflammatory Microglia", "Acute-response Macrophage", "Clec7a+ lipid-handling Microglia", 
            "Gpnmb+ tissue-repair Macrophage", "MHC-II+ antigen-presenting Macrophage", "Proliferating Microglia", 
            "Tnf+ inflammatory Microglia", "Ccl12+ acute Microglia", "Neutrophil"), colnames(rctd))
        fibrotic_cols <- intersect(c("Fibroblast/Meningeal"), colnames(rctd))
        if (length(normal_cols) > 0L) {
            rctd$rctd_normal_weight <- rowSums(rctd[, normal_cols, drop = FALSE], na.rm = TRUE)
        }
        else {
            rctd$rctd_normal_weight <- NA_real_
        }
        if (length(injury_cols) > 0L) {
            rctd$rctd_injury_myeloid_weight <- rowSums(rctd[, injury_cols, drop = FALSE], na.rm = TRUE)
        }
        else {
            rctd$rctd_injury_myeloid_weight <- NA_real_
        }
        if (length(fibrotic_cols) > 0L) {
            rctd$rctd_fibrotic_weight <- rowSums(rctd[, fibrotic_cols, drop = FALSE], na.rm = TRUE)
        }
        else {
            rctd$rctd_fibrotic_weight <- NA_real_
        }
        rctd_sections <- unique(rctd$section_uid_internal)
        rctd_section_summary <- do.call(rbind, lapply(rctd_sections, function(uid) {
            d <- rctd[rctd$section_uid_internal == uid, , drop = FALSE]
            data.frame(section_uid_internal = uid, sample_label = as.character(d$sample_label[1]), GSM = as.character(d$GSM[1]), 
                n_RCTD_spots = nrow(d), mean_RCTD_normal_weight = safe_mean(d$rctd_normal_weight), mean_RCTD_injury_myeloid_weight = safe_mean(d$rctd_injury_myeloid_weight), 
                mean_RCTD_fibrotic_weight = safe_mean(d$rctd_fibrotic_weight), stringsAsFactors = FALSE)
        }))
        if (all(c("section_uid_internal", "spatial_region_v03") %in% colnames(rctd))) {
            region_tab <- as.data.frame(table(section_uid_internal = rctd$section_uid_internal, spatial_region_v03 = rctd$spatial_region_v03), 
                stringsAsFactors = FALSE)
            region_tab <- region_tab[region_tab$Freq > 0, , drop = FALSE]
            totals <- aggregate(Freq ~ section_uid_internal, data = region_tab, FUN = sum)
            region_tab <- merge(region_tab, totals, by = "section_uid_internal", suffixes = c("", "_total"))
            region_tab$fraction <- region_tab$Freq/region_tab$Freq_total
            rctd_region_proportions <- region_tab
        }
    }
    else {
        warning("Optional RCTD file exists but required section-aware columns are missing; ", "RCTD cross-check skipped.")
        rctd_available <- FALSE
    }
} else {
}

if (nrow(rctd_section_summary) > 0L) {
    section_summary <- merge(section_summary, rctd_section_summary, by = c("section_uid_internal", "sample_label", 
        "GSM"), all.x = TRUE, sort = FALSE)
    candidate_validation <- section_summary[section_summary$putative_uninjured_candidate, , drop = FALSE]
}

write_csv_utf8(matrix_audit, file.path(OUT_DIR, paste0("00_114_matrix_input_audit_", VERSION, ".csv")))

write_csv_utf8(gene_coverage, file.path(OUT_DIR, paste0("01_114_gene_set_coverage_", VERSION, ".csv")))

write_csv_gz(spot_scores_all, file.path(OUT_DIR, paste0("02_114_spot_signature_scores_", VERSION, ".csv.gz")))

write_csv_utf8(section_summary, file.path(OUT_DIR, paste0("03_114_section_signature_summary_", VERSION, 
    ".csv")))

write_csv_utf8(candidate_validation, file.path(OUT_DIR, paste0("04_114_candidate_section_validation_", 
    VERSION, ".csv")))

if (nrow(rctd_section_summary) > 0L) {
    write_csv_utf8(rctd_section_summary, file.path(OUT_DIR, paste0("05_114_RCTD_section_crosscheck_", 
        VERSION, ".csv")))
}

if (nrow(rctd_region_proportions) > 0L) {
    write_csv_utf8(rctd_region_proportions, file.path(OUT_DIR, paste0("06_114_RCTD_region_composition_", 
        VERSION, ".csv")))
}

p_rank <- ggplot(section_summary, aes(x = section_position_internal, y = normal_minus_injury_median_z, 
    color = putative_uninjured_candidate, shape = putative_uninjured_candidate)) + geom_hline(yintercept = 0, 
    linewidth = 0.35, color = "grey65") + geom_point(size = 2.8) + facet_wrap(~sample_label, ncol = 4) + 
    scale_color_manual(values = c(`FALSE` = "grey35", `TRUE` = "#D55E00"), labels = c(`FALSE` = "Other section", 
        `TRUE` = "Putative uninjured")) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "Other section", 
    `TRUE` = "Putative uninjured")) + labs(title = "Normal-versus-injury signature balance by coordinate-defined section", 
    subtitle = "Within-capture-area gene-wise z-score; candidate sections are hypotheses", x = "Internal coordinate position", 
    y = "Median normal-minus-injury score", color = NULL, shape = NULL) + theme_classic(base_size = 9) + 
    theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(size = 8), legend.position = "bottom", 
        plot.title = element_text(size = 11))

ggsave(file.path(OUT_DIR, paste0("07_114_section_signature_rank_", VERSION, ".pdf")), p_rank, width = 12, 
    height = 7, units = "in", device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, paste0("07_114_section_signature_rank_", VERSION, ".png")), p_rank, width = 12, 
    height = 7, units = "in", dpi = 300, bg = "white")

p_spatial <- ggplot(spot_scores_all, aes(x = array_x, y = array_y, color = injury_combined_score_z)) + 
    geom_point(size = 0.65, alpha = 0.95) + coord_equal() + scale_y_reverse() + scale_color_gradient2(low = "#2166AC", 
    mid = "grey95", high = "#B2182B", midpoint = 0, na.value = "grey80") + facet_wrap(~section_uid_internal, 
    ncol = 4) + labs(title = "Spatial distribution of injury-associated transcriptional scores", subtitle = "Spot-level targeted score; Visium array coordinates only, no H&E registration", 
    x = "Visium array X", y = "Visium array Y", color = "Injury score") + theme_classic(base_size = 7) + 
    theme(strip.text = element_text(size = 6.5), legend.position = "right", plot.title = element_text(size = 11))

ggsave(file.path(OUT_DIR, paste0("08_114_spatial_injury_score_maps_", VERSION, ".pdf")), p_spatial, width = 14, 
    height = 18, units = "in", device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, paste0("08_114_spatial_injury_score_maps_", VERSION, ".png")), p_spatial, width = 14, 
    height = 18, units = "in", dpi = 300, bg = "white")

heat_metrics <- c("normal_neuronal_median_z", "normal_myelin_median_z", "injury_myeloid_median_z", "injury_ECM_median_z", 
    "reactive_glial_median_z")

heat_labels <- c(normal_neuronal_median_z = "Normal neuronal", normal_myelin_median_z = "Normal myelin/oligodendrocyte", 
    injury_myeloid_median_z = "Injury myeloid", injury_ECM_median_z = "Injury fibrotic/ECM", reactive_glial_median_z = "Reactive glial")

heat_list <- lapply(heat_metrics, function(m) {
    data.frame(sample_label = section_summary$sample_label, section_uid_internal = section_summary$section_uid_internal, 
        section_position_internal = section_summary$section_position_internal, putative_uninjured_candidate = section_summary$putative_uninjured_candidate, 
        signature = unname(heat_labels[m]), score = section_summary[[m]], stringsAsFactors = FALSE)
})
heat_df <- do.call(rbind, heat_list)

heat_df$score_display_z <- ave(heat_df$score, heat_df$signature, FUN = function(x) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) 
        return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE))/sx
})

p_heat <- ggplot(heat_df, aes(x = section_position_internal, y = signature, fill = score_display_z)) + 
    geom_tile(color = "white", linewidth = 0.25) + facet_wrap(~sample_label, ncol = 4) + scale_fill_gradient2(low = "#2166AC", 
    mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") + labs(title = "Section-level transcriptional signature overview", 
    subtitle = "Display z-scores standardized across all 32 coordinate-defined sections", x = "Internal coordinate position", 
    y = NULL, fill = "Display z") + theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 35, 
    hjust = 1), strip.text = element_text(size = 8), plot.title = element_text(size = 11))

ggsave(file.path(OUT_DIR, paste0("09_114_section_signature_heatmap_", VERSION, ".pdf")), p_heat, width = 12, 
    height = 8, units = "in", device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, paste0("09_114_section_signature_heatmap_", VERSION, ".png")), p_heat, width = 12, 
    height = 8, units = "in", dpi = 300, bg = "white")

qc <- data.frame(check = c("Eight expected samples are present in Step 60 manifest", "Thirty-two unique coordinate-defined sections are present", 
    "Four coordinate-defined sections are present per capture area", "All five user-proposed candidate section IDs are present", 
    "All eight extracted count matrices were processed", "At least 95 percent of each matrix barcode set matched Step 68", 
    "All 32 sections have matched raw-count spots", "All four primary gene signatures pass feature-table coverage in every sample", 
    "No RDS/RData input was specified or loaded", "No ZIP/TAR archive was used", "No clustering, differential expression, or RCTD was rerun", 
    "Optional RCTD cross-check was available"), pass = c(nrow(sample_manifest) == 8L && setequal(sample_manifest$sample_label, 
    EXPECTED_SAMPLES), nrow(section_manifest) == 32L && length(unique(section_manifest$section_uid_internal)) == 
    32L, all(table(section_manifest$sample_label)[EXPECTED_SAMPLES] == 4L), all(CANDIDATE_SECTIONS %in% 
    section_manifest$section_uid_internal), nrow(matrix_audit) == 8L, all(matrix_audit$barcode_match_rate >= 
    0.95), all(section_summary$n_spots > 0L), all(gene_coverage$coverage_pass[gene_coverage$gene_set %in% 
    primary_signature_sets]), TRUE, TRUE, TRUE, rctd_available), stringsAsFactors = FALSE)

write_csv_utf8(qc, file.path(OUT_DIR, paste0("10_114_FINAL_QC_", VERSION, ".csv")))

candidate_lines <- character(0)

for (i in seq_len(nrow(candidate_validation))) {
    d <- candidate_validation[i, , drop = FALSE]
    candidate_lines <- c(candidate_lines, paste0("  ", d$section_uid_internal, " | normal-minus-injury median=", 
        round(d$normal_minus_injury_median_z, 3), " | injury median=", round(d$injury_combined_median_z, 
            3), " | delta normal vs other sections=", round(d$delta_normal_vs_other_sections, 3), " | delta injury vs other sections=", 
        round(d$delta_injury_vs_other_sections, 3), " | injury-high fraction delta=", round(d$delta_injury_high_fraction_q75_vs_others, 
            3), " | support=", d$candidate_support_pattern, " | rank normal-minus-injury=", d$rank_normal_minus_injury_within_capture, 
        "/4"))
}

coverage_lines <- apply(gene_coverage, 1, function(x) {
    paste0("  ", x[["sample_label"]], " | ", x[["gene_set"]], ": ", x[["n_genes_in_feature_table"]], 
        "/", x[["n_genes_requested"]], " genes; pass=", x[["coverage_pass"]])
})

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = file.path(OUT_DIR, paste0("11_114_REPORT_", VERSION, ".txt")), useBytes = TRUE)

print(candidate_validation[, c("section_uid_internal", "n_spots", "normal_combined_median_z", "injury_combined_median_z", 
    "normal_minus_injury_median_z", "injury_high_spot_fraction_q75", "candidate_support_pattern", "rank_normal_minus_injury_within_capture"), 
    drop = FALSE], row.names = FALSE)

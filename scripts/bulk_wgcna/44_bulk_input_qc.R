source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

DATA_DIR <- file.path(ROOT, "data", "17702045")

EXPR_FILE <- file.path(DATA_DIR, "sci_bulkrna_expression.csv")

META_FILE <- file.path(DATA_DIR, "sci_bulkrna_samples.csv")

SCRNA_FILE <- file.path(DATA_DIR, "sci_scrna.rds")

OUT_DIR <- file.path(ROOT, "results", "44_Bulk_17702045_input_QC")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "44_Bulk_17702045_input_QC_v04"

EXPECTED_TIME_LEVELS <- c("Uninjured", "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi")

SEED <- 1234
set.seed(SEED)

required_pkgs <- c("ggplot2")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall before rerunning, e.g. install.packages('ggplot2').")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

fmt <- format_number

clean_chr <- function(x) {
  trimws(as.character(x))
}

standardize_timepoint <- function(x) {
    z <- tolower(trimws(as.character(x)))
    z <- gsub("[_\\.-]+", " ", z)
    z <- gsub("\\s+", " ", z)
    out <- rep(NA_character_, length(z))
    out[grepl("uninjured|control|sham|^0$", z)] <- "Uninjured"
    out[grepl("15\\s*min|15min|0\\.25\\s*h", z)] <- "15 min"
    out[grepl("(^|[^0-9])1\\s*dpi|day\\s*1$|^1\\s*d$|^d1$", z)] <- "1 dpi"
    out[grepl("(^|[^0-9])3\\s*dpi|day\\s*3$|^3\\s*d$|^d3$", z)] <- "3 dpi"
    out[grepl("(^|[^0-9])7\\s*dpi|day\\s*7$|^7\\s*d$|^d7$", z)] <- "7 dpi"
    out[grepl("14\\s*dpi|day\\s*14$|^14\\s*d$|^d14$", z)] <- "14 dpi"
    out[grepl("28\\s*dpi|day\\s*28$|^28\\s*d$|^d28$", z)] <- "28 dpi"
    out[grepl("42\\s*dpi|day\\s*42$|^42\\s*d$|^d42$", z)] <- "42 dpi"
    out
}

integer_like_fraction <- function(x, tol = 1e-08) {
  x <- x[is.finite(x)]
  if (length(x) == 0)
    return(NA_real_)
  mean(abs(x - round(x)) < tol)
}

input_tbl <- data.frame(file = c(EXPR_FILE, META_FILE, SCRNA_FILE), exists = file.exists(c(EXPR_FILE,
  META_FILE, SCRNA_FILE)), size_MB = sapply(c(EXPR_FILE, META_FILE, SCRNA_FILE), function(f) {
  if (file.exists(f)) {
    round(file.info(f)$size/1024^2, 3)
  } else {
    NA_real_
  }
}), stringsAsFactors = FALSE)

write.csv(input_tbl, file.path(OUT_DIR, paste0(PREFIX, "_input_files.csv")), row.names = FALSE)

if (!file.exists(EXPR_FILE)) {
  stop("Expression file not found:\n", EXPR_FILE)
}

if (!file.exists(META_FILE)) {
  stop("Sample metadata file not found:\n", META_FILE)
}

print(input_tbl)

expr_raw <- read.csv(EXPR_FILE, check.names = FALSE, stringsAsFactors = FALSE)

meta_raw <- read.csv(META_FILE, check.names = FALSE, stringsAsFactors = FALSE)

meta_names_original <- colnames(meta_raw)

blank_meta_names <- is.na(meta_names_original) | trimws(meta_names_original) == ""

meta_names_safe <- meta_names_original
meta_names_safe[blank_meta_names] <- paste0("unnamed_col_", which(blank_meta_names))

meta_names_safe <- make.unique(meta_names_safe, sep = "_dup")

colnames(meta_raw) <- meta_names_safe

write.csv(data.frame(column_index = seq_along(meta_names_original), original_name = meta_names_original,
  safe_name = meta_names_safe, stringsAsFactors = FALSE), file.path(OUT_DIR, paste0(PREFIX, "_metadata_column_name_map.csv")),
  row.names = FALSE, fileEncoding = "UTF-8")

if (nrow(expr_raw) == 0 || ncol(expr_raw) < 2) {
  stop("Expression CSV is empty or has fewer than 2 columns.")
}

if (nrow(meta_raw) == 0 || ncol(meta_raw) == 0) {
  stop("Sample metadata CSV is empty.")
}

writeLines(c("Expression columns:", paste(colnames(expr_raw), collapse = "\n"), "", "Metadata columns:",
  paste(colnames(meta_raw), collapse = "\n")), con = file.path(OUT_DIR, paste0(PREFIX, "_column_names.txt")))

expr_col_samples <- colnames(expr_raw)[-1]

meta_overlap <- data.frame(metadata_column = colnames(meta_raw), n_overlap_with_expr_columns = NA_integer_,
  overlap_fraction_of_metadata = NA_real_, stringsAsFactors = FALSE)

for (j in seq_along(meta_raw)) {
    vals <- unique(clean_chr(meta_raw[[j]]))
    vals <- vals[!is.na(vals) & vals != ""]
    n_hit <- sum(vals %in% expr_col_samples)
    meta_overlap$n_overlap_with_expr_columns[j] <- n_hit
    meta_overlap$overlap_fraction_of_metadata[j] <- if (length(vals) > 0) {
        n_hit/length(vals)
    }
    else {
        0
    }
}

write.csv(meta_overlap, file.path(OUT_DIR, paste0(PREFIX, "_metadata_sampleID_overlap.csv")), row.names = FALSE)

best_j <- which.max(meta_overlap$n_overlap_with_expr_columns)

best_overlap <- meta_overlap$n_overlap_with_expr_columns[best_j]

sample_id_idx <- best_j

orientation <- "genes_by_samples"

if (length(best_overlap) == 0 || is.na(best_overlap) || best_overlap == 0) {

  expr_first_col <- clean_chr(expr_raw[[1]])

  reverse_overlap <- data.frame(metadata_column = colnames(meta_raw), n_overlap_with_expr_first_column = NA_integer_,
    stringsAsFactors = FALSE)

  for (j in seq_along(meta_raw)) {

    vals <- unique(clean_chr(meta_raw[[j]]))

    vals <- vals[!is.na(vals) & vals != ""]

    reverse_overlap$n_overlap_with_expr_first_column[j] <- sum(vals %in% expr_first_col)
  }

  best_rev_j <- which.max(reverse_overlap$n_overlap_with_expr_first_column)

  best_rev <- reverse_overlap$n_overlap_with_expr_first_column[best_rev_j]

  write.csv(reverse_overlap, file.path(OUT_DIR, paste0(PREFIX, "_metadata_reverse_overlap.csv")), row.names = FALSE)

  if (length(best_rev) > 0 && is.finite(best_rev) && best_rev > 0) {
    orientation <- "samples_by_genes"
    best_j <- best_rev_j
    sample_id_idx <- best_rev_j
    best_overlap <- best_rev
  } else {
    stop(paste0("Could not match any metadata column to the expression matrix.\n", "Inspect:\n",
      file.path(OUT_DIR, paste0(PREFIX, "_column_names.txt")), "\nand the overlap CSV files."))
  }
}

sample_id_col <- colnames(meta_raw)[sample_id_idx]

sample_id_col_original <- meta_names_original[sample_id_idx]

cat("Detected metadata sample-ID column: ", sample_id_col, ifelse(is.na(sample_id_col_original) || trimws(sample_id_col_original) ==
  "", " [original header was blank]", paste0(" [original header: ", sample_id_col_original, "]")),
  "\n", sep = "")

if (orientation == "genes_by_samples") {
    gene_id <- clean_chr(expr_raw[[1]])
    expr_df <- expr_raw[, -1, drop = FALSE]
    sample_names <- colnames(expr_df)
    expr_num <- as.data.frame(lapply(expr_df, function(x) suppressWarnings(as.numeric(x))), check.names = FALSE)
    expr_mat <- as.matrix(expr_num)
    rownames(expr_mat) <- make.unique(gene_id)
} else {
    sample_names <- clean_chr(expr_raw[[1]])
    gene_names <- colnames(expr_raw)[-1]
    tmp_df <- expr_raw[, -1, drop = FALSE]
    tmp_num <- as.data.frame(lapply(tmp_df, function(x) suppressWarnings(as.numeric(x))), check.names = FALSE)
    tmp_mat <- as.matrix(tmp_num)
    rownames(tmp_mat) <- make.unique(sample_names)
    expr_mat <- t(tmp_mat)
    rownames(expr_mat) <- make.unique(gene_names)
}

storage.mode(expr_mat) <- "numeric"

keep_gene <- rowSums(is.finite(expr_mat)) > 0

expr_mat <- expr_mat[keep_gene, , drop = FALSE]

if (anyNA(expr_mat)) {
    na_n <- sum(is.na(expr_mat))
    warning("Expression matrix contains ", fmt(na_n), " NA values after numeric conversion.")
}

meta_raw[[sample_id_idx]] <- clean_chr(meta_raw[[sample_id_idx]])

expr_samples <- colnames(expr_mat)

if (!all(expr_samples %in% meta_raw[[sample_id_idx]])) {
    missing_meta_samples <- setdiff(expr_samples, meta_raw[[sample_id_idx]])
    writeLines(missing_meta_samples, con = file.path(OUT_DIR, paste0(PREFIX, "_expression_samples_missing_metadata.txt")))
    stop("Some expression-matrix samples are absent from metadata. ", "See *_expression_samples_missing_metadata.txt")
}

meta <- meta_raw[match(expr_samples, meta_raw[[sample_id_idx]]), , drop = FALSE]

rownames(meta) <- expr_samples

if (anyDuplicated(meta[[sample_id_idx]]) > 0) {
  stop("Duplicate sample IDs detected in metadata column: ", sample_id_col)
}

sample_match_audit <- data.frame(expression_sample_id = expr_samples, metadata_sample_id = meta[[sample_id_idx]],
  exact_match = expr_samples == meta[[sample_id_idx]], stringsAsFactors = FALSE)

write.csv(sample_match_audit, file.path(OUT_DIR, paste0(PREFIX, "_sampleID_match_audit.csv")), row.names = FALSE,
  fileEncoding = "UTF-8")

time_scan <- data.frame(metadata_column = colnames(meta), n_standardized_timepoints = 0L, stringsAsFactors = FALSE)

time_standardized_candidates <- vector("list", ncol(meta))

for (j in seq_along(meta)) {
    xx <- standardize_timepoint(meta[[j]])
    time_standardized_candidates[[j]] <- xx
    time_scan$n_standardized_timepoints[j] <- sum(!is.na(xx))
}

write.csv(time_scan, file.path(OUT_DIR, paste0(PREFIX, "_metadata_timepoint_scan.csv")), row.names = FALSE)

time_candidate_idx <- setdiff(seq_len(ncol(meta)), sample_id_idx)

if (length(time_candidate_idx) > 0) {
  best_time_j <- time_candidate_idx[which.max(time_scan$n_standardized_timepoints[time_candidate_idx])]
} else {
  best_time_j <- which.max(time_scan$n_standardized_timepoints)
}

n_time_hits <- time_scan$n_standardized_timepoints[best_time_j]

if (length(n_time_hits) == 0 || is.na(n_time_hits) || n_time_hits == 0) {
    meta$timepoint_std <- NA_character_
    time_col <- NA_character_
    warning("No metadata column could be confidently standardized to expected time points.")
} else {
    time_col <- colnames(meta)[best_time_j]
    meta$timepoint_original <- meta[[time_col]]
    meta$timepoint_std <- time_standardized_candidates[[best_time_j]]
    meta$timepoint_std <- factor(meta$timepoint_std, levels = EXPECTED_TIME_LEVELS)
}

finite_vals <- expr_mat[is.finite(expr_mat)]

if (length(finite_vals) == 0) {
  stop("No finite numeric values found in expression matrix.")
}

min_expr <- min(finite_vals)

max_expr <- max(finite_vals)

pct_negative <- 100 * mean(finite_vals < 0)

integer_frac <- integer_like_fraction(finite_vals)

raw_count_like <- (min_expr >= 0 && is.finite(integer_frac) && integer_frac >= 0.995)

matrix_qc <- data.frame(n_genes = nrow(expr_mat), n_samples = ncol(expr_mat), min_value = min_expr, max_value = max_expr,
  percent_negative_values = pct_negative, integer_like_fraction = integer_frac, raw_count_like = raw_count_like,
  sample_id_column = sample_id_col, sample_id_column_original = ifelse(is.na(sample_id_col_original),
    "", sample_id_col_original), timepoint_column = ifelse(is.na(time_col), "", time_col), orientation_detected = orientation,
  stringsAsFactors = FALSE)

write.csv(matrix_qc, file.path(OUT_DIR, paste0(PREFIX, "_matrix_QC.csv")), row.names = FALSE)

print(matrix_qc)

sample_qc <- data.frame(sample_id = colnames(expr_mat), library_size = colSums(expr_mat, na.rm = TRUE),
  detected_genes = colSums(expr_mat > 0, na.rm = TRUE), stringsAsFactors = FALSE)

if ("timepoint_std" %in% colnames(meta)) {
  sample_qc$timepoint <- as.character(meta$timepoint_std)
} else {
  sample_qc$timepoint <- NA_character_
}

write.csv(sample_qc, file.path(OUT_DIR, paste0(PREFIX, "_sample_QC.csv")), row.names = FALSE)

if ("timepoint_std" %in% colnames(meta) && any(!is.na(meta$timepoint_std))) {
    time_summary <- as.data.frame(table(factor(meta$timepoint_std, levels = EXPECTED_TIME_LEVELS), useNA = "ifany"))
    colnames(time_summary) <- c("timepoint", "n_samples")
} else {
    time_summary <- data.frame(timepoint = character(), n_samples = integer())
}

write.csv(time_summary, file.path(OUT_DIR, paste0(PREFIX, "_timepoint_sample_counts.csv")), row.names = FALSE)

print(time_summary)

qc_mat <- log2(pmax(expr_mat, 0) + 1)

gene_var <- apply(qc_mat, 1, var, na.rm = TRUE)

gene_var[!is.finite(gene_var)] <- 0

n_top <- min(3000, sum(gene_var > 0))

if (n_top < 100) {
  stop("Fewer than 100 variable genes available for QC PCA.")
}

top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_top)]

pca_input <- t(qc_mat[top_genes, , drop = FALSE])

pca_col_sd <- apply(pca_input, 2, sd, na.rm = TRUE)

pca_input <- pca_input[, is.finite(pca_col_sd) & pca_col_sd > 0, drop = FALSE]

pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)

var_exp <- (pca$sdev^2/sum(pca$sdev^2)) * 100

pca_df <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)

if ("timepoint_std" %in% colnames(meta)) {
  pca_df$timepoint <- as.character(meta[pca_df$sample_id, "timepoint_std"])
} else {
  pca_df$timepoint <- NA_character_
}

write.csv(pca_df, file.path(OUT_DIR, paste0(PREFIX, "_PCA_coordinates.csv")), row.names = FALSE)

sample_cor <- cor(qc_mat, use = "pairwise.complete.obs", method = "pearson")

write.csv(sample_cor, file.path(OUT_DIR, paste0(PREFIX, "_sample_correlation.csv")), row.names = TRUE)

pdf(file.path(OUT_DIR, paste0(PREFIX, "_diagnostics.pdf")), width = 10, height = 7.5, onefile = TRUE)

p1 <- ggplot(sample_qc, aes(x = reorder(sample_id, library_size), y = library_size, fill = timepoint)) +
  geom_col(width = 0.8) + coord_flip() + labs(title = "Bulk RNA-seq sample library size", x = NULL,
  y = "Column sum", fill = "Time") + theme_classic(base_size = 10)

print(p1)

p2 <- ggplot(sample_qc, aes(x = reorder(sample_id, detected_genes), y = detected_genes, fill = timepoint)) +
  geom_col(width = 0.8) + coord_flip() + labs(title = "Detected genes per bulk sample", x = NULL, y = "Genes with expression > 0",
  fill = "Time") + theme_classic(base_size = 10)

print(p2)

p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint, label = sample_id)) + geom_point(size = 3) +
  geom_text(size = 2.5, vjust = -0.7, check_overlap = TRUE) + labs(title = "Bulk RNA-seq PCA (QC)",
  x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2]), color = "Time") +
  theme_classic(base_size = 10)

print(p3)

op <- par(mar = c(7, 7, 3, 2))

heatmap(sample_cor, Rowv = NA, Colv = NA, scale = "none", margins = c(8, 8), main = "Sample-sample Pearson correlation")

par(op)

dev.off()

write.csv(meta, file.path(OUT_DIR, paste0(PREFIX, "_metadata_clean.csv")), row.names = TRUE)

saveRDS(expr_mat, file.path(OUT_DIR, paste0(PREFIX, "_expression_matrix.rds")), compress = FALSE)

saveRDS(meta, file.path(OUT_DIR, paste0(PREFIX, "_metadata_clean.rds")), compress = FALSE)

saveRDS(list(expression = expr_mat, metadata = meta, sample_qc = sample_qc, matrix_qc = matrix_qc, expected_time_levels = EXPECTED_TIME_LEVELS),
  file.path(OUT_DIR, paste0(PREFIX, "_bulk_QC_bundle.rds")), compress = FALSE)

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = file.path(OUT_DIR, paste0(PREFIX, "_REPORT.txt")))

capture.output(sessionInfo(), file = file.path(OUT_DIR, paste0(PREFIX, "_sessionInfo.txt")))

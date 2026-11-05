options(stringsAsFactors = FALSE, width = 160, scipen = 30)

ROOT <- "D:/zz"
OUT <- file.path(ROOT, "results", "141_Figure6E_TF_gene_correlation_robustness_v03")
VERSION <- "v03"
ORIGINAL_RHO_TOL <- 1e-6
RESIDUAL_RELATIVE_NORM_TOL <- 1e-8
RUN_SELF_TESTS <- TRUE
EXP_SCALE <- 10000
TARGET <- "Clec7a+ lipid-handling Microglia"
TARGET_LEGACY <- "Clec7a+ lipid-associated Microglia"
PLOT_FIGURE <- TRUE

P <- list(original_object = file.path(ROOT, "results/20_AUCell_lipid_metabolism_scoring/20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds"), 
    metadata = file.path(ROOT, "results/84_pySCENIC_input_preparation/09_84_pySCENIC_cell_metadata_v01.csv.gz"), 
    tf_roster = file.path(ROOT, "results/90_pySCENIC_top10_TF_network_signature_correlation/02_90_top10_RSS_TFs_v05.csv"), 
    gene_panel = file.path(ROOT, "results/90_pySCENIC_top10_TF_network_signature_correlation/07_90_Figure6E_key_functional_gene_panel_v05.csv"), 
    original_target = file.path(ROOT, "results/90_pySCENIC_top10_TF_network_signature_correlation/08_90_Figure6E_TF_functional_gene_correlation_target_cells_v05.csv"))

if (!dir.exists(ROOT)) stop("Project directory not found: ", ROOT, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT, "141_Figure6E_TF_gene_correlation_robustness_v03.log")

log_msg <- function(...) {
    s <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
    cat(s, "\n")
    cat(s, "\n", file = LOG, append = TRUE)
    invisible(s)
}

fail <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

write_csv <- function(x, filename) {
  path <- file.path(OUT, filename)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

write_txt <- function(x, filename) {
  path <- file.path(OUT, filename)
  writeLines(as.character(x), path, useBytes = TRUE)
  invisible(path)
}

read_csv <- function(path, select = NULL) {
    if (!file.exists(path)) 
        fail("Required file missing: ", path)
    if (requireNamespace("data.table", quietly = TRUE)) {
        z <- tryCatch(data.table::fread(path, select = select, data.table = FALSE, check.names = FALSE, 
            showProgress = FALSE), error = function(e) NULL)
        if (!is.null(z)) 
            return(z)
    }
    z <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        utils::read.csv(gzfile(path), check.names = FALSE)
    }
    else {
        utils::read.csv(path, check.names = FALSE)
    }
    if (!is.null(select)) {
        absent <- setdiff(select, names(z))
        if (length(absent)) {
            fail("Missing selected columns in ", path, ": ", paste(absent, collapse = ", "))
        }
        z <- z[, select, drop = FALSE]
    }
    z
}

require_fields <- function(z, keys, label) {
  missing <- setdiff(keys, names(z))
  if (length(missing)) {
    fail(label, " missing columns: ", paste(missing, collapse = ", "))
  }
}

num <- function(x) suppressWarnings(as.numeric(x))

calc_p <- function(r, n, n_cov) {
    if (!is.finite(r) || n <= n_cov + 2L) 
        return(NA_real_)
    r <- min(max(r, -0.999999999), 0.999999999)
    2 * stats::pt(-abs(r) * sqrt((n - n_cov - 2)/pmax(1 - r * r, .Machine$double.eps)), df = n - n_cov - 
        2)
}

blank_pair <- function(tfs, genes, reason, n = 0L) {
  o <- expand.grid(TF = tfs, key_gene = genes, stringsAsFactors = FALSE)
  o$n_cells <- as.integer(n)
  o$rho_raw <- NA_real_
  o$rho_adj_technical <- NA_real_
  o$p_raw_cell_approx <- NA_real_
  o$p_adj_technical_cell_approx <- NA_real_
  o$FDR_raw_cell_approx <- NA_real_
  o$FDR_adj_technical_cell_approx <- NA_real_
  o$n_covariates_technical <- NA_integer_
  o$raw_TF_nonconstant <- NA
  o$raw_gene_nonconstant <- NA
  o$technical_TF_residual_valid <- NA
  o$technical_gene_residual_valid <- NA
  o$raw_estimate_status <- reason
  o$technical_estimate_status <- reason
  o$estimate_status <- reason
  o
}

residual_variation_flags <- function(rank_matrix, residual_matrix, originally_nonconstant, tolerance = RESIDUAL_RELATIVE_NORM_TOL) {
    rank_matrix <- as.matrix(rank_matrix)
    residual_matrix <- as.matrix(residual_matrix)
    centered <- sweep(rank_matrix, 2, colMeans(rank_matrix), FUN = "-")
    original_ss <- colSums(centered * centered)
    residual_ss <- colSums(residual_matrix * residual_matrix)
    ratio <- rep(NA_real_, ncol(rank_matrix))
    valid_base <- is.finite(original_ss) & original_ss > 0 & originally_nonconstant
    ratio[valid_base] <- sqrt(pmax(residual_ss[valid_base], 0)/original_ss[valid_base])
    ok <- valid_base & is.finite(ratio) & ratio > tolerance
    names(ok) <- colnames(rank_matrix)
    list(ok = ok, ratio = setNames(ratio, colnames(rank_matrix)))
}

pair_cor <- function(expr, tfs, genes, depth, features) {
  n <- ncol(expr)
  if (n < 4L) return(blank_pair(tfs, genes, "TOO_FEW_CELLS", n))

  xx <- t(expr[unique(c(tfs, genes)), , drop = FALSE])
  xnames <- colnames(xx)

  ranks <- apply(xx, 2, rank, ties.method = "average")
  if (is.null(dim(ranks))) {
    ranks <- matrix(ranks, ncol = 1, dimnames = list(NULL, xnames))
  }
  colnames(ranks) <- xnames

  raw_nonconstant <- vapply(
    seq_len(ncol(xx)),
    function(k) {
      v <- xx[, k]
      all(is.finite(v)) && any(v != v[1L])
    },
    logical(1)
  )
  names(raw_nonconstant) <- xnames

  tech <- data.frame(
    log1p_UMI = rank(log1p(pmax(depth, 0))),
    log1p_detected = rank(log1p(pmax(features, 0)))
  )

  qr_names <- character(0)
  design <- matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "intercept"))

  for (nm in names(tech)) {
    v <- num(tech[[nm]])
    if (
      all(is.finite(v)) &&
      stats::sd(v) > 0 &&
      qr(cbind(design, v))$rank > ncol(design)
    ) {
      design <- cbind(design, v)
      colnames(design)[ncol(design)] <- nm
      qr_names <- c(qr_names, nm)
    }
  }

  qtech <- qr(design)
  residual_tech <- qr.resid(qtech, ranks)
  tech_flags <- residual_variation_flags(
    ranks,
    residual_tech,
    raw_nonconstant
  )

  corr_raw <- suppressWarnings(stats::cor(ranks, method = "pearson"))
  corr_tech <- suppressWarnings(stats::cor(residual_tech, method = "pearson"))

  o <- expand.grid(TF = tfs, key_gene = genes, stringsAsFactors = FALSE)
  o$n_cells <- n
  o$rho_raw <- mapply(function(a, b) corr_raw[a, b], o$TF, o$key_gene)
  o$rho_adj_technical <- mapply(
    function(a, b) corr_tech[a, b],
    o$TF,
    o$key_gene
  )
  o$n_covariates_technical <- length(qr_names)
  o$raw_TF_nonconstant <- unname(raw_nonconstant[o$TF])
  o$raw_gene_nonconstant <- unname(raw_nonconstant[o$key_gene])
  o$technical_TF_residual_valid <- unname(tech_flags$ok[o$TF])
  o$technical_gene_residual_valid <- unname(tech_flags$ok[o$key_gene])

  raw_ok <- o$raw_TF_nonconstant &
    o$raw_gene_nonconstant &
    is.finite(o$rho_raw)
  o$rho_raw[!raw_ok] <- NA_real_

  tech_ok <- raw_ok &
    o$technical_TF_residual_valid &
    o$technical_gene_residual_valid &
    is.finite(o$rho_adj_technical)
  tech_ok[is.na(tech_ok)] <- FALSE
  o$rho_adj_technical[!tech_ok] <- NA_real_

  o$raw_estimate_status <- ifelse(
    raw_ok,
    "ESTIMATED",
    ifelse(
      !o$raw_TF_nonconstant & !o$raw_gene_nonconstant,
      "BOTH_ORIGINAL_EXPRESSION_CONSTANT",
      ifelse(
        !o$raw_TF_nonconstant,
        "TF_ORIGINAL_EXPRESSION_CONSTANT",
        ifelse(
          !o$raw_gene_nonconstant,
          "GENE_ORIGINAL_EXPRESSION_CONSTANT",
          "RAW_CORRELATION_UNDEFINED"
        )
      )
    )
  )

  o$technical_estimate_status <- ifelse(
    tech_ok,
    "ESTIMATED",
    ifelse(
      !raw_ok,
      "NOT_ESTIMABLE_ORIGINAL_EXPRESSION_CONSTANT_OR_UNDEFINED",
      ifelse(
        !o$technical_TF_residual_valid,
        "TF_TECHNICAL_RESIDUAL_DEGENERATE",
        ifelse(
          !o$technical_gene_residual_valid,
          "GENE_TECHNICAL_RESIDUAL_DEGENERATE",
          "TECHNICAL_CORRELATION_UNDEFINED"
        )
      )
    )
  )

  o$p_raw_cell_approx <- vapply(
    o$rho_raw,
    calc_p,
    numeric(1),
    n = n,
    n_cov = 0L
  )
  o$p_adj_technical_cell_approx <- vapply(
    o$rho_adj_technical,
    calc_p,
    numeric(1),
    n = n,
    n_cov = length(qr_names)
  )
  o$FDR_raw_cell_approx <- stats::p.adjust(o$p_raw_cell_approx, method = "BH")
  o$FDR_adj_technical_cell_approx <- stats::p.adjust(
    o$p_adj_technical_cell_approx,
    method = "BH"
  )

  o$estimate_status <- o$technical_estimate_status
  if (length(qr_names) < 2L) {
    o$estimate_status[o$estimate_status == "ESTIMATED"] <-
      "ESTIMATED_WITH_FEWER_TECHNICAL_COVARIATES"
  }

  orphan <- !is.finite(o$rho_raw) &
    (is.finite(o$rho_adj_technical) | is.finite(o$p_adj_technical_cell_approx))
  if (any(orphan)) {
    fail(
      "Internal v03 validation failed: adjusted orphan estimates: ",
      sum(orphan)
    )
  }

  o
}

run_pair_cor_self_tests <- function() {
  n <- 80L
  set.seed(14103)

  xx <- rbind(
    TFvary = as.numeric(seq_len(n)),
    TFconst = rep(2, n),
    GeneVary = as.numeric(sample(seq_len(n))),
    GeneConst = rep(0, n)
  )

  d <- as.numeric(seq_len(n) * 10)
  f <- as.numeric(sample(seq_len(n)))

  a <- pair_cor(
    xx,
    c("TFvary", "TFconst"),
    c("GeneVary", "GeneConst"),
    d,
    f
  )

  g <- function(tf, gene) {
    a[a$TF == tf & a$key_gene == gene, , drop = FALSE]
  }

  p_const <- g("TFvary", "GeneConst")
  tf_const <- g("TFconst", "GeneVary")
  valid_pair <- g("TFvary", "GeneVary")

  stopifnot(
    nrow(p_const) == 1L,
    is.na(p_const$rho_raw),
    is.na(p_const$rho_adj_technical),
    is.na(p_const$p_adj_technical_cell_approx),
    is.na(tf_const$rho_raw),
    is.na(tf_const$rho_adj_technical),
    is.finite(valid_pair$rho_raw),
    is.na(valid_pair$rho_adj_technical),
    valid_pair$technical_estimate_status == "TF_TECHNICAL_RESIDUAL_DEGENERATE"
  )

  b <- pair_cor(
    xx,
    "TFvary",
    "GeneVary",
    sample(seq_len(n)) * 12,
    sample(seq_len(n)) * 3
  )
  stopifnot(is.finite(b$rho_raw[1]), is.finite(b$rho_adj_technical[1]))

  invisible(TRUE)
}

if (RUN_SELF_TESTS) {
  run_pair_cor_self_tests()
}

audit <- data.frame(key = names(P), path = unlist(P, use.names = FALSE), exists = file.exists(unlist(P, 
    use.names = FALSE)), stringsAsFactors = FALSE)
audit$bytes <- ifelse(audit$exists, as.numeric(file.info(audit$path)$size), NA_real_)
write_csv(audit, "00_141_INPUT_AUDIT_v03.csv")

if (any(!audit$exists)) {
  fail("Required Step90/84/20 input missing; see 00_141_INPUT_AUDIT_v03.csv")
}

if (!requireNamespace("SeuratObject", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
    fail("Required R packages SeuratObject and Matrix must be installed.")
}

roster <- read_csv(P$tf_roster)
panel <- read_csv(P$gene_panel)
original <- read_csv(P$original_target)

require_fields(roster, c("TF", "RSS_rank"), "Step90 roster")
require_fields(panel, c("selection_rank", "gene"), "Step90 functional panel")
require_fields(original, c("TF", "key_gene", "rho", "n_cells"), "Step90 original Figure6E")

roster <- roster[order(num(roster$RSS_rank)), , drop = FALSE]
panel <- panel[order(num(panel$selection_rank)), , drop = FALSE]

tfs <- as.character(roster$TF)
genes <- as.character(panel$gene)

if (length(tfs) != 10L || length(genes) != 20L || anyDuplicated(tfs) || anyDuplicated(genes)) {
    fail("Figure6E roster must have 10 unique TFs and 20 unique genes")
}

if (nrow(original) != 200L || anyDuplicated(paste(original$TF, original$key_gene))) {
 fail("Figure6E v05 table must contain exactly 200 unique TF-gene pairs")
}

if (!setequal(tfs, original$TF) || !setequal(genes, original$key_gene)) {
  fail("TF-gene panel mismatches original Figure6E; STOP to protect analysis comparability")
}

write_csv(data.frame(type = c(rep("TF", length(tfs)), rep("Functional_gene", length(genes))), name = c(tfs, 
    genes), order = c(seq_along(tfs), seq_along(genes))), "01_141_frozen_Figure6E_TF_gene_panel_v03.csv")

meta_keys <- c("cell", "myeloid_state_frozen", "nCount_RNA", "nFeature_RNA")
meta <- read_csv(P$metadata, select = meta_keys)
require_fields(meta, meta_keys, "Step84 metadata")

if (anyDuplicated(meta$cell)) {
  fail("Step84 metadata contains duplicate cell identifiers")
}

meta$cell <- as.character(meta$cell)
meta$myeloid_state_frozen <- as.character(meta$myeloid_state_frozen)
meta$myeloid_state_frozen[meta$myeloid_state_frozen == TARGET_LEGACY] <- TARGET

obj <- readRDS(P$original_object)

if (!"RNA" %in% names(obj@assays)) {
  fail("RNA assay missing from Step20 RDS")
}

counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
if (!inherits(counts, "Matrix")) {
  fail("Step20 RNA counts is not a Matrix; inspect source object")
}
if (is.null(colnames(counts)) || is.null(rownames(counts))) {
  fail("Step20 RNA counts has no cell/gene names")
}

all_cells <- colnames(counts)
if (length(all_cells) != 58533L) {
 fail("Step20 source cell count differs from Step90: ", length(all_cells))
}

idx_meta <- match(all_cells, meta$cell)
if (anyNA(idx_meta)) {
  fail("Step84 metadata does not cover all original Step20 RNA cells")
}
meta <- meta[idx_meta, , drop = FALSE]

if (anyDuplicated(meta$cell)) {
  fail("Step20 RNA metadata has duplicated cells")
}

state_rds <- as.character(obj@meta.data[all_cells, "myeloid_state"])
state_rds[state_rds == TARGET_LEGACY] <- TARGET
if (any(is.na(state_rds)) || !identical(state_rds, meta$myeloid_state_frozen)) {
    fail("Step84 state labels differ from original Step20 object; do NOT silently continue")
}

missing_genes <- setdiff(unique(c(tfs, genes)), rownames(counts))
if (length(missing_genes)) {
  fail("Missing required RNA genes: ", paste(missing_genes, collapse = ", "))
}

lib <- as.numeric(Matrix::colSums(counts))
features <- as.numeric(Matrix::colSums(counts != 0))
valid <- is.finite(lib) & lib > 0 & is.finite(features) & features > 0

meta$RNA_library_from_counts <- lib
meta$RNA_features_from_counts <- features
meta$RNA_valid <- valid
meta$RNA_count_metadata_abs_difference <- abs(num(meta$nCount_RNA) - lib)
meta$RNA_feature_metadata_abs_difference <- abs(num(meta$nFeature_RNA) - features)
meta$is_target <- meta$myeloid_state_frozen == TARGET

expected_target_all <- sum(meta$is_target)
if (expected_target_all != 17713L) {
 fail("target state size changed: ", expected_target_all)
}

valid_target <- which(meta$is_target & valid)
if (length(valid_target) != 17709L) {
    fail("Expected 17,709 valid Figure6E target cells (Step90 v05), got ", length(valid_target), ". Check RNA counts and zero-library exclusions.")
}

needed <- unique(c(tfs, genes))
small <- as.matrix(counts[needed, which(valid), drop = FALSE])
expr <- log1p(sweep(small, 2, lib[valid], "/") * EXP_SCALE)
rownames(expr) <- needed
colnames(expr) <- all_cells[valid]
valid_meta <- meta[valid, , drop = FALSE]

rm(small, counts, obj)
invisible(gc(verbose = FALSE))

if (!identical(colnames(expr), valid_meta$cell)) {
  fail("Expression column order / metadata alignment failed")
}

target_idx <- which(valid_meta$is_target)

primary <- pair_cor(expr[, target_idx, drop = FALSE], tfs, genes, valid_meta$RNA_library_from_counts[target_idx], 
    valid_meta$RNA_features_from_counts[target_idx])

if (!all(is.finite(primary$rho_raw))) {
 fail("A Figure6E raw correlation is non-finite")
}

original$key <- paste(original$TF, original$key_gene, sep = "||")
primary$key <- paste(primary$TF, primary$key_gene, sep = "||")
ord <- match(primary$key, original$key)
if (anyNA(ord)) fail("Original Figure6E 200-pair matching failed")

primary$rho_step90 <- num(original$rho[ord])
primary$delta_vs_step90 <- primary$rho_raw - primary$rho_step90
primary$original_n_cells <- num(original$n_cells[ord])
primary$reproduction_pass <- abs(primary$delta_vs_step90) <= ORIGINAL_RHO_TOL & primary$n_cells == primary$original_n_cells
primary$key <- NULL

write_csv(primary[, c("TF", "key_gene", "n_cells", "original_n_cells", "rho_step90", "rho_raw", "delta_vs_step90", 
    "reproduction_pass")], "02_141_Figure6E_original_vs_reproduced_v03.csv")

if (!all(primary$reproduction_pass)) {
    fail("Original Figure6E reproduction FAILED; check 02_141 table. Maximum absolute error: ", signif(max(abs(primary$delta_vs_step90)), 
        5))
}

log_msg("Original Figure6E reproduced: maximum |delta rho| = ", max(abs(primary$delta_vs_step90)))

result <- primary[, c("TF", "key_gene", "n_cells", "rho_step90", "rho_raw", "rho_adj_technical", "p_raw_cell_approx", 
    "p_adj_technical_cell_approx", "FDR_raw_cell_approx", "FDR_adj_technical_cell_approx", "n_covariates_technical", 
    "raw_estimate_status", "technical_estimate_status", "estimate_status")]

result$delta_rho_after_technical_adjustment <- result$rho_adj_technical - result$rho_raw
result$abs_rho_attenuated <- abs(result$rho_adj_technical) < abs(result$rho_raw)
result$sign_changed_after_technical_adjustment <- sign(result$rho_adj_technical) != sign(result$rho_raw)

write_csv(result, "03_141_Figure6E_pair_level_raw_vs_technical_adjusted_v03.csv")

tf_summary <- do.call(rbind, lapply(tfs, function(tf) {
    z <- result[result$TF == tf, , drop = FALSE]
    data.frame(TF = tf, n_functional_genes = nrow(z), n_estimable_raw = sum(is.finite(z$rho_raw)), n_estimable_technical_adjusted = sum(is.finite(z$rho_adj_technical)), 
        n_positive_raw = sum(z$rho_raw > 0, na.rm = TRUE), n_positive_technical_adjusted = sum(z$rho_adj_technical > 
            0, na.rm = TRUE), mean_rho_raw = mean(z$rho_raw, na.rm = TRUE), mean_rho_technical_adjusted = mean(z$rho_adj_technical, 
            na.rm = TRUE), median_rho_raw = stats::median(z$rho_raw, na.rm = TRUE), median_rho_technical_adjusted = stats::median(z$rho_adj_technical, 
            na.rm = TRUE), median_abs_rho_raw = stats::median(abs(z$rho_raw), na.rm = TRUE), median_abs_rho_technical_adjusted = stats::median(abs(z$rho_adj_technical), 
            na.rm = TRUE), n_abs_rho_attenuated = sum(z$abs_rho_attenuated, na.rm = TRUE), n_sign_changed = sum(z$sign_changed_after_technical_adjustment, 
            na.rm = TRUE), stringsAsFactors = FALSE)
}))

write_csv(tf_summary, "04_141_TF_level_raw_vs_technical_adjusted_summary_v03.csv")

technical_qc <- data.frame(item = c("n_total_cells", "n_valid_RNA_cells", "n_target_cells_before_zero_library_filter", 
    "n_target_cells_used", "median_RNA_library_size_target", "median_detected_genes_target", "median_abs_metadata_vs_count_library_difference", 
    "median_abs_metadata_vs_count_feature_difference", "n_technical_covariates_retained"), value = c(nrow(meta), 
    sum(valid), sum(meta$is_target), length(target_idx), stats::median(valid_meta$RNA_library_from_counts[target_idx]), 
    stats::median(valid_meta$RNA_features_from_counts[target_idx]), stats::median(meta$RNA_count_metadata_abs_difference, 
        na.rm = TRUE), stats::median(meta$RNA_feature_metadata_abs_difference, na.rm = TRUE), unique(primary$n_covariates_technical)[1]), 
    stringsAsFactors = FALSE)

write_csv(technical_qc, "05_141_technical_complexity_QC_v03.csv")

if (PLOT_FIGURE) {
  make_heatmap <- function() {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    par(
      mfrow = c(1, 2),
      mar = c(8.2, 5.6, 4.1, 2.4),
      oma = c(0.1, 0, 1.4, 0)
    )

    lim <- max(
      0.1,
      max(
        abs(c(result$rho_raw, result$rho_adj_technical)),
        na.rm = TRUE
      )
    )

    pal <- grDevices::colorRampPalette(
      c("#276FBF", "white", "#CE4F3D")
    )(201)

    for (field in c("rho_raw", "rho_adj_technical")) {
      m <- matrix(
        NA_real_,
        length(tfs),
        length(genes),
        dimnames = list(tfs, genes)
      )

      for (k in seq_len(nrow(result))) {
        m[result$TF[k], result$key_gene[k]] <- result[[field]][k]
      }

      image(
        x = seq_len(ncol(m)),
        y = seq_len(nrow(m)),
        z = t(m),
        col = pal,
        zlim = c(-lim, lim),
        axes = FALSE,
        xlab = "",
        ylab = "",
        main = switch(
          field,
          rho_raw = "Original pooled Spearman",
          rho_adj_technical = "Technical-adjusted partial Spearman"
        )
      )

      axis(
        1,
        at = seq_len(ncol(m)),
        labels = colnames(m),
        las = 2,
        cex.axis = 0.47,
        tick = FALSE
      )
      axis(
        2,
        at = seq_len(nrow(m)),
        labels = rownames(m),
        las = 2,
        cex.axis = 0.75,
        tick = FALSE
      )
      box()
    }

    mtext(
      "Figure 6E correlations: sensitivity to RNA technical complexity",
      outer = TRUE,
      side = 3,
      cex = 0.84
    )
  }

  grDevices::pdf(
    file.path(OUT, "06_141_Figure6E_technical_complexity_sensitivity_v03.pdf"),
    width = 12,
    height = 8.5
  )
  make_heatmap()
  grDevices::dev.off()

  grDevices::png(
    file.path(OUT, "06_141_Figure6E_technical_complexity_sensitivity_v03.png"),
    width = 3400,
    height = 2200,
    res = 260
  )
  make_heatmap()
  grDevices::dev.off()

}

orphan_n <- sum(!is.finite(result$rho_raw) & is.finite(result$rho_adj_technical))

qc <- data.frame(check = c("required_input_files", "frozen_Figure6E_10_TF_20_gene_panel", "original_Figure6E_200_coefficients_reproduced", 
    "zero_library_exclusion_reproduced_17709", "both_technical_covariates_available", "no_adjusted_coefficient_when_raw_undefined", 
    "sample_or_unit_covariates_used", "source_stage_treatment_covariates_used"), status = c(if (all(audit$exists)) "PASS" else "FAIL", 
    if (length(tfs) == 10L && length(genes) == 20L) "PASS" else "FAIL", if (all(primary$reproduction_pass)) "PASS" else "FAIL", 
    if (length(target_idx) == 17709L) "PASS" else "FAIL", if (all(primary$n_covariates_technical == 2L)) "PASS" else "REVIEW", 
    if (orphan_n == 0L) "PASS" else "FAIL", "NO", "NO"), details = c("All original D:/zz files available; no ZIP used", 
    "Step90 v05 identifiers used", paste0("max_abs_rho_delta=", format(max(abs(primary$delta_vs_step90)), 
        digits = 8)), "Target atlas state 17713 before zero-library filtering", "Ranked log1p total RNA UMI and ranked log1p detected genes", 
    paste0("finite adjusted with undefined raw n=", orphan_n), "v03 intentionally performs no sample/unit fixed-effect adjustment", 
    "v03 intentionally performs no source/stage/treatment adjustment"), stringsAsFactors = FALSE)

write_csv(qc, "07_141_FINAL_QC_v03.csv")

write_txt("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    "08_141_REPORT_v03.txt")
write_txt(capture.output(sessionInfo()), "09_141_sessionInfo_v03.txt")

produced <- list.files(OUT, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
index <- data.frame(file = basename(produced), relative_path = substring(produced, nchar(OUT) + 2L), 
    size_bytes = as.numeric(file.info(produced)$size), stringsAsFactors = FALSE)
write_csv(index, "00_141_OUTPUT_INDEX_v03.csv")

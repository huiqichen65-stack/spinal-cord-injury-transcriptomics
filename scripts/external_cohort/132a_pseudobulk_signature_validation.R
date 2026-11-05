source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

ROOT <- PROJECT_ROOT
VERSION <- "v01"
IN119 <- file.path(ROOT, "results", "119_GSE281205_scDblFinder", "RDS")
IN131 <- file.path(ROOT, "results", "131_GSE281205_myeloid_subtype_annotation_and_microglia_groups",
  "13_131_READY_FOR_STEP132_original_UMI_cell_manifest_v01.csv")
SIGNATURE_FILE <- file.path(ROOT, "results", "21_Clec7a_DEG_volcano_myeloid", "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv")
OUT <- file.path(ROOT, "results", "132_GSE281205_original_UMI_pseudobulk_signature_validation")
RDS_DIR <- file.path(OUT, "RDS")
FIG_DIR <- file.path(OUT, "figures")
LOG_DIR <- file.path(OUT, "logs")

GSMS <- paste0("GSM", 8614662:8614667)
EXPECTED_DPI <- setNames(rep(c(14L, 28L), each = 3L), GSMS)
MIN_CELLS_PER_GSM_GROUP <- 20L
EXPECTED_N_SIGNATURE <- 1320L
EXPECTED_N_MANIFEST <- 28716L
PRIOR_COUNT <- 2

IDENTITY_MARKERS <- c("Clec7a", "Spp1", "Lpl", "Gpnmb", "Cst7", "Apoe", "Trem2", "Tyrobp", "P2ry12",
  "Tmem119", "Fcrls", "Hexb", "Cx3cr1", "Lyz2", "Ms4a7", "Mrc1", "Cd74", "H2-Aa", "Flt3", "Cd209a",
  "Pclaf", "Birc5", "Mki67", "Top2a")
CONTRASTS <- c("contrast_A_target_vs_other_microglia_131", "contrast_B_target_vs_all_other_highconf_myeloid_131",
  "contrast_C_target_vs_nonmicroglial_myeloid_131", "sensitivity_A_broad_microglia_lt50", "sensitivity_B_broad_all_myeloid_lt50",
  "sensitivity_A_partition_consistent")

stop_if <- function(cond, ...) {
  if (isTRUE(cond))
    stop(paste0(...), call. = FALSE)
  invisible(TRUE)
}
need <- function(cond, ...) stop_if(!isTRUE(cond), ...)
need(file.exists(IN131), "Missing Step131 input: ", IN131)
need(file.exists(SIGNATURE_FILE), "Missing frozen Step21 signature: ", SIGNATURE_FILE, "\nDo not substitute a reselected or new gene list.")
need(dir.exists(IN119), "Missing Step119 original UMI RDS directory: ", IN119)
for (pkg in c("Matrix", "edgeR", "limma")) {
  need(requireNamespace(pkg, quietly = TRUE), "Missing R package: ", pkg)
}
for (d in c(OUT, RDS_DIR, FIG_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(LOG_DIR, "00_132_run_log_v01.log")
cat("", file = LOG_FILE)
log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = LOG_FILE, append = TRUE, sep = "")
  invisible(msg)
}
write_csv <- function(d, filename) {
  need(is.data.frame(d), "Not a data.frame: ", filename)
  f <- file.path(OUT, filename)
  utils::write.csv(d, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(f)
}
write_txt <- function(v, filename) {
  writeLines(as.character(v), con = file.path(OUT, filename), useBytes = TRUE)
}
check_unique <- function(v, label) {
  need(!anyNA(v) && !any(!nzchar(as.character(v))), "Missing ", label)
  need(!anyDuplicated(v), "Duplicated ", label)
}
mean_or_na <- function(x) if (length(x) && any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
med_or_na <- function(x) if (length(x) && any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
safe_sd <- function(x) if (length(x) >= 2L) stats::sd(x) else NA_real_

for (p in c("Matrix", "edgeR", "limma")) log_msg("Package ", p, ": ", as.character(utils::packageVersion(p)))

signature <- utils::read.csv(SIGNATURE_FILE, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
need(all(c("gene", "avg_log2FC") %in% names(signature)), "Step21 signature needs gene and avg_log2FC fields")
signature$gene <- trimws(as.character(signature$gene))
signature$avg_log2FC <- suppressWarnings(as.numeric(signature$avg_log2FC))
need(nrow(signature) == EXPECTED_N_SIGNATURE, "Frozen input must have exactly 1,320 rows; found ", nrow(signature))
check_unique(signature$gene, "frozen signature gene symbols")
need(all(is.finite(signature$avg_log2FC)) && all(signature$avg_log2FC > 0), "Frozen signature must be an upregulated 1,320-gene list")
SIG_GENES <- signature$gene
SIG_MARKER_OVERLAP <- intersect(SIG_GENES, IDENTITY_MARKERS)
log_msg("Frozen signature n=", length(SIG_GENES), "; marker overlap n=", length(SIG_MARKER_OVERLAP),
  "; Step21 MD5=", unname(tools::md5sum(SIGNATURE_FILE)))

manifest <- utils::read.csv(IN131, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
mandatory <- c("cell_id", "gsm", "dpi", "primary_group_131", CONTRASTS)
need(all(mandatory %in% names(manifest)), "Step131 manifest missing columns: ", paste(setdiff(mandatory,
  names(manifest)), collapse = ", "))
need(nrow(manifest) == EXPECTED_N_MANIFEST, "Step131 expected 28,716 cells; got ", nrow(manifest))
check_unique(manifest$cell_id, "Step131 cell_id")
need(identical(sort(unique(manifest$gsm)), GSMS), "Step131 GSM set/order mismatch")
manifest$dpi <- as.integer(manifest$dpi)
need(all(manifest$dpi == EXPECTED_DPI[manifest$gsm]), "Step131 GSM vs dpi mismatch")
for (nm in CONTRASTS) {
  manifest[[nm]] <- as.character(manifest[[nm]])
  need(all(!is.na(manifest[[nm]]) & manifest[[nm]] %in% c("TARGET", "OTHER", "NOT_IN_CONTRAST")), "Invalid contrast values in ",
    nm)

  for (g in GSMS) {
    z <- manifest[manifest$gsm == g, nm]
    need(sum(z == "TARGET") >= MIN_CELLS_PER_GSM_GROUP && sum(z == "OTHER") >= MIN_CELLS_PER_GSM_GROUP,
      nm, "/", g, ": fewer than ", MIN_CELLS_PER_GSM_GROUP, " cells in TARGET or OTHER")
  }
}

for (g in GSMS) {
  rows <- manifest$gsm == g
  aa <- manifest[[CONTRASTS[1]]][rows] == "TARGET"
  bb <- manifest[[CONTRASTS[2]]][rows] == "TARGET"
  cc <- manifest[[CONTRASTS[3]]][rows] == "TARGET"
  need(identical(aa, bb) && identical(bb, cc), "Main target population mismatch for ", g)
}

write_csv(data.frame(input_role = c("Step119_RDS_dir", "Step131_cell_manifest", "Step21_frozen_signature"),
  input_path = c(IN119, IN131, SIGNATURE_FILE), md5 = c(NA_character_, unname(tools::md5sum(IN131)),
    unname(tools::md5sum(SIGNATURE_FILE))), stringsAsFactors = FALSE), "01_132_input_provenance_v01.csv")
write_csv(data.frame(gene = SIG_GENES, discovered_avg_log2FC = signature$avg_log2FC, marker_overlap = SIG_GENES %in%
  IDENTITY_MARKERS, stringsAsFactors = FALSE), "02_132_frozen_signature_metadata_v01.csv")

PB_CHECKPOINT <- file.path(RDS_DIR, "132_all_comparisons_per_GSM_raw_UMI_pseudobulk_checkpoint_v01.rds")
manifest_md5 <- unname(tools::md5sum(IN131))

if (file.exists(PB_CHECKPOINT)) {
  checkpoint <- readRDS(PB_CHECKPOINT)
  need(is.list(checkpoint) && all(c("manifest_md5", "pb", "feature_metadata", "cell_sizes") %in% names(checkpoint)),
    "Malformed checkpoint. Remove it and rerun: ", PB_CHECKPOINT)
  need(identical(checkpoint$manifest_md5, manifest_md5), "Manifest has changed since checkpoint; remove old checkpoint and rerun")
  need(identical(names(checkpoint$pb), CONTRASTS), "Checkpoint comparison definitions changed; remove checkpoint")
  pb <- checkpoint$pb
  feature_md <- checkpoint$feature_metadata
  cell_sizes <- checkpoint$cell_sizes
  rm(checkpoint)
} else {
  pb <- setNames(lapply(CONTRASTS, function(x) list()), CONTRASTS)
  size_rows <- list()
  feature_md <- NULL
  gene_ids <- NULL
  k <- 0L
  for (gsm in GSMS) {
    raw_file <- file.path(IN119, paste0("119_", gsm, "_singlet_raw_counts_v02.rds"))
    need(file.exists(raw_file), "Missing original INTEGER UMI file: ", raw_file)
    obj <- readRDS(raw_file)
    need(is.list(obj) && all(c("counts", "cell_metadata", "feature_metadata") %in% names(obj)), "Wrong Step119 list format: ",
      gsm)
    x <- obj$counts
    need(inherits(x, "sparseMatrix"), "Raw counts not sparse: ", gsm)
    x <- as(x, "dgCMatrix")
    need(!is.null(rownames(x)) && !is.null(colnames(x)), "Gene/cell dimnames missing: ", gsm)
    check_unique(rownames(x), paste(gsm, "raw gene IDs"))
    check_unique(colnames(x), paste(gsm, "raw cell IDs"))
    need(all(is.finite(x@x)) && all(x@x >= 0) && all(abs(x@x - round(x@x)) < 1e-08), "Noninteger or invalid raw UMI counts: ",
      gsm)
    md <- obj$cell_metadata
    feat <- obj$feature_metadata
    need(is.data.frame(md) && is.data.frame(feat), "Raw metadata not data.frames: ", gsm)
    need("cell_id" %in% names(md) && identical(as.character(md$cell_id), colnames(x)), "Cell metadata mismatch: ",
      gsm)
    need(all(c("gene_id", "gene_symbol") %in% names(feat)), "Feature metadata missing gene_id/gene_symbol: ",
      gsm)
    need(identical(as.character(feat$gene_id), rownames(x)), "Gene mapping mismatch: ", gsm)
    need(all(as.character(md$gsm) == gsm), "Raw GSM metadata error: ", gsm)
    if ("scDblFinder_class" %in% names(md))
      need(all(as.character(md$scDblFinder_class) == "singlet"), "Step119 file contains non-singlet cells: ",
        gsm)
    if (is.null(feature_md)) {
      feature_md <- feat[, c("gene_id", "gene_symbol"), drop = FALSE]
      gene_ids <- rownames(x)
    } else {
      need(identical(gene_ids, rownames(x)) && identical(as.character(feature_md$gene_symbol),
        as.character(feat$gene_symbol)), "Feature IDs/order/symbol mapping inconsistent across GSMs: ",
        gsm, ". Resolve mapping upstream, do not silently intersect genes.")
    }
    msub <- manifest[manifest$gsm == gsm, , drop = FALSE]
    pos <- match(msub$cell_id, colnames(x))
    need(!anyNA(pos), "Step131 cells missing from raw singlet UMI for ", gsm, " (count=", sum(is.na(pos)),
      ")")
    need(length(pos) == nrow(msub), "Internal cell matching error: ", gsm)
    for (ct in CONTRASTS) {
      for (group in c("OTHER", "TARGET")) {
        local <- pos[msub[[ct]] == group]
        need(length(local) >= MIN_CELLS_PER_GSM_GROUP, "Insufficient cells ", gsm, " ", ct, " ",
          group)
        vec <- Matrix::rowSums(x[, local, drop = FALSE])
        need(length(vec) == nrow(x) && all(is.finite(vec)) && all(vec >= 0) && all(abs(vec -
          round(vec)) < 1e-08), "Aggregation lost integer-valued original UMIs")
        col_id <- paste(gsm, group, sep = "__")
        pb[[ct]][[col_id]] <- vec
        k <- k + 1L
        size_rows[[k]] <- data.frame(comparison = ct, gsm = gsm, dpi = as.integer(EXPECTED_DPI[[gsm]]),
          group = group, n_cells = length(local), sum_original_UMI = sum(vec), n_genes_detected = sum(vec >
          0), source_raw_RDS = raw_file, stringsAsFactors = FALSE)
      }
    }
    rm(obj, x, md, feat, msub, pos)
    invisible(gc(verbose = FALSE))
  }
  cell_sizes <- do.call(rbind, size_rows)

  checkpoint <- list(manifest_md5 = manifest_md5, feature_metadata = feature_md, pb = pb, cell_sizes = cell_sizes,
    note = "Exact integer-valued original UMI sums from Step119, no DecontX")
  saveRDS(checkpoint, PB_CHECKPOINT, compress = "gzip")
  rm(checkpoint, size_rows)
}
need(is.data.frame(feature_md) && nrow(feature_md) > 10000L, "Bad feature mapping in checkpoint")
need(nrow(cell_sizes) == length(CONTRASTS) * 12L, "Bad count of pseudobulk groups")
need(all(cell_sizes$n_cells >= MIN_CELLS_PER_GSM_GROUP) && all(cell_sizes$sum_original_UMI > 0), "Low cells / zero original-UMI pseudobulk detected")
write_csv(cell_sizes, "03_132_pseudobulk_library_cell_and_UMI_QC_v01.csv")

gene_symbol <- trimws(as.character(feature_md$gene_symbol))
gene_id <- as.character(feature_md$gene_id)
bad_symbols <- is.na(gene_symbol) | !nzchar(gene_symbol)
gene_symbol[bad_symbols] <- gene_id[bad_symbols]
all_symbol_names <- unique(gene_symbol)
sig_map <- data.frame(gene = SIG_GENES, discovered_avg_log2FC = signature$avg_log2FC, mapped_in_raw_features = SIG_GENES %in%
  all_symbol_names, n_raw_feature_ids = as.integer(table(factor(gene_symbol, levels = SIG_GENES))),
  marker_overlap = SIG_GENES %in% IDENTITY_MARKERS, stringsAsFactors = FALSE)
write_csv(sig_map, "04_132_frozen_signature_gene_mapping_audit_v01.csv")
need(sum(sig_map$mapped_in_raw_features) >= 700L, "Less than 700/1320 signature genes mapped; inspect feature mappings")

summary_rows <- list()
per_gsm_rows <- list()
per_gene_rows <- list()
model_rows <- list()
camera_rows <- list()
model_num <- 0L

for (ct in CONTRASTS) {
  vecs <- pb[[ct]]
  need(length(vecs) == 12L, "Expected 12 GSM/group vectors: ", ct)
  for (dpi in c(14L, 28L)) {
    stage_gsms <- GSMS[EXPECTED_DPI[GSMS] == dpi]
    wanted <- as.vector(rbind(paste0(stage_gsms, "__OTHER"), paste0(stage_gsms, "__TARGET")))
    need(all(wanted %in% names(vecs)), "Missing some GSM x group pseudobulk vectors: ", ct, "/",
      dpi)
    original_counts <- do.call(cbind, vecs[wanted])
    rownames(original_counts) <- gene_id
    colnames(original_counts) <- wanted
    need(all(is.finite(original_counts)) && all(original_counts >= 0) && all(abs(original_counts -
      round(original_counts)) < 1e-08), "Raw pseudobulk is not integer-valued")

    counts <- rowsum(original_counts, group = gene_symbol, reorder = FALSE)
    need(!anyDuplicated(rownames(counts)), "Duplicated gene symbol after rowsum")
    need(identical(colnames(counts), wanted), "Pseudobulk column order changed")
    raw_sum_before <- colSums(original_counts)
    need(isTRUE(all.equal(as.numeric(colSums(counts)), as.numeric(raw_sum_before), tolerance = 1e-09)),
      "Raw UMI count conservation failed during symbol collapse")

    meta <- data.frame(library_id = wanted, gsm = factor(rep(stage_gsms, each = 2L), levels = stage_gsms),
      group = factor(rep(c("OTHER", "TARGET"), times = 3L), levels = c("OTHER", "TARGET")), dpi = dpi,
      stringsAsFactors = FALSE)
    rownames(meta) <- meta$library_id
    design <- stats::model.matrix(~gsm + group, data = meta)
    need(qr(design)$rank == ncol(design), "Rank-deficient paired design: ", ct, " / ", dpi)
    residual_df <- ncol(counts) - ncol(design)
    need(residual_df >= 2L, "Too few residual degrees of freedom")
    y_all <- edgeR::DGEList(counts = counts)
    keep <- edgeR::filterByExpr(y_all, design = design)
    need(sum(keep) > 500L, "Too few genes pass edgeR filterByExpr: ", ct, " / ", dpi)
    y <- y_all[keep, , keep.lib.sizes = TRUE]
    y <- edgeR::calcNormFactors(y, method = "TMM")
    need(all(is.finite(y$samples$norm.factors) & y$samples$norm.factors > 0), "Bad TMM factors")
    y <- edgeR::estimateDisp(y, design = design, robust = FALSE)
    fit <- edgeR::glmQLFit(y, design = design, robust = FALSE)
    coef_index <- match("groupTARGET", colnames(design))
    need(!is.na(coef_index), "Missing TARGET coefficient in paired design")
    qlf <- edgeR::glmQLFTest(fit, coef = coef_index)
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    tt$gene <- rownames(tt)
    need(all(c("gene", "logFC", "logCPM", "F", "PValue", "FDR") %in% names(tt)), "Unexpected edgeR topTags columns")
    tt <- tt[match(rownames(y$counts), tt$gene), , drop = FALSE]
    need(!anyNA(tt$gene), "edgeR gene ordering mismatch")

    lib_eff <- y$samples$lib.size * y$samples$norm.factors
    lcpm <- edgeR::cpm(counts, lib.size = lib_eff, log = TRUE, prior.count = PRIOR_COUNT)
    need(all(is.finite(lcpm)), "Nonfinite logCPM result")
    mapped_sig <- intersect(SIG_GENES, rownames(lcpm))
    test_sig <- intersect(SIG_GENES, tt$gene)
    deltas <- vapply(seq_along(stage_gsms), function(i) {
      lcpm[, 2L * i] - lcpm[, 2L * i - 1L]
    }, numeric(nrow(lcpm)))
    rownames(deltas) <- rownames(lcpm)
    colnames(deltas) <- stage_gsms
    need(all(is.finite(deltas)), "Invalid paired logCPM delta")

    for (i in seq_along(stage_gsms)) {
      d_sig <- deltas[mapped_sig, i]
      d_clean <- deltas[setdiff(mapped_sig, IDENTITY_MARKERS), i]
      d_rest <- deltas[setdiff(rownames(deltas), SIG_GENES), i]
      per_gsm_rows[[length(per_gsm_rows) + 1L]] <- data.frame(comparison = ct, dpi = dpi, gsm = stage_gsms[i],
        n_signature_genes_mapped = length(d_sig), mean_signature_target_minus_other_logCPM = mean_or_na(d_sig),
        median_signature_target_minus_other_logCPM = med_or_na(d_sig), fraction_signature_genes_higher_in_target = mean(d_sig >
          0), n_signature_without_predefined_identity_markers = length(d_clean), mean_marker_excluded_signature_delta = mean_or_na(d_clean),
        fraction_marker_excluded_signature_higher = mean(d_clean > 0), mean_nonsignature_delta = mean_or_na(d_rest),
        n_target_cells = cell_sizes$n_cells[cell_sizes$comparison == ct & cell_sizes$gsm == stage_gsms[i] &
          cell_sizes$group == "TARGET"], n_other_cells = cell_sizes$n_cells[cell_sizes$comparison ==
          ct & cell_sizes$gsm == stage_gsms[i] & cell_sizes$group == "OTHER"], stringsAsFactors = FALSE)
    }

    res <- sig_map
    res$comparison <- ct
    res$dpi <- dpi
    res$pass_edgeR_expression_filter <- res$gene %in% tt$gene
    m <- match(res$gene, tt$gene)
    for (field in c("logFC", "logCPM", "F", "PValue", "FDR")) res[[field]] <- tt[[field]][m]
    dmatch <- match(res$gene, rownames(deltas))
    ds <- deltas[dmatch, , drop = FALSE]
    res$GSM1_delta <- ds[, 1]
    res$GSM2_delta <- ds[, 2]
    res$GSM3_delta <- ds[, 3]
    res$mean_paired_logCPM_delta <- rowMeans(ds, na.rm = TRUE)
    res$n_GSM_with_target_higher <- rowSums(ds > 0, na.rm = TRUE)
    res$n_GSM_with_nonmissing_delta <- rowSums(is.finite(ds))
    res$mean_paired_logCPM_delta[res$n_GSM_with_nonmissing_delta == 0] <- NA_real_
    res$GSM1 <- stage_gsms[1]
    res$GSM2 <- stage_gsms[2]
    res$GSM3 <- stage_gsms[3]
    res <- res[, c("comparison", "dpi", "gene", "discovered_avg_log2FC", "mapped_in_raw_features",
      "n_raw_feature_ids", "marker_overlap", "pass_edgeR_expression_filter", "logFC", "logCPM",
      "F", "PValue", "FDR", "GSM1", "GSM2", "GSM3", "GSM1_delta", "GSM2_delta", "GSM3_delta", "mean_paired_logCPM_delta",
      "n_GSM_with_target_higher", "n_GSM_with_nonmissing_delta"), drop = FALSE]
    per_gene_rows[[length(per_gene_rows) + 1L]] <- res
    positive_tested <- sum(tt$gene %in% SIG_GENES & tt$logFC > 0)
    neg_tested <- sum(tt$gene %in% SIG_GENES & tt$logFC < 0)
    sig_test_res <- tt[tt$gene %in% SIG_GENES, , drop = FALSE]
    d <- do.call(rbind, per_gsm_rows[(length(per_gsm_rows) - 2L):length(per_gsm_rows)])
    del_mean <- d$mean_signature_target_minus_other_logCPM
    del_clean <- d$mean_marker_excluded_signature_delta
    sem <- safe_sd(del_mean)/sqrt(length(del_mean))

    tcrit <- stats::qt(0.975, df = length(del_mean) - 1L)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(comparison = ct, dpi = dpi, n_GSM_pairs = length(stage_gsms),
      residual_df_edgeR_QL = residual_df, frozen_signature_size = length(SIG_GENES), n_signature_mapped = length(mapped_sig),
      n_signature_unmapped = length(SIG_GENES) - length(mapped_sig), n_signature_passed_filter = length(test_sig),
      n_signature_tested_logFC_positive = positive_tested, n_signature_tested_logFC_negative = neg_tested,
      fraction_tested_signature_logFC_positive = if (length(test_sig))
        positive_tested/length(test_sig) else NA_real_, n_signature_FDR_lt_005_positive = sum(sig_test_res$FDR < 0.05 & sig_test_res$logFC >
        0), n_signature_FDR_lt_005_negative = sum(sig_test_res$FDR < 0.05 & sig_test_res$logFC <
        0), n_signature_higher_in_all_three_GSM = sum(res$n_GSM_with_target_higher == 3L & res$n_GSM_with_nonmissing_delta ==
        3L), mean_of_three_GSM_signature_deltas = mean(del_mean), median_of_three_GSM_signature_deltas = stats::median(del_mean),
      exploratory_95CI_lower_three_GSM = mean(del_mean) - tcrit * sem, exploratory_95CI_upper_three_GSM = mean(del_mean) +
        tcrit * sem, min_GSM_signature_delta = min(del_mean), max_GSM_signature_delta = max(del_mean),
      n_GSM_with_positive_signature_delta = sum(del_mean > 0), mean_GSM_delta_excluding_identity_markers = mean(del_clean),
      n_GSM_positive_excluding_identity_markers = sum(del_clean > 0), n_genes_after_edgeR_filter = sum(keep),
      stringsAsFactors = FALSE)
    model_num <- model_num + 1L
    model_rows[[model_num]] <- data.frame(comparison = ct, dpi = dpi, n_GSM_pairs = 3L, n_pseudobulk = ncol(counts),
      n_design_parameters = ncol(design), residual_df = residual_df, model = "~ gsm + group (OTHER baseline; TARGET positive)",
      n_genes_total = nrow(counts), n_genes_tested = sum(keep), TMM_effective_library_sizes = paste(signif(lib_eff,
        5), collapse = ";"), TMM_norm_factors = paste(signif(y$samples$norm.factors, 5), collapse = ";"),
      dispersion_common = y$common.dispersion, stringsAsFactors = FALSE)

    camera_out <- tryCatch({
      idx <- match(test_sig, rownames(y$counts))
      idx <- idx[!is.na(idx)]
      need(length(idx) > 20L, "Too few signature genes for camera")
      xnorm <- edgeR::cpm(y, log = TRUE, prior.count = PRIOR_COUNT)
      cam <- limma::camera(xnorm, index = list(Frozen_1320_gene_signature = idx), design = design,
        contrast = as.numeric(colnames(design) == "groupTARGET"))
      data.frame(comparison = ct, dpi = dpi, gene_set = rownames(cam), NGenes = cam$NGenes, Direction = cam$Direction,
        PValue = cam$PValue, FDR = cam$FDR, method = "limma camera; exploratory n_GSM=3", notes = "",
        stringsAsFactors = FALSE)
    }, error = function(e) {
      log_msg("CAMERA not available for ", ct, "/", dpi, ": ", conditionMessage(e))
      data.frame(comparison = ct, dpi = dpi, gene_set = "Frozen_1320_gene_signature", NGenes = NA_integer_,
        Direction = NA_character_, PValue = NA_real_, FDR = NA_real_, method = "limma camera exploratory (failed)",
        notes = conditionMessage(e), stringsAsFactors = FALSE)
    })
    camera_rows[[length(camera_rows) + 1L]] <- camera_out

    saveRDS(list(original_integer_UMI_by_gene_symbol = counts, library_metadata = meta, design_matrix = design,
      filterByExpr = keep, tmm_norm_factors = y$samples$norm.factors, effective_library_sizes = lib_eff,
      residual_df = residual_df, original_data_source = "Step119 original integer UMI; not DecontX",
      frozen_signature_MD5 = unname(tools::md5sum(SIGNATURE_FILE))), file.path(RDS_DIR, sprintf("132_pseudobulk_%s_%ddpi_v01.rds",
      ct, dpi)), compress = "gzip")
    rm(original_counts, counts, y_all, y, fit, qlf, tt, lcpm, deltas, res, ds, meta, camera_out)
    invisible(gc(verbose = FALSE))
  }
}

summary_tab <- do.call(rbind, summary_rows)
per_gsm_tab <- do.call(rbind, per_gsm_rows)
per_gene_tab <- do.call(rbind, per_gene_rows)
model_tab <- do.call(rbind, model_rows)
camera_tab <- do.call(rbind, camera_rows)
need(nrow(summary_tab) == 12L && nrow(per_gsm_tab) == 36L && nrow(per_gene_tab) == 12L * EXPECTED_N_SIGNATURE,
  "Unexpected summary output row counts")
write_csv(model_tab, "05_132_stage_paired_edgeR_design_and_normalization_v01.csv")
write_csv(per_gsm_tab, "06_132_per_GSM_frozen_signature_paired_deltas_v01.csv")
write_csv(summary_tab, "07_132_stage_signature_validation_summary_v01.csv")
write_csv(per_gene_tab, "08_132_all_1320_signature_gene_results_v01.csv")
write_csv(camera_tab, "09_132_exploratory_camera_gene_set_test_v01.csv")

extract_summary <- function(ct) {
  d <- summary_tab[summary_tab$comparison == ct, c("dpi", "mean_of_three_GSM_signature_deltas", "fraction_tested_signature_logFC_positive",
    "n_GSM_with_positive_signature_delta"), drop = FALSE]
  names(d)[-1] <- paste0(ct, "__", names(d)[-1])
  d
}
robustness <- Reduce(function(a, b) merge(a, b, by = "dpi", all = TRUE, sort = TRUE), lapply(CONTRASTS,
  extract_summary))
write_csv(robustness, "10_132_comparison_sensitivity_matrix_v01.csv")

tryCatch({
  f <- file.path(FIG_DIR, "Fig_132_signature_paired_delta_by_GSM_v01.png")
  grDevices::png(f, width = 1850, height = 900, res = 145)
  graphics::par(mfrow = c(1, 3), mar = c(7, 4.8, 3.5, 1), oma = c(0, 0, 2, 0))
  for (ct in CONTRASTS[1:3]) {
    d <- per_gsm_tab[per_gsm_tab$comparison == ct, , drop = FALSE]
    vals <- d$mean_signature_target_minus_other_logCPM
    cols <- ifelse(d$dpi == 14, "#317A9B", "#BD7447")
    yl <- range(c(0, vals), finite = TRUE)
    pad <- max(0.05, diff(yl) * 0.2)
    graphics::barplot(vals, col = cols, border = NA, names.arg = paste0(substr(d$gsm, 8, 10), "\n",
      d$dpi, "d"), las = 1, ylim = c(yl[1] - pad, yl[2] + pad), main = sub("_131$", "", sub("^contrast_",
      "", ct)), ylab = "Mean frozen-signature TARGET - OTHER logCPM")
    graphics::abline(h = 0, lty = 2)
  }
  graphics::mtext("Original UMI pseudobulk | each bar = one GSM pair", outer = TRUE, line = 0.2, cex = 1.2,
    font = 2)
  grDevices::dev.off()
}, error = function(e) {
  try(grDevices::dev.off(), silent = TRUE)
  log_msg("WARNING: figure 1 failed: ", conditionMessage(e))
})
tryCatch({
  f <- file.path(FIG_DIR, "Fig_132_signature_direction_14_vs_28_v01.png")
  grDevices::png(f, width = 1500, height = 850, res = 145)
  graphics::par(mfrow = c(1, 2), mar = c(9, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  for (dpi in c(14L, 28L)) {
    d <- summary_tab[summary_tab$dpi == dpi, ]
    bp <- graphics::barplot(100 * d$fraction_tested_signature_logFC_positive, col = "#4D8F95", border = NA,
      ylim = c(0, 100), las = 2, names.arg = c("A", "B", "C", "A broad", "B broad", "A stable"),
      ylab = "Frozen signature genes with positive edgeR logFC (%)", main = paste0(dpi, " dpi"))
    graphics::abline(h = 50, lty = 2, col = "gray50")
    graphics::text(bp, 100 * d$fraction_tested_signature_logFC_positive + 3, labels = sprintf("%.0f%%",
      100 * d$fraction_tested_signature_logFC_positive), cex = 0.75)
  }
  graphics::mtext("Positive = Clec7a-high TARGET enriched | 3 GSM pairs/stage", outer = TRUE, line = 0.2,
    cex = 1.05)
  grDevices::dev.off()
}, error = function(e) {
  try(grDevices::dev.off(), silent = TRUE)
  log_msg("WARNING: figure 2 failed: ", conditionMessage(e))
})

write_txt("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    "11_132_REPORT_v01.txt")
write_txt(capture.output(sessionInfo()), "12_132_sessionInfo_v01.txt")

created <- list.files(OUT, full.names = FALSE, recursive = TRUE)
write_csv(data.frame(file = created, stringsAsFactors = FALSE), "13_132_OUTPUT_INDEX_v01.csv")

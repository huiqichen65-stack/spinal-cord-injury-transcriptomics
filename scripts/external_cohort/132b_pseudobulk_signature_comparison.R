source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

ROOT <- PROJECT_ROOT
OUT <- file.path(ROOT, "results", "132_GSE281205_original_UMI_pseudobulk_signature_validation")
RDS_DIR <- file.path(OUT, "RDS")
FIG_DIR <- file.path(OUT, "figures")
LOG_DIR <- file.path(OUT, "logs")
SIG_FILE <- file.path(ROOT, "results", "21_Clec7a_DEG_volcano_myeloid", "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv")
FEATURE_MAP_FILE <- file.path(ROOT, "results", "124_GSE281205_lineage_review_myeloid_extraction", "12_124_gene_feature_mapping_v01.csv")
MANIFEST_FILE <- file.path(ROOT, "results", "131_GSE281205_myeloid_subtype_annotation_and_microglia_groups",
  "13_131_READY_FOR_STEP132_original_UMI_cell_manifest_v01.csv")
OLD_SUMMARY_FILE <- file.path(OUT, "07_132_stage_signature_validation_summary_v01.csv")
OLD_MAPPING_FILE <- file.path(OUT, "04_132_frozen_signature_gene_mapping_audit_v01.csv")
OLD_CAMERA_FILE <- file.path(OUT, "09_132_exploratory_camera_gene_set_test_v01.csv")
GSMS <- paste0("GSM", 8614662:8614667)
EXPECTED_DPI <- setNames(rep(c(14L, 28L), each = 3L), GSMS)
EXPECTED_GENES <- 1320L
PRIOR_COUNT <- 2
CONTRASTS <- c("contrast_A_target_vs_other_microglia_131", "contrast_B_target_vs_all_other_highconf_myeloid_131",
  "contrast_C_target_vs_nonmicroglial_myeloid_131", "sensitivity_A_broad_microglia_lt50", "sensitivity_B_broad_all_myeloid_lt50",
  "sensitivity_A_partition_consistent")
IDENTITY_MARKERS <- c("Clec7a", "Spp1", "Lpl", "Gpnmb", "Cst7", "Apoe", "Trem2", "Tyrobp", "P2ry12",
  "Tmem119", "Fcrls", "Hexb", "Cx3cr1", "Lyz2", "Ms4a7", "Mrc1", "Cd74", "H2-Aa", "Flt3", "Cd209a",
  "Pclaf", "Birc5", "Mki67", "Top2a")

need <- function(ok, ...) {
  if (!isTRUE(ok))
    stop(paste(..., collapse = ""), call. = FALSE)
  invisible(TRUE)
}
need(file.exists(SIG_FILE), "Step21 gene list missing: ", SIG_FILE)
need(file.exists(FEATURE_MAP_FILE), "Step124 gene ID mapping missing: ", FEATURE_MAP_FILE)
need(file.exists(MANIFEST_FILE), "Step131 cell manifest missing: ", MANIFEST_FILE)
for (p in c("edgeR", "limma", "AnnotationDbi", "org.Mm.eg.db")) {
  need(requireNamespace(p, quietly = TRUE), "Required package not installed: ", p, "\nInstall with BiocManager::install(c(\"edgeR\",\"limma\",\"AnnotationDbi\",\"org.Mm.eg.db\"))")
}
for (d in c(OUT, RDS_DIR, FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(LOG_DIR, "00_132_run_log_v02.log")
cat("", file = LOG)
log_msg <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n", sep = "")
  cat(s, "\n", sep = "", file = LOG, append = TRUE)
  invisible(s)
}
write_csv <- function(d, name) {
  need(is.data.frame(d), "Internal error: CSV is not data.frame: ", name)
  path <- file.path(OUT, name)
  utils::write.csv(d, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}
write_txt <- function(lines, name) {
  writeLines(as.character(lines), file.path(OUT, name), useBytes = TRUE)
}
read_csv <- function(path) utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
mean_ok <- function(x) {
  x <- x[is.finite(x)]
  if (length(x))
    mean(x) else NA_real_
}
median_ok <- function(x) {
  x <- x[is.finite(x)]
  if (length(x))
    stats::median(x) else NA_real_
}
md5 <- function(p) unname(tools::md5sum(p))

for (p in c("edgeR", "limma", "AnnotationDbi", "org.Mm.eg.db")) log_msg("PACKAGE ", p, ": ", as.character(utils::packageVersion(p)))

sig <- read_csv(SIG_FILE)
need(all(c("gene", "avg_log2FC") %in% names(sig)), "Step21 frozen signature requires gene, avg_log2FC columns")
sig$gene <- trimws(as.character(sig$gene))
sig$avg_log2FC <- suppressWarnings(as.numeric(sig$avg_log2FC))
need(nrow(sig) == EXPECTED_GENES, "Frozen signature is not 1,320 rows: ", nrow(sig))
need(!anyNA(sig$gene) && !any(!nzchar(sig$gene)) && !anyDuplicated(sig$gene), "Frozen signature symbols missing or duplicated")
need(all(is.finite(sig$avg_log2FC) & sig$avg_log2FC > 0), "Step21 signature not a valid positive-logFC 1,320-gene signature")
SIG_GENES <- sig$gene
SIG_MD5 <- md5(SIG_FILE)
MANIFEST_MD5 <- md5(MANIFEST_FILE)

manifest <- read_csv(MANIFEST_FILE)
need(nrow(manifest) == 28716L && all(c("cell_id", "gsm", "dpi", CONTRASTS) %in% names(manifest)), "Step131 manifest invalid or different from v01 contract")
need(!anyDuplicated(manifest$cell_id) && !anyNA(manifest$cell_id), "Step131 duplicated/missing cell IDs")
need(identical(sort(unique(as.character(manifest$gsm))), GSMS), "Unexpected GSM set")
need(all(as.integer(manifest$dpi) == EXPECTED_DPI[manifest$gsm]), "Manifest GSM and stage conflict")
for (ct in CONTRASTS) {
  need(all(as.character(manifest[[ct]]) %in% c("TARGET", "OTHER", "NOT_IN_CONTRAST")), "Invalid groups for ",
    ct)
}

rds_path <- function(ct, dpi) file.path(RDS_DIR, sprintf("132_pseudobulk_%s_%ddpi_v01.rds", ct, dpi))
all_rds_paths <- unlist(lapply(CONTRASTS, function(ct) vapply(c(14L, 28L), function(dpi) rds_path(ct,
  dpi), character(1))))
missing_rds <- all_rds_paths[!file.exists(all_rds_paths)]
need(!length(missing_rds), "Step132 v01 original-integer-UMI pseudobulk RDS missing:\n", paste(missing_rds,
  collapse = "\n"), "\nRun the original Step132 v01 once to create these; v02 NEVER requires ZIP inputs.")

first <- readRDS(all_rds_paths[1])
need(is.list(first) && all(c("original_integer_UMI_by_gene_symbol", "library_metadata", "design_matrix") %in%
  names(first)), "Step132 v01 pseudobulk RDS has an unrecognized structure")
raw_symbols <- rownames(first$original_integer_UMI_by_gene_symbol)
need(length(raw_symbols) > 10000L && !anyDuplicated(raw_symbols), "Invalid v01 raw feature-symbol rows")
if (!is.null(first$frozen_signature_MD5)) need(identical(as.character(first$frozen_signature_MD5), SIG_MD5),
  "Step21 frozen signature file changed since Step132 v01. Stop.")
rm(first)

feat <- read_csv(FEATURE_MAP_FILE)
need(all(c("gene_id", "gene_symbol") %in% names(feat)), "Step124 feature mapping missing gene_id or gene_symbol")
feat$gene_id <- as.character(feat$gene_id)
feat$gene_symbol <- trimws(as.character(feat$gene_symbol))
feat$stable_ensembl <- sub("\\.[0-9]+$", "", feat$gene_id)
feat <- feat[!is.na(feat$gene_symbol) & nzchar(feat$gene_symbol) & feat$gene_symbol %in% raw_symbols,
  , drop = FALSE]
need(nrow(feat) > 10000L, "Too few Step124 mappings overlap Step132 raw features")

raw_symbols <- unique(as.character(raw_symbols))
stopifnot(!anyDuplicated(raw_symbols))

org_db <- org.Mm.eg.db::org.Mm.eg.db
need(all(c("ALIAS", "SYMBOL", "ENSEMBL") %in% AnnotationDbi::keytypes(org_db)), "Unexpected org.Mm.eg.db key types")
annotation_version <- as.character(utils::packageVersion("org.Mm.eg.db"))

mapped <- data.frame(gene = SIG_GENES, discovered_avg_log2FC = sig$avg_log2FC, mapped_in_raw_features_v01 = SIG_GENES %in%
  raw_symbols, resolved_symbol = rep(NA_character_, length(SIG_GENES)), mapping_status = rep("unresolved_no_official_match",
  length(SIG_GENES)), mapping_method = rep(NA_character_, length(SIG_GENES)), n_candidate_raw_symbols = integer(length(SIG_GENES)),
  candidate_raw_symbols = rep("", length(SIG_GENES)), annotation_source = rep(paste0("org.Mm.eg.db ",
    annotation_version), length(SIG_GENES)), stringsAsFactors = FALSE)
exact <- mapped$mapped_in_raw_features_v01
mapped$resolved_symbol[exact] <- SIG_GENES[exact]
mapped$mapping_status[exact] <- "exact_symbol"
mapped$mapping_method[exact] <- "original_symbol_exact"
mapped$n_candidate_raw_symbols[exact] <- 1L
mapped$candidate_raw_symbols[exact] <- SIG_GENES[exact]

not_direct <- SIG_GENES[!exact]
valid_alias_keys <- intersect(not_direct, AnnotationDbi::keys(org_db, keytype = "ALIAS"))
raw_alias <- data.frame(ALIAS = character(), SYMBOL = character(), ENSEMBL = character(), stringsAsFactors = FALSE)
if (length(valid_alias_keys)) {
  raw_alias <- suppressMessages(AnnotationDbi::select(org_db, keys = valid_alias_keys, keytype = "ALIAS",
    columns = c("SYMBOL", "ENSEMBL")))
  raw_alias <- as.data.frame(raw_alias, stringsAsFactors = FALSE)
  need(all(c("ALIAS", "SYMBOL", "ENSEMBL") %in% names(raw_alias)), "Official alias selection did not return required fields")
  raw_alias$ALIAS <- as.character(raw_alias$ALIAS)
  raw_alias$SYMBOL <- as.character(raw_alias$SYMBOL)
  raw_alias$ENSEMBL <- sub("\\.[0-9]+$", "", as.character(raw_alias$ENSEMBL))
  raw_alias <- unique(raw_alias[, c("ALIAS", "SYMBOL", "ENSEMBL"), drop = FALSE])
}

ens2sym <- split(as.character(feat$gene_symbol), feat$stable_ensembl)
ens2sym <- lapply(ens2sym, unique)
alias_details <- list()
for (i in which(!exact)) {
  old <- SIG_GENES[i]
  anno <- raw_alias[!is.na(raw_alias$ALIAS) & raw_alias$ALIAS == old, , drop = FALSE]
  candidate <- character()
  direct_symbol_hits <- character()
  ensemble_hits <- character()
  if (nrow(anno)) {
    direct_symbol_hits <- unique(anno$SYMBOL[!is.na(anno$SYMBOL) & anno$SYMBOL %in% raw_symbols])
    ids <- unique(anno$ENSEMBL[!is.na(anno$ENSEMBL) & nzchar(anno$ENSEMBL)])
    ensemble_hits <- unique(unlist(ens2sym[ids], use.names = FALSE))
    candidate <- sort(unique(c(direct_symbol_hits, ensemble_hits)))
  }
  mapped$n_candidate_raw_symbols[i] <- length(candidate)
  mapped$candidate_raw_symbols[i] <- paste(candidate, collapse = ";")
  if (length(candidate) == 1L) {
    mapped$resolved_symbol[i] <- candidate
    mapped$mapping_status[i] <- "official_alias_unique"
    if (candidate %in% direct_symbol_hits && candidate %in% ensemble_hits) {
      mapped$mapping_method[i] <- "org.Mm.eg.db:ALIAS_to_SYMBOL_and_ENSEMBL"
    } else if (candidate %in% direct_symbol_hits) {
      mapped$mapping_method[i] <- "org.Mm.eg.db:ALIAS_to_SYMBOL"
    } else {
      mapped$mapping_method[i] <- "org.Mm.eg.db:ALIAS_to_ENSEMBL"
    }
  } else if (length(candidate) > 1L) {
    mapped$mapping_status[i] <- "ambiguous_official_alias"
  }
  alias_details[[length(alias_details) + 1L]] <- data.frame(frozen_gene = old, is_valid_org_alias = old %in%
    valid_alias_keys, n_raw_matching_candidates = length(candidate), raw_candidate_symbols = paste(candidate,
    collapse = ";"), annotation_SYMBOL_matches = paste(sort(direct_symbol_hits), collapse = ";"),
    annotation_ENSEMBL_matches = paste(sort(ensemble_hits), collapse = ";"), candidate_annotation_symbols = paste(sort(unique(anno$SYMBOL[!is.na(anno$SYMBOL)])),
      collapse = ";"), candidate_annotation_ensembl = paste(sort(unique(anno$ENSEMBL[!is.na(anno$ENSEMBL)])),
      collapse = ";"), decision = mapped$mapping_status[i], stringsAsFactors = FALSE)
}
alias_audit <- if (length(alias_details)) do.call(rbind, alias_details) else data.frame(frozen_gene = character(),
  stringsAsFactors = FALSE)

known <- which(!is.na(mapped$resolved_symbol))
duplicated_targets <- unique(mapped$resolved_symbol[known][duplicated(mapped$resolved_symbol[known])])
for (sym in duplicated_targets) {
  ii <- which(!is.na(mapped$resolved_symbol) & mapped$resolved_symbol == sym)
  direct_ii <- ii[mapped$mapping_status[ii] == "exact_symbol"]
  reject <- if (length(direct_ii) == 1L)
    setdiff(ii, direct_ii) else ii
  mapped$mapping_status[reject] <- "collision_with_other_frozen_gene"
  mapped$mapping_method[reject] <- NA_character_
  mapped$resolved_symbol[reject] <- NA_character_
}
need(!anyDuplicated(stats::na.omit(mapped$resolved_symbol)), "Two symbols still resolve to the same raw feature")

mapped$marker_overlap <- mapped$gene %in% IDENTITY_MARKERS | (!is.na(mapped$resolved_symbol) & mapped$resolved_symbol %in%
  IDENTITY_MARKERS)
mapped$mapped_in_raw_features <- !is.na(mapped$resolved_symbol)
mapped$n_raw_feature_ids <- vapply(mapped$resolved_symbol, function(g) {
  if (is.na(g))
    return(0L)
  sum(feat$gene_symbol == g)
}, integer(1))
mapped$unresolved_reason <- ifelse(mapped$mapped_in_raw_features, "", mapped$mapping_status)
need(sum(mapped$mapped_in_raw_features) >= sum(exact), "Reconciliation unexpectedly reduced exact-match coverage")

write_csv(data.frame(input_role = c("Step21_frozen_signature", "Step131_cell_manifest", "Step124_feature_mapping",
  "Step132_v01_pseudobulk_RDS"), path = c(SIG_FILE, MANIFEST_FILE, FEATURE_MAP_FILE, RDS_DIR), md5 = c(SIG_MD5,
  MANIFEST_MD5, md5(FEATURE_MAP_FILE), NA_character_)), "01_132_input_provenance_v02.csv")
write_csv(data.frame(gene = SIG_GENES, discovered_avg_log2FC = sig$avg_log2FC, marker_overlap = mapped$marker_overlap),
  "02_132_frozen_signature_metadata_v02.csv")
write_csv(mapped, "04_132_frozen_signature_gene_mapping_audit_v02.csv")
write_csv(alias_audit, "14_132_official_alias_candidate_audit_v02.csv")

coverage_cutoffs <- c(10L, 50L, 100L, 250L, 500L, EXPECTED_GENES)
coverage <- do.call(rbind, lapply(coverage_cutoffs, function(k) {
  ii <- order(-sig$avg_log2FC)[seq_len(k)]
  data.frame(top_ranked_by_discovery_log2FC = k, n_exact_v01 = sum(exact[ii]), n_official_alias_recovered = sum(mapped$mapping_status[ii] ==
    "official_alias_unique"), n_mapped_v02 = sum(mapped$mapped_in_raw_features[ii]), n_unresolved_v02 = sum(!mapped$mapped_in_raw_features[ii]),
    pct_mapped_v02 = 100 * mean(mapped$mapped_in_raw_features[ii]), stringsAsFactors = FALSE)
}))
write_csv(coverage, "15_132_top_signature_gene_mapping_before_after_v02.csv")

summary_rows <- list()
per_gsm_rows <- list()
per_gene_rows <- list()
model_rows <- list()
camera_rows <- list()
qc_rows <- list()

for (ct in CONTRASTS) for (dpi in c(14L, 28L)) {
  f <- rds_path(ct, dpi)
  obj <- readRDS(f)
  need(is.list(obj) && all(c("original_integer_UMI_by_gene_symbol", "library_metadata", "design_matrix") %in%
    names(obj)), "Malformed RDS: ", f)
  if (!is.null(obj$frozen_signature_MD5))
    need(identical(as.character(obj$frozen_signature_MD5), SIG_MD5), "Frozen signature MD5 mismatch: ",
      f)
  counts <- as.matrix(obj$original_integer_UMI_by_gene_symbol)
  need(identical(rownames(counts), raw_symbols), "Gene order differs across original-UMI pseudobulk RDS: ",
    f)
  need(ncol(counts) == 6L && all(is.finite(counts)) && all(counts >= 0) && all(abs(counts - round(counts)) <
    1e-08), "Expected six original INTEGER UMI pseudobulks: ", f)
  stage_gsms <- GSMS[EXPECTED_DPI[GSMS] == dpi]
  wanted <- as.vector(rbind(paste0(stage_gsms, "__OTHER"), paste0(stage_gsms, "__TARGET")))
  need(identical(colnames(counts), wanted), "Pseudobulk GSM/group ordering differs from input: ", f)
  meta <- as.data.frame(obj$library_metadata)
  need(identical(as.character(meta$library_id), wanted), "Original pseudobulk library metadata mismatch: ",
    f)
  meta$gsm <- factor(as.character(meta$gsm), levels = stage_gsms)
  meta$group <- factor(as.character(meta$group), levels = c("OTHER", "TARGET"))
  need(identical(as.character(meta$group), rep(c("OTHER", "TARGET"), 3L)), "Target/control order mismatch: ",
    f)
  design <- stats::model.matrix(~gsm + group, data = meta)
  need(qr(design)$rank == ncol(design) && ncol(counts) - ncol(design) == 2L, "Paired model rank or residual df invalid: ",
    f)
  need(identical(unname(obj$design_matrix), unname(design)), "New paired design differs from Step132 v01: ",
    f)
  for (g in stage_gsms) for (grp in c("OTHER", "TARGET")) {
    n_from_manifest <- sum(manifest$gsm == g & manifest[[ct]] == grp)
    need(n_from_manifest >= 20L, "Unexpected small group ", ct, "/", g, "/", grp)
    qc_rows[[length(qc_rows) + 1L]] <- data.frame(comparison = ct, dpi = dpi, gsm = g, group = grp,
      n_cells_from_frozen_manifest = n_from_manifest, original_integer_UMI_library_sum = sum(counts[,
        paste0(g, "__", grp)]), n_genes_with_UMI = sum(counts[, paste0(g, "__", grp)] > 0), input_RDS = basename(f),
      stringsAsFactors = FALSE)
  }
  y0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(y0, design = design)
  need(sum(keep) > 500L, "Too few passing genes: ", f)
  y <- edgeR::calcNormFactors(y0[keep, , keep.lib.sizes = TRUE], method = "TMM")
  y <- edgeR::estimateDisp(y, design = design, robust = FALSE)
  fit <- edgeR::glmQLFit(y, design = design, robust = FALSE)
  coef_idx <- match("groupTARGET", colnames(design))
  need(!is.na(coef_idx), "groupTARGET absent in design")
  test <- edgeR::glmQLFTest(fit, coef = coef_idx)
  tt <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  tt$resolved_symbol <- rownames(tt)
  tt <- tt[match(rownames(y$counts), tt$resolved_symbol), , drop = FALSE]
  need(!anyNA(tt$resolved_symbol), "edgeR result gene ordering failed: ", f)
  lib_eff <- y$samples$lib.size * y$samples$norm.factors
  lcpm <- edgeR::cpm(counts, lib.size = lib_eff, log = TRUE, prior.count = PRIOR_COUNT)
  need(all(is.finite(lcpm)), "logCPM failed: ", f)
  deltas <- sapply(seq_along(stage_gsms), function(i) {
    lcpm[, 2L * i] - lcpm[, 2L * i - 1L]
  })
  rownames(deltas) <- rownames(lcpm)
  colnames(deltas) <- stage_gsms
  res <- mapped
  res$comparison <- ct
  res$dpi <- dpi
  res$pass_edgeR_expression_filter <- res$mapped_in_raw_features & res$resolved_symbol %in% tt$resolved_symbol
  m <- match(res$resolved_symbol, tt$resolved_symbol)
  for (nm in c("logFC", "logCPM", "F", "PValue", "FDR")) res[[nm]] <- tt[[nm]][m]
  m_all <- match(res$resolved_symbol, rownames(deltas))
  ds <- deltas[m_all, , drop = FALSE]
  res$GSM1 <- stage_gsms[1]
  res$GSM2 <- stage_gsms[2]
  res$GSM3 <- stage_gsms[3]
  res$GSM1_delta <- ds[, 1]
  res$GSM2_delta <- ds[, 2]
  res$GSM3_delta <- ds[, 3]
  res$n_GSM_with_nonmissing_delta <- rowSums(is.finite(ds))
  res$n_GSM_with_target_higher <- rowSums(ds > 0, na.rm = TRUE)
  res$mean_paired_logCPM_delta <- rowMeans(ds, na.rm = TRUE)
  res$mean_paired_logCPM_delta[res$n_GSM_with_nonmissing_delta == 0] <- NA_real_
  per_gene_rows[[length(per_gene_rows) + 1L]] <- res

  all_resolved <- mapped$resolved_symbol[mapped$mapped_in_raw_features]
  all_clean <- mapped$resolved_symbol[mapped$mapped_in_raw_features & !mapped$marker_overlap]
  other_genes <- setdiff(rownames(deltas), all_resolved)
  stage_gsm_means <- numeric(length(stage_gsms))
  stage_clean_means <- numeric(length(stage_gsms))
  for (j in seq_along(stage_gsms)) {
    x <- deltas[all_resolved, j]
    xx <- deltas[all_clean, j]
    stage_gsm_means[j] <- mean_ok(x)
    stage_clean_means[j] <- mean_ok(xx)
    per_gsm_rows[[length(per_gsm_rows) + 1L]] <- data.frame(comparison = ct, dpi = dpi, gsm = stage_gsms[j],
      n_signature_genes_mapped = length(x), mean_signature_target_minus_other_logCPM = mean_ok(x),
      median_signature_target_minus_other_logCPM = median_ok(x), fraction_signature_genes_higher_in_target = mean(x >
        0), n_signature_without_predefined_identity_markers = length(xx), mean_marker_excluded_signature_delta = mean_ok(xx),
      fraction_marker_excluded_signature_higher = mean(xx > 0), mean_nonsignature_delta = mean_ok(deltas[other_genes,
        j]), stringsAsFactors = FALSE)
  }
  test_symbols <- tt$resolved_symbol
  tested <- mapped$resolved_symbol[mapped$mapped_in_raw_features & mapped$resolved_symbol %in% test_symbols]
  matched_tt <- tt[match(tested, tt$resolved_symbol), , drop = FALSE]
  sem <- stats::sd(stage_gsm_means)/sqrt(3)
  ci_half <- stats::qt(0.975, df = 2) * sem
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(comparison = ct, dpi = dpi, n_GSM_pairs = 3L,
    residual_df_edgeR_QL = 2L, frozen_signature_size = EXPECTED_GENES, n_signature_exact_v01 = sum(exact),
    n_signature_official_alias_recovered = sum(mapped$mapping_status == "official_alias_unique"),
    n_signature_mapped = length(all_resolved), n_signature_unmapped = EXPECTED_GENES - length(all_resolved),
    n_signature_passed_filter = length(tested), n_signature_tested_logFC_positive = sum(matched_tt$logFC >
      0), n_signature_tested_logFC_negative = sum(matched_tt$logFC < 0), fraction_tested_signature_logFC_positive = mean(matched_tt$logFC >
      0), n_signature_FDR_lt_005_positive = sum(matched_tt$FDR < 0.05 & matched_tt$logFC > 0),
    n_signature_FDR_lt_005_negative = sum(matched_tt$FDR < 0.05 & matched_tt$logFC < 0), n_signature_higher_in_all_three_GSM = sum(res$n_GSM_with_target_higher ==
      3 & res$n_GSM_with_nonmissing_delta == 3), mean_of_three_GSM_signature_deltas = mean(stage_gsm_means),
    median_of_three_GSM_signature_deltas = stats::median(stage_gsm_means), exploratory_95CI_lower_three_GSM = mean(stage_gsm_means) -
      ci_half, exploratory_95CI_upper_three_GSM = mean(stage_gsm_means) + ci_half, min_GSM_signature_delta = min(stage_gsm_means),
    max_GSM_signature_delta = max(stage_gsm_means), n_GSM_with_positive_signature_delta = sum(stage_gsm_means >
      0), mean_GSM_delta_excluding_identity_markers = mean(stage_clean_means), n_GSM_positive_excluding_identity_markers = sum(stage_clean_means >
      0), n_genes_after_edgeR_filter = sum(keep), stringsAsFactors = FALSE)
  model_rows[[length(model_rows) + 1L]] <- data.frame(comparison = ct, dpi = dpi, n_GSM_pairs = 3L,
    n_pseudobulk = ncol(counts), n_design_parameters = ncol(design), residual_df = 2L, model = "~ gsm + group; paired GSM (OTHER baseline, TARGET positive)",
    n_genes_total = nrow(counts), n_genes_tested = sum(keep), TMM_effective_library_sizes = paste(signif(lib_eff,
      5), collapse = ";"), TMM_norm_factors = paste(signif(y$samples$norm.factors, 5), collapse = ";"),
    dispersion_common = y$common.dispersion, stringsAsFactors = FALSE)

  lcpm_tested <- edgeR::cpm(y, log = TRUE, prior.count = PRIOR_COUNT)
  camera_sets <- list(frozen_full_1320 = all_resolved, frozen_without_identity_markers = all_clean)
  for (set_name in names(camera_sets)) {
    idx <- match(intersect(camera_sets[[set_name]], rownames(lcpm_tested)), rownames(lcpm_tested))
    idx <- idx[!is.na(idx)]
    cam_one <- tryCatch({
      need(length(idx) >= 20L, "Fewer than 20 tested genes: CAMERA skipped")
      cam <- limma::camera(lcpm_tested, index = setNames(list(idx), set_name), design = design,
        contrast = as.numeric(colnames(design) == "groupTARGET"))
      need(is.data.frame(cam) && nrow(cam) == 1L, "Unexpected single-set CAMERA result shape")
      need(all(c("NGenes", "Direction", "PValue") %in% names(cam)), "CAMERA missing mandatory fields: ",
        paste(names(cam), collapse = ","))
      p <- as.numeric(cam[["PValue"]][1])
      need(length(p) == 1L && is.finite(p) && p >= 0 && p <= 1, "CAMERA returned non-finite PValue")
      data.frame(comparison = ct, dpi = dpi, gene_set = set_name, n_input_resolved_members = length(camera_sets[[set_name]]),
        NGenes = as.integer(cam[["NGenes"]][1]), Direction = as.character(cam[["Direction"]][1]),
        PValue = p, status = "success", error_message = "", stringsAsFactors = FALSE)
    }, error = function(e) data.frame(comparison = ct, dpi = dpi, gene_set = set_name, n_input_resolved_members = length(camera_sets[[set_name]]),
      NGenes = NA_integer_, Direction = NA_character_, PValue = NA_real_, status = "failed", error_message = conditionMessage(e),
      stringsAsFactors = FALSE))
    camera_rows[[length(camera_rows) + 1L]] <- cam_one
  }
  rm(obj, counts, y0, y, fit, test, tt, lcpm, deltas, ds, res, lcpm_tested, cam_one, matched_tt)
  invisible(gc(verbose = FALSE))
}

summary_tab <- do.call(rbind, summary_rows)
per_gsm_tab <- do.call(rbind, per_gsm_rows)
per_gene_tab <- do.call(rbind, per_gene_rows)
model_tab <- do.call(rbind, model_rows)
camera_tab <- do.call(rbind, camera_rows)
qc_tab <- do.call(rbind, qc_rows)
need(nrow(summary_tab) == 12L && nrow(per_gsm_tab) == 36L && nrow(per_gene_tab) == 12L * EXPECTED_GENES &&
  nrow(camera_tab) == 24L, "Final analysis dimensions do not match the intended 12 comparisons")

camera_tab$FDR_BH_across_12_stage_comparisons <- NA_real_
for (gene_set in unique(camera_tab$gene_set)) {
  ii <- which(camera_tab$gene_set == gene_set & is.finite(camera_tab$PValue))
  if (length(ii))
    camera_tab$FDR_BH_across_12_stage_comparisons[ii] <- stats::p.adjust(camera_tab$PValue[ii], method = "BH",
      n = 12L)
}
camera_tab$method <- "limma camera(log2-TMM CPM); GSM-paired ~gsm+group"
camera_tab$interpretation <- "exploratory competitive gene-set test; 3 GSM pairs/stage"

write_csv(qc_tab, "03_132_pseudobulk_library_cell_and_UMI_QC_v02.csv")
write_csv(model_tab, "05_132_stage_paired_edgeR_design_and_normalization_v02.csv")
write_csv(per_gsm_tab, "06_132_per_GSM_frozen_signature_paired_deltas_v02.csv")
write_csv(summary_tab, "07_132_stage_signature_validation_summary_v02.csv")
write_csv(per_gene_tab, "08_132_all_1320_signature_gene_results_v02.csv")
write_csv(camera_tab, "09_132_exploratory_camera_gene_set_test_v02.csv")

extract_summary <- function(ct) {
  d <- summary_tab[summary_tab$comparison == ct, c("dpi", "mean_of_three_GSM_signature_deltas", "fraction_tested_signature_logFC_positive",
    "n_GSM_with_positive_signature_delta"), drop = FALSE]
  names(d)[-1] <- paste0(ct, "__", names(d)[-1])
  d
}
sensitivity <- Reduce(function(a, b) merge(a, b, by = "dpi", all = TRUE, sort = TRUE), lapply(CONTRASTS,
  extract_summary))
write_csv(sensitivity, "10_132_comparison_sensitivity_matrix_v02.csv")

if (file.exists(OLD_SUMMARY_FILE)) {
  old <- read_csv(OLD_SUMMARY_FILE)
  cols <- c("comparison", "dpi", "n_signature_mapped", "n_signature_passed_filter", "fraction_tested_signature_logFC_positive",
    "mean_of_three_GSM_signature_deltas", "n_signature_FDR_lt_005_positive", "n_signature_FDR_lt_005_negative")
  need(all(cols %in% names(old)), "Step132 v01 summary schema changed")
  old <- old[, cols, drop = FALSE]
  new <- summary_tab[, cols, drop = FALSE]
  names(old)[-(1:2)] <- paste0(names(old)[-(1:2)], "_v01")
  names(new)[-(1:2)] <- paste0(names(new)[-(1:2)], "_v02")
  before_after <- merge(old, new, by = c("comparison", "dpi"), sort = FALSE)
  for (nm in cols[-(1:2)]) before_after[[paste0(nm, "_change")]] <- before_after[[paste0(nm, "_v02")]] -
    before_after[[paste0(nm, "_v01")]]
} else {
  log_msg("NOTE Step132 v01 summary CSV missing, before/after summary uses v02 only")
  before_after <- summary_tab
}
write_csv(before_after, "16_132_stage_signature_before_after_v02.csv")

if (file.exists(OLD_MAPPING_FILE)) {
  om <- read_csv(OLD_MAPPING_FILE)
  need(nrow(om) == EXPECTED_GENES && identical(as.character(om$gene), SIG_GENES), "v01 mapping gene order inconsistent with frozen signature")
  need(identical(as.logical(om$mapped_in_raw_features), as.logical(exact)), "v01 mapping disagrees with raw-symbol exact match; check provenance")
}
if (file.exists(OLD_CAMERA_FILE)) {
  oc <- read_csv(OLD_CAMERA_FILE)
  log_msg("Previous Step132 v01 CAMERA failed rows: ", sum(is.na(oc$PValue)), "/", nrow(oc))
}

tryCatch({
  path <- file.path(FIG_DIR, "Fig_132_v02_signature_mean_delta_by_GSM.png")
  grDevices::png(path, width = 1600, height = 850, res = 135)
  graphics::par(mfrow = c(1, 3), mar = c(6.5, 4.5, 3, 0.6), oma = c(0, 0, 2, 0))
  for (ct in CONTRASTS[1:3]) {
    d <- per_gsm_tab[per_gsm_tab$comparison == ct, , drop = FALSE]
    z <- d$mean_signature_target_minus_other_logCPM
    ylims <- range(c(0, z))
    padding <- max(0.04, diff(ylims) * 0.2)
    graphics::barplot(z, names.arg = paste0(substr(d$gsm, 8, 10), "\n", d$dpi, "d"), col = ifelse(d$dpi ==
      14, "#427A9F", "#D1904E"), border = NA, ylim = c(ylims[1] - padding, ylims[2] + padding),
      las = 2, ylab = "Mean TARGET - OTHER logCPM", main = sub("^contrast_", "", ct))
    graphics::abline(h = 0, lty = 2)
  }
  graphics::mtext("132 v02 | official gene alias reconciliation | each bar: one GSM", outer = TRUE,
    line = 0.2, cex = 1.15, font = 2)
  grDevices::dev.off()
}, error = function(e) {
  try(grDevices::dev.off(), silent = TRUE)
  log_msg("WARNING plot failed: ", conditionMessage(e))
})

tryCatch({
  path <- file.path(FIG_DIR, "Fig_132_v02_camera_competitive_results.png")
  grDevices::png(path, width = 1600, height = 1050, res = 135)
  graphics::par(mfrow = c(2, 1), mar = c(6, 9, 3, 1))
  for (gs in unique(camera_tab$gene_set)) {
    d <- camera_tab[camera_tab$gene_set == gs, , drop = FALSE]
    d <- d[match(as.vector(rbind(paste(CONTRASTS, "14", sep = "__"), paste(CONTRASTS, "28", sep = "__"))),
      paste(d$comparison, d$dpi, sep = "__")), , drop = FALSE]
    z <- -log10(pmax(d$PValue, 1e-300))
    z[!is.finite(z)] <- 0
    lab <- paste0(ifelse(d$dpi == 14, "14d", "28d"), " ", sub("^contrast_", "", sub("_131$", "",
      d$comparison)))
    bar_col <- ifelse(d$Direction == "Up", "#3E8993", "#BD7D5A")
    bar_col[is.na(bar_col)] <- "gray65"
    graphics::barplot(z, horiz = TRUE, las = 1, names.arg = lab, col = bar_col, border = NA, cex.names = 0.67,
      main = gs, xlab = "-log10 CAMERA nominal P")
    graphics::abline(v = -log10(0.05), lty = 2, col = "gray45")
  }
  grDevices::dev.off()
}, error = function(e) {
  try(grDevices::dev.off(), silent = TRUE)
  log_msg("WARNING plot failed: ", conditionMessage(e))
})

write_txt("Numerical results are provided in the accompanying output tables.", 
    "11_132_REPORT_v02.txt")
write_txt(capture.output(sessionInfo()), "12_132_sessionInfo_v02.txt")
created <- list.files(OUT, recursive = TRUE, full.names = FALSE)
write_csv(data.frame(file = created, stringsAsFactors = FALSE), "13_132_OUTPUT_INDEX_v02.csv")

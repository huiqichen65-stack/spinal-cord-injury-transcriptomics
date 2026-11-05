source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
ROOT <- PROJECT_ROOT
STEP <- "142_GSE281205_lipid_composition_specimen_robustness"
OUT <- file.path(ROOT, "results", STEP)
FIG <- file.path(OUT, "figures")
LOGDIR <- file.path(OUT, "logs")
VER <- "v01"
EXPECTED_GSM <- paste0("GSM", 8614662:8614667)
EXPECTED_DPI <- setNames(rep(c(14L, 28L), each = 3L), EXPECTED_GSM)
MIN_SET_GENES <- 5L
PRIOR_COUNT <- 2
PRIMARY_CT <- "contrast_A_target_vs_other_microglia_131"
LIPID_IDENTITY <- c("Clec7a", "Apoe", "Trem2", "Lpl", "Spp1", "Gpnmb", "Cst7", "Tyrobp", "P2ry12", "Tmem119",
  "Hexb", "Fcrls", "Cx3cr1")

D20 <- file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring")
D131 <- file.path(ROOT, "results", "131_GSE281205_myeloid_subtype_annotation_and_microglia_groups")
D132 <- file.path(ROOT, "results", "132_GSE281205_original_UMI_pseudobulk_signature_validation")
FILE_SETS <- file.path(D20, "20_AUCell_lipid_gene_set_manifest_v03.csv")
FILE_COMP <- file.path(D131, "05_131_GSM_primary_group_composition_v01.csv")
FILE_PB <- setNames(vapply(c(14L, 28L), function(t) {
  file.path(D132, "RDS", sprintf("132_pseudobulk_%s_%ddpi_v01.rds", PRIMARY_CT, t))
}, character(1)), c("14", "28"))
INPUTS <- c(FILE_SETS, FILE_COMP, unname(FILE_PB))

must <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
must(dir.exists(ROOT), paste("Missing root:", ROOT))
missing <- INPUTS[!file.exists(INPUTS)]
must(length(missing) == 0L, paste0("Required input file(s) missing. Run Step20/131/132 beforehand:\n",
  paste(missing, collapse = "\n")))
must(requireNamespace("edgeR", quietly = TRUE), "Missing Bioconductor edgeR: BiocManager::install(\"edgeR\")")
for (p in c(OUT, FIG, LOGDIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(LOGDIR, "00_142_run_log_v01.log")
cat("", file = LOG)
logmsg <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n", sep = "")
  cat(s, "\n", sep = "", file = LOG, append = TRUE)
  invisible(s)
}
logmsg("BEGIN Step142 v01 ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logmsg("ROOT=", ROOT, " | edgeR=", as.character(utils::packageVersion("edgeR")))
write_csv <- function(x, name) {
  must(is.data.frame(x), paste("Not a data frame:", name))
  utils::write.csv(x, file.path(OUT, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  logmsg("SAVED ", name, " | rows=", nrow(x))
  invisible(file.path(OUT, name))
}
read_csv <- function(p) utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
fmt <- function(x, n = 4) ifelse(is.finite(x), sprintf(paste0("%.", n, "f"), x), "NA")
mean_ok <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
min_ok <- function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
max_ok <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_

sets <- read_csv(FILE_SETS)
must(all(c("gene_set", "gene") %in% names(sets)), "Step20 manifest requires gene_set and gene columns")
sets$gene_set <- trimws(as.character(sets$gene_set))
sets$gene <- trimws(as.character(sets$gene))
sets <- unique(sets[nzchar(sets$gene_set) & nzchar(sets$gene), c("gene_set", "gene")])
set_order <- unique(sets$gene_set)
must(length(set_order) == 11L, paste("Expected 11 Step20 lipid programs, found", length(set_order)))

comp <- read_csv(FILE_COMP)
cols_comp <- c("gsm", "group_131", "n_cells", "dpi", "total_candidate_myeloid_in_GSM")
must(all(cols_comp %in% names(comp)), paste("Step131 composition missing:", paste(setdiff(cols_comp,
  names(comp)), collapse = ", ")))
comp$gsm <- as.character(comp$gsm)
comp$group_131 <- as.character(comp$group_131)
comp$n_cells <- suppressWarnings(as.integer(comp$n_cells))
comp$dpi <- suppressWarnings(as.integer(comp$dpi))
comp$total_candidate_myeloid_in_GSM <- suppressWarnings(as.integer(comp$total_candidate_myeloid_in_GSM))
must(identical(sort(unique(comp$gsm)), EXPECTED_GSM), "Unexpected GSM identity or count")
must(all(is.finite(comp$n_cells) & comp$n_cells >= 0), "Invalid cell count")
must(!anyDuplicated(paste(comp$gsm, comp$group_131, sep = "|")), "Duplicate GSM by Step131 state rows")
must(all(comp$dpi == EXPECTED_DPI[comp$gsm]), "GSM to stage disagreement")
needed_states <- c("TARGET_CLEC7A_HIGH_MICROGLIA", "REFERENCE_P2RY12_HIGH_MICROGLIA", "PROLIFERATING_MICROGLIA",
  "EXCLUDED_AMBIGUOUS_OR_QC")
must(all(needed_states %in% comp$group_131), "Missing required Step131 states for denominator-sensitivity analysis")

provenance <- data.frame(input_type = c("frozen_11_lipid_sets", "GSM_state_composition", "pseudobulk_original_UMI_14dpi",
  "pseudobulk_original_UMI_28dpi"), file = normalizePath(INPUTS, winslash = "/", mustWork = TRUE),
  bytes = file.info(INPUTS)$size, md5 = unname(tools::md5sum(INPUTS)), stringsAsFactors = FALSE)
write_csv(provenance, "01_142_input_provenance_v01.csv")

comp_rows <- list()
for (g in EXPECTED_GSM) {
  d <- comp[comp$gsm == g, , drop = FALSE]
  n_by_state <- setNames(d$n_cells, d$group_131)
  get_n <- function(nm) if (nm %in% names(n_by_state))
    n_by_state[[nm]] else 0L
  target <- get_n("TARGET_CLEC7A_HIGH_MICROGLIA")
  reference <- get_n("REFERENCE_P2RY12_HIGH_MICROGLIA")
  proliferating <- get_n("PROLIFERATING_MICROGLIA")
  excluded <- get_n("EXCLUDED_AMBIGUOUS_OR_QC")
  denom_all <- unique(d$total_candidate_myeloid_in_GSM)
  must(length(denom_all) == 1L && is.finite(denom_all) && denom_all > 0, paste("Bad candidate-myeloid denominator:",
    g))
  must(sum(d$n_cells) == denom_all, paste("Step131 group counts do not sum to denominator:", g))
  d_highconf <- denom_all - excluded
  d_all_micro <- target + reference + proliferating
  d_core_micro <- target + reference
  denom <- c(candidate_myeloid_all = denom_all, candidate_myeloid_QC_included = d_highconf, all_curated_microglia = d_all_micro,
    target_plus_P2ry12_microglia = d_core_micro)
  must(all(denom >= target & denom > 0), paste("Invalid denominator:", g))
  for (k in names(denom)) comp_rows[[length(comp_rows) + 1L]] <- data.frame(gsm = g, dpi = EXPECTED_DPI[[g]],
    denominator_type = k, n_Clec7a_high = target, n_denominator = as.integer(denom[[k]]), n_reference_microglia = reference,
    n_proliferating_microglia = proliferating, n_QC_excluded = excluded, proportion = target/denom[[k]],
    stringsAsFactors = FALSE)
}
composition <- do.call(rbind, comp_rows)
write_csv(composition, "06_142_per_GSM_composition_denominator_sensitivity_v01.csv")

exact_permutation <- function(v, groups) {
  must(length(v) == 6L && all(is.finite(v)) && identical(sort(as.integer(table(groups))), c(3L, 3L)),
    "Exact label permutation requires 3+3 independent GSM-level values")
  idxs <- utils::combn(seq_along(v), 3L)
  obs <- mean(v[groups == 28L]) - mean(v[groups == 14L])
  permutations <- apply(idxs, 2L, function(sel) mean(v[sel]) - mean(v[-sel]))
  p <- mean(abs(permutations) >= abs(obs) - 1e-12)
  c(observed_diff = obs, p_two_sided = p, n_permutations = length(permutations))
}
co_summary <- list()
co_loo <- list()
for (den in unique(composition$denominator_type)) {
  z <- composition[composition$denominator_type == den, ]
  z <- z[match(EXPECTED_GSM, z$gsm), ]
  must(identical(z$gsm, EXPECTED_GSM), "Lost GSM in composition")
  a <- z$proportion[z$dpi == 14L]
  b <- z$proportion[z$dpi == 28L]
  test <- exact_permutation(z$proportion, z$dpi)
  loo_diff <- c()
  for (i in 1:3) for (j in 1:3) {
    diff <- mean(b[-j]) - mean(a[-i])
    co_loo[[length(co_loo) + 1L]] <- data.frame(denominator_type = den, omitted_14dpi = EXPECTED_GSM[i],
      omitted_28dpi = EXPECTED_GSM[3L + j], n_per_stage_after_leave_one_out = 2L, delta_28_minus_14 = diff,
      direction_positive = diff > 0, stringsAsFactors = FALSE)
    loo_diff <- c(loo_diff, diff)
  }
  pooled_a <- sum(z$n_Clec7a_high[z$dpi == 14L])/sum(z$n_denominator[z$dpi == 14L])
  pooled_b <- sum(z$n_Clec7a_high[z$dpi == 28L])/sum(z$n_denominator[z$dpi == 28L])
  co_summary[[length(co_summary) + 1L]] <- data.frame(denominator_type = den, n_GSM_14 = 3L, n_GSM_28 = 3L,
    mean_fraction_14dpi = mean(a), sd_fraction_14dpi = stats::sd(a), range_14dpi_min = min(a), range_14dpi_max = max(a),
    mean_fraction_28dpi = mean(b), sd_fraction_28dpi = stats::sd(b), range_28dpi_min = min(b), range_28dpi_max = max(b),
    delta_28_minus_14_GSM_equally_weighted = test["observed_diff"], exact_permutation_p_two_sided = test["p_two_sided"],
    n_exact_label_assignments = as.integer(test["n_permutations"]), all_28dpi_higher_than_all_14dpi = min(b) >
      max(a), n_28dpi_GSM_above_max_14dpi = sum(b > max(a)), leave_one_each_stage_min_delta = min(loo_diff),
    leave_one_each_stage_max_delta = max(loo_diff), leave_one_each_stage_all_positive = all(loo_diff >
      0), pooled_cell_fraction_14dpi_DESCRIPTIVE_ONLY = pooled_a, pooled_cell_fraction_28dpi_DESCRIPTIVE_ONLY = pooled_b,
    pooled_cell_fraction_delta_DESCRIPTIVE_ONLY = pooled_b - pooled_a, unit = "GSM; pooled cells are not independent replicates",
    inference = "exploratory; label exchangeability/batch confounding unverified", stringsAsFactors = FALSE)
}
co_summary <- do.call(rbind, co_summary)
co_loo <- do.call(rbind, co_loo)
write_csv(co_summary, "07_142_stage_composition_exact_permutation_v01.csv")
write_csv(co_loo, "08_142_stage_composition_leave_one_out_v01.csv")
logmsg("Composition: primary 14->28 GSM-equal delta=", fmt(co_summary$delta_28_minus_14_GSM_equally_weighted[1]),
  " | exact p=", fmt(co_summary$exact_permutation_p_two_sided[1]))

exact_signflip <- function(x) {
  x <- as.numeric(x)
  must(length(x) >= 2L && all(is.finite(x)), "Insufficient data for sign flip")
  masks <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(x))))
  val <- as.numeric(masks %*% x)/length(x)
  obs <- mean(x)
  mean(abs(val) >= abs(obs) - 1e-12)
}
lipid_per_sample <- list()
lipid_coverage <- list()
input_qc <- list()
for (stage in c(14L, 28L)) {
  f <- FILE_PB[[as.character(stage)]]
  logmsg("Reading cached original-UMI pseudobulk (not single-cell): ", f)
  o <- readRDS(f)
  must(is.list(o) && all(c("original_integer_UMI_by_gene_symbol", "library_metadata", "effective_library_sizes") %in%
    names(o)), paste("Unrecognized Step132 v01 RDS structure:", f))
  mat <- o$original_integer_UMI_by_gene_symbol
  must(is.matrix(mat) && is.numeric(mat), paste("Pseudobulk counts not a numeric matrix:", f))
  ids <- EXPECTED_GSM[EXPECTED_DPI == stage]
  wanted <- as.vector(rbind(paste0(ids, "__OTHER"), paste0(ids, "__TARGET")))
  must(identical(colnames(mat), wanted), paste("Unexpected Step132 paired library order:", f))
  must(nrow(mat) > 10000L && ncol(mat) == 6L && !anyDuplicated(rownames(mat)) && all(is.finite(mat)) &&
    all(mat >= 0) && all(abs(mat - round(mat)) < 1e-08), paste("Step132 pseudobulk gene/count integrity failed:",
    f))
  m <- as.data.frame(o$library_metadata)
  must(all(c("library_id", "gsm", "group") %in% names(m)) && identical(as.character(m$library_id),
    wanted), paste("Missing/misaligned Step132 library metadata:", f))
  eff <- as.numeric(o$effective_library_sizes)
  must(length(eff) == 6L && all(is.finite(eff) & eff > 0), paste("Missing Step132 effective TMM library sizes:",
    f))

  lcpm <- edgeR::cpm(mat, lib.size = eff, log = TRUE, prior.count = PRIOR_COUNT)
  must(identical(rownames(lcpm), rownames(mat)) && all(is.finite(lcpm)), "Normalized pseudobulk logCPM invalid")
  delta <- sapply(seq_along(ids), function(i) lcpm[, 2L * i] - lcpm[, 2L * i - 1L])
  if (is.null(dim(delta)))
    delta <- matrix(delta, ncol = length(ids))
  rownames(delta) <- rownames(mat)
  colnames(delta) <- ids

  input_qc[[length(input_qc) + 1L]] <- data.frame(dpi = stage, gsm = rep(ids, each = 2L), group = rep(c("OTHER",
    "TARGET"), length(ids)), library_id = wanted, library_total_raw_UMI = as.numeric(colSums(mat)),
    library_effective_TMM_size = eff, n_genes_with_UMI = as.integer(colSums(mat > 0)), n_total_genes = nrow(mat),
    stringsAsFactors = FALSE)

  for (program in set_order) {
    genes <- unique(sets$gene[sets$gene_set == program])
    mapped <- intersect(genes, rownames(mat))
    absent <- setdiff(genes, mapped)
    no_marker <- setdiff(mapped, LIPID_IDENTITY)
    lipid_coverage[[length(lipid_coverage) + 1L]] <- data.frame(dpi = stage, gene_set = program,
      n_frozen_genes = length(genes), n_mapped_in_GSE281205 = length(mapped), n_marker_excluded_mapped = length(no_marker),
      n_missing = length(absent), missing_genes = paste(absent, collapse = ";"), gene_list_frozen_source = "Step20 v03 (original integrated atlas)",
      can_score = length(mapped) >= MIN_SET_GENES, stringsAsFactors = FALSE)
    for (g in ids) {
      if (length(mapped) >= MIN_SET_GENES) {
        x <- delta[mapped, g]
        xx <- delta[no_marker, g]
        res <- c(mean_delta = mean(x), median_delta = stats::median(x), frac_up = mean(x > 0),
          marker_excluded_mean = if (length(xx) >= MIN_SET_GENES) mean(xx) else NA_real_, marker_excluded_fraction = if (length(xx) >=
          MIN_SET_GENES) mean(xx > 0) else NA_real_)
      } else res <- c(mean_delta = NA_real_, median_delta = NA_real_, frac_up = NA_real_, marker_excluded_mean = NA_real_,
        marker_excluded_fraction = NA_real_)
      lipid_per_sample[[length(lipid_per_sample) + 1L]] <- data.frame(gene_set = program, dpi = stage,
        gsm = g, comparison = PRIMARY_CT, n_frozen_genes = length(genes), n_mapped_genes = length(mapped),
        n_marker_excluded_genes = length(no_marker), mean_TARGET_minus_OTHER_logCPM = unname(res["mean_delta"]),
        median_TARGET_minus_OTHER_logCPM = unname(res["median_delta"]), fraction_genes_TARGET_higher = unname(res["frac_up"]),
        mean_delta_excluding_identity_markers = unname(res["marker_excluded_mean"]), fraction_genes_higher_excluding_markers = unname(res["marker_excluded_fraction"]),
        metric = "paired mean across gene set; log2(TMM-CPM); NOT AUCell", stringsAsFactors = FALSE)
    }
  }
  rm(o, mat, m, lcpm, delta)
  invisible(gc(verbose = FALSE))
}
coverage <- do.call(rbind, lipid_coverage)
per_gsm <- do.call(rbind, lipid_per_sample)
qc <- do.call(rbind, input_qc)
write_csv(coverage, "02_142_lipid_gene_set_coverage_v01.csv")
write_csv(per_gsm, "03_142_per_GSM_lipid_program_effects_v01.csv")
write_csv(qc, "09_142_pseudobulk_library_integrity_v01.csv")

lipid_summ <- list()
lipid_loo <- list()
for (program in set_order) for (stage in c(14L, 28L)) {
  x <- per_gsm[per_gsm$gene_set == program & per_gsm$dpi == stage, ]
  x <- x[match(EXPECTED_GSM[EXPECTED_DPI == stage], x$gsm), ]
  vals <- as.numeric(x$mean_TARGET_minus_OTHER_logCPM)
  marker_vals <- as.numeric(x$mean_delta_excluding_identity_markers)
  good <- all(is.finite(vals))
  if (good) {
    sf_p <- exact_signflip(vals)
    loo <- vapply(seq_along(vals), function(i) mean(vals[-i]), numeric(1))
    se <- stats::sd(vals)/sqrt(3)
    tcrit <- stats::qt(0.975, df = 2)
    ci <- mean(vals) + c(-1, 1) * tcrit * se
    for (i in 1:3) lipid_loo[[length(lipid_loo) + 1L]] <- data.frame(gene_set = program, dpi = stage,
      omitted_gsm = x$gsm[i], n_GSM_remaining = 2L, leave_one_out_mean_delta = loo[i], same_direction_as_full = sign(loo[i]) ==
        sign(mean(vals)), stringsAsFactors = FALSE)
  } else {
    sf_p <- NA_real_
    loo <- rep(NA_real_, 3L)
    ci <- rep(NA_real_, 2L)
  }
  lipid_summ[[length(lipid_summ) + 1L]] <- data.frame(gene_set = program, dpi = stage, n_GSM = 3L,
    n_mapped_genes = unique(x$n_mapped_genes)[1], mean_paired_logCPM_delta = mean_ok(vals), sd_across_GSM = if (good)
      stats::sd(vals) else NA_real_, median_GSM_delta = if (good)
      stats::median(vals) else NA_real_, n_GSM_positive = sum(vals > 0, na.rm = TRUE), n_GSM_negative = sum(vals < 0, na.rm = TRUE),
    min_GSM_delta = min_ok(vals), max_GSM_delta = max_ok(vals), signflip_exact_two_sided_p = sf_p,
    exploratory_t_CI95_lower = ci[1], exploratory_t_CI95_upper = ci[2], n_GSM_positive_marker_excluded = sum(marker_vals >
      0, na.rm = TRUE), mean_marker_excluded_delta = mean_ok(marker_vals), LOO_mean_min = min_ok(loo),
    LOO_mean_max = max_ok(loo), LOO_all_keep_full_direction = if (good)
      all(sign(loo) == sign(mean(vals))) else NA, statistical_unit = "GSM; within-GSM target vs other microglia", status = if (good)
      "complete_exploratory" else "insufficient_gene_set_coverage", stringsAsFactors = FALSE)
}
lipid_summ <- do.call(rbind, lipid_summ)
lipid_loo <- if (length(lipid_loo)) do.call(rbind, lipid_loo) else data.frame()
lipid_summ$BH_FDR_across_22_program_stage_tests <- stats::p.adjust(lipid_summ$signflip_exact_two_sided_p,
  method = "BH")
write_csv(lipid_summ, "04_142_lipid_program_stage_robustness_v01.csv")
write_csv(lipid_loo, "05_142_lipid_program_leave_one_GSM_out_v01.csv")

pooled <- lapply(set_order, function(program) {
  y <- per_gsm[per_gsm$gene_set == program, ]
  v <- y$mean_TARGET_minus_OTHER_logCPM
  data.frame(gene_set = program, n_GSM = sum(is.finite(v)), n_GSM_positive = sum(v > 0, na.rm = TRUE),
    mean_GSM_paired_delta = mean_ok(v), exploratory_6_GSM_signflip_p = if (length(v) == 6L && all(is.finite(v)))
      exact_signflip(v) else NA_real_, limitation = "6 GSM pooled from two stages; not equivalent to six stage-matched animals",
    stringsAsFactors = FALSE)
})
pooled <- do.call(rbind, pooled)
pooled$BH_FDR_across_11_sets <- stats::p.adjust(pooled$exploratory_6_GSM_signflip_p, method = "BH")
write_csv(pooled, "10_142_lipid_program_six_GSM_descriptive_consistency_v01.csv")

plot_to_files <- function(stem, fun, width = 11, height = 7) {
  f_pdf <- file.path(FIG, paste0(stem, ".pdf"))
  f_png <- file.path(FIG, paste0(stem, ".png"))
  grDevices::pdf(f_pdf, width = width, height = height, onefile = FALSE, family = "Helvetica")
  tryCatch(fun(), finally = grDevices::dev.off())
  grDevices::png(f_png, width = round(width * 250), height = round(height * 250), res = 250)
  tryCatch(fun(), finally = grDevices::dev.off())

  f_tif <- file.path(FIG, paste0(stem, ".tiff"))
  tryCatch({
    grDevices::tiff(f_tif, width = round(width * 300), height = round(height * 300), res = 300, compression = "lzw")
    tryCatch(fun(), finally = grDevices::dev.off())
  }, error = function(e) logmsg("TIFF skipped: ", conditionMessage(e)))
  logmsg("SAVED figure ", stem, " [PDF, PNG", if (file.exists(f_tif))
    ", TIFF" else "", "]")
}

plot_to_files("Fig_142_GSE281205_lipid_program_GSM_effects_v01", function() {
  op <- graphics::par(mar = c(5, 12, 3.5, 2), las = 1)
  on.exit(graphics::par(op), add = TRUE)

  z <- matrix(NA_real_, nrow = length(set_order), ncol = 6L, dimnames = list(set_order, EXPECTED_GSM))
  for (i in seq_len(nrow(per_gsm))) {
    r <- match(per_gsm$gene_set[i], set_order)
    c <- match(per_gsm$gsm[i], EXPECTED_GSM)
    z[r, c] <- per_gsm$mean_TARGET_minus_OTHER_logCPM[i]
  }
  lim <- max(abs(z), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0)
    lim <- 0.01

  graphics::plot(NA, xlim = c(0.5, 6.5), ylim = c(0.5, length(set_order) + 0.5), xaxt = "n", yaxt = "n",
    xlab = "", ylab = "", main = "Per-GSM lipid-program paired effects (not AUCell)")
  graphics::axis(1, at = 1:6, labels = paste0(sub("GSM861", "", EXPECTED_GSM), "\n", rep(c("14 dpi",
    "28 dpi"), each = 3L)), las = 1, cex.axis = 0.85)
  graphics::axis(2, at = seq_along(set_order), labels = set_order, las = 1, cex.axis = 0.84)
  graphics::abline(v = 3.5, lty = 2, col = "gray65")
  for (i in seq_len(nrow(z))) for (j in seq_len(ncol(z))) {
    v <- z[i, j]
    if (!is.finite(v))
      next
    graphics::points(j, i, pch = if (v >= 0)
      16 else 1, cex = 0.6 + 1.4 * sqrt(abs(v)/lim), col = if (v >= 0)
      "#176ea6" else "#bb6224", lwd = 1.3)
  }
  graphics::mtext("TARGET minus OTHER microglia | circle size = magnitude; filled blue = positive",
    side = 1, line = 3.2, cex = 0.75)
}, width = 12.8, height = 8)

plot_to_files("Fig_142_GSE281205_composition_by_GSM_denominator_v01", function() {
  den_order <- c("candidate_myeloid_all", "candidate_myeloid_QC_included", "all_curated_microglia",
    "target_plus_P2ry12_microglia")
  pretty_names <- c("All candidate myeloid", "QC-included candidate myeloid", "Curated microglia",
    "Target + P2ry12-high microglia")
  op <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3.8, 1), oma = c(0, 0, 1, 0))
  on.exit(graphics::par(op), add = TRUE)
  for (k in seq_along(den_order)) {
    d <- composition[composition$denominator_type == den_order[k], ]
    a <- d$proportion[d$dpi == 14L]
    b <- d$proportion[d$dpi == 28L]
    graphics::plot(c(rep(1, 3), rep(2, 3)), 100 * c(a, b), xlim = c(0.6, 2.4), ylim = c(0, 100),
      xaxt = "n", xlab = "", ylab = "Clec7a-high (%)", pch = 16, cex = 1.25, col = c(rep("#257eaf",
        3), rep("#d17f30", 3)), main = pretty_names[k])
    graphics::axis(1, at = 1:2, labels = c("14 dpi", "28 dpi"))
    graphics::segments(0.8, 100 * mean(a), 1.2, 100 * mean(a), col = "#257eaf", lwd = 3)
    graphics::segments(1.8, 100 * mean(b), 2.2, 100 * mean(b), col = "#d17f30", lwd = 3)
    graphics::text(c(rep(1, 3), rep(2, 3)) + c(-0.14, 0, 0.14, -0.14, 0, 0.14), 100 * c(a, b) + 4,
      labels = sub("GSM861", "", d$gsm), cex = 0.63)
    graphics::mtext(sprintf("GSM-equal difference = %+.1f pp", 100 * (mean(b) - mean(a))), side = 3,
      line = 0.15, cex = 0.72)
  }
  graphics::mtext("Each point is one GSM. Mean bars are specimen-weighted, not cell-weighted.", outer = TRUE,
    side = 3, line = -0.4, cex = 0.8)
}, width = 11, height = 9)

for (i in seq_len(nrow(co_summary))) {
  a <- co_summary[i, ]
}
for (stage in c(14L, 28L)) {
  x <- lipid_summ[lipid_summ$dpi == stage & lipid_summ$status == "complete_exploratory", ]
}
writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT, "11_142_REPORT_v01.txt"), useBytes = TRUE)
limits <- c("LIMITATIONS AND METHODS | Step142 v01", "Biological unit: GEO GSM/sample-library, not individual cells.",
  "Independent animals: NOT verified from the identifiers by this script.", "Groups are from Step131; no new clustering or identity redefinition.",
  "Lipid sets are Step20 v03 sets; existing source present_in_object flags", "refer to the OLD integrated atlas, not GSE281205. Mapping here uses symbols",
  "present in the Step132 GSE281205 original-count pseudobulk.", "Lipid-program metric: within GSM, mean over available genes of log2",
  "TMM CPM (Clec7a-high - other microglia), Step132 prior.count=2.", "This metric is not AUCell; shifts can reflect broad transcript differences.",
  "No claim that all genes or all 11 lipid programs replicate.", "Identity-gene exclusion sensitivity is reported per program, with coverage.",
  "Paired sign flips provide exploratory exact 2-sided sample-level P (n=3).", "A sign-flip test assumes sample-level symmetric null errors; it is exploratory.",
  "Stage permutation provides exploratory exact 2-sided P (3+3 GSM).", "Stage labels might not be exchangeable due to experimental/batch conditions.",
  "95% t intervals for n=3, when given, are descriptive and unstable.", "No cell-wise binomial P; no independence assigned to thousands of cells.",
  "Four denominators test how QC and state definitions affect proportions.", "Temporal data comprise 14 and 28 dpi only, not full injury progression.",
  "28 dpi GEO age metadata discrepancy should be resolved or disclosed.", "Potential animal pooling and sample identity must be checked externally.",
  "The prior original Figure 2D 88-90% is from DIFFERENT datasets/denominators", "and cannot be validated directly by these six GSM proportions.",
  "No direct functional lipid assay, and no causal/experimental validation.")
writeLines(limits, file.path(OUT, "12_142_methods_limitations_v01.txt"), useBytes = TRUE)
writeLines(capture.output(utils::sessionInfo()), file.path(OUT, "13_142_sessionInfo_v01.txt"), useBytes = TRUE)

all_output <- list.files(OUT, recursive = TRUE, full.names = TRUE)
all_output <- all_output[file.info(all_output)$isdir == FALSE]
idx <- data.frame(relative_path = substring(normalizePath(all_output, winslash = "/", mustWork = TRUE),
  nchar(normalizePath(OUT, winslash = "/", mustWork = TRUE)) + 2L), bytes = file.info(all_output)$size,
  extension = tolower(tools::file_ext(all_output)), stringsAsFactors = FALSE)
write_csv(idx, "00_142_OUTPUT_INDEX_v01.csv")
logmsg("DONE Step142 v01 | Output: ", OUT)

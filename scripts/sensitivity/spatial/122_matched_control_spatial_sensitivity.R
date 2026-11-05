options(stringsAsFactors = FALSE, warn = 1)
ROOT <- "D:/zz"
VERSION <- "v02"
STEP <- "122_matched_controls_spatial_dependence"
OUT <- file.path(ROOT, "results", STEP)
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

SEED <- 12209L
N_CONTROL_PROGRAMS <- 500L
MATCH_TOP_K <- 30L
MIN_SPOTS_PER_SECTION <- 100L
MIN_SPOTS_PER_REGION_WITHIN_SECTION <- 30L
MIN_SECTIONS_PER_REGION_SUMMARY <- 3L
SPATIAL_RADIUS_MULTIPLIER <- 1.55
MAX_FROZEN_SCORE_DIFFERENCE <- 1e-05
EXPECTED_SECTIONS <- 32L
EXPECTED_AREAS <- 8L
STAGES <- c("3dpi", "7dpi", "14dpi", "28dpi")
SAMPLES <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1",
  "28dpi_rep2")
SCENARIOS <- c("all32_GEO_labels", "inferred_injury_only27")
TARGET <- "Clec7a+ lipid-handling Microglia"
set.seed(SEED)

IN <- list(spatial = file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation",
  "12_64_GSE195783_spatial_region_annotated_v03.rds"), rctd = file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping",
  "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv.gz"), lipid_manifest = file.path(ROOT,
  "results", "73_GSE195783_Clec7a_lipid_spatial_coupling", "01_73_lipid_program_gene_matching_v03.csv"),
  frozen_scores = file.path(ROOT, "results", "73_GSE195783_Clec7a_lipid_spatial_coupling", "02_73_Clec7a_lipid_spot_scores_section_aware_v03.csv.gz"),
  portions117 = file.path(ROOT, "results", "117_capture_area_heterogeneity", "tables", "01_117_frozen_32piece_status_and_provenance_v02.csv"),
  rctd_features120 = file.path(ROOT, "results", "120_RCTD_lipid_feature_overlap", "tables", "03_120_reconstructed_RCTD_selected_feature_lists_v02.csv"),
  spotQC60 = file.path(ROOT, "results", "60_GSE195783_Visium_import_coordinate_reconstruction", "03_60_GSE195783_final_spot_metadata_v05.csv.gz"))

LOG <- file.path(OUT, paste0(STEP, "_", VERSION, ".log"))
logcon <- file(LOG, open = "wt", encoding = "UTF-8")

sink(logcon, type = "output", split = TRUE)
step_log <- function(...) cat(paste0(..., collapse = ""), "\n", sep = "")
tryCatch({

  need_pkg <- c("Matrix", "SeuratObject")
  missing_pkgs <- need_pkg[!vapply(need_pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs))
    stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(SeuratObject))

  read_csv <- function(p) {
    if (!file.exists(p))
      stop("Missing required input: ", p)
    read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  }
  save_csv <- function(x, name) {
    write.csv(x, file.path(TAB, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
  }
  must <- function(x, cols, label) {
    z <- setdiff(cols, names(x))
    if (length(z))
      stop(label, " missing columns: ", paste(z, collapse = ", "))
  }
  key <- function(df) paste(as.character(df$GSM), as.character(df$array_x), as.character(df$array_y),
    sep = "||")
  fz <- function(x) atanh(pmax(-0.999999, pmin(0.999999, x)))
  good_rho <- function(x, y, nmin = 10L) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < nmin || sd(x[ok]) < 1e-12 || sd(y[ok]) < 1e-12)
      return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  }
  medfin <- function(x) if (any(is.finite(x)))
    median(x[is.finite(x)]) else NA_real_
  meanfin <- function(x) if (any(is.finite(x)))
    mean(x[is.finite(x)]) else NA_real_
  zsd <- function(x) {
    q <- sd(x)
    if (!is.finite(q) || q < 1e-08)
      return(rep(0, length(x)))
    as.numeric((x - mean(x))/q)
  }
  partial_spearman <- function(x, y, covariates, min_n = 30L) {
    cc <- is.finite(x) & is.finite(y) & apply(is.finite(covariates), 1, all)
    if (sum(cc) < min_n)
      return(NA_real_)
    z <- cbind(Intercept = 1, covariates[cc, , drop = FALSE])
    if (nrow(z) <= ncol(z) + 4L)
      return(NA_real_)
    r <- qr(z)
    xx <- qr.resid(r, rank(x[cc], ties.method = "average"))
    yy <- qr.resid(r, rank(y[cc], ties.method = "average"))
    if (sd(xx) < 1e-12 || sd(yy) < 1e-12)
      return(NA_real_)
    as.numeric(cor(xx, yy))
  }
  spatial_covariates <- function(xx, yy) {
    a <- zsd(xx)
    b <- zsd(yy)
    cbind(x = a, y = b, x2 = a^2, xy = a * b, y2 = b^2)
  }
  spatial_graph <- function(x, y) {
    p <- cbind(as.numeric(x), as.numeric(y))
    d <- as.matrix(dist(p))
    diag(d) <- Inf
    mednn <- median(apply(d, 1L, min), na.rm = TRUE)
    if (!is.finite(mednn) || mednn <= 0)
      stop("Invalid nearest-neighbor distance")
    rad <- SPATIAL_RADIUS_MULTIPLIER * mednn
    a <- (d <= rad)
    diag(a) <- FALSE
    storage.mode(a) <- "double"
    deg <- rowSums(a)
    list(a = a, deg = deg, radius = rad, nn_median = mednn, isolated = sum(deg == 0))
  }
  moran_descriptive <- function(x, graph) {
    z <- x - mean(x, na.rm = TRUE)
    if (any(!is.finite(z)) || sum(z * z) < 1e-10)
      return(NA_real_)
    idx <- which(graph$deg > 0)
    if (!length(idx))
      return(NA_real_)
    lag <- as.numeric(graph$a %*% z)
    lag[idx] <- lag[idx]/graph$deg[idx]
    (length(z)/length(idx)) * sum(z[idx] * lag[idx])/sum(z * z)
  }
  exact_signflip <- function(area_effect_z) {
    v <- area_effect_z[is.finite(area_effect_z)]
    if (length(v) < 3L || length(v) > 15L)
      return(NA_real_)
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(v))))
    vals <- as.numeric(signs %*% v)/length(v)
    mean(vals >= mean(v) - 1e-12)
  }

  inputs <- data.frame(input_name = names(IN), path = unlist(IN, use.names = FALSE), required = names(IN) !=
    "spotQC60", stringsAsFactors = FALSE)
  inputs$exists <- file.exists(inputs$path)
  inputs$bytes <- ifelse(inputs$exists, file.info(inputs$path)$size, NA_real_)
  save_csv(inputs, "00_122_INPUT_INDEX_v02.csv")
  if (any(inputs$required & !inputs$exists)) {
    stop("Missing mandatory inputs:\n", paste(inputs$path[inputs$required & !inputs$exists], collapse = "\n"))
  }

  spot_scores <- read_csv(IN$frozen_scores)
  rctd <- read_csv(IN$rctd)
  gene_map <- read_csv(IN$lipid_manifest)
  p117 <- read_csv(IN$portions117)
  features <- read_csv(IN$rctd_features120)
  must(spot_scores, c("GSM", "array_x", "array_y", "sample_label", "timepoint", "section_uid_internal",
    "score_Core_lipid_handling", "Clec7a_RCTDnorm"), "Step73")
  must(rctd, c("GSM", "array_x", "array_y", "section_uid_internal", "spatial_region_v03", TARGET),
    "Step65")
  must(gene_map, c("module", "requested_gene", "matched_gene", "detected"), "Step73 genes")
  must(p117, c("section_uid_internal", "inferred_uninjured"), "Step117")
  must(features, c("sample_label", "gene", "RCTD_regression_fit", "RCTD_platform_normalization"), "Step120")

  if (anyDuplicated(key(spot_scores)) || anyDuplicated(key(rctd)))
    stop("Duplicate spatial keys")
  j <- match(key(spot_scores), key(rctd))
  if (anyNA(j))
    stop("Step73 spot coordinates not fully present in Step65 RCTD")
  if (any(abs(as.numeric(spot_scores$Clec7a_RCTDnorm) - as.numeric(rctd[[TARGET]][j])) > 1e-07)) {
    stop("Step73 RCTD signal differs from Step65 values")
  }
  spot_scores$spatial_region_v03 <- as.character(rctd$spatial_region_v03[j])
  myeloid_cols <- grep("Microglia|Macrophage|Neutrophil", names(rctd), value = TRUE)
  if (!length(myeloid_cols))
    stop("No myeloid RCTD columns found")
  spot_scores$total_myeloid_weight <- rowSums(sapply(myeloid_cols, function(g) as.numeric(rctd[[g]][j])))
  spot_scores$Clec7a_RCTDnorm <- as.numeric(spot_scores$Clec7a_RCTDnorm)
  spot_scores$score_Core_lipid_handling <- as.numeric(spot_scores$score_Core_lipid_handling)
  spot_scores$array_x <- as.numeric(spot_scores$array_x)
  spot_scores$array_y <- as.numeric(spot_scores$array_y)
  if (anyNA(spot_scores$section_uid_internal) || any(!nzchar(spot_scores$section_uid_internal)))
    stop("Incomplete section IDs")
  if (!setequal(unique(spot_scores$sample_label), SAMPLES))
    stop("Expected 8 Step73 capture areas")
  if (length(unique(spot_scores$section_uid_internal)) != EXPECTED_SECTIONS)
    stop("Expected 32 Step73 sections")
  if (anyNA(spot_scores$spatial_region_v03))
    stop("Missing Step65 spatial regions")
  q117 <- match(spot_scores$section_uid_internal, p117$section_uid_internal)
  if (anyNA(q117))
    stop("Step117 has no status for one or more sections")
  spot_scores$inferred_uninjured <- as.logical(p117$inferred_uninjured[q117])
  if (anyNA(spot_scores$inferred_uninjured))
    stop("Invalid Step117 inferred_uninjured flags")
  if (sum(tapply(spot_scores$inferred_uninjured, spot_scores$section_uid_internal, function(x) unique(x)[1])) !=
    5L)
    stop("Step117 sensitivity expects 5 inferred uninjured portions")

  obj <- readRDS(IN$spatial)
  if (!("Spatial" %in% names(obj@assays)))
    stop("Spatial assay absent")
  expr <- tryCatch(SeuratObject::LayerData(obj, assay = "Spatial", layer = "data"), error = function(e) NULL)
  if (is.null(expr))
    expr <- tryCatch(SeuratObject::GetAssayData(obj, assay = "Spatial", slot = "data"), error = function(e) NULL)
  if (is.null(expr) || !nrow(expr) || !ncol(expr))
    stop("Cannot load Spatial log-normalized data layer")
  if (!inherits(expr, "Matrix"))
    expr <- Matrix::Matrix(expr, sparse = TRUE)
  md <- obj@meta.data
  if (!all(colnames(expr) %in% rownames(md)))
    stop("Step64 metadata not aligned to expression columns")
  md <- md[colnames(expr), , drop = FALSE]
  must(md, c("GSM", "array_x", "array_y", "sample_label"), "Step64 metadata")
  kk <- key(md)
  if (anyDuplicated(kk))
    stop("Step64 duplicate spatial coordinate keys")
  pos <- match(key(spot_scores), kk)
  if (anyNA(pos))
    stop("Missing Step73 spots in Step64 matrix; coordinate join invalid")
  if (!all(as.character(md$sample_label[pos]) == as.character(spot_scores$sample_label)))
    stop("Step64/73 sample labels do not agree")

  full_sg <- split(seq_len(ncol(expr)), factor(as.character(md$sample_label), levels = SAMPLES), drop = TRUE)
  if (length(full_sg) != EXPECTED_AREAS)
    stop("Step64 missing full capture-area groups")
  expr_full <- expr
  expr <- expr_full[, pos, drop = FALSE]
  rm(obj, md)
  invisible(gc())
  stopifnot(ncol(expr) == nrow(spot_scores), identical(key(spot_scores), kk[pos]))

  qc <- if (file.exists(IN$spotQC60))
    read_csv(IN$spotQC60) else NULL
  spot_scores$nCount_Spatial <- NA_real_
  spot_scores$nFeature_Spatial <- NA_real_
  if (!is.null(qc)) {
    must(qc, c("GSM", "array_x", "array_y", "nCount_Spatial", "nFeature_Spatial"), "Step60 QC")
    if (anyDuplicated(key(qc)))
      stop("Duplicate Step60 QC coordinate keys")
    jq <- match(key(spot_scores), key(qc))
    if (anyNA(jq))
      stop("Incomplete Step60 QC join")
    spot_scores$nCount_Spatial <- as.numeric(qc$nCount_Spatial[jq])
    spot_scores$nFeature_Spatial <- as.numeric(qc$nFeature_Spatial[jq])
  } else {
    step_log("WARN: optional Step60 spot-QC file missing: combined QC-adjusted effect will be NA")
  }

  orig_tab <- gene_map[gene_map$module == "Core_lipid_handling_union" & as.logical(gene_map$detected),
    , drop = FALSE]
  orig <- unique(orig_tab$matched_gene[nzchar(orig_tab$matched_gene) & !is.na(orig_tab$matched_gene)])
  if (length(orig) != 37L)
    stop("Expected 37 scored core lipid genes; found ", length(orig))
  if (!all(tolower(orig) %in% tolower(rownames(expr))))
    stop("A scored gene is absent from Step64")
  orig <- rownames(expr)[match(tolower(orig), tolower(rownames(expr)))]
  if ("Clec7a" %in% orig)
    stop("Clec7a must not be scored")

  n <- nrow(expr)
  smp <- as.character(spot_scores$sample_label)
  sg <- split(seq_len(ncol(expr)), factor(smp, levels = SAMPLES), drop = TRUE)
  if (length(sg) != EXPECTED_AREAS)
    stop("Missing capture-area spots")
  gene_stats <- vector("list", length(sg))
  names(gene_stats) <- names(sg)
  for (s in names(sg)) {
    idx <- full_sg[[s]]
    m <- expr_full[, idx, drop = FALSE]
    mu <- as.numeric(Matrix::rowMeans(m))
    second <- as.numeric(Matrix::rowSums(m^2))
    sdv <- sqrt(pmax(0, (second - length(idx) * mu^2)/(length(idx) - 1)))
    det <- as.numeric(Matrix::rowSums(m > 0))/length(idx)
    gene_stats[[s]] <- data.frame(gene = rownames(expr), sample_label = s, mean_logexpr = mu, detection = det,
      sd_logexpr = sdv)
  }
  rm(expr_full)
  invisible(gc())
  all_stats <- do.call(rbind, gene_stats)
  save_csv(all_stats, "01_122_expression_detection_by_capture_v02.csv")

  metric_mean <- apply(sapply(gene_stats, function(x) x$mean_logexpr), 1, median)
  metric_det <- apply(sapply(gene_stats, function(x) x$detection), 1, median)
  metric_sdmin <- apply(sapply(gene_stats, function(x) x$sd_logexpr), 1, min)
  candidate <- data.frame(gene = rownames(expr), mean_logexpr = metric_mean, detection = metric_det,
    min_sd = metric_sdmin)

  feature_key <- tolower(features$gene)
  if (!setequal(unique(features$sample_label), SAMPLES))
    stop("Step120 area labels differ from Step73")
  features$selected_either <- as.logical(features$RCTD_regression_fit) | as.logical(features$RCTD_platform_normalization)
  if (anyNA(features$selected_either))
    stop("Bad Step120 feature status values")
  pattern_mat <- matrix(FALSE, nrow = nrow(candidate), ncol = length(SAMPLES), dimnames = list(tolower(candidate$gene),
    SAMPLES))
  for (s in SAMPLES) {
    ft <- features[features$sample_label == s & features$selected_either, , drop = FALSE]
    ft_gene <- intersect(tolower(ft$gene), rownames(pattern_mat))
    pattern_mat[ft_gene, s] <- TRUE
  }
  candidate$RCTD_pattern <- apply(pattern_mat, 1, function(x) paste(as.integer(x), collapse = ""))
  orig_idx <- match(tolower(orig), tolower(candidate$gene))
  candidate$original_lipid <- seq_len(nrow(candidate)) %in% orig_idx
  candidate$eligible_control <- !candidate$original_lipid & !(tolower(candidate$gene) %in% c("clec7a")) &
    !grepl("^(mt-|Rpl|Rps)", candidate$gene, ignore.case = TRUE) & candidate$min_sd > 1e-08 & candidate$detection >=
    0.005
  save_csv(candidate, "02_122_candidate_gene_pool_and_RCTD_pattern_v02.csv")
  if (sum(candidate$eligible_control) < 200L)
    stop("Insufficient eligible candidate genes")

  score_panel <- function(genes) {
    out <- numeric(ncol(expr))
    for (s in names(sg)) {
      idx <- sg[[s]]
      st <- gene_stats[[s]]
      gidx <- match(genes, st$gene)
      if (anyNA(gidx))
        stop("Unmatched gene in panel")
      keep <- which(is.finite(st$sd_logexpr[gidx]) & st$sd_logexpr[gidx] > 1e-08)
      if (length(keep) < 5L)
        stop("Fewer than 5 variable genes in scoring panel")

      g_use <- genes[keep]
      gidx <- gidx[keep]
      w <- 1/st$sd_logexpr[gidx]

      sub <- expr[g_use, idx, drop = FALSE]
      out[idx] <- (as.numeric(Matrix::crossprod(w, sub)) - sum(st$mean_logexpr[gidx] * w))/length(g_use)
    }
    out
  }
  computed_original <- score_panel(orig)
  score_err <- max(abs(computed_original - spot_scores$score_Core_lipid_handling), na.rm = TRUE)
  if (!is.finite(score_err) || score_err > MAX_FROZEN_SCORE_DIFFERENCE) {
    stop("Score method did not reproduce Step73. STOP instead of comparing non-equivalent controls.")
  }

  eligible <- which(candidate$eligible_control)
  scl_m <- sd(candidate$mean_logexpr[eligible])
  scl_d <- sd(candidate$detection[eligible])
  if (scl_m < 1e-10 || scl_d < 1e-10)
    stop("Matching metrics have no variation")
  opts <- vector("list", length(orig))
  for (i in seq_along(orig)) {
    x <- candidate[orig_idx[i], ]
    pool <- eligible[candidate$RCTD_pattern[eligible] == x$RCTD_pattern]
    if (length(pool) < 30L)
      stop("Fewer than 30 RCTD-pattern-matched control genes for ", orig[i], ". Do not relax overlap matching silently.")
    dist2 <- ((candidate$mean_logexpr[pool] - x$mean_logexpr)/scl_m)^2 + ((candidate$detection[pool] -
      x$detection)/scl_d)^2
    opts[[i]] <- pool[order(dist2)[seq_len(min(length(pool), max(60L, MATCH_TOP_K * 4L)))]]
  }

  target_order <- order(vapply(opts, length, integer(1)), decreasing = FALSE)
  controls <- vector("list", N_CONTROL_PROGRAMS)
  control_match_records <- vector("list", N_CONTROL_PROGRAMS)
  for (b in seq_len(N_CONTROL_PROGRAMS)) {
    used <- integer(0)
    selected <- integer(length(orig))
    for (ii in target_order) {
      cand <- setdiff(opts[[ii]], used)
      if (!length(cand))
        stop("Cannot form distinct size-matched controls; program ", b)
      cand <- cand[seq_len(min(MATCH_TOP_K, length(cand)))]

      sampled <- cand[sample.int(length(cand), 1L, prob = rev(seq_along(cand)))]
      selected[ii] <- sampled
      used <- c(used, sampled)
    }
    controls[[b]] <- candidate$gene[selected]
    control_match_records[[b]] <- data.frame(control_id = b, original_gene = orig, control_gene = candidate$gene[selected],
      original_mean = candidate$mean_logexpr[orig_idx], control_mean = candidate$mean_logexpr[selected],
      original_detection = candidate$detection[orig_idx], control_detection = candidate$detection[selected],
      original_pattern = candidate$RCTD_pattern[orig_idx], control_pattern = candidate$RCTD_pattern[selected])
  }
  control_matches <- do.call(rbind, control_match_records)
  save_csv(control_matches, "03_122_matched_control_gene_manifest_v02.csv")
  control_qc <- do.call(rbind, lapply(split(control_matches, control_matches$control_id), function(x) {
    data.frame(control_id = x$control_id[1], n_genes = nrow(x), n_unique = length(unique(x$control_gene)),
      rctd_pattern_matches = all(x$original_pattern == x$control_pattern), mean_absolute_logexpr_difference = mean(abs(x$original_mean -
        x$control_mean)), mean_absolute_detection_difference = mean(abs(x$original_detection -
        x$control_detection)))
  }))
  save_csv(control_qc, "04_122_matching_quality_by_control_v02.csv")
  if (any(!control_qc$rctd_pattern_matches | control_qc$n_unique != length(orig)))
    stop("Bad control matching")

  groups <- split(seq_len(nrow(spot_scores)), spot_scores$section_uid_internal)
  section_info <- do.call(rbind, lapply(names(groups), function(u) {
    k <- groups[[u]]
    data.frame(section_uid_internal = u, sample_label = unique(spot_scores$sample_label[k]), timepoint = unique(spot_scores$timepoint[k]),
      n_spots = length(k), inferred_uninjured = unique(spot_scores$inferred_uninjured[k]), stringsAsFactors = FALSE)
  }))
  if (any(section_info$n_spots < MIN_SPOTS_PER_SECTION))
    stop("One or more sections below minimum count")
  save_csv(section_info, "05_122_frozen_portion_manifest_v02.csv")

  graphs <- lapply(groups, function(idx) spatial_graph(spot_scores$array_x[idx], spot_scores$array_y[idx]))
  graph_qc <- do.call(rbind, lapply(names(graphs), function(u) {
    g <- graphs[[u]]
    data.frame(section_uid_internal = u, n_spots = length(groups[[u]]), radius = g$radius, median_nearest_neighbor = g$nn_median,
      isolated_spots = g$isolated, mean_neighbor_count = mean(g$deg))
  }))
  save_csv(graph_qc, "06_122_spatial_graph_QC_v02.csv")

  calc_section_rhos <- function(score, spatial = FALSE) {
    ret <- lapply(names(groups), function(u) {
      idx <- groups[[u]]
      x <- spot_scores$Clec7a_RCTDnorm[idx]
      y <- score[idx]
      out <- data.frame(section_uid_internal = u, rho = good_rho(x, y))
      if (spatial) {
        loc <- spatial_covariates(spot_scores$array_x[idx], spot_scores$array_y[idx])
        out$rho_spatial_detrended <- partial_spearman(x, y, loc)
        covs <- cbind(loc, log1p(spot_scores$nCount_Spatial[idx]), spot_scores$nFeature_Spatial[idx],
          spot_scores$total_myeloid_weight[idx])
        out$rho_spatial_QC_myeloid_adjusted <- partial_spearman(x, y, covs)
        g <- graphs[[u]]
        out$moran_target <- moran_descriptive(x, g)
        out$moran_score <- moran_descriptive(y, g)

        nei <- as.numeric(g$a %*% y)
        nei[g$deg > 0] <- nei[g$deg > 0]/g$deg[g$deg > 0]
        nei[g$deg == 0] <- NA_real_
        out$rho_neighbor <- good_rho(x, nei)
      }
      out
    })
    z <- do.call(rbind, ret)
    merge(section_info, z, by = "section_uid_internal", sort = FALSE)
  }
  obs_sec <- calc_section_rhos(spot_scores$score_Core_lipid_handling, spatial = TRUE)
  save_csv(obs_sec, "07_122_observed_per_portion_coupling_spatial_sensitivity_v02.csv")

  region_grps <- split(seq_len(nrow(spot_scores)), paste(spot_scores$section_uid_internal, spot_scores$spatial_region_v03,
    sep = "@@@"))
  region_spots <- do.call(rbind, lapply(region_grps, function(idx) {
    u <- unique(spot_scores$section_uid_internal[idx])
    r <- unique(spot_scores$spatial_region_v03[idx])
    if (length(u) != 1L || length(r) != 1L)
      stop("Region group mixes sections/labels")
    x <- spot_scores$Clec7a_RCTDnorm[idx]
    y <- spot_scores$score_Core_lipid_handling[idx]
    ok <- length(idx) >= MIN_SPOTS_PER_REGION_WITHIN_SECTION
    cc <- spatial_covariates(spot_scores$array_x[idx], spot_scores$array_y[idx])
    data.frame(section_uid_internal = u, sample_label = unique(spot_scores$sample_label[idx]), timepoint = unique(spot_scores$timepoint[idx]),
      region = r, n_spots = length(idx), eligible = ok, rho = if (ok)
        good_rho(x, y) else NA_real_, rho_spatial_detrended = if (ok)
        partial_spearman(x, y, cc) else NA_real_, mean_target = mean(x), mean_lipid = mean(y), inferred_uninjured = unique(spot_scores$inferred_uninjured[idx]))
  }))
  save_csv(region_spots, "08_122_observed_within_portion_region_coupling_v02.csv")

  summarise_effect <- function(sec, effect_col, scenario) {
    z <- if (scenario == SCENARIOS[1])
      sec else sec[!sec$inferred_uninjured, , drop = FALSE]
    values <- z[[effect_col]]
    data.frame(scenario = scenario, metric = effect_col, n_portions = sum(is.finite(values)), n_capture_areas = length(unique(z$sample_label[is.finite(values)])),
      median_rho = medfin(values), median_fisher_z = medfin(fz(values[is.finite(values)])), mean_fisher_z = meanfin(fz(values[is.finite(values)])),
      proportion_positive = if (any(is.finite(values)))
        mean(values[is.finite(values)] > 0) else NA_real_)
  }
  obs_overall <- do.call(rbind, lapply(SCENARIOS, function(sc) do.call(rbind, lapply(c("rho", "rho_neighbor",
    "rho_spatial_detrended", "rho_spatial_QC_myeloid_adjusted"), function(v) summarise_effect(obs_sec,
    v, sc)))))
  save_csv(obs_overall, "09_122_observed_scenario_effect_summaries_v02.csv")

  cap_summary <- do.call(rbind, lapply(SCENARIOS, function(sc) {
    dat <- if (sc == SCENARIOS[1])
      obs_sec else obs_sec[!obs_sec$inferred_uninjured, , drop = FALSE]
    do.call(rbind, lapply(split(dat, dat$sample_label), function(s) {
      data.frame(scenario = sc, sample_label = s$sample_label[1], timepoint = s$timepoint[1], n_portions = nrow(s),
        n_spots = sum(s$n_spots), median_direct_rho = medfin(s$rho), median_direct_fisher_z = medfin(fz(s$rho)),
        median_neighbor_rho = medfin(s$rho_neighbor), median_spatial_detrended_rho = medfin(s$rho_spatial_detrended),
        median_spatial_QC_myeloid_adjusted_rho = medfin(s$rho_spatial_QC_myeloid_adjusted))
    }))
  }))
  save_csv(cap_summary, "10_122_observed_capture_area_effects_v02.csv")
  cap_signflip <- do.call(rbind, lapply(SCENARIOS, function(sc) {
    d <- cap_summary[cap_summary$scenario == sc, ]
    data.frame(scenario = sc, n_capture_areas = nrow(d), one_sided_exact_signflip_P_direct = exact_signflip(d$median_direct_fisher_z),
      one_sided_exact_signflip_P_spatial_detrended = exact_signflip(fz(d$median_spatial_detrended_rho)),
      note = "8-area sign-flip sensitivity; captures are technical grouping, not mouse counts")
  }))
  save_csv(cap_signflip, "11_122_capture_area_signflip_sensitivity_v02.csv")

  reg_summary <- do.call(rbind, lapply(SCENARIOS, function(sc) {
    z <- if (sc == SCENARIOS[1])
      region_spots else region_spots[!region_spots$inferred_uninjured, ]
    do.call(rbind, lapply(split(z, z$region), function(g) {
      h <- g[g$eligible & is.finite(g$rho), , drop = FALSE]
      data.frame(scenario = sc, region = g$region[1], n_portions = nrow(h), n_capture_areas = length(unique(h$sample_label)),
        total_spots = sum(h$n_spots), region_summary_eligible = nrow(h) >= MIN_SECTIONS_PER_REGION_SUMMARY,
        median_within_portion_rho = medfin(h$rho), median_spatial_detrended_rho = medfin(h$rho_spatial_detrended),
        n_positive = sum(h$rho > 0), proportion_positive = if (nrow(h))
          mean(h$rho > 0) else NA_real_)
    }))
  }))
  save_csv(reg_summary, "12_122_region_resolved_summary_across_portions_v02.csv")

  cap_region <- do.call(rbind, lapply(SCENARIOS, function(sc) {
    z <- if (sc == SCENARIOS[1])
      region_spots else region_spots[!region_spots$inferred_uninjured, ]
    z <- z[z$eligible & is.finite(z$rho), ]
    do.call(rbind, lapply(split(z, paste(z$sample_label, z$region, sep = "@@")), function(g) data.frame(scenario = sc,
      sample_label = g$sample_label[1], timepoint = g$timepoint[1], region = g$region[1], n_portions = nrow(g),
      n_spots = sum(g$n_spots), median_within_portion_rho = medfin(g$rho), median_spatial_detrended_rho = medfin(g$rho_spatial_detrended))))
  }))
  save_csv(cap_region, "13_122_capture_area_by_region_resolved_effects_v02.csv")

  null_portion <- vector("list", N_CONTROL_PROGRAMS)
  null_capture <- vector("list", N_CONTROL_PROGRAMS)
  null_overall <- vector("list", N_CONTROL_PROGRAMS)
  for (b in seq_len(N_CONTROL_PROGRAMS)) {
    sc <- score_panel(controls[[b]])
    sec <- calc_section_rhos(sc, spatial = FALSE)
    null_portion[[b]] <- data.frame(control_id = b, section_uid_internal = sec$section_uid_internal,
      sample_label = sec$sample_label, timepoint = sec$timepoint, inferred_uninjured = sec$inferred_uninjured,
      rho = sec$rho)
    cap <- do.call(rbind, lapply(SCENARIOS, function(scname) {
      zz <- if (scname == SCENARIOS[1])
        sec else sec[!sec$inferred_uninjured, ]
      do.call(rbind, lapply(split(zz, zz$sample_label), function(g) data.frame(control_id = b,
        scenario = scname, sample_label = g$sample_label[1], n_portions = nrow(g), median_direct_rho = medfin(g$rho),
        median_direct_fisher_z = medfin(fz(g$rho)))))
    }))
    null_capture[[b]] <- cap
    null_overall[[b]] <- do.call(rbind, lapply(SCENARIOS, function(scname) {
      z <- if (scname == SCENARIOS[1])
        sec else sec[!sec$inferred_uninjured, ]
      data.frame(control_id = b, scenario = scname, n_portions = nrow(z), median_direct_rho = medfin(z$rho),
        median_direct_fisher_z = medfin(fz(z$rho)), median_capture_fisher_z = medfin(cap$median_direct_fisher_z[cap$scenario ==
          scname]))
    }))
    if (b%%25L == 0L)
      step_log("  finished control ", b, "/", N_CONTROL_PROGRAMS)
  }
  null_portion <- do.call(rbind, null_portion)
  null_capture <- do.call(rbind, null_capture)
  null_overall <- do.call(rbind, null_overall)
  save_csv(null_portion, "14_122_matched_control_portion_correlations_v02.csv")
  save_csv(null_capture, "15_122_matched_control_capture_area_summaries_v02.csv")
  save_csv(null_overall, "16_122_matched_control_null_distribution_v02.csv")

  comparison <- do.call(rbind, lapply(SCENARIOS, function(sc) {
    z <- null_overall[null_overall$scenario == sc, ]
    oo <- obs_sec[if (sc == SCENARIOS[1])
      rep(TRUE, nrow(obs_sec)) else !obs_sec$inferred_uninjured, ]
    cobs <- cap_summary[cap_summary$scenario == sc, ]
    observed <- c(medfin(fz(oo$rho)), medfin(cobs$median_direct_fisher_z))
    names(observed) <- c("portion_median_fisher_z", "capture_median_fisher_z")
    do.call(rbind, lapply(names(observed), function(stat) {
      nul <- if (stat == "portion_median_fisher_z")
        z$median_direct_fisher_z else z$median_capture_fisher_z
      data.frame(scenario = sc, statistic = stat, n_control_programs = length(nul), observed_fisher_z = observed[[stat]],
        null_median_fisher_z = median(nul), null_p025 = unname(quantile(nul, 0.025)), null_p975 = unname(quantile(nul,
          0.975)), empirical_one_sided_P = (1 + sum(nul >= observed[[stat]]))/(1 + length(nul)),
        empirical_two_sided_P = (1 + sum(abs(nul - median(nul)) >= abs(observed[[stat]] - median(nul))))/(1 +
          length(nul)))
    }))
  }))
  save_csv(comparison, "17_122_matched_control_empirical_comparison_v02.csv")

  capture_null_compare <- do.call(rbind, lapply(SCENARIOS, function(sc) {
    cobs <- cap_summary[cap_summary$scenario == sc, ]
    do.call(rbind, lapply(seq_len(nrow(cobs)), function(k) {
      row <- cobs[k, ]
      z <- null_capture[null_capture$scenario == sc & null_capture$sample_label == row$sample_label,
        ]
      data.frame(scenario = sc, sample_label = row$sample_label, n_portions = row$n_portions, observed_median_fisher_z = row$median_direct_fisher_z,
        null_median_fisher_z = median(z$median_direct_fisher_z), null_p025 = unname(quantile(z$median_direct_fisher_z,
          0.025)), null_p975 = unname(quantile(z$median_direct_fisher_z, 0.975)), empirical_one_sided_P = (1 +
          sum(z$median_direct_fisher_z >= row$median_direct_fisher_z))/(1 + nrow(z)))
    }))
  }))
  save_csv(capture_null_compare, "18_122_capture_area_control_comparison_exploratory_v02.csv")

  qc <- data.frame(check = c("original_gene_count_37", "eight_capture_areas", "32_major_sections",
    "27_inferred_injured_sections", "frozen_score_reproduced", "all_control_panels_size37", "all_control_panels_unique",
    "rctd_overlap_pattern_matched", "all_spatial_joins_complete", "control_empirical_P_not_spot_level"),
    pass = c(length(orig) == 37L, length(sg) == 8L, nrow(section_info) == 32L, sum(!section_info$inferred_uninjured) ==
      27L, score_err <= MAX_FROZEN_SCORE_DIFFERENCE, all(control_qc$n_genes == 37L), all(control_qc$n_unique ==
      37L), all(control_qc$rctd_pattern_matches), !anyNA(pos) && !anyNA(j) && !anyNA(q117), TRUE))
  save_csv(qc, "19_122_FINAL_QC_v02.csv")
  if (!all(qc$pass))
    stop("One or more final QC checks failed; inspect 19_122_FINAL_QC_v02.csv")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    grDevices::pdf(file.path(FIG, "20_122_spatial_robustness_diagnostics_v02.pdf"), width = 10, height = 6.5,
      onefile = TRUE)
    p1 <- ggplot2::ggplot(obs_sec, ggplot2::aes(x = timepoint, y = rho, color = sample_label)) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0), size = 2) +
      ggplot2::facet_wrap(~inferred_uninjured, labeller = ggplot2::label_both) + ggplot2::geom_hline(yintercept = 0,
      linetype = 2) + ggplot2::labs(title = "Per-portion Clec7a–lipid correlation (putative mouse units)",
      y = "Within-portion Spearman rho") + ggplot2::theme_bw()
    print(p1)
    observed_vline_df <- do.call(rbind, lapply(SCENARIOS, function(s) {
      v <- if (s == SCENARIOS[1])
        obs_sec$rho else obs_sec$rho[!obs_sec$inferred_uninjured]
      data.frame(scenario = s, observed = medfin(fz(v)))
    }))
    p2 <- ggplot2::ggplot(null_overall, ggplot2::aes(median_direct_fisher_z)) + ggplot2::geom_histogram(bins = 35,
      fill = "grey70", color = "white") + ggplot2::geom_vline(data = observed_vline_df, ggplot2::aes(xintercept = observed),
      color = "#B45309", linewidth = 1) + ggplot2::facet_wrap(~scenario, scales = "free_y") + ggplot2::labs(title = "Expression-, size- and RCTD-pattern-matched negative controls",
      x = "Median per-portion Fisher z", y = "Control-program count") + ggplot2::theme_bw()
    print(p2)
    p3 <- ggplot2::ggplot(obs_sec, ggplot2::aes(x = rho, y = rho_spatial_detrended, color = timepoint)) +
      ggplot2::geom_point(size = 2) + ggplot2::geom_abline(linetype = 2) + ggplot2::geom_hline(yintercept = 0,
      linetype = 3) + ggplot2::labs(title = "Spatial detrending sensitivity within each portion",
      x = "Original Spearman rho", y = "Correlation of rank residuals after quadratic x/y trend") +
      ggplot2::theme_bw()
    print(p3)
    grDevices::dev.off()
  } else step_log("ggplot2 not installed: numerical tables generated; PDF figures skipped")

  sink_file <- file.path(OUT, "21_122_sessionInfo_v02.txt")
  capture.output(sessionInfo(), file = sink_file)
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "22_122_REPORT_v02.txt"), useBytes = TRUE)
}, error = function(e) {
  step_log("STEP122 FAILED: ", conditionMessage(e))
  capture.output(traceback(2), file = file.path(OUT, "99_122_ERROR_TRACE_v02.txt"))
  stop(e)
}, finally = {
  if (sink.number(type = "output") > 0L)
    sink(type = "output")
  try(close(logcon), silent = TRUE)
})

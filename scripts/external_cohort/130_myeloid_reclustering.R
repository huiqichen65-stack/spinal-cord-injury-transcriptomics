source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- PROJECT_ROOT
STEM <- "GSE281205"
DIR124 <- file.path(ROOT, "results", paste0("124_", STEM, "_lineage_review_myeloid_extraction"))
DIR128 <- file.path(ROOT, "results", paste0("128_", STEM, "_reviewed_lineage_finalization"))
OUT <- file.path(ROOT, "results", paste0("130_", STEM, "_myeloid_reclustering"))
FIG_DIR <- file.path(OUT, "figures")
RDS_DIR <- file.path(OUT, "RDS")
LOG_DIR <- file.path(OUT, "logs")

IN_RDS <- file.path(DIR124, "RDS", "124_GSE281205_extended_myeloid_DecontX_counts_v01.rds")
IN_IDS <- file.path(DIR128, "06_128_myeloid_include_cell_ids_v02.csv")
IN_CURATED <- file.path(DIR128, "02_128_curated_24cluster_annotation_v02.csv")
IN_FEATURES <- file.path(DIR124, "12_124_gene_feature_mapping_v01.csv")
OUT_RDS <- file.path(RDS_DIR, "130_GSE281205_myeloid_LogNorm_PCA_Harmony_clusters_v02.rds")
LOG_FILE <- file.path(LOG_DIR, "00_130_run_log_v02.log")

EXPECTED_GSM <- paste0("GSM", 8614662:8614667)
EXPECTED_CELLS <- 28716L
N_HVG <- 3000L
N_PCS <- 30L
K_NEIGHBORS <- 20L
MAIN_RES <- 0.1
HARMONY_RES <- c(0.1, 0.15, 0.2)
RANDOM_SEED <- 13001L
HARMONY_THETA <- 2
RESUME_FROM_SAVED_OBJECT <- TRUE
RDS_COMPRESS <- FALSE
RUN_PCA_CLUSTER_SENSITIVITY <- TRUE

required <- c("Seurat", "SeuratObject", "Matrix", "harmony", "ggplot2")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop(paste0("Missing R packages: ", paste(missing_packages, collapse = ", "),
  ". Install missing packages before running Step130."), call. = FALSE)
for (d in c(OUT, FIG_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
if (!file.exists(LOG_FILE)) file.create(LOG_FILE)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG_FILE, append = TRUE)
  invisible(msg)
}
assert <- function(ok, message) {
  if (!isTRUE(ok))
    stop(message, call. = FALSE)
}
read_tab <- function(path) {
  assert(file.exists(path), paste("Missing input file:", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}
write_tab <- function(x, basename) {
  assert(is.data.frame(x), paste("Not a data.frame:", basename))
  path <- file.path(OUT, basename)
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}
need_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  assert(!length(miss), paste(label, "missing columns:", paste(miss, collapse = ", ")))
}
median_safe <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x)))
    return(NA_real_)
  stats::median(x[is.finite(x)])
}
safe_plot <- function(p, filename, width = 11, height = 7, dpi = 150) {
  dest <- file.path(FIG_DIR, filename)
  tryCatch({
    ggplot2::ggsave(filename = dest, plot = p, width = width, height = height, units = "in", dpi = dpi,
      limitsize = FALSE)
  }, error = function(e) {
    log_msg("FIGURE WARNING (other outputs retained): ", filename, " | ", conditionMessage(e))
  })
  invisible(dest)
}
ari_from_labels <- function(a, b) {
  assert(length(a) == length(b), "ARI label vectors are different lengths")
  tab <- table(as.character(a), as.character(b))
  choose2 <- function(x) x * (x - 1)/2
  n <- sum(tab)
  if (n <= 1)
    return(NA_real_)
  same <- sum(choose2(as.vector(tab)))
  left <- sum(choose2(rowSums(tab)))
  right <- sum(choose2(colSums(tab)))
  whole <- choose2(n)
  expected <- left * right/whole
  limit <- (left + right)/2
  if (abs(limit - expected) < 1e-10)
    return(NA_real_)
  as.numeric((same - expected)/(limit - expected))
}

variance_explained <- function(embedding, group) {
  g <- as.character(group)
  assert(nrow(embedding) == length(g), "Embedding and group lengths differ")
  assert(!anyNA(g), "Group contains missing labels")
  vapply(seq_len(ncol(embedding)), function(j) {
    x <- as.numeric(embedding[, j])
    ss_total <- sum((x - mean(x))^2)
    if (ss_total == 0)
      return(NA_real_)
    gmean <- tapply(x, g, mean)
    ss_between <- sum(as.numeric(table(g)) * (as.numeric(gmean[names(table(g))]) - mean(x))^2)
    max(0, min(1, ss_between/ss_total))
  }, numeric(1))
}
normalize_dpi <- function(x) {
  ans <- suppressWarnings(as.integer(as.character(x)))
  assert(!anyNA(ans) && all(ans %in% c(14L, 28L)), "dpi metadata must be 14/28 (numeric or numeric character)")
  ans
}
validate_input <- function(seu, include, curated) {
  assert(inherits(seu, "Seurat"), "Step124 RDS is not a Seurat object")
  assert("DecontX" %in% SeuratObject::Assays(seu), "Step124 Seurat object lacks the DecontX assay")
  assert(ncol(seu) == EXPECTED_CELLS, paste("Step124 cell count is", ncol(seu), "instead of", EXPECTED_CELLS))
  need_cols(include, c("cell_id", "gsm", "dpi", "cluster_harmony_r04", "curated_major_lineage_128"),
    "Step128 included IDs")
  assert(nrow(include) == EXPECTED_CELLS && !anyNA(include$cell_id) && !anyDuplicated(include$cell_id) &&
    setequal(as.character(include$cell_id), colnames(seu)), "Step124 RDS cell IDs do not match Step128 INCLUDE exactly")
  assert(all(as.character(include$curated_major_lineage_128) %in% c("Microglia", "Monocyte_macrophage",
    "APC_dendritic")), "Step128 INCLUDE contains a non-myeloid lineage")
  assert(setequal(as.character(include$gsm), EXPECTED_GSM), "Step128 INCLUDE does not contain the expected six GSMs")
  assert(!anyNA(include$dpi), "Step128 INCLUDE has missing dpi")
  include$dpi <- normalize_dpi(include$dpi)
  for (gsm in EXPECTED_GSM) {
    observed <- unique(include$dpi[as.character(include$gsm) == gsm])
    intended <- if (gsm %in% EXPECTED_GSM[1:3])
      14L else 28L
    assert(length(observed) == 1 && observed == intended, paste("GSM/dpi mapping inconsistent for",
      gsm))
  }
  need_cols(curated, c("cluster", "curated_major_lineage", "myeloid_gate"), "Step128 curated clusters")
  expected_includes <- as.character(c(0, 1, 2, 3, 6, 18, 22))
  got_includes <- as.character(curated$cluster[curated$myeloid_gate == "INCLUDE"])
  assert(setequal(got_includes, expected_includes), "Step128 INCLUDE cluster map differs from the adopted seven clusters")
  inc <- include[match(colnames(seu), as.character(include$cell_id)), , drop = FALSE]
  assert(identical(as.character(inc$cell_id), colnames(seu)), "Could not align Step128 ID rows to Step124 Seurat cells")
  md <- seu[[]]
  need_cols(md, c("gsm", "dpi", "cluster_harmony_r04", "decontX_contamination_percent"), "Step124 Seurat metadata")
  assert(identical(as.character(md$gsm), as.character(inc$gsm)) && identical(normalize_dpi(md$dpi),
    inc$dpi) && identical(as.character(md$cluster_harmony_r04), as.character(inc$cluster_harmony_r04)),
    "Step124 metadata conflicts with Step128 INCLUDE GSM/stage/parent-cluster")
  inc
}

run <- function() {
  start <- Sys.time()
  for (p in c(IN_RDS, IN_IDS, IN_CURATED, IN_FEATURES)) {
    assert(file.exists(p), paste("Required input is missing:", p))
  }
  write_tab(data.frame(input_role = c("Step124_DecontX_myeloid_Seurat", "Step128_INCLUDES", "Step128_curated_clusters",
    "Step124_gene_ID_symbol_mapping"), path = c(IN_RDS, IN_IDS, IN_CURATED, IN_FEATURES), bytes = as.numeric(file.info(c(IN_RDS,
    IN_IDS, IN_CURATED, IN_FEATURES))$size), stringsAsFactors = FALSE), "01_130_input_provenance_v02.csv")
  include <- read_tab(IN_IDS)
  curated <- read_tab(IN_CURATED)
  features <- read_tab(IN_FEATURES)
  need_cols(features, c("gene_id", "gene_symbol"), "Step124 feature map")
  assert(!anyDuplicated(as.character(features$gene_id)), "Feature map has duplicate gene IDs; fix upstream before analysis")

  set.seed(RANDOM_SEED)
  resumed <- FALSE
  if (isTRUE(RESUME_FROM_SAVED_OBJECT) && file.exists(OUT_RDS) && is.finite(file.info(OUT_RDS)$size) &&
    file.info(OUT_RDS)$size > 1e+05) {
    candidate <- tryCatch(readRDS(OUT_RDS), error = function(e) {
      log_msg("Saved RDS read failed; will recompute: ", conditionMessage(e))
      NULL
    })
    if (!is.null(candidate) && inherits(candidate, "Seurat") && !is.null(candidate@misc$step130) &&
      identical(candidate@misc$step130$version, "v02") && identical(candidate@misc$step130$seed,
      RANDOM_SEED) && all(c("pca", "harmony", "umap_myeloid_pca", "umap_myeloid_harmony") %in%
      SeuratObject::Reductions(candidate)) && all(c("myeloid_harmony_r010", "myeloid_harmony_r015",
      "myeloid_harmony_r020", "myeloid_pca_r010") %in% colnames(candidate[[]]))) {
      seu <- candidate
      rm(candidate)
      validate_input(seu, include, curated)
      resumed <- TRUE
    } else {
      if (!is.null(candidate))
        rm(candidate)
    }
  }
  if (!resumed) {
    seu <- readRDS(IN_RDS)
    inc <- validate_input(seu, include, curated)
    md <- seu[[]]
    md$curated_major_lineage_128 <- as.character(inc$curated_major_lineage_128)
    md$myeloid_gate_128 <- "INCLUDE"
    md$input_parent_cluster_123 <- as.character(inc$cluster_harmony_r04)
    md$stage_label_130 <- paste0(inc$dpi, " dpi")
    md$high_contam_ge50_130 <- as.numeric(md$decontX_contamination_percent) >= 50
    md$note_28dpi_age_discordance_130 <- inc$dpi == 28L
    seu <- SeuratObject::AddMetaData(seu, metadata = md[, c("curated_major_lineage_128", "myeloid_gate_128",
      "input_parent_cluster_123", "stage_label_130", "high_contam_ge50_130", "note_28dpi_age_discordance_130"),
      drop = FALSE])
    rm(md, inc)
    SeuratObject::DefaultAssay(seu) <- "DecontX"
    seu <- Seurat::NormalizeData(seu, assay = "DecontX", normalization.method = "LogNormalize", scale.factor = 10000,
      verbose = TRUE)
    seu <- Seurat::FindVariableFeatures(seu, assay = "DecontX", selection.method = "vst", nfeatures = N_HVG,
      verbose = TRUE)
    hvgs <- SeuratObject::VariableFeatures(seu)
    assert(length(hvgs) >= N_PCS + 1L, "Too few variable genes to perform 30-PC PCA")
    seu <- Seurat::ScaleData(seu, assay = "DecontX", features = hvgs, verbose = TRUE)
    seu <- Seurat::RunPCA(seu, assay = "DecontX", features = hvgs, npcs = N_PCS, seed.use = RANDOM_SEED,
      verbose = TRUE)
    seu <- harmony::RunHarmony(seu, group.by.vars = "gsm", reduction.use = "pca", dims.use = seq_len(N_PCS),
      reduction.save = "harmony", project.dim = FALSE, theta = HARMONY_THETA, verbose = TRUE)
    assert("harmony" %in% SeuratObject::Reductions(seu), "Harmony reduction was not created")
    seu <- Seurat::RunUMAP(seu, reduction = "pca", dims = seq_len(N_PCS), reduction.name = "umap_myeloid_pca",
      reduction.key = "UMP_", n.neighbors = 30, min.dist = 0.3, seed.use = RANDOM_SEED, verbose = TRUE)
    seu <- Seurat::RunUMAP(seu, reduction = "harmony", dims = seq_len(N_PCS), reduction.name = "umap_myeloid_harmony",
      reduction.key = "UMH_", n.neighbors = 30, min.dist = 0.3, seed.use = RANDOM_SEED, verbose = TRUE)
    seu <- Seurat::FindNeighbors(seu, reduction = "harmony", dims = seq_len(N_PCS), k.param = K_NEIGHBORS,
      graph.name = c("myeloid_harmony_nn", "myeloid_harmony_snn"), verbose = TRUE)
    field_names <- c("myeloid_harmony_r010", "myeloid_harmony_r015", "myeloid_harmony_r020")
    for (i in seq_along(HARMONY_RES)) {
      seu <- Seurat::FindClusters(seu, graph.name = "myeloid_harmony_snn", resolution = HARMONY_RES[i],
        cluster.name = field_names[i], algorithm = 1, random.seed = RANDOM_SEED, verbose = TRUE)
      assert(field_names[i] %in% colnames(seu[[]]), paste("Missing expected cluster field", field_names[i]))
    }
    if (isTRUE(RUN_PCA_CLUSTER_SENSITIVITY)) {
      seu <- Seurat::FindNeighbors(seu, reduction = "pca", dims = seq_len(N_PCS), k.param = K_NEIGHBORS,
        graph.name = c("myeloid_pca_nn", "myeloid_pca_snn"), verbose = TRUE)
      seu <- Seurat::FindClusters(seu, graph.name = "myeloid_pca_snn", resolution = MAIN_RES, cluster.name = "myeloid_pca_r010",
        algorithm = 1, random.seed = RANDOM_SEED, verbose = TRUE)
    } else {
      seu$myeloid_pca_r010 <- NA_character_
    }
    SeuratObject::Idents(seu) <- "myeloid_harmony_r010"
    seu@misc$step130 <- list(version = "v02", seed = RANDOM_SEED, input_124 = IN_RDS, include_128 = IN_IDS,
      original_integer_UMI = "Step119 singlet_raw_counts_v02.rds (not stored here)", corrected_estimated_expression = TRUE,
      signature_1320_loaded = FALSE, target_state_annotation_performed = FALSE, r010_main_cluster = "myeloid_harmony_r010",
      batch_group = "gsm", stage_batch_confounded = TRUE, age_28dpi_metadata_discordance = TRUE)

    tmp <- paste0(OUT_RDS, ".partial")
    if (file.exists(tmp))
      unlink(tmp)
    saveRDS(seu, file = tmp, compress = RDS_COMPRESS)
    assert(file.exists(tmp) && file.info(tmp)$size > 1e+05, "Saved RDS is missing or unexpectedly small")
    if (file.exists(OUT_RDS))
      unlink(OUT_RDS)
    assert(isTRUE(file.rename(tmp, OUT_RDS)), "Cannot finalize saved Seurat RDS; .partial file is retained")
  }

  md <- seu[[]]
  idx <- match(colnames(seu), as.character(include$cell_id))
  assert(!anyNA(idx), "Input manifest no longer aligns with the saved Seurat object")
  if (!"curated_major_lineage_128" %in% names(md)) {
    md$curated_major_lineage_128 <- include$curated_major_lineage_128[idx]
  }
  assert(nrow(md) == EXPECTED_CELLS && identical(as.character(md$gsm), as.character(include$gsm[idx])) &&
    identical(normalize_dpi(md$dpi), normalize_dpi(include$dpi[idx])), "Saved Seurat cell metadata inconsistent with Step128 manifest")
  fields <- c("myeloid_harmony_r010", "myeloid_harmony_r015", "myeloid_harmony_r020", "myeloid_pca_r010")
  assert(all(fields %in% colnames(md)), "Expected cluster fields missing")
  assert(!anyNA(md[, fields[1:3], drop = FALSE]), "Harmony cluster fields have missing values")
  hvg <- SeuratObject::VariableFeatures(seu)
  write_tab(data.frame(rank = seq_along(hvg), gene_id = hvg, gene_symbol = features$gene_symbol[match(hvg,
    features$gene_id)]), "02_130_HVG_3000_v02.csv")

  pca_std <- seu[["pca"]]@stdev
  if (length(pca_std) < N_PCS)
    stop("Seurat PCA object has <30 PCs", call. = FALSE)
  write_tab(data.frame(PC = seq_len(N_PCS), stddev = pca_std[seq_len(N_PCS)], variance = pca_std[seq_len(N_PCS)]^2,
    pct_among_first_30 = 100 * pca_std[seq_len(N_PCS)]^2/sum(pca_std[seq_len(N_PCS)]^2)), "03_130_PCA_variance_v02.csv")
  cl_h10 <- as.character(md$myeloid_harmony_r010)
  cl_h15 <- as.character(md$myeloid_harmony_r015)
  cl_h20 <- as.character(md$myeloid_harmony_r020)
  cl_p10 <- as.character(md$myeloid_pca_r010)
  distinct_clusters <- function(x, source) {
    tab <- table(as.character(x))
    data.frame(method = source, cluster = names(tab), n_cells = as.integer(tab), stringsAsFactors = FALSE)
  }
  write_tab(do.call(rbind, list(distinct_clusters(cl_h10, "Harmony_r0.10"), distinct_clusters(cl_h15,
    "Harmony_r0.15"), distinct_clusters(cl_h20, "Harmony_r0.20"), distinct_clusters(cl_p10, "PCA_r0.10"))),
    "04_130_cluster_counts_by_resolution_v02.csv")

  membership <- data.frame(cell_id = rownames(md), gsm = as.character(md$gsm), dpi = normalize_dpi(md$dpi),
    replicate = as.integer(md$replicate), parent_cluster_123 = as.character(md$cluster_harmony_r04),
    curated_major_lineage_128 = as.character(md$curated_major_lineage_128), myeloid_harmony_r010 = cl_h10,
    myeloid_harmony_r015 = cl_h15, myeloid_harmony_r020 = cl_h20, myeloid_pca_r010 = cl_p10, contamination_percent = as.numeric(md$decontX_contamination_percent),
    high_contam_ge50 = as.numeric(md$decontX_contamination_percent) >= 50, stringsAsFactors = FALSE)
  assert(!anyDuplicated(membership$cell_id) && identical(as.character(membership$cell_id), colnames(seu)),
    "Per-cell manifest not aligned with Seurat object")
  write_tab(membership, "05_130_per_cell_myeloid_cluster_membership_v02.csv")

  cluster_levels <- as.character(sort(as.integer(unique(cl_h10))))
  gsm_cluster <- as.data.frame(table(gsm = factor(membership$gsm, levels = EXPECTED_GSM), cluster = factor(cl_h10,
    levels = cluster_levels)), stringsAsFactors = FALSE)
  names(gsm_cluster)[3] <- "n_cells"
  gsm_cluster$dpi <- ifelse(gsm_cluster$gsm %in% EXPECTED_GSM[1:3], 14L, 28L)
  gsm_n <- as.integer(table(factor(membership$gsm, levels = EXPECTED_GSM)))
  gsm_cluster$denominator_GSM <- gsm_n[match(as.character(gsm_cluster$gsm), EXPECTED_GSM)]
  gsm_cluster$fraction_within_GSM <- gsm_cluster$n_cells/gsm_cluster$denominator_GSM
  write_tab(gsm_cluster[, c("gsm", "dpi", "cluster", "n_cells", "denominator_GSM", "fraction_within_GSM")],
    "06_130_sample_by_myeloid_cluster_v02.csv")
  stage_cluster <- stats::aggregate(n_cells ~ dpi + cluster, data = gsm_cluster, FUN = sum)
  denom <- tapply(gsm_cluster$denominator_GSM[!duplicated(gsm_cluster$gsm)], as.character(gsm_cluster$dpi[!duplicated(gsm_cluster$gsm)]),
    sum)
  stage_cluster$denominator_stage <- as.integer(denom[as.character(stage_cluster$dpi)])
  stage_cluster$fraction_within_stage <- stage_cluster$n_cells/stage_cluster$denominator_stage
  write_tab(stage_cluster, "07_130_stage_myeloid_cluster_fractions_v02.csv")

  qc <- do.call(rbind, lapply(as.character(sort(as.integer(unique(cl_h10)))), function(k) {
    sub <- md[cl_h10 == k, , drop = FALSE]
    data.frame(cluster = k, n_cells = nrow(sub), n_libraries = length(unique(as.character(sub$gsm))),
      n_14dpi = sum(normalize_dpi(sub$dpi) == 14L), n_28dpi = sum(normalize_dpi(sub$dpi) == 28L),
      median_DecontX_contamination_percent = median_safe(sub$decontX_contamination_percent), pct_DecontX_contamination_ge50 = 100 *
        mean(as.numeric(sub$decontX_contamination_percent) >= 50), median_raw_nFeature_RNA = if ("nFeature_RNA" %in%
        names(sub))
        median_safe(sub$nFeature_RNA) else NA_real_, median_raw_nCount_RNA = if ("nCount_RNA" %in% names(sub))
        median_safe(sub$nCount_RNA) else NA_real_, stringsAsFactors = FALSE)
  }))
  write_tab(qc, "08_130_cluster_QC_and_contamination_v02.csv")

  canonical <- list(microglia_identity = c("P2ry12", "Tmem119", "Sall1", "Fcrls", "Hexb", "Cx3cr1"),
    injury_activation_lipid = c("Clec7a", "Apoe", "Lpl", "Trem2", "Gpnmb", "Spp1", "Cst7", "Itgax"),
    monocytes = c("Ccr2", "Ly6c2", "Lyz2", "Ms4a7", "S100a8", "S100a9"), macrophage_related = c("Mrc1",
      "Pf4", "F13a1", "Cd68", "C1qa"), APC_DC = c("Flt3", "Cd209a", "Cd74", "H2-Aa", "H2-Ab1"),
    interferon = c("Oasl1", "Ifit2", "Isg15", "Cxcl10"), proliferation = c("Mki67", "Pclaf", "Birc5"),
    non_myeloid_exclusion = c("Pecam1", "Cldn5", "Rgs5", "Pdgfrb", "Aqp4", "Mbp", "Plp1", "Cd3d",
      "Nkg7", "Col1a1"))
  symbols <- unique(unlist(canonical, use.names = FALSE))
  mapping <- do.call(rbind, lapply(symbols, function(sym) {
    ix <- which(as.character(features$gene_symbol) == sym)
    ix <- ix[features$gene_id[ix] %in% rownames(seu)]
    if (!length(ix))
      return(data.frame(gene_symbol = sym, gene_id = NA_character_, is_present = FALSE, n_geneID_matches = 0L))
    data.frame(gene_symbol = sym, gene_id = as.character(features$gene_id[ix[1]]), is_present = TRUE,
      n_geneID_matches = length(ix))
  }))
  mapping$marker_category <- vapply(mapping$gene_symbol, function(s) {
    names(canonical)[which(vapply(canonical, function(v) s %in% v, logical(1)))[1]]
  }, character(1))
  write_tab(mapping, "09_130_canonical_marker_mapping_v02.csv")
  if (sum(mapping$is_present) < 25) {
    log_msg("WARNING: only ", sum(mapping$is_present), " canonical marker genes available")
  }
  present <- mapping[mapping$is_present, , drop = FALSE]
  SeuratObject::DefaultAssay(seu) <- "DecontX"
  cnt <- SeuratObject::LayerData(seu[["DecontX"]], layer = "counts")
  dat <- SeuratObject::LayerData(seu[["DecontX"]], layer = "data")
  assert(identical(colnames(cnt), colnames(seu)) && identical(colnames(dat), colnames(seu)), "Expression matrix columns do not match Seurat cells")
  marker_rows <- vector("list", length(unique(cl_h10)))
  for (kidx in seq_along(marker_rows)) {
    k <- as.character(sort(as.integer(unique(cl_h10))))[kidx]
    ci <- which(cl_h10 == k)
    count_sub <- cnt[present$gene_id, ci, drop = FALSE]
    norm_sub <- dat[present$gene_id, ci, drop = FALSE]
    marker_rows[[kidx]] <- data.frame(cluster = k, n_cells = length(ci), gene_symbol = present$gene_symbol,
      gene_id = present$gene_id, marker_category = present$marker_category, pct_cells_corrected_counts_ge1 = 100 *
        as.numeric(Matrix::rowMeans(count_sub >= 1)), pct_cells_corrected_counts_ge5 = 100 *
        as.numeric(Matrix::rowMeans(count_sub >= 5)), mean_log_normalized_expr = as.numeric(Matrix::rowMeans(norm_sub)),
      stringsAsFactors = FALSE)
  }
  mark_table <- do.call(rbind, marker_rows)
  write_tab(mark_table, "10_130_canonical_marker_evidence_v02.csv")

  hvdata <- dat[hvg, , drop = FALSE]
  total_expr <- Matrix::rowSums(hvdata)
  exploratory <- vector("list", length(unique(cl_h10)))
  for (kidx in seq_along(exploratory)) {
    k <- as.character(sort(as.integer(unique(cl_h10))))[kidx]
    ii <- which(cl_h10 == k)
    assert(length(ii) < ncol(hvdata), "Cannot compute rest-of-myeloid mean")
    inside <- Matrix::rowSums(hvdata[, ii, drop = FALSE])
    difference <- inside/length(ii) - (total_expr - inside)/(ncol(hvdata) - length(ii))
    order_idx <- head(order(difference, decreasing = TRUE, na.last = NA), 20)
    exploratory[[kidx]] <- data.frame(cluster = k, rank = seq_along(order_idx), gene_id = rownames(hvdata)[order_idx],
      gene_symbol = features$gene_symbol[match(rownames(hvdata)[order_idx], features$gene_id)],
      mean_logexpr_difference_vs_other_myeloid = as.numeric(difference[order_idx]), note = "Descriptive HVG mean difference; NOT specimen-level DEG / signature validation",
      stringsAsFactors = FALSE)
  }
  write_tab(do.call(rbind, exploratory), "11_130_exploratory_HVG_cluster_markers_v02.csv")
  rm(cnt, dat, hvdata)
  invisible(gc(verbose = FALSE))

  ari <- data.frame(comparison = c("Harmony_0.10_vs_0.15", "Harmony_0.10_vs_0.20", "Harmony_0.10_vs_PCA_0.10"),
    ARI = c(ari_from_labels(cl_h10, cl_h15), ari_from_labels(cl_h10, cl_h20), ari_from_labels(cl_h10,
      cl_p10)), stringsAsFactors = FALSE)
  write_tab(ari, "12_130_cluster_stability_ARI_v02.csv")

  pca_embedding <- SeuratObject::Embeddings(seu[["pca"]])[, seq_len(N_PCS), drop = FALSE]
  harmony_embedding <- SeuratObject::Embeddings(seu[["harmony"]])[, seq_len(N_PCS), drop = FALSE]
  assert(identical(rownames(pca_embedding), rownames(md)) && identical(rownames(harmony_embedding),
    rownames(md)), "Reduction cell order differs from metadata")
  sep_for <- function(emb, method) {
    dd <- data.frame(dim = seq_len(N_PCS), r2_GSM = variance_explained(emb, md$gsm), r2_stage = variance_explained(emb,
      normalize_dpi(md$dpi)))
    dd$method <- method
    dd
  }
  batch <- rbind(sep_for(pca_embedding, "PCA_uncorrected"), sep_for(harmony_embedding, "Harmony_GSM_corrected"))
  write_tab(batch[, c("method", "dim", "r2_GSM", "r2_stage")], "13_130_GSM_stage_descriptive_variance_v02.csv")

  per_gsm <- as.data.frame(table(gsm = factor(as.character(md$gsm), levels = EXPECTED_GSM)), stringsAsFactors = FALSE)
  names(per_gsm)[2] <- "n_myeloid_cells"
  per_gsm$dpi <- ifelse(as.character(per_gsm$gsm) %in% EXPECTED_GSM[1:3], 14L, 28L)
  per_gsm$n_harmony_clusters_r010_present <- vapply(as.character(per_gsm$gsm), function(g) length(unique(cl_h10[membership$gsm ==
    g])), integer(1))
  per_gsm$n_high_contamination_ge50 <- vapply(as.character(per_gsm$gsm), function(g) sum(membership$high_contam_ge50[membership$gsm ==
    g]), integer(1))
  per_gsm$metadata_note <- ifelse(per_gsm$dpi == 28L, "GEO Young/Old treatment field conflict; BioSample supports Young; preserve warning",
    "14dpi metadata consistent")
  write_tab(per_gsm, "14_130_per_GSM_myeloid_feasibility_v02.csv")

  e_h <- SeuratObject::Embeddings(seu[["umap_myeloid_harmony"]])
  e_p <- SeuratObject::Embeddings(seu[["umap_myeloid_pca"]])
  df_h <- data.frame(x = e_h[, 1], y = e_h[, 2], membership, stringsAsFactors = FALSE)
  df_p <- data.frame(x = e_p[, 1], y = e_p[, 2], membership, stringsAsFactors = FALSE)
  plot_umap <- function(df, attr, title, subtitle = "") {
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = .data[[attr]])) + ggplot2::geom_point(size = 0.22,
      alpha = 0.52, stroke = 0) + ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank()) +
      ggplot2::labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2", color = attr) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
  }
  safe_plot(plot_umap(df_h, "myeloid_harmony_r010", "GSE281205 | myeloid Harmony clusters (r=0.10)",
    "Working myeloid gate; not final cell-state annotation"), "Fig_130_Harmony_myeloid_clusters_r010_v02.png")
  safe_plot(plot_umap(df_h, "myeloid_harmony_r015", "Myeloid Harmony clusters (r=0.15)"), "Fig_130_Harmony_myeloid_clusters_r015_v02.png")
  safe_plot(plot_umap(df_h, "myeloid_harmony_r020", "Myeloid Harmony clusters (r=0.20)"), "Fig_130_Harmony_myeloid_clusters_r020_v02.png")
  safe_plot(plot_umap(df_p, "myeloid_pca_r010", "Uncorrected PCA-based myeloid clusters (r=0.10)",
    "Sensitivity check, same retained cells"), "Fig_130_PCA_myeloid_clusters_r010_v02.png")
  safe_plot(plot_umap(df_h, "gsm", "Myeloid Harmony UMAP by GSM"), "Fig_130_Harmony_myeloid_by_GSM_v02.png")
  df_h$stage <- factor(paste0(df_h$dpi, " dpi"), levels = c("14 dpi", "28 dpi"))
  safe_plot(plot_umap(df_h, "stage", "Myeloid Harmony UMAP by injury stage"), "Fig_130_Harmony_myeloid_by_stage_v02.png")
  safe_plot(plot_umap(df_h, "curated_major_lineage_128", "Myeloid Harmony UMAP by parent broad lineage"),
    "Fig_130_Harmony_myeloid_by_parent_lineage_v02.png")
  df_h$high_contamination <- ifelse(df_h$high_contam_ge50, ">=50% contamination estimate", "<50% contamination estimate")
  safe_plot(plot_umap(df_h, "high_contamination", "High-DecontX-contamination candidate cells"), "Fig_130_Harmony_myeloid_high_contamination_v02.png")
  plot_stage <- ggplot2::ggplot(stage_cluster, ggplot2::aes(x = factor(dpi), y = fraction_within_stage,
    fill = factor(cluster))) + ggplot2::geom_col() + ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 *
    v), "%")) + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = "Working myeloid cluster composition by injury stage",
    subtitle = "Cell fractions descriptive only; library-level replication not yet modeled", x = "Days post injury",
    y = "Fraction within sampled myeloid", fill = "Cluster")
  safe_plot(plot_stage, "Fig_130_myeloid_cluster_stage_composition_v02.png", 9, 6)
  if (nrow(mark_table) > 0) {
    plot_mark <- ggplot2::ggplot(mark_table, ggplot2::aes(x = factor(cluster), y = gene_symbol, color = mean_log_normalized_expr,
      size = pct_cells_corrected_counts_ge1)) + ggplot2::geom_point() + ggplot2::facet_grid(marker_category ~
      ., scales = "free_y", space = "free_y") + ggplot2::scale_color_viridis_c() + ggplot2::theme_bw(base_size = 10) +
      ggplot2::labs(title = "Canonical markers across myeloid subclusters", subtitle = "DecontX-corrected estimated expression; detection >=1 estimated count",
        x = "Myeloid Harmony r0.10 cluster", y = NULL, color = "Mean log-normalized", size = "% >=1") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    safe_plot(plot_mark, "Fig_130_myeloid_canonical_marker_dotplot_v02.png", 12, 13)
  }
  plot_batch <- stats::aggregate(cbind(r2_GSM, r2_stage) ~ method, data = batch, FUN = mean)
  plot_batch_long <- rbind(data.frame(method = plot_batch$method, component = "GSM", r2 = plot_batch$r2_GSM),
    data.frame(method = plot_batch$method, component = "dpi", r2 = plot_batch$r2_stage))
  safe_plot(ggplot2::ggplot(plot_batch_long, ggplot2::aes(x = method, y = r2, fill = component)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = "PCA / Harmony sample and stage separation",
    subtitle = "Mean coordinate-level descriptive R2; GSM nested in dpi", x = NULL, y = "Mean descriptive R2",
    fill = "Grouping"), "Fig_130_PCA_Harmony_GSM_stage_variance_v02.png", 9, 6)

  cl_counts <- vapply(fields, function(f) length(unique(as.character(md[[f]]))), integer(1))
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "15_130_REPORT_v02.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(OUT, "16_130_sessionInfo_v02.txt"), useBytes = TRUE)
  output_files <- list.files(OUT, full.names = TRUE, recursive = TRUE, include.dirs = FALSE)
  output_index <- data.frame(relative_path = substring(normalizePath(output_files, winslash = "/",
    mustWork = FALSE), nchar(normalizePath(OUT, winslash = "/", mustWork = FALSE)) + 2L), bytes = as.numeric(file.info(output_files)$size))
  write_tab(output_index, "17_130_OUTPUT_INDEX_v02.csv")
  invisible(TRUE)
}

tryCatch(run(), error = function(e) {
  log_msg("ERROR | Step130 v02 | ", conditionMessage(e))
  stop(e)
})

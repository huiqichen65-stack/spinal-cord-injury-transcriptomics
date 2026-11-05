options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

ROOT <- "D:/zz"
IN_DIR <- file.path(ROOT, "results", "119_GSE281205_scDblFinder")
IN_RDS <- file.path(IN_DIR, "RDS")
OUT_DIR <- file.path(ROOT, "results", "121_GSE281205_normalize_PCA_batch_audit")
RDS_DIR <- file.path(OUT_DIR, "RDS")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
VERSION <- "v02"
GSM_IDS <- paste0("GSM", 8614662:8614667)
EXPECTED_DPI <- rep(c(14L, 28L), each = 3L)
EXPECTED_REP <- rep(1:3, times = 2L)
SEED <- 12101L
N_VARIABLE_FEATURES <- 3000L
N_PCS <- 30L
KNN_K <- 30L
RUN_UMAP <- TRUE
RUN_KNN_SAMPLE_MIXING <- TRUE
HARMONY_GROUP <- "gsm"
HARMONY_THETA <- 2
SAVE_SEURAT_RDS <- TRUE

assert <- function(cond, msg) {
  if (!isTRUE(cond))
    stop(msg, call. = FALSE)
}
required <- c("Matrix", "Seurat", "SeuratObject", "ggplot2", "harmony")
if (RUN_KNN_SAMPLE_MIXING) required <- c(required, "RANN")
for (pkg in required) {
  assert(requireNamespace(pkg, quietly = TRUE), paste0("Missing R package '", pkg, "'. Install it before running Step 121."))
}
assert(dir.exists(IN_RDS), paste("Missing Step 119 RDS directory:", IN_RDS))
for (p in c(OUT_DIR, RDS_DIR, FIG_DIR, LOG_DIR)) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
}

LOG_FILE <- file.path(OUT_DIR, "00_121_run_log_v02.log")
cat("", file = LOG_FILE)
log_line <- function(...) {
  x <- paste0(..., collapse = "")
  cat(x, "\n", sep = "")
  cat(x, "\n", sep = "", file = LOG_FILE, append = TRUE)
  invisible(x)
}
write_csv <- function(x, fname) {
  assert(is.data.frame(x), paste("Output is not a data.frame:", fname))
  utils::write.csv(x, file.path(OUT_DIR, fname), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
med <- function(x) {
  if (length(x) == 0L || all(is.na(x)))
    return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}
q95 <- function(x) {
  if (length(x) == 0L || all(is.na(x)))
    return(NA_real_)
  as.numeric(stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE))
}

group_r2 <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- as.numeric(x[ok])
  g <- as.character(g[ok])
  if (length(x) < 2L || length(unique(g)) < 2L)
    return(NA_real_)
  grand <- mean(x)
  denom <- sum((x - grand)^2)
  if (!is.finite(denom) || denom <= 0)
    return(NA_real_)
  means <- tapply(x, g, mean)
  ns <- table(g)
  between <- sum(as.numeric(ns[names(means)]) * (as.numeric(means) - grand)^2)
  pmin(1, pmax(0, between/denom))
}

read_119_rds <- function(file, gsm, corrected = FALSE) {
  assert(file.exists(file), paste("Input RDS missing:", file))
  x <- readRDS(file)
  name <- if (corrected)
    "counts_decontX" else "counts"
  needed <- c(name, "cell_metadata", "feature_metadata", "provenance")
  assert(is.list(x) && all(needed %in% names(x)), paste("Unexpected 119 RDS structure:", file))
  m <- x[[name]]
  assert(inherits(m, "sparseMatrix"), paste("Not sparse:", file))
  m <- as(m, "dgCMatrix")
  md <- x$cell_metadata
  fm <- x$feature_metadata
  assert(is.data.frame(md) && is.data.frame(fm), paste("Metadata missing:", file))
  assert(all(c("cell_id", "gsm", "dpi", "replicate", "age_assignment_status", "nCount_RNA", "nFeature_RNA",
    "percent.mt", "decontX_contamination_percent", "scDblFinder_class") %in% names(md)), paste("Mandatory cell metadata missing:",
    file))
  assert(all(c("gene_id", "gene_symbol") %in% names(fm)), paste("Gene mapping incomplete:", file))
  assert(nrow(m) == nrow(fm) && ncol(m) == nrow(md), paste("Matrix/metadata dimension mismatch:", file))
  assert(identical(rownames(m), as.character(fm$gene_id)) && identical(colnames(m), as.character(md$cell_id)),
    paste("Matrix row/column IDs not aligned:", file))
  assert(!anyDuplicated(rownames(m)) && !anyDuplicated(colnames(m)), paste("Duplicate gene IDs or cell IDs:",
    file))
  assert(all(as.character(md$gsm) == gsm) && all(as.character(md$scDblFinder_class) == "singlet"),
    paste("Wrong GSM or non-singlets found:", file))
  assert(all(is.finite(m@x)) && all(m@x >= 0), paste("Invalid expression matrix:", file))
  if (!corrected) {
    assert(all(abs(m@x - round(m@x)) < 1e-08), paste("Raw UMI matrix is not integer-valued:", file))
  }
  x[[name]] <- m
  x
}

plot_embedding <- function(dat, x, y, colour, path, title_text, colour_type = c("discrete", "numeric")) {
  colour_type <- match.arg(colour_type)
  pp <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[colour]])) +
    ggplot2::geom_point(size = 0.21, alpha = 0.38) + ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    legend.position = "right") + ggplot2::labs(title = title_text, x = x, y = y, subtitle = "Exploratory cell-level embedding; not independent biological replicates")
  if (colour_type == "discrete") {
    pp <- pp + ggplot2::scale_color_discrete(drop = FALSE)
  } else {
    pp <- pp + ggplot2::scale_color_gradient(low = "gray75", high = "firebrick3")
  }
  ggplot2::ggsave(filename = path, plot = pp, width = 10.5, height = 7.3, dpi = 170, limitsize = FALSE)
  invisible(path)
}

audit_mixing <- function(pc, metadata, k = 30L) {
  out <- list()
  for (stage in c(14L, 28L)) {
    cells <- which(as.integer(metadata$dpi) == stage)
    X <- pc[cells, , drop = FALSE]
    gsm <- as.character(metadata$gsm[cells])
    n <- nrow(X)
    assert(n > k + 1L, paste("Too few cells for kNN mixing:", stage))
    nn <- RANN::nn2(data = X, query = X, k = k + 1L, searchtype = "standard")$nn.idx
    frac <- numeric(n)
    for (i in seq_len(n)) {
      ix <- nn[i, ]
      ix <- ix[ix != i]
      ix <- ix[seq_len(min(k, length(ix)))]
      frac[[i]] <- if (length(ix))
        mean(gsm[ix] == gsm[[i]]) else NA_real_
    }
    ttab <- table(gsm)
    for (key in GSM_IDS[GSM_IDS %in% names(ttab)]) {
      j <- which(gsm == key)
      obs <- mean(frac[j], na.rm = TRUE)
      expected <- (as.numeric(ttab[key]) - 1)/(n - 1)
      out[[length(out) + 1L]] <- data.frame(dpi = stage, gsm = key, n_cells = length(j), knn_k = k,
        observed_same_gsm_percent = 100 * obs, random_label_expected_same_gsm_percent = 100 *
          expected, observed_minus_expected_pp = 100 * (obs - expected), interpretation = "Descriptive only; biology, technical variation and sampling can all affect mixing",
        stringsAsFactors = FALSE)
    }
    rm(X, nn)
    invisible(gc(verbose = FALSE))
  }
  do.call(rbind, out)
}

main <- function() {
  log_line("R version: ", getRversion())
  log_line("Seurat version: ", utils::packageVersion("Seurat"))
  log_line("Harmony version: ", utils::packageVersion("harmony"))

  set.seed(SEED)
  sum119_file <- file.path(IN_DIR, "03_119_per_sample_doublet_summary_v02.csv")
  assert(file.exists(sum119_file), paste("Missing 119 per-sample summary:", sum119_file))
  s119 <- utils::read.csv(sum119_file, stringsAsFactors = FALSE)
  assert(all(c("gsm", "dpi", "n_predicted_singlets") %in% names(s119)), "119 summary lacks expected columns")
  assert(nrow(s119) == length(GSM_IDS) && setequal(as.character(s119$gsm), GSM_IDS), "119 summary GSMs do not match all six expected samples")

  corrected_by_gsm <- vector("list", length(GSM_IDS))
  md_by_gsm <- vector("list", length(GSM_IDS))
  reference_fm <- NULL
  manifests <- vector("list", length(GSM_IDS))

  for (i in seq_along(GSM_IDS)) {
    gsm <- GSM_IDS[[i]]
    dpi <- EXPECTED_DPI[[i]]
    raw_path <- file.path(IN_RDS, paste0("119_", gsm, "_singlet_raw_counts_v02.rds"))
    corrected_path <- file.path(IN_RDS, paste0("119_", gsm, "_singlet_DecontX_counts_v02.rds"))
    raw_obj <- read_119_rds(raw_path, gsm, corrected = FALSE)
    cor_obj <- read_119_rds(corrected_path, gsm, corrected = TRUE)
    rawx <- raw_obj$counts
    corx <- cor_obj$counts_decontX
    mdr <- raw_obj$cell_metadata
    mdc <- cor_obj$cell_metadata
    fm <- cor_obj$feature_metadata

    assert(identical(dimnames(rawx), dimnames(corx)), paste("Raw and corrected gene/cell IDs differ:",
      gsm))
    assert(identical(as.character(mdr$cell_id), as.character(mdc$cell_id)) && identical(as.character(raw_obj$feature_metadata$gene_id),
      as.character(fm$gene_id)), paste("Raw/corrected RDS annotation mismatch:", gsm))
    assert(length(unique(mdc$dpi)) == 1L && as.integer(unique(mdc$dpi)) == dpi, paste("Wrong injury day:",
      gsm))
    assert(length(unique(mdc$replicate)) == 1L && as.integer(unique(mdc$replicate)) == EXPECTED_REP[[i]],
      paste("Wrong replicate label:", gsm))
    assert(length(unique(mdc$age_assignment_status)) == 1L, paste("Age metadata mismatch within GSM:",
      gsm))
    assert(ncol(corx) == as.integer(s119$n_predicted_singlets[match(gsm, s119$gsm)]), paste("Singlet count disagrees with 119 summary:",
      gsm))
    assert(all(is.finite(mdc$decontX_contamination_percent)) && all(mdc$decontX_contamination_percent >=
      0 & mdc$decontX_contamination_percent <= 100), paste("Invalid DecontX contamination field:",
      gsm))
    assert(!anyNA(mdc$cell_id) && !anyDuplicated(mdc$cell_id), paste("Invalid cell IDs:", gsm))
    if (is.null(reference_fm)) {
      reference_fm <- fm
    } else {
      assert(identical(as.character(reference_fm$gene_id), as.character(fm$gene_id)) && identical(as.character(reference_fm$gene_symbol),
        as.character(fm$gene_symbol)), paste("Gene IDs or symbols not aligned between GSMs:",
        gsm))
    }
    if (dpi == 28L) {
      log_line("WARNING: ", gsm, " GEO title/treatment young-old conflict; age status retained: ",
        unique(mdc$age_assignment_status))
    }
    manifests[[i]] <- data.frame(gsm = gsm, dpi = dpi, replicate = EXPECTED_REP[[i]], age_assignment_status = as.character(unique(mdc$age_assignment_status)),
      singlets = ncol(corx), genes = nrow(corx), raw_umi_total = sum(rawx@x), corrected_estimated_umi_total = sum(corx@x),
      n_contamination_ge50 = sum(mdc$decontX_contamination_percent >= 50), source_raw_RDS = raw_path,
      source_corrected_RDS = corrected_path, stringsAsFactors = FALSE)
    corrected_by_gsm[[i]] <- corx
    md_by_gsm[[i]] <- mdc
    rm(raw_obj, cor_obj, rawx, corx, mdr, mdc, fm)
    invisible(gc(verbose = FALSE))
  }
  n_expected <- sum(vapply(corrected_by_gsm, ncol, integer(1)))
  assert(n_expected > 30000L, "Unexpectedly few post-doublet cells")
  corrected_mat <- do.call(cbind, corrected_by_gsm)
  md <- do.call(rbind, md_by_gsm)
  rownames(md) <- as.character(md$cell_id)
  assert(identical(colnames(corrected_mat), rownames(md)), "Combined corrected count columns and cell metadata differ")
  assert(identical(rownames(corrected_mat), as.character(reference_fm$gene_id)), "Combined gene ID reference mismatch")
  assert(!anyDuplicated(colnames(corrected_mat)), "Duplicated cell IDs in combined object")
  assert(sum(table(md$gsm)) == n_expected, "Combined metadata has wrong size")
  rm(corrected_by_gsm, md_by_gsm)
  invisible(gc(verbose = FALSE))

  write_csv(do.call(rbind, manifests), "01_121_input_sample_manifest_v02.csv")
  write_csv(reference_fm, "02_121_gene_feature_mapping_v02.csv")
  md$gsm <- factor(as.character(md$gsm), levels = GSM_IDS)
  md$dpi <- as.integer(md$dpi)
  md$dpi_label <- factor(paste0(md$dpi, " dpi"), levels = c("14 dpi", "28 dpi"))
  md$high_contamination_ge50 <- factor(ifelse(md$decontX_contamination_percent >= 50, ">=50% estimated contamination",
    "<50% estimated contamination"), levels = c("<50% estimated contamination", ">=50% estimated contamination"))
  md$known_28dpi_age_discordance <- (md$dpi == 28L)
  md$original_UMI_source <- "119 per-GSM singlet_raw_counts_v02.rds (not this assay)"

  seu <- SeuratObject::CreateSeuratObject(counts = corrected_mat, assay = "DecontX", project = "GSE281205",
    meta.data = md, min.cells = 0, min.features = 0)
  SeuratObject::DefaultAssay(seu) <- "DecontX"
  assert(ncol(seu) == n_expected && identical(colnames(seu), colnames(corrected_mat)), "Seurat creation changed cell set/order")
  rm(corrected_mat, md)
  invisible(gc(verbose = FALSE))

  seu <- Seurat::NormalizeData(object = seu, assay = "DecontX", normalization.method = "LogNormalize",
    scale.factor = 10000, verbose = TRUE)
  seu <- Seurat::FindVariableFeatures(object = seu, assay = "DecontX", selection.method = "vst", nfeatures = N_VARIABLE_FEATURES,
    verbose = TRUE)
  hvgs <- SeuratObject::VariableFeatures(seu)
  assert(length(hvgs) >= N_PCS + 1L, "Not enough variable genes for PCA")
  write_csv(data.frame(rank = seq_along(hvgs), gene_id = hvgs, gene_symbol = as.character(reference_fm$gene_symbol[match(hvgs,
    reference_fm$gene_id)]), stringsAsFactors = FALSE), "03_121_HVG_3000_v02.csv")
  seu <- Seurat::ScaleData(object = seu, assay = "DecontX", features = hvgs, verbose = TRUE)
  seu <- Seurat::RunPCA(object = seu, assay = "DecontX", features = hvgs, npcs = N_PCS, seed.use = SEED,
    verbose = TRUE)
  pcs <- SeuratObject::Embeddings(seu[["pca"]])
  assert(nrow(pcs) == ncol(seu) && ncol(pcs) >= N_PCS, "PCA embeddings have unexpected dimensions")
  assert(identical(rownames(pcs), colnames(seu)), "PCA embedding order not aligned to Seurat object")
  pc_names <- colnames(pcs)[seq_len(N_PCS)]
  pc_stdev <- seu[["pca"]]@stdev[seq_len(N_PCS)]

  pc_fraction <- 100 * pc_stdev^2/sum(pc_stdev^2)
  write_csv(data.frame(PC = pc_names, stdev = pc_stdev, variance = pc_stdev^2, variance_percent_among_first_30 = pc_fraction,
    cumulative_percent_among_first_30 = cumsum(pc_fraction), stringsAsFactors = FALSE), "04_121_PCA_eigenvalue_summary_v02.csv")

  stopifnot(identical(HARMONY_GROUP, "gsm"))
  seu <- harmony::RunHarmony(object = seu, group.by.vars = HARMONY_GROUP, reduction.use = "pca", dims.use = seq_len(N_PCS),
    reduction.save = "harmony", project.dim = FALSE, theta = HARMONY_THETA, verbose = TRUE)
  assert("harmony" %in% Seurat::Reductions(seu), "Harmony reduction missing after RunHarmony")
  har <- SeuratObject::Embeddings(seu[["harmony"]])
  assert(nrow(har) == nrow(pcs) && ncol(har) >= N_PCS && identical(rownames(har), rownames(pcs)) &&
    all(is.finite(har)), "Harmony embeddings are not finite or cells are misaligned")
  harm_stdev <- apply(har[, seq_len(N_PCS), drop = FALSE], 2, stats::sd)
  assert(all(is.finite(harm_stdev)) && all(harm_stdev > 0), "Harmony embedding variance is invalid")
  harm_fraction <- 100 * harm_stdev^2/sum(harm_stdev^2)
  write_csv(data.frame(dimension = seq_len(N_PCS), stdev = harm_stdev, variance = harm_stdev^2, variance_percent_among_first_30 = harm_fraction,
    stringsAsFactors = FALSE), "12_121_Harmony_dimension_variance_v02.csv")

  emb <- data.frame(cell_id = colnames(seu), gsm = as.character(seu$gsm), dpi = as.integer(seu$dpi),
    dpi_label = as.character(seu$dpi_label), replicate = as.integer(seu$replicate), age_assignment_status = as.character(seu$age_assignment_status),
    nCount_RNA = as.numeric(seu$nCount_RNA), nFeature_RNA = as.numeric(seu$nFeature_RNA), percent.mt = as.numeric(seu$percent.mt),
    decontX_contamination_percent = as.numeric(seu$decontX_contamination_percent), high_contamination_ge50 = as.character(seu$high_contamination_ge50),
    stringsAsFactors = FALSE)
  for (j in seq_len(min(10L, ncol(pcs)))) {
    emb[[paste0("PC", j)]] <- as.numeric(pcs[, j])
  }

  for (j in seq_len(min(10L, ncol(har)))) {
    emb[[paste0("HARM", j)]] <- as.numeric(har[, j])
  }

  if (RUN_UMAP) {
    seu <- Seurat::RunUMAP(object = seu, reduction = "pca", dims = seq_len(N_PCS), reduction.name = "umap_pca",
      reduction.key = "PCAUMAP_", seed.use = SEED, n.neighbors = 30L, min.dist = 0.3, return.model = FALSE,
      verbose = TRUE)
    um <- SeuratObject::Embeddings(seu[["umap_pca"]])
    assert(identical(rownames(um), emb$cell_id), "Pre-Harmony UMAP cell IDs not aligned")
    emb$UMAP1 <- as.numeric(um[, 1L])
    emb$UMAP2 <- as.numeric(um[, 2L])
    seu <- Seurat::RunUMAP(object = seu, reduction = "harmony", dims = seq_len(N_PCS), reduction.name = "umap_harmony",
      reduction.key = "HARMUMAP_", seed.use = SEED, n.neighbors = 30L, min.dist = 0.3, return.model = FALSE,
      verbose = TRUE)
    hum <- SeuratObject::Embeddings(seu[["umap_harmony"]])
    assert(identical(rownames(hum), emb$cell_id), "Harmony UMAP cell IDs not aligned")
    emb$HARMUMAP1 <- as.numeric(hum[, 1L])
    emb$HARMUMAP2 <- as.numeric(hum[, 2L])
  }

  write_csv(emb, "05_121_per_cell_embedding_metadata_v02.csv")
  sample_stats <- do.call(rbind, lapply(GSM_IDS, function(gsm) {
    z <- emb[emb$gsm == gsm, , drop = FALSE]
    assert(nrow(z) > 0L, paste("Missing GSM in merged metadata:", gsm))
    data.frame(gsm = gsm, dpi = unique(z$dpi)[[1L]], n_cells = nrow(z), original_UMI_median = med(z$nCount_RNA),
      original_features_median = med(z$nFeature_RNA), original_percent_mt_median = med(z$percent.mt),
      contamination_median_percent = med(z$decontX_contamination_percent), contamination_p95_percent = q95(z$decontX_contamination_percent),
      n_contamination_ge50 = sum(z$decontX_contamination_percent >= 50), PC1_centroid = mean(z$PC1),
      PC2_centroid = mean(z$PC2), PC1_sd = stats::sd(z$PC1), PC2_sd = stats::sd(z$PC2), Harmony1_centroid = mean(z$HARM1),
      Harmony2_centroid = mean(z$HARM2), Harmony1_sd = stats::sd(z$HARM1), Harmony2_sd = stats::sd(z$HARM2),
      stringsAsFactors = FALSE)
  }))
  write_csv(sample_stats, "06_121_per_sample_PCA_QC_summary_v02.csv")

  r2 <- do.call(rbind, lapply(seq_len(N_PCS), function(j) {
    x <- as.numeric(pcs[, j])
    data.frame(PC = paste0("PC", j), variance_percent_among_first_30 = pc_fraction[[j]], stage_14_vs_28_R2 = group_r2(x,
      emb$dpi), GSM_all_six_R2 = group_r2(x, emb$gsm), within_14dpi_GSM_R2 = group_r2(x[emb$dpi ==
      14L], emb$gsm[emb$dpi == 14L]), within_28dpi_GSM_R2 = group_r2(x[emb$dpi == 28L], emb$gsm[emb$dpi ==
      28L]), interpretation = "Cell-level descriptive variance only; no replicate-based inference",
      stringsAsFactors = FALSE)
  }))
  write_csv(r2, "07_121_PCA_stage_and_GSM_variance_audit_v02.csv")
  harm_r2 <- do.call(rbind, lapply(seq_len(N_PCS), function(j) {
    x <- as.numeric(har[, j])
    data.frame(dimension = paste0("Harmony", j), variance_percent_among_first_30 = harm_fraction[[j]],
      stage_14_vs_28_R2 = group_r2(x, emb$dpi), GSM_all_six_R2 = group_r2(x, emb$gsm), within_14dpi_GSM_R2 = group_r2(x[emb$dpi ==
        14L], emb$gsm[emb$dpi == 14L]), within_28dpi_GSM_R2 = group_r2(x[emb$dpi == 28L], emb$gsm[emb$dpi ==
        28L]), interpretation = "Descriptive cell-level effect; stage can be attenuated by GSM correction",
      stringsAsFactors = FALSE)
  }))
  write_csv(harm_r2, "13_121_Harmony_stage_and_GSM_variance_audit_v02.csv")
  weighted <- function(tab, field) {
    w <- as.numeric(tab$variance_percent_among_first_30)/100
    sum(w * as.numeric(tab[[field]]), na.rm = TRUE)
  }
  stage_before <- weighted(r2, "stage_14_vs_28_R2")
  stage_after <- weighted(harm_r2, "stage_14_vs_28_R2")
  stage_ratio <- if (stage_before > 0)
    stage_after/stage_before else NA_real_
  comparison <- data.frame(reduction = c("uncorrected_pca", "gsm_harmony"), total_cells = nrow(emb),
    stage_R2_variance_weighted = c(stage_before, stage_after), GSM_all_six_R2_variance_weighted = c(weighted(r2,
      "GSM_all_six_R2"), weighted(harm_r2, "GSM_all_six_R2")), within_14dpi_GSM_R2_variance_weighted = c(weighted(r2,
      "within_14dpi_GSM_R2"), weighted(harm_r2, "within_14dpi_GSM_R2")), within_28dpi_GSM_R2_variance_weighted = c(weighted(r2,
      "within_28dpi_GSM_R2"), weighted(harm_r2, "within_28dpi_GSM_R2")), note = "Descriptive PC-coordinate R2; weights normalized within reduction; not a batch-free biological validation",
    stringsAsFactors = FALSE)
  write_csv(comparison, "14_121_PCA_Harmony_comparison_v02.csv")
  if (is.finite(stage_ratio) && stage_ratio < 0.5) {
    log_line("WARNING: Apparent stage separation was attenuated >50% by Harmony.")
    log_line("WARNING: GSM nested within dpi; strongly evaluate overcorrection.")
  }

  mixing <- NULL
  if (RUN_KNN_SAMPLE_MIXING) {
    mixing <- rbind(transform(audit_mixing(pcs[, seq_len(N_PCS), drop = FALSE], emb, k = KNN_K),
      reduction = "uncorrected_pca"), transform(audit_mixing(har[, seq_len(N_PCS), drop = FALSE],
      emb, k = KNN_K), reduction = "gsm_harmony"))
    write_csv(mixing, "08_121_within_stage_sample_kNN_mixing_v02.csv")
  }

  plot_embedding(emb, "PC1", "PC2", "gsm", file.path(FIG_DIR, "Fig_121_PCA_by_GSM_v02.png"), "PCA before batch correction | 6 GSM libraries")
  plot_embedding(emb, "PC1", "PC2", "dpi_label", file.path(FIG_DIR, "Fig_121_PCA_by_dpi_v02.png"),
    "PCA before batch correction | 14 vs 28 dpi")
  plot_embedding(emb, "PC1", "PC2", "decontX_contamination_percent", file.path(FIG_DIR, "Fig_121_PCA_by_contamination_v02.png"),
    "PCA | DecontX estimated contamination (%)", "numeric")
  if (RUN_UMAP) {
    plot_embedding(emb, "UMAP1", "UMAP2", "gsm", file.path(FIG_DIR, "Fig_121_UMAP_by_GSM_v02.png"),
      "Unintegrated UMAP | 6 GSM libraries")
    plot_embedding(emb, "UMAP1", "UMAP2", "dpi_label", file.path(FIG_DIR, "Fig_121_UMAP_by_dpi_v02.png"),
      "Unintegrated UMAP | 14 vs 28 dpi")
    plot_embedding(emb, "UMAP1", "UMAP2", "high_contamination_ge50", file.path(FIG_DIR, "Fig_121_UMAP_by_high_contamination_v02.png"),
      "Unintegrated UMAP | high contamination flags (NOT filtered)")
    plot_embedding(emb, "UMAP1", "UMAP2", "nFeature_RNA", file.path(FIG_DIR, "Fig_121_UMAP_by_original_nFeature_v02.png"),
      "Unintegrated UMAP | original genes detected", "numeric")
    plot_embedding(emb, "HARMUMAP1", "HARMUMAP2", "gsm", file.path(FIG_DIR, "Fig_121_Harmony_UMAP_by_GSM_v02.png"),
      "Harmony-adjusted UMAP | 6 GSM libraries")
    plot_embedding(emb, "HARMUMAP1", "HARMUMAP2", "dpi_label", file.path(FIG_DIR, "Fig_121_Harmony_UMAP_by_dpi_v02.png"),
      "Harmony-adjusted UMAP | 14 vs 28 dpi")
    plot_embedding(emb, "HARMUMAP1", "HARMUMAP2", "high_contamination_ge50", file.path(FIG_DIR, "Fig_121_Harmony_UMAP_high_contamination_v02.png"),
      "Harmony-adjusted UMAP | high contamination (not removed)")
    plot_embedding(emb, "HARMUMAP1", "HARMUMAP2", "nFeature_RNA", file.path(FIG_DIR, "Fig_121_Harmony_UMAP_original_nFeature_v02.png"),
      "Harmony-adjusted UMAP | original detected genes", "numeric")

    for (groupcol in c("gsm", "dpi_label")) {
      co <- rbind(data.frame(x = emb$UMAP1, y = emb$UMAP2, group = emb[[groupcol]], embedding = "Before Harmony (PCA-UMAP)"),
        data.frame(x = emb$HARMUMAP1, y = emb$HARMUMAP2, group = emb[[groupcol]], embedding = "After Harmony (Harmony-UMAP)"))
      cmp <- ggplot2::ggplot(co, ggplot2::aes(x = x, y = y, color = group)) + ggplot2::geom_point(size = 0.18,
        alpha = 0.3) + ggplot2::facet_wrap(~embedding, ncol = 2, scales = "free") + ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::labs(title = paste0("Before/after Harmony | ",
        groupcol), subtitle = "GSM nested within stage; Harmony may reduce true stage differences",
        color = groupcol, x = "UMAP 1", y = "UMAP 2")
      ggplot2::ggsave(file.path(FIG_DIR, paste0("Fig_121_before_after_Harmony_by_", groupcol, "_v02.png")),
        cmp, width = 15, height = 7, dpi = 170)
    }
  }
  elbow <- ggplot2::ggplot(data.frame(PC = seq_len(N_PCS), stdev = pc_stdev), ggplot2::aes(x = PC,
    y = stdev)) + ggplot2::geom_line(size = 0.65) + ggplot2::geom_point(size = 1.8) + ggplot2::scale_x_continuous(breaks = seq(1,
    N_PCS, by = 2)) + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = "PCA elbow: first 30 PCs",
    y = "PC standard deviation", subtitle = "DecontX LogNormalize, 3000 HVGs, no batch correction")
  ggplot2::ggsave(file.path(FIG_DIR, "Fig_121_PCA_elbow_v02.png"), elbow, width = 9, height = 6, dpi = 170)

  r2long <- rbind(data.frame(PC = seq_len(N_PCS), R2 = r2$stage_14_vs_28_R2, variable = "Injury stage (biology)"),
    data.frame(PC = seq_len(N_PCS), R2 = r2$within_14dpi_GSM_R2, variable = "GSM within 14 dpi"),
    data.frame(PC = seq_len(N_PCS), R2 = r2$within_28dpi_GSM_R2, variable = "GSM within 28 dpi"))
  p_r2 <- ggplot2::ggplot(r2long, ggplot2::aes(x = PC, y = R2, color = variable)) + ggplot2::geom_line(size = 0.8) +
    ggplot2::geom_point(size = 1) + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = "Descriptive PCA variance attribution",
    subtitle = "Cell-level R-squared only; no P-values or batch correction", y = "Fraction of PC score variance (R-squared)",
    color = "Grouping")
  ggplot2::ggsave(file.path(FIG_DIR, "Fig_121_PCA_stage_GSM_R2_v02.png"), p_r2, width = 11, height = 6.5,
    dpi = 170)
  r2_compare <- rbind(data.frame(dimension = seq_len(N_PCS), R2 = r2$stage_14_vs_28_R2, reduction = "PCA"),
    data.frame(dimension = seq_len(N_PCS), R2 = harm_r2$stage_14_vs_28_R2, reduction = "Harmony"))
  p_compare <- ggplot2::ggplot(r2_compare, ggplot2::aes(x = dimension, y = R2, color = reduction)) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1) + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Stage (14 vs 28 dpi) separation before / after Harmony", subtitle = "Per-coordinate cell-level descriptive R2; not biological replicate inference",
      x = "Dimension index", y = "Stage R2", color = "Reduction")
  ggplot2::ggsave(file.path(FIG_DIR, "Fig_121_stage_preservation_PCA_Harmony_v02.png"), p_compare,
    width = 10, height = 6.2, dpi = 170)

  if (SAVE_SEURAT_RDS) {
    sout <- file.path(RDS_DIR, "121_GSE281205_singlets_DecontX_LogNorm_PCA_Harmony_v02.rds")
    seu@misc$step121 <- list(input_step = "119_v02", input_assay = "DecontX corrected estimated counts",
      original_UMI_source_directory = IN_RDS, original_UMI_filename_template = "119_{GSM}_singlet_raw_counts_v02.rds",
      corrected_RDS_filename_template = "119_{GSM}_singlet_DecontX_counts_v02.rds", normalization = "LogNormalize, scale.factor=10000",
      HVGs = N_VARIABLE_FEATURES, PCA_dims = N_PCS, harmony = TRUE, harmony_group_by = HARMONY_GROUP,
      harmony_theta = HARMONY_THETA, harmony_reduction = "harmony", uncorrected_reduction = "pca",
      uncorrected_umap = if (RUN_UMAP) "umap_pca" else NA_character_, corrected_umap = if (RUN_UMAP) "umap_harmony" else NA_character_,
      rpca = FALSE, sctransform = FALSE, correction_of_stage = FALSE, clustering = FALSE, age_28dpi_warning = "GEO title/treatment field discordance retained",
      high_contamination_cells_removed = FALSE, signature_1320_used = FALSE, seed = SEED)
    saveRDS(seu, sout, compress = "gzip")
    assert(file.exists(sout) && file.info(sout)$size > 1000, "Output Seurat RDS was not saved successfully")
  }

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, "09_121_REPORT_v02.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(OUT_DIR, "10_121_sessionInfo_v02.txt"),
    useBytes = TRUE)
  file_names <- c("00_121_run_log_v02.log", "01_121_input_sample_manifest_v02.csv", "02_121_gene_feature_mapping_v02.csv",
    "03_121_HVG_3000_v02.csv", "04_121_PCA_eigenvalue_summary_v02.csv", "05_121_per_cell_embedding_metadata_v02.csv",
    "06_121_per_sample_PCA_QC_summary_v02.csv", "07_121_PCA_stage_and_GSM_variance_audit_v02.csv",
    if (RUN_KNN_SAMPLE_MIXING) "08_121_within_stage_sample_kNN_mixing_v02.csv", "09_121_REPORT_v02.txt",
    "10_121_sessionInfo_v02.txt", "12_121_Harmony_dimension_variance_v02.csv", "13_121_Harmony_stage_and_GSM_variance_audit_v02.csv",
    "14_121_PCA_Harmony_comparison_v02.csv", if (SAVE_SEURAT_RDS) file.path("RDS", "121_GSE281205_singlets_DecontX_LogNorm_PCA_Harmony_v02.rds"),
    file.path("figures", list.files(FIG_DIR, pattern = "^Fig_121_.*_v02\\.png$")))
  index <- data.frame(relative_path = file_names, exists = file.exists(file.path(OUT_DIR, file_names)),
    stringsAsFactors = FALSE)
  write_csv(index, "11_121_OUTPUT_INDEX_v02.csv")
  assert(all(index$exists), "Some expected Step 121 outputs were not generated")
  invisible(list(n_cells = ncol(seu), n_hvg = length(hvgs), n_samples = nrow(sample_stats), output_directory = OUT_DIR))
}

tryCatch(main(), error = function(e) {
  log_line("FAILED | Step121 v02: ", conditionMessage(e))
  stop(e)
})

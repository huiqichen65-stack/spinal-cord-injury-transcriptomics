options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- "D:/zz"
STEP121 <- file.path(ROOT, "results", "121_GSE281205_normalize_PCA_batch_audit")
INPUT_RDS <- file.path(STEP121, "RDS", "121_GSE281205_singlets_DecontX_LogNorm_PCA_Harmony_v02.rds")
FEATURE_CSV <- file.path(STEP121, "02_121_gene_feature_mapping_v02.csv")
OUT_DIR <- file.path(ROOT, "results", "123_GSE281205_major_lineage_clustering")
FIG_DIR <- file.path(OUT_DIR, "figures")
RDS_DIR <- file.path(OUT_DIR, "RDS")
LOG_DIR <- file.path(OUT_DIR, "logs")
VERSION <- "v03"
SEED <- 12301L

RESUME_FROM_COMPLETE_V01_CSV <- TRUE
DIMS_USE <- 1:30
K_PARAM <- 20L
HARMONY_RESOLUTIONS <- c(0.4, 0.6, 0.8)
MAIN_CLUSTER_COL <- "cluster_harmony_r04"
PCA_CLUSTER_COL <- "cluster_pca_r04"
RUN_PCA_SENSITIVITY <- TRUE
RUN_EXPLORATORY_HVG_MARKERS <- TRUE
MARKER_MAX_CELLS_PER_IDENT <- 400L
SAVE_FULL_SEURAT_RDS <- TRUE

RDS_COMPRESS <- FALSE

stopif <- function(cond, message) {
  if (!isTRUE(cond))
    stop(message, call. = FALSE)
}
for (pkg in c("Seurat", "SeuratObject", "Matrix", "ggplot2")) {
  stopif(requireNamespace(pkg, quietly = TRUE), paste0("Missing package: ", pkg, ". Install it and rerun Step 123."))
}
stopif(file.exists(INPUT_RDS), paste("Missing Step121 RDS:", INPUT_RDS))
stopif(file.exists(FEATURE_CSV), paste("Missing Step121 mapping CSV:", FEATURE_CSV))
for (p in c(OUT_DIR, FIG_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
}
LOG_FILE <- file.path(OUT_DIR, "00_123_run_log_v03.log")
cat("", file = LOG_FILE)
log_line <- function(...) {
  v <- paste0(..., collapse = "")
  cat(v, "\n", sep = "")
  cat(v, "\n", sep = "", file = LOG_FILE, append = TRUE)
  invisible(v)
}
write_csv <- function(x, filename) {
  stopif(is.data.frame(x), paste("Not a data.frame:", filename))
  utils::write.csv(x, file.path(OUT_DIR, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
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
choose2 <- function(x) x * (x - 1)/2
adjusted_rand_index <- function(a, b) {
  tab <- table(as.character(a), as.character(b))
  n <- sum(tab)
  if (n < 2L)
    return(NA_real_)
  s <- sum(choose2(as.numeric(tab)))
  aa <- sum(choose2(rowSums(tab)))
  bb <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- aa * bb/total
  maximum <- (aa + bb)/2
  denom <- maximum - expected
  if (abs(denom) < 1e-12)
    return(if (abs(s - expected) < 1e-12) 1 else NA_real_)
  (s - expected)/denom
}

LINEAGE_MARKERS <- list(Microglia = c("P2ry12", "Tmem119", "Sall1", "Hexb", "Fcrls"), Monocyte_macrophage = c("Ccr2",
  "Ly6c2", "Lyz2", "Ms4a7", "Fcgr1"), APC_dendritic = c("Flt3", "Itgax", "Cd74", "H2-Ab1", "H2-Aa"),
  Neutrophil = c("S100a8", "S100a9", "Retnlg", "Ly6g"), T_cell = c("Cd3d", "Cd3e", "Trbc2"), NK_cytotoxic = c("Nkg7",
    "Klrb1c", "Prf1", "Gzmb"), B_cell = c("Cd79a", "Cd79b", "Ms4a1"), Astrocyte = c("Aqp4", "Slc1a3",
    "Aldh1l1", "Gfap"), Oligodendrocyte = c("Mbp", "Mog", "Plp1", "Cnp"), OPC = c("Pdgfra", "Cspg4",
    "Olig1", "Olig2"), Neuron = c("Snap25", "Rbfox3", "Syt1", "Tubb3"), Ependymal = c("Foxj1", "Dnah5",
    "Tubb4b"), Endothelial = c("Pecam1", "Kdr", "Cldn5", "Vwf"), Pericyte_vSMC = c("Rgs5", "Pdgfrb",
    "Acta2", "Tagln"), Fibroblast_meningeal = c("Col1a1", "Col1a2", "Dcn", "Lum"))

DISPLAY_ONLY_MARKERS <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Gpnmb", "Spp1", "Abca1", "Cd68", "Tyrobp",
  "Fabp5")

save_plot <- function(g, filename, width = 11, height = 7) {
  ggplot2::ggsave(file.path(FIG_DIR, filename), plot = g, width = width, height = height, dpi = 160,
    limitsize = FALSE)
}
plot_embedding <- function(seu, reduction, colour_col, title, filename, label_clusters = FALSE) {
  xy <- SeuratObject::Embeddings(seu[[reduction]])
  stopif(ncol(xy) >= 2, paste("Missing 2D coordinates:", reduction))
  md <- seu[[]]
  stopif(identical(rownames(xy), rownames(md)), paste("Cell order mismatch in embedding:", reduction))
  dat <- data.frame(x = xy[, 1], y = xy[, 2], group = as.character(md[[colour_col]]))
  dat$group <- factor(dat$group, levels = unique(dat$group))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, colour = group)) + ggplot2::geom_point(size = 0.19,
    alpha = 0.48) + ggplot2::labs(title = title, x = paste(reduction, "1"), y = paste(reduction,
    "2"), colour = colour_col, subtitle = "Exploratory embedding; plotted cells are NOT animal replicates") +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (label_clusters) {
    cent <- stats::aggregate(cbind(x, y) ~ group, dat, stats::median)
    p <- p + ggplot2::geom_text(data = cent, mapping = ggplot2::aes(x = x, y = y, label = group),
      inherit.aes = FALSE, colour = "black", size = 3.4, fontface = "bold", check_overlap = TRUE)
  }
  save_plot(p, filename, 12, 8)
}

resume_complete_v01 <- function() {
  required <- c("01_123_input_sample_audit_v01.csv", "02_123_cluster_counts_by_resolution_v01.csv",
    "03_123_per_cell_cluster_membership_v01.csv", "04_123_sample_by_cluster_counts_v01.csv", "05_123_stage_cluster_fractions_v01.csv",
    "06_123_cluster_QC_contamination_v01.csv", "07_123_marker_gene_mapping_v01.csv", "08_123_main_cluster_marker_evidence_v01.csv",
    "09_123_provisional_major_lineage_scores_v01.csv", "10_123_manual_lineage_review_template_v01.csv",
    "11_123_clustering_stability_ARI_v01.csv", "12_123_harmony_vs_pca_crosstab_v01.csv", "13_123_exploratory_HVG_cluster_markers_v01.csv",
    "13a_123_top25_exploratory_cluster_markers_v01.csv")
  required_paths <- file.path(OUT_DIR, required)
  if (!all(file.exists(required_paths))) {
    return(FALSE)
  }
  read_old <- function(file) utils::read.csv(file.path(OUT_DIR, file), stringsAsFactors = FALSE, check.names = FALSE)
  membership <- read_old("03_123_per_cell_cluster_membership_v01.csv")
  expected_columns <- c("cell_id", "gsm", "dpi", "replicate", "contamination_percent", "cluster_harmony_r04",
    "cluster_harmony_r06", "cluster_harmony_r08", "cluster_pca_r04")
  stopif(all(expected_columns %in% names(membership)), "v01 cell membership CSV does not contain all four cluster labels")
  stopif(!anyNA(membership[, expected_columns]), "v01 cell membership contains NA values")
  stopif(!anyDuplicated(membership$cell_id), "v01 cell membership contains duplicate cell IDs")
  seu <- readRDS(INPUT_RDS)
  stopif(inherits(seu, "Seurat"), "Input Step121 RDS is not a Seurat object")
  stopif(ncol(seu) == 42361L, "The Step121 singlet count changed; unsafe to resume v01")
  stopif(identical(as.character(membership$cell_id), colnames(seu)), "Step121 cell IDs/order differ from the v01 cluster CSV")
  md <- seu[[]]
  stopif(identical(as.character(membership$gsm), as.character(md$gsm)), "GSM ID mismatch in the existing clustering CSV")
  stopif(identical(as.integer(membership$dpi), as.integer(md$dpi)), "Stage mismatch in the existing clustering CSV")
  stopif(identical(as.integer(membership$replicate), as.integer(md$replicate)), "Replicate mismatch in the existing clustering CSV")
  stopif(isTRUE(all.equal(as.numeric(membership$contamination_percent), as.numeric(md$decontX_contamination_percent),
    tolerance = 1e-06, check.attributes = FALSE)), "Contamination metadata mismatch in the existing clustering CSV")
  stopif(all(c("pca", "harmony", "umap_harmony", "umap_pca") %in% Seurat::Reductions(seu)), "Step121 embeddings missing; unsafe to resume")
  stopif("DecontX" %in% SeuratObject::Assays(seu), "Expected DecontX assay is missing")
  SeuratObject::DefaultAssay(seu) <- "DecontX"

  seu <- Seurat::FindNeighbors(seu, reduction = "harmony", dims = DIMS_USE, k.param = K_PARAM, graph.name = c("harmony_nn_123",
    "harmony_snn_123"), verbose = FALSE)
  if (RUN_PCA_SENSITIVITY) {
    seu <- Seurat::FindNeighbors(seu, reduction = "pca", dims = DIMS_USE, k.param = K_PARAM, graph.name = c("pca_nn_123",
      "pca_snn_123"), verbose = FALSE)
  }
  for (nm in expected_columns[6:9]) {
    seu[[nm]] <- as.character(membership[[nm]])
  }
  primary <- as.character(membership[[MAIN_CLUSTER_COL]])
  stopif(all(nzchar(primary)), "Missing main cluster values in cache")
  seu$seurat_clusters <- factor(primary, levels = sort(unique(primary)))
  seu <- SeuratObject::`Idents<-`(seu, value = seu$seurat_clusters)
  counts_res <- read_old("02_123_cluster_counts_by_resolution_v01.csv")
  stopif(all(c("cluster_column", "cluster", "n_cells") %in% names(counts_res)), "Cluster count CSV malformed")
  for (nm in expected_columns[6:9]) {
    recorded <- counts_res[counts_res$cluster_column == nm, c("cluster", "n_cells"), drop = FALSE]
    actual <- table(as.character(membership[[nm]]))
    stopif(nrow(recorded) == length(actual) && all(as.integer(recorded$n_cells[match(names(actual),
      recorded$cluster)]) == as.integer(actual)), paste("Cached membership/counts disagree:", nm))
  }
  review <- read_old("10_123_manual_lineage_review_template_v01.csv")
  stopif(all(c("cluster", "suggested_lineage") %in% names(review)) && !anyDuplicated(review$cluster),
    "Review template lacks unique provisional lineages")
  stopif(setequal(review$cluster, unique(primary)), "Cached review clusters differ from cached main clusters")
  suggested <- setNames(as.character(review$suggested_lineage), as.character(review$cluster))
  seu$provisional_major_lineage_123 <- unname(suggested[primary])
  seu$dpi_label <- paste0(as.integer(md$dpi), " dpi")
  seu$high_contam_label_123 <- ifelse(as.numeric(md$decontX_contamination_percent) >= 50, ">=50% model-estimated contamination",
    "<50% model-estimated contamination")
  stage_tab <- read_old("05_123_stage_cluster_fractions_v01.csv")
  stopif(all(c("dpi", "cluster", "percent_within_stage") %in% names(stage_tab)), "v01 stage composition table is malformed")
  stopif(setequal(as.character(stage_tab$dpi), c("14", "28")), "Stage table must include 14 and 28 dpi")
  stopif(all(is.finite(stage_tab$percent_within_stage)), "Nonfinite stage composition percent")
  for (d in c("14", "28")) {
    subtotal <- sum(stage_tab$percent_within_stage[as.character(stage_tab$dpi) == d])
    stopif(abs(subtotal - 100) < 1e-05, paste("Stage composition does not sum to 100:", d))
  }
  stage_tab$dpi <- factor(as.character(stage_tab$dpi), levels = c("14", "28"))
  stage_tab$cluster <- factor(stage_tab$cluster, levels = sort(unique(primary)))
  pstage <- ggplot2::ggplot(stage_tab, ggplot2::aes(x = dpi, y = percent_within_stage, fill = cluster)) +
    ggplot2::geom_col(position = "stack") + ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Sampled cells by stage and Harmony cluster", subtitle = paste("Cell-level proportions only; 3 libraries per stage;",
      "not an animal-level test"), x = "Days post injury", y = "Percent of recovered singlets")
  save_plot(pstage, "Fig_123_stage_cluster_composition_v03.png", 12, 7)

  for (old in required_paths) {
    dest <- file.path(OUT_DIR, sub("_v01", "_v03", basename(old), fixed = TRUE))
    stopif(isTRUE(file.copy(old, dest, overwrite = TRUE)) && file.exists(dest), paste("Cannot copy v01 table to v03:",
      dest))
  }

  oldfig <- list.files(FIG_DIR, pattern = "_v01[.]png$", full.names = TRUE)
  for (old in oldfig) {
    dest <- file.path(FIG_DIR, sub("_v01", "_v03", basename(old), fixed = TRUE))
    if (basename(dest) == "Fig_123_stage_cluster_composition_v03.png") {
      next
    }
    stopif(isTRUE(file.copy(old, dest, overwrite = TRUE)) && file.exists(dest), paste("Cannot copy v01 figure to v03:",
      dest))
  }
  md_save <- seu[[]]
  saveRDS(md_save, file.path(RDS_DIR, "123_GSE281205_cell_cluster_metadata_v03.rds"), compress = "gzip")
  seu@misc$step123 <- list(input = INPUT_RDS, input_step = "121_v02", assay = "DecontX corrected normalized estimates; NOT raw UMI",
    raw_integer_UMI_source = paste0("D:/zz/results/", "119_GSE281205_scDblFinder/RDS/", "119_{GSM}_singlet_raw_counts_v02.rds"),
    reduction_primary = "harmony", reduction_sensitivity = "pca", dims = DIMS_USE, k_param = K_PARAM,
    resolutions = HARMONY_RESOLUTIONS, primary_cluster_column = MAIN_CLUSTER_COL, other_cluster_columns = c("cluster_harmony_r04",
      "cluster_harmony_r06", "cluster_harmony_r08", PCA_CLUSTER_COL), provisional_major_lineage_column = "provisional_major_lineage_123",
    manual_review_completed = FALSE, high_contamination_cells_removed = FALSE, signature_1320_used = FALSE,
    pseudobulk_validation_performed = FALSE, recovered_from_verified_v01_cell_assignment_CSV = TRUE,
    warning = paste("GSM nested within dpi; stage biology may be overcorrected;", "28dpi GEO title/treatment age discordance requires disclosure"),
    seed = SEED)
  final <- file.path(RDS_DIR, "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
  stopif(!file.exists(final), paste("Final v03 RDS already exists; do not overwrite:", final))
  partial <- paste0(final, ".partial")
  stopif(!file.exists(partial), paste("Existing partial v03 RDS must be inspected first:", partial))
  saveRDS(seu, partial, compress = RDS_COMPRESS)
  stopif(file.exists(partial) && is.finite(file.info(partial)$size) && file.info(partial)$size > 1e+06,
    "Final v03 RDS save incomplete")
  stopif(isTRUE(file.rename(partial, final)), paste("Cannot rename partial file to final RDS:", final))
  qc_table <- read_old("06_123_cluster_QC_contamination_v01.csv")
  compare <- read_old("11_123_clustering_stability_ARI_v01.csv")
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, "14_123_REPORT_v03.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(OUT_DIR, "15_123_sessionInfo_v03.txt"),
    useBytes = TRUE)
  outputs <- list.files(OUT_DIR, recursive = TRUE, full.names = FALSE, include.dirs = FALSE)
  indx <- data.frame(relative_path = outputs, file_size_bytes = as.numeric(file.info(file.path(OUT_DIR,
    outputs))$size), stringsAsFactors = FALSE)
  write_csv(indx, "16_123_OUTPUT_INDEX_v03.csv")
  TRUE
}

main <- function() {
  set.seed(SEED)
  log_line("Seurat: ", as.character(utils::packageVersion("Seurat")))
  if (RESUME_FROM_COMPLETE_V01_CSV) {
    if (resume_complete_v01())
      return(invisible(TRUE))
  }
  final_rds <- file.path(RDS_DIR, "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
  if (SAVE_FULL_SEURAT_RDS && file.exists(final_rds)) {
    stop("The final Step123 Seurat RDS already exists: ", final_rds, ". Move/rename it before an intentional rerun.",
      call. = FALSE)
  }

  seu <- readRDS(INPUT_RDS)
  stopif(inherits(seu, "Seurat"), "Step121 input is not a Seurat object")
  stopif(all(c("DecontX") %in% SeuratObject::Assays(seu)), "Expected DecontX assay missing")
  SeuratObject::DefaultAssay(seu) <- "DecontX"
  reds <- Seurat::Reductions(seu)
  stopif(all(c("pca", "harmony", "umap_harmony", "umap_pca") %in% reds), paste("Step121 reductions incomplete:",
    paste(reds, collapse = ", ")))
  for (red in c("pca", "harmony")) {
    m <- SeuratObject::Embeddings(seu[[red]])
    stopif(nrow(m) == ncol(seu) && ncol(m) >= max(DIMS_USE) && identical(rownames(m), colnames(seu)) &&
      all(is.finite(m[, DIMS_USE, drop = FALSE])), paste("Invalid embeddings:", red))
  }
  md <- seu[[]]
  req <- c("gsm", "dpi", "replicate", "age_assignment_status", "nFeature_RNA", "nCount_RNA", "percent.mt",
    "decontX_contamination_percent")
  stopif(all(req %in% names(md)), "Missing required Step121 cell metadata")
  stopif(identical(rownames(md), colnames(seu)), "Cell metadata not aligned")
  gsm_ids <- paste0("GSM", 8614662:8614667)
  stopif(setequal(unique(as.character(md$gsm)), gsm_ids), "Unexpected GSMs: expected six GSM8614662..667")
  stopif(ncol(seu) == 42361L, paste("Expected 42,361 Step121 singlets but found", ncol(seu)))
  sample_stage <- unique(data.frame(gsm = as.character(md$gsm), dpi = as.integer(md$dpi)))
  stopif(nrow(sample_stage) == 6L && all(sample_stage$dpi[match(gsm_ids, sample_stage$gsm)] == rep(c(14L,
    28L), each = 3L)), "Incorrect GSM-stage assignment")
  stopif(all(is.finite(md$decontX_contamination_percent)) && all(md$decontX_contamination_percent >=
    0 & md$decontX_contamination_percent <= 100), "Invalid DecontX contamination metadata")
  stopif(!is.null(seu@misc$step121) && isTRUE(seu@misc$step121$harmony) && identical(seu@misc$step121$harmony_group_by,
    "gsm"), "Input is not the expected Step121-v02 GSM-Harmony object")
  audit <- data.frame(gsm = gsm_ids, dpi = rep(c(14L, 28L), each = 3L), n_cells = as.integer(table(factor(as.character(md$gsm),
    gsm_ids))), n_contamination_ge50 = vapply(gsm_ids, function(g) sum(as.character(md$gsm) == g &
    md$decontX_contamination_percent >= 50), numeric(1)), median_nFeature_RNA = vapply(gsm_ids, function(g) med(md$nFeature_RNA[as.character(md$gsm) ==
    g]), numeric(1)), median_percent_mt = vapply(gsm_ids, function(g) med(md$percent.mt[as.character(md$gsm) ==
    g]), numeric(1)), age_assignment_status = vapply(gsm_ids, function(g) paste(unique(as.character(md$age_assignment_status[as.character(md$gsm) ==
    g])), collapse = ";"), character(1)), stringsAsFactors = FALSE)
  write_csv(audit, "01_123_input_sample_audit_v03.csv")

  seu <- Seurat::FindNeighbors(object = seu, reduction = "harmony", dims = DIMS_USE, k.param = K_PARAM,
    graph.name = c("harmony_nn_123", "harmony_snn_123"), verbose = TRUE)
  res_cols <- c("cluster_harmony_r04", "cluster_harmony_r06", "cluster_harmony_r08")
  for (i in seq_along(HARMONY_RESOLUTIONS)) {
    r <- HARMONY_RESOLUTIONS[[i]]
    target <- res_cols[[i]]
    seu <- Seurat::FindClusters(object = seu, graph.name = "harmony_snn_123", resolution = r, cluster.name = target,
      algorithm = 1, random.seed = SEED, n.start = 10, n.iter = 10, verbose = TRUE)
    stopif(target %in% names(seu[[]]), paste("FindClusters did not create:", target))
  }
  if (RUN_PCA_SENSITIVITY) {
    seu <- Seurat::FindNeighbors(object = seu, reduction = "pca", dims = DIMS_USE, k.param = K_PARAM,
      graph.name = c("pca_nn_123", "pca_snn_123"), verbose = TRUE)
    seu <- Seurat::FindClusters(object = seu, graph.name = "pca_snn_123", resolution = 0.4, cluster.name = PCA_CLUSTER_COL,
      algorithm = 1, random.seed = SEED, n.start = 10, n.iter = 10, verbose = TRUE)
  }
  md <- seu[[]]
  for (nm in c(res_cols, if (RUN_PCA_SENSITIVITY) PCA_CLUSTER_COL)) {
    stopif(all(!is.na(md[[nm]])), paste("Missing cluster assignments:", nm))
  }

  seu$seurat_clusters <- factor(as.character(md[[MAIN_CLUSTER_COL]]), levels = sort(unique(as.character(md[[MAIN_CLUSTER_COL]]))))
  seu <- SeuratObject::`Idents<-`(seu, value = seu$seurat_clusters)
  md <- seu[[]]

  counts_res <- do.call(rbind, lapply(res_cols, function(x) {
    tab <- table(as.character(md[[x]]))
    data.frame(reduction = "harmony", resolution = as.numeric(HARMONY_RESOLUTIONS[match(x, res_cols)]),
      cluster_column = x, cluster = names(tab), n_cells = as.integer(tab), stringsAsFactors = FALSE)
  }))
  if (RUN_PCA_SENSITIVITY) {
    tab <- table(as.character(md[[PCA_CLUSTER_COL]]))
    counts_res <- rbind(counts_res, data.frame(reduction = "pca", resolution = 0.4, cluster_column = PCA_CLUSTER_COL,
      cluster = names(tab), n_cells = as.integer(tab)))
  }
  write_csv(counts_res, "02_123_cluster_counts_by_resolution_v03.csv")
  membership <- data.frame(cell_id = colnames(seu), gsm = as.character(md$gsm), dpi = as.integer(md$dpi),
    replicate = as.integer(md$replicate), contamination_percent = as.numeric(md$decontX_contamination_percent),
    cluster_harmony_r04 = as.character(md$cluster_harmony_r04), cluster_harmony_r06 = as.character(md$cluster_harmony_r06),
    cluster_harmony_r08 = as.character(md$cluster_harmony_r08), cluster_pca_r04 = if (RUN_PCA_SENSITIVITY)
      as.character(md[[PCA_CLUSTER_COL]]) else NA_character_, stringsAsFactors = FALSE)
  write_csv(membership, "03_123_per_cell_cluster_membership_v03.csv")
  t_gsm <- as.data.frame(table(gsm = factor(as.character(md$gsm), gsm_ids), cluster = as.character(md[[MAIN_CLUSTER_COL]])),
    responseName = "n_cells", stringsAsFactors = FALSE)
  t_gsm$dpi <- ifelse(t_gsm$gsm %in% gsm_ids[1:3], 14L, 28L)
  t_gsm <- t_gsm[, c("dpi", "gsm", "cluster", "n_cells")]
  t_gsm$percent_within_sample <- 100 * t_gsm$n_cells/as.numeric(table(factor(as.character(md$gsm),
    gsm_ids))[match(t_gsm$gsm, gsm_ids)])
  write_csv(t_gsm, "04_123_sample_by_cluster_counts_v03.csv")
  stage_tab <- as.data.frame(table(dpi = as.integer(md$dpi), cluster = as.character(md[[MAIN_CLUSTER_COL]])),
    responseName = "n_cells", stringsAsFactors = FALSE)
  stage_tab$percent_within_stage <- 100 * stage_tab$n_cells/as.numeric(table(as.integer(md$dpi))[match(stage_tab$dpi,
    names(table(as.integer(md$dpi))))])
  stage_tab$note <- "Cell-level composition; not animal-level proportions"
  write_csv(stage_tab, "05_123_stage_cluster_fractions_v03.csv")

  clusters <- sort(unique(as.character(md[[MAIN_CLUSTER_COL]])))
  clu_idx <- lapply(clusters, function(k) which(as.character(md[[MAIN_CLUSTER_COL]]) == k))
  names(clu_idx) <- clusters
  qc_table <- do.call(rbind, lapply(clusters, function(k) {
    idx <- clu_idx[[k]]
    contamination <- as.numeric(md$decontX_contamination_percent[idx])
    data.frame(cluster = k, n_cells = length(idx), median_nFeature_RNA = med(md$nFeature_RNA[idx]),
      median_nCount_RNA = med(md$nCount_RNA[idx]), median_percent_mt = med(md$percent.mt[idx]),
      median_contamination_percent = med(contamination), p95_contamination_percent = q95(contamination),
      n_contamination_ge50 = sum(contamination >= 50), pct_contamination_ge50 = 100 * mean(contamination >=
        50), n_contamination_ge80 = sum(contamination >= 80), pct_contamination_ge80 = 100 *
        mean(contamination >= 80), largest_gsm_fraction = 100 * max(table(as.character(md$gsm[idx])))/length(idx),
      high_contamination_review = (mean(contamination >= 50) >= 0.2 || med(contamination) >= 30),
      stringsAsFactors = FALSE)
  }))
  write_csv(qc_table, "06_123_cluster_QC_contamination_v03.csv")

  f <- utils::read.csv(FEATURE_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  stopif(all(c("gene_id", "gene_symbol") %in% names(f)), "Step121 gene mapping lacks gene_id / gene_symbol")
  stopif(!anyDuplicated(f$gene_id), "Duplicate Ensembl IDs in feature CSV")
  stopif(setequal(f$gene_id, rownames(seu)), "Step121 mapping IDs do not match Seurat features")
  all_lineage <- unique(unlist(LINEAGE_MARKERS, use.names = FALSE))
  all_symbols <- unique(c(all_lineage, DISPLAY_ONLY_MARKERS))
  map_rows <- lapply(all_symbols, function(sym) {
    ids <- as.character(f$gene_id[f$gene_symbol == sym & f$gene_id %in% rownames(seu)])
    if (length(ids) == 0L) {
      return(data.frame(gene_symbol = sym, gene_id = NA_character_, mapping_status = "not_found",
        duplicate_ids = 0L, role = ifelse(sym %in% DISPLAY_ONLY_MARKERS && !sym %in% all_lineage,
          "display_only", "broad_lineage")))
    }

    if (length(ids) > 1L) {
      mm <- SeuratObject::LayerData(seu, assay = "DecontX", layer = "data", features = ids)
      nn <- Matrix::rowSums(mm > 0)
      best <- rownames(mm)[which.max(as.numeric(nn))]
      status <- paste0("duplicate_symbol_selected_by_detection;alternatives=", paste(ids, collapse = ";"))
    } else {
      best <- ids[[1]]
      status <- "unique"
    }
    data.frame(gene_symbol = sym, gene_id = best, mapping_status = status, duplicate_ids = length(ids),
      role = ifelse(sym %in% all_lineage, "broad_lineage", "display_only"), stringsAsFactors = FALSE)
  })
  marker_map <- do.call(rbind, map_rows)
  write_csv(marker_map, "07_123_marker_gene_mapping_v03.csv")
  mapped <- marker_map[!is.na(marker_map$gene_id), , drop = FALSE]
  stopif(nrow(mapped) >= 20L, "Too few canonical markers mapped: verify gene ID/symbol mapping")
  log_line("Mapped canonical/display symbols: ", nrow(mapped), "/", nrow(marker_map), "; missing: ",
    paste(marker_map$gene_symbol[is.na(marker_map$gene_id)], collapse = ", "))
  marker_ids <- unique(as.character(mapped$gene_id))
  expr <- SeuratObject::LayerData(seu, assay = "DecontX", layer = "data", features = marker_ids)
  stopif(setequal(rownames(expr), marker_ids) && setequal(colnames(expr), colnames(seu)), "Marker data layer not aligned with Seurat cells")
  expr <- expr[marker_ids, colnames(seu), drop = FALSE]
  marker_data <- do.call(rbind, lapply(clusters, function(k) {
    idx <- clu_idx[[k]]
    mm <- expr[, idx, drop = FALSE]
    data.frame(cluster = k, gene_id = marker_ids, mean_log_normalized_expr = as.numeric(Matrix::rowMeans(mm)),
      pct_cells_detected = 100 * as.numeric(Matrix::rowMeans(mm > 0)), n_cluster_cells = length(idx),
      stringsAsFactors = FALSE)
  }))
  marker_data$gene_symbol <- marker_map$gene_symbol[match(marker_data$gene_id, marker_map$gene_id)]
  marker_data$role <- marker_map$role[match(marker_data$gene_id, marker_map$gene_id)]
  marker_data <- marker_data[, c("cluster", "gene_symbol", "gene_id", "role", "n_cluster_cells", "mean_log_normalized_expr",
    "pct_cells_detected")]
  write_csv(marker_data, "08_123_main_cluster_marker_evidence_v03.csv")

  expr_mean <- xtabs(mean_log_normalized_expr ~ gene_symbol + cluster, marker_data)
  expr_z <- t(apply(expr_mean, 1L, function(x) {
    s <- stats::sd(x)
    if (!is.finite(s) || s < 1e-10)
      return(rep(0, length(x)))
    (x - mean(x))/s
  }))
  rownames(expr_z) <- rownames(expr_mean)
  colnames(expr_z) <- colnames(expr_mean)
  expr_pct <- xtabs(pct_cells_detected ~ gene_symbol + cluster, marker_data)
  lineage_rows <- list()
  for (k in clusters) {
    for (lin in names(LINEAGE_MARKERS)) {
      syms <- intersect(LINEAGE_MARKERS[[lin]], rownames(expr_z))
      n_hit <- if (length(syms))
        sum(expr_pct[syms, k] >= 10) else 0L
      lineage_rows[[length(lineage_rows) + 1L]] <- data.frame(cluster = k, lineage = lin, n_panel_genes = length(LINEAGE_MARKERS[[lin]]),
        n_mapped_genes = length(syms), n_markers_detected_in_ge10pct = n_hit, mean_z_across_markers = if (length(syms))
          mean(expr_z[syms, k]) else NA_real_, stringsAsFactors = FALSE)
    }
  }
  lineage_score <- do.call(rbind, lineage_rows)
  write_csv(lineage_score, "09_123_provisional_major_lineage_scores_v03.csv")

  review <- do.call(rbind, lapply(clusters, function(k) {
    sub <- lineage_score[lineage_score$cluster == k & lineage_score$n_mapped_genes >= 2L, , drop = FALSE]
    sub <- sub[order(-sub$mean_z_across_markers), , drop = FALSE]
    best <- if (nrow(sub))
      sub[1, , drop = FALSE] else NULL
    runner <- if (nrow(sub) >= 2L)
      sub[2, , drop = FALSE] else NULL
    gap <- if (nrow(sub) >= 2L)
      as.numeric(best$mean_z_across_markers - runner$mean_z_across_markers) else NA_real_
    qrow <- qc_table[qc_table$cluster == k, , drop = FALSE]

    possible <- !is.null(best) && best$n_markers_detected_in_ge10pct >= 2L && best$mean_z_across_markers >=
      0.35 && is.finite(gap) && gap >= 0.15 && !qrow$high_contamination_review && qrow$n_cells >=
      30L
    data.frame(cluster = k, n_cells = qrow$n_cells, suggested_lineage = if (possible)
      as.character(best$lineage) else "Unresolved", best_lineage_candidate = if (!is.null(best))
      as.character(best$lineage) else NA_character_, best_score = if (!is.null(best))
      as.numeric(best$mean_z_across_markers) else NA_real_, second_candidate = if (!is.null(runner))
      as.character(runner$lineage) else NA_character_, score_margin = gap, supporting_marker_count_ge10pct = if (!is.null(best))
      as.integer(best$n_markers_detected_in_ge10pct) else NA_integer_, pct_contamination_ge50 = qrow$pct_contamination_ge50, high_contamination_review = qrow$high_contamination_review,
      needs_manual_review = TRUE, manual_confirmed_lineage = "", manual_review_notes = "", stringsAsFactors = FALSE)
  }))
  write_csv(review, "10_123_manual_lineage_review_template_v03.csv")

  suggestion <- setNames(as.character(review$suggested_lineage), review$cluster)
  seu$provisional_major_lineage_123 <- unname(suggestion[as.character(seu[[MAIN_CLUSTER_COL, drop = TRUE]])])

  compare <- data.frame(clustering_pair = c("Harmony r0.4 vs r0.6", "Harmony r0.4 vs r0.8", if (RUN_PCA_SENSITIVITY) "Harmony r0.4 vs PCA r0.4"),
    adjusted_rand_index = c(adjusted_rand_index(md[[res_cols[1]]], md[[res_cols[2]]]), adjusted_rand_index(md[[res_cols[1]]],
      md[[res_cols[3]]]), if (RUN_PCA_SENSITIVITY) adjusted_rand_index(md[[res_cols[1]]], md[[PCA_CLUSTER_COL]])),
    n_cells = ncol(seu), note = "Cell-wise cluster stability, descriptive; not biological replication",
    stringsAsFactors = FALSE)
  write_csv(compare, "11_123_clustering_stability_ARI_v03.csv")
  if (RUN_PCA_SENSITIVITY) {
    ct <- as.data.frame(table(harmony_cluster = as.character(md[[MAIN_CLUSTER_COL]]), pca_cluster = as.character(md[[PCA_CLUSTER_COL]])),
      responseName = "n_cells", stringsAsFactors = FALSE)
    write_csv(ct, "12_123_harmony_vs_pca_crosstab_v03.csv")
  }
  rm(expr, expr_mean, expr_z, expr_pct)
  invisible(gc(verbose = FALSE))

  if (RUN_EXPLORATORY_HVG_MARKERS) {
    tryCatch({
      seu <- SeuratObject::`Idents<-`(seu, value = factor(as.character(seu[[MAIN_CLUSTER_COL, drop = TRUE]])))
      panel <- unique(c(SeuratObject::VariableFeatures(seu), marker_ids))
      markers <- Seurat::FindAllMarkers(object = seu, assay = "DecontX", features = panel, only.pos = TRUE,
        test.use = "wilcox", min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = MARKER_MAX_CELLS_PER_IDENT,
        return.thresh = 1, random.seed = SEED, verbose = FALSE)
      if (!is.null(markers) && nrow(markers) > 0L) {
        fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))
        stopif(length(fc_col) > 0L, "FindAllMarkers has no fold-change column")
        markers$gene_id <- as.character(markers$gene)
        markers$gene_symbol <- f$gene_symbol[match(markers$gene_id, f$gene_id)]
        markers <- markers[order(as.character(markers$cluster), -markers[[fc_col[[1]]]]), , drop = FALSE]
        write_csv(markers, "13_123_exploratory_HVG_cluster_markers_v03.csv")
        top_n <- do.call(rbind, lapply(split(markers, as.character(markers$cluster)), function(d) head(d,
          25)))
        write_csv(top_n, "13a_123_top25_exploratory_cluster_markers_v03.csv")
      } else {
        log_line("WARNING: FindAllMarkers returned no rows; use marker panel evidence.")
      }
    }, error = function(e) {
      log_line("WARNING: exploratory FindAllMarkers failed: ", conditionMessage(e), "; core cluster outputs remain valid")
    })
  }

  plot_embedding(seu, "umap_harmony", MAIN_CLUSTER_COL, "Whole-cell Harmony SNN clusters | resolution 0.4",
    "Fig_123_Harmony_clusters_r04_v03.png", TRUE)
  for (j in 2:3) {
    plot_embedding(seu, "umap_harmony", res_cols[[j]], paste0("Whole-cell Harmony SNN clusters | resolution ",
      HARMONY_RESOLUTIONS[[j]]), paste0("Fig_123_Harmony_clusters_r0", c(6, 8)[j - 1], "_v03.png"),
      TRUE)
  }
  if (RUN_PCA_SENSITIVITY) {
    plot_embedding(seu, "umap_pca", PCA_CLUSTER_COL, "Uncorrected PCA SNN clusters | resolution 0.4",
      "Fig_123_PCA_clusters_r04_v03.png", TRUE)
    plot_embedding(seu, "umap_harmony", PCA_CLUSTER_COL, "PCA clusters mapped onto Harmony UMAP",
      "Fig_123_PCA_clusters_on_Harmony_UMAP_v03.png", FALSE)
  }
  plot_embedding(seu, "umap_harmony", "gsm", "Harmony UMAP | six independent GSM libraries", "Fig_123_Harmony_by_GSM_v03.png")
  if (!"dpi_label" %in% colnames(seu[[]])) {
    seu$dpi_label <- paste0(seu$dpi, " dpi")
  }
  plot_embedding(seu, "umap_harmony", "dpi_label", "Harmony UMAP | 14 dpi vs 28 dpi", "Fig_123_Harmony_by_stage_v03.png")
  plot_embedding(seu, "umap_harmony", "provisional_major_lineage_123", "Broad-lineage suggestions ONLY (awaiting manual review)",
    "Fig_123_provisional_major_lineage_v03.png")
  md2 <- seu[[]]
  md2$high_contam_label <- ifelse(md2$decontX_contamination_percent >= 50, ">=50% model-estimated contamination",
    "<50% model-estimated contamination")
  seu$high_contam_label_123 <- md2$high_contam_label
  plot_embedding(seu, "umap_harmony", "high_contam_label_123", "DecontX >=50% flags: visible, NOT automatically removed",
    "Fig_123_Harmony_high_contamination_v03.png")

  tile <- marker_data[marker_data$role == "broad_lineage", , drop = FALSE]
  tile$cluster <- factor(tile$cluster, levels = clusters)
  tile$gene_symbol <- factor(tile$gene_symbol, levels = rev(unique(unlist(LINEAGE_MARKERS, use.names = FALSE))))
  pheat <- ggplot2::ggplot(tile, ggplot2::aes(x = cluster, y = gene_symbol, fill = mean_log_normalized_expr)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient(low = "white", high = "steelblue4") + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7), axis.text.x = ggplot2::element_text(size = 9),
      panel.grid = ggplot2::element_blank()) + ggplot2::labs(title = "Canonical major-lineage markers by Harmony cluster",
    subtitle = "Mean DecontX LogNormalize expression; not a definitive cell-type assignment", x = "Harmony cluster (resolution 0.4)",
    y = "Gene symbol", fill = "Mean log norm.")
  save_plot(pheat, "Fig_123_major_lineage_marker_heatmap_v03.png", max(12, 0.42 * length(clusters) +
    6), 16)

  pc <- ggplot2::ggplot(qc_table, ggplot2::aes(x = factor(cluster, levels = clusters), y = pct_contamination_ge50)) +
    ggplot2::geom_col(fill = "firebrick3") + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = "High-contamination estimates across main clusters",
    subtitle = "Flags only; >=50% DecontX estimates do not automatically imply doublets", x = "Harmony cluster 0.4",
    y = "Percent with model-estimated contamination >=50%")
  save_plot(pc, "Fig_123_cluster_high_contamination_v03.png", 12, 6)

  stage_tab$dpi <- factor(as.character(stage_tab$dpi), levels = c("14", "28"))
  pstage <- ggplot2::ggplot(stage_tab, ggplot2::aes(x = dpi, y = percent_within_stage, fill = cluster)) +
    ggplot2::geom_col(position = "stack") + ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Sampled cells by stage and Harmony cluster", subtitle = "Cell-level proportions only; 3 libraries per stage; not an animal-level test",
      x = "Days post injury", y = "Percent of recovered singlets")
  save_plot(pstage, "Fig_123_stage_cluster_composition_v03.png", 12, 7)

  md_save <- seu[[]]
  saveRDS(md_save, file.path(RDS_DIR, "123_GSE281205_cell_cluster_metadata_v03.rds"), compress = "gzip")
  if (SAVE_FULL_SEURAT_RDS) {
    seu@misc$step123 <- list(input = INPUT_RDS, input_step = "121_v02", assay = "DecontX corrected normalized estimates; NOT raw UMI",
      raw_integer_UMI_source = paste0("D:/zz/results/", "119_GSE281205_scDblFinder/RDS/", "119_{GSM}_singlet_raw_counts_v02.rds"),
      reduction_primary = "harmony", reduction_sensitivity = "pca", dims = DIMS_USE, k_param = K_PARAM,
      resolutions = HARMONY_RESOLUTIONS, primary_cluster_column = MAIN_CLUSTER_COL, other_cluster_columns = c(res_cols,
        PCA_CLUSTER_COL), provisional_major_lineage_column = "provisional_major_lineage_123",
      manual_review_completed = FALSE, high_contamination_cells_removed = FALSE, signature_1320_used = FALSE,
      pseudobulk_validation_performed = FALSE, warning = paste("GSM nested within dpi; stage biology may be overcorrected;",
        "28dpi GEO title/treatment age discordance requires disclosure"), seed = SEED)
    outfile <- file.path(RDS_DIR, "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
    if (file.exists(outfile)) {
      stop("Output full Seurat RDS already exists; move it aside before rerunning: ", outfile,
        call. = FALSE)
    }
    tmpfile <- paste0(outfile, ".partial")
    if (file.exists(tmpfile))
      unlink(tmpfile)
    saveRDS(seu, tmpfile, compress = RDS_COMPRESS)
    stopif(file.exists(tmpfile) && file.info(tmpfile)$size > 1e+06, "Full Seurat RDS save incomplete or unexpectedly small")
    stopif(isTRUE(file.rename(tmpfile, outfile)), paste("Could not finalize RDS; partial file remains:",
      tmpfile))
  }

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, "14_123_REPORT_v03.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(OUT_DIR, "15_123_sessionInfo_v03.txt"),
    useBytes = TRUE)
  results <- list.files(OUT_DIR, recursive = TRUE, full.names = FALSE, include.dirs = FALSE)
  index <- data.frame(relative_path = results, file_size_bytes = as.numeric(file.info(file.path(OUT_DIR,
    results))$size), stringsAsFactors = FALSE)
  write_csv(index, "16_123_OUTPUT_INDEX_v03.csv")
  invisible(list(n_cells = ncol(seu), n_main_clusters = length(clusters), output_directory = OUT_DIR))
}

tryCatch(main(), error = function(e) {
  log_line("ERROR | Step123 v03 | ", conditionMessage(e))
  stop(e)
})

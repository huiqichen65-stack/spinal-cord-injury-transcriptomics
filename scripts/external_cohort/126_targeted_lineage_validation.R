source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
ROOT <- PROJECT_ROOT
SRC123 <- file.path(ROOT, "results", "123_GSE281205_major_lineage_clustering")
SRC125 <- file.path(ROOT, "results", "125_GSE281205_major_lineage_annotation")
SRC121 <- file.path(ROOT, "results", "121_GSE281205_normalize_PCA_batch_audit")
OUT_DIR <- file.path(ROOT, "results", "126_GSE281205_targeted_lineage_review")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
VERSION <- "v01"
EXPECTED_N <- 42361L
GSM_IDS <- paste0("GSM", 8614662:8614667)
EXPECTED_CLUSTERS <- as.character(0:23)
FOCUS_CLUSTERS <- c("0", "3", "8", "10", "13", "18", "19", "22")
REFERENCE_CLUSTERS <- c("1", "2", "6")
DISPLAY_CLUSTERS <- sort(unique(c(FOCUS_CLUSTERS, REFERENCE_CLUSTERS)))
MAIN_COL <- "cluster_harmony_r04"

IN_RDS <- file.path(SRC123, "RDS", "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
IN_MEM <- file.path(SRC123, "03_123_per_cell_cluster_membership_v03.csv")
IN_MARKER <- file.path(SRC123, "07_123_marker_gene_mapping_v03.csv")
IN_FM <- file.path(SRC121, "02_121_gene_feature_mapping_v02.csv")
IN_EVID <- file.path(SRC125, "02_125_24cluster_evidence_summary_v01.csv")
IN_DEC <- file.path(SRC125, "04_125_manual_decisions_v01.csv")
LOG_FILE <- file.path(LOG_DIR, "00_126_run_log_v01.log")

for (d in c(OUT_DIR, FIG_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

assert <- function(ok, message) {
  if (!isTRUE(ok))
    stop(message, call. = FALSE)
}
for (p in c(IN_RDS, IN_MEM, IN_MARKER, IN_FM, IN_EVID, IN_DEC)) assert(file.exists(p), paste("Required source file not found:",
  p))
for (pkg in c("SeuratObject", "Matrix", "ggplot2")) assert(requireNamespace(pkg, quietly = TRUE), paste0("Package ",
  pkg, " missing; install.packages('", pkg, "')"))

cat("", file = LOG_FILE)
logline <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG_FILE, append = TRUE)
  invisible(msg)
}
read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}
write_csv <- function(x, name) {
  assert(is.data.frame(x), paste("Invalid table:", name))
  utils::write.csv(x, file.path(OUT_DIR, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  logline("Wrote ", name, " [", nrow(x), " rows]")
}
save_plot <- function(p, name, width = 11, height = 7) {
  ggplot2::ggsave(filename = file.path(FIG_DIR, name), plot = p, width = width, height = height, dpi = 160,
    limitsize = FALSE)
  logline("Saved figure: ", name)
}
median_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x))
    return(NA_real_)
  stats::median(x)
}
prc <- function(x) {
  if (!length(x))
    return(NA_real_)
  100 * mean(x, na.rm = TRUE)
}

MARKER_PANELS <- list(microglia_identity = c("P2ry12", "Tmem119", "Sall1", "Hexb", "Fcrls"), infiltrating_myeloid = c("Ccr2",
  "Ly6c2", "Ms4a7", "Lyz2", "Fcgr1"), myeloid_activation = c("Clec7a", "Trem2", "Apoe", "Lpl", "Gpnmb",
  "Spp1"), APC_DC = c("Flt3", "Cd209a", "Ccl17", "Itgax", "Cd74", "H2-Aa", "H2-Ab1"), neutrophil = c("Ly6g",
  "S100a8", "S100a9", "Retnlg"), endothelial = c("Pecam1", "Cldn5", "Kdr", "Vwf"), mural = c("Rgs5",
  "Pdgfrb", "Acta2", "Tagln", "Myh11"), proliferation = c("Mki67", "Pclaf", "Birc5", "Top2a", "Ube2c"),
  stromal = c("Col1a1", "Dcn", "Lum"), nonmyeloid_control = c("Cd3e", "Nkg7", "Aqp4", "Mbp", "Pdgfra",
    "Snap25"))
ALL_SYMS <- unique(unlist(MARKER_PANELS, use.names = FALSE))
PAIR_DEFS <- data.frame(pair = c("resident_P2ry12_Fcrls", "resident_Tmem119_Fcrls", "resident_Fcrls_Clec7a",
  "resident_P2ry12_Lpl", "infiltrating_Ccr2_Ly6c2", "myeloid_Lyz2_Fcgr1", "activation_Apoe_Lpl", "activation_Lpl_Gpnmb",
  "endothelial_Pecam1_Cldn5", "mural_Rgs5_Pdgfrb", "vascular_cross_Pecam1_Rgs5", "vascular_cross_Cldn5_Pdgfrb",
  "APC_Cd74_H2_Aa", "DC_Flt3_Cd209a", "neutrophil_Ly6g_S100a8", "cycling_Pclaf_Birc5", "cycling_Mki67_P2ry12"),
  marker_a = c("P2ry12", "Tmem119", "Fcrls", "P2ry12", "Ccr2", "Lyz2", "Apoe", "Lpl", "Pecam1", "Rgs5",
    "Pecam1", "Cldn5", "Cd74", "Flt3", "Ly6g", "Pclaf", "Mki67"), marker_b = c("Fcrls", "Fcrls",
    "Clec7a", "Lpl", "Ly6c2", "Fcgr1", "Lpl", "Gpnmb", "Cldn5", "Pdgfrb", "Rgs5", "Pdgfrb", "H2-Aa",
    "Cd209a", "S100a8", "Birc5", "P2ry12"), stringsAsFactors = FALSE)

REVIEW_QUESTIONS <- c(`0` = "Resident microglial markers versus infiltrating macrophage markers; check each GSM",
  `3` = "Lipid/activation signal with resident microglial identity versus infiltrating/other myeloid",
  `8` = "Very high estimated contamination and low depth: is there coherent native lineage?", `10` = "Mural markers and estimated contamination: check whether cells retain coherent vascular identity",
  `13` = "Endothelial and mural gene co-detection at CELL LEVEL versus mixed cluster of distinct cells",
  `18` = "APC/DC versus neutrophil markers in individual cells, accounting for contamination", `19` = "Low-depth/high-contamination cells: unresolved lineage or coherent identifiable subpopulation?",
  `22` = "Proliferation with resident microglia co-detection; verify source identity")

run_step126 <- function() {
  logline("==========================================================")
  logline("STEP126 v01 | GSE281205 targeted major-lineage adjudication")
  logline("Start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  logline("Inputs: direct Step121/123/125 files, NO ZIP")
  logline("Focus: ", paste(FOCUS_CLUSTERS, collapse = ","))
  logline("No 1320-gene signature; NO DE/pseudobulk; NO automatic approval")
  logline("==========================================================")

  mem <- read_csv(IN_MEM)
  ev <- read_csv(IN_EVID)
  dec <- read_csv(IN_DEC)
  fm <- read_csv(IN_FM)
  known <- read_csv(IN_MARKER)
  assert(all(c("cell_id", "gsm", "dpi", "replicate", MAIN_COL, "cluster_harmony_r06", "cluster_harmony_r08",
    "cluster_pca_r04", "contamination_percent") %in% names(mem)), "Step123 cell membership columns mismatch")
  assert(all(c("cluster", "n_cells", "draft_lineage") %in% names(ev)), "Step125 evidence columns mismatch")
  assert(all(c("cluster", "manual_lineage", "review_status", "myeloid_gate") %in% names(dec)), "Step125 decision template columns mismatch")
  assert(all(c("gene_id", "gene_symbol") %in% names(fm)), "Step121 feature mapping columns mismatch")
  assert(all(c("gene_symbol", "gene_id") %in% names(known)), "Step123 marker mapping columns mismatch")
  assert(nrow(mem) == EXPECTED_N && !anyDuplicated(mem$cell_id), "Step123 expected 42,361 unique cells")
  assert(setequal(as.character(mem$gsm), GSM_IDS), "Unexpected GSM identifiers")
  assert(all(as.integer(mem$dpi) %in% c(14L, 28L)), "Expected only 14/28 dpi")
  assert(!anyNA(mem$contamination_percent), "Missing DecontX contamination")
  assert(setequal(as.character(mem[[MAIN_COL]]), EXPECTED_CLUSTERS), "Unexpected cluster ID set (expected 0..23)")
  assert(!anyDuplicated(ev$cluster) && setequal(as.character(ev$cluster), EXPECTED_CLUSTERS), "Step125 evidence not 24 unique clusters")
  assert(!anyDuplicated(dec$cluster) && setequal(as.character(dec$cluster), EXPECTED_CLUSTERS), "Step125 decisions not 24 unique clusters")
  for (g in GSM_IDS) {
    stages <- unique(as.integer(mem$dpi[mem$gsm == g]))
    e <- if (g %in% GSM_IDS[1:3])
      14L else 28L
    assert(length(stages) == 1L && stages == e, paste("Unexpected stage for", g))
  }
  mcnt <- table(as.character(mem[[MAIN_COL]]))
  assert(all(as.integer(mcnt[as.character(ev$cluster)]) == as.integer(ev$n_cells)), "Step125 cluster counts mismatch")

  audit <- data.frame(source_file = c(IN_RDS, IN_MEM, IN_MARKER, IN_FM, IN_EVID, IN_DEC), available = file.exists(c(IN_RDS,
    IN_MEM, IN_MARKER, IN_FM, IN_EVID, IN_DEC)), purpose = c("expression/UMAP", "cell IDs and clusters",
    "canonical mappings", "all gene ID/symbol mappings", "cluster marker/QC summary", "human-editable decisions (READ ONLY)"),
    stringsAsFactors = FALSE)
  write_csv(audit, "01_126_input_audit_v01.csv")

  fm$gene_id <- as.character(fm$gene_id)
  fm$gene_symbol <- as.character(fm$gene_symbol)
  mapping <- do.call(rbind, lapply(ALL_SYMS, function(sym) {
    ids <- unique(fm$gene_id[!is.na(fm$gene_symbol) & fm$gene_symbol == sym])
    ids <- ids[!is.na(ids) & nzchar(ids)]
    data.frame(gene_symbol = sym, gene_id = if (length(ids) == 1L)
      ids else NA_character_, n_matching_ids = length(ids), status = if (length(ids) == 1L)
      "unique" else if (!length(ids))
      "not_found" else "ambiguous_multiple_gene_ids", stringsAsFactors = FALSE)
  }))
  write_csv(mapping, "02_126_marker_mapping_audit_v01.csv")
  mapped <- mapping[mapping$status == "unique", , drop = FALSE]
  assert(nrow(mapped) >= 35L, "Too few unique canonical markers: cannot safely interpret clusters")

  logline("Loading large Step123 Seurat RDS once (several GiB RAM): ", IN_RDS)
  seu <- readRDS(IN_RDS)
  assert(inherits(seu, "Seurat"), "Expected Seurat RDS")
  assert(ncol(seu) == EXPECTED_N, "Seurat cell count mismatch")
  assert(identical(as.character(colnames(seu)), as.character(mem$cell_id)), "RDS cell order/IDs differ from Step123 manifest")
  md <- seu[[]]
  assert(all(c("gsm", "dpi", MAIN_COL, "nCount_RNA", "nFeature_RNA", "percent.mt", "decontX_contamination_percent") %in%
    names(md)), "Seurat metadata missing fields")
  assert(identical(as.character(md$gsm), as.character(mem$gsm)) && identical(as.integer(md$dpi), as.integer(mem$dpi)) &&
    identical(as.character(md[[MAIN_COL]]), as.character(mem[[MAIN_COL]])), "GSM/stage/cluster mismatch between Step123 RDS and CSV")
  assert(isTRUE(all.equal(as.numeric(md$decontX_contamination_percent), as.numeric(mem$contamination_percent),
    tolerance = 1e-06, check.attributes = FALSE)), "DecontX contamination mismatch between CSV and RDS")
  assert("DecontX" %in% SeuratObject::Assays(seu), "Seurat object lacks DecontX assay")

  valid_map <- mapped[mapped$gene_id %in% rownames(seu), , drop = FALSE]
  absent_ids <- mapped[!(mapped$gene_id %in% rownames(seu)), , drop = FALSE]
  if (nrow(absent_ids))
    logline("Mapped genes not in Seurat object: ", paste(absent_ids$gene_symbol, collapse = ", "))
  assert(nrow(valid_map) >= 35L, "Too few marker features actually present in Seurat")
  ids <- as.character(valid_map$gene_id)
  mat <- SeuratObject::LayerData(seu, assay = "DecontX", layer = "data", features = ids)
  assert(setequal(rownames(mat), ids) && setequal(colnames(mat), colnames(seu)), "DecontX normalized marker matrix IDs do not match source")
  mat <- mat[ids, colnames(seu), drop = FALSE]
  rownames(mat) <- valid_map$gene_symbol[match(rownames(mat), valid_map$gene_id)]
  assert(!anyDuplicated(rownames(mat)), "Duplicate mapped marker symbols")
  logline("Loaded marker matrix: ", nrow(mat), " genes x ", ncol(mat), " cells. Detection threshold: DecontX LogNormalized > 0")

  cl <- as.character(mem[[MAIN_COL]])
  ci <- split(seq_len(ncol(seu)), factor(cl, levels = EXPECTED_CLUSTERS), drop = TRUE)
  cluster_markers <- do.call(rbind, lapply(EXPECTED_CLUSTERS, function(k) {
    ix <- ci[[k]]
    mm <- mat[, ix, drop = FALSE]
    data.frame(cluster = k, gene_symbol = rownames(mat), panel = vapply(rownames(mat), function(s) {
      p <- names(MARKER_PANELS)[vapply(MARKER_PANELS, function(g) s %in% g, logical(1))]
      paste(p, collapse = ";")
    }, character(1)), n_cells = length(ix), detection_percent = 100 * as.numeric(Matrix::rowMeans(mm >
      0)), mean_log_normalized = as.numeric(Matrix::rowMeans(mm)), stringsAsFactors = FALSE)
  }))
  write_csv(cluster_markers, "03_126_cluster_marker_detection_v01.csv")

  pair_rows <- list()
  append_pair <- function(k, g, ix, pd) {
    if (!length(ix))
      return(NULL)
    a <- pd$marker_a
    b <- pd$marker_b
    if (!all(c(a, b) %in% rownames(mat)))
      return(NULL)
    present_a <- as.vector(mat[a, ix, drop = FALSE] > 0)
    present_b <- as.vector(mat[b, ix, drop = FALSE] > 0)
    ca <- as.numeric(mem$contamination_percent[ix])
    data.frame(cluster = k, gsm = g, dpi = if (g == "ALL")
      NA_integer_ else as.integer(unique(mem$dpi[ix])), pair = pd$pair, marker_a = a, marker_b = b, n_cells = length(ix),
      n_a_detected = sum(present_a), n_b_detected = sum(present_b), n_both_detected = sum(present_a &
        present_b), pct_a_detected = prc(present_a), pct_b_detected = prc(present_b), pct_both_detected = prc(present_a &
        present_b), pct_a_only = prc(present_a & !present_b), pct_b_only = prc(!present_a & present_b),
      pct_neither = prc(!present_a & !present_b), pct_both_and_contamination_ge50 = prc(present_a &
        present_b & ca >= 50), stringsAsFactors = FALSE)
  }
  for (k in DISPLAY_CLUSTERS) {
    ix_all <- ci[[k]]
    for (j in seq_len(nrow(PAIR_DEFS))) {
      pr <- PAIR_DEFS[j, , drop = FALSE]
      pair_rows[[length(pair_rows) + 1L]] <- append_pair(k, "ALL", ix_all, pr)
      for (g in GSM_IDS) {
        ix <- ix_all[mem$gsm[ix_all] == g]
        item <- append_pair(k, g, ix, pr)
        if (!is.null(item))
          pair_rows[[length(pair_rows) + 1L]] <- item
      }
    }
  }
  pair_rows <- Filter(Negate(is.null), pair_rows)
  assert(length(pair_rows) > 0, "No co-detection pairs could be evaluated")
  pairs <- do.call(rbind, pair_rows)
  write_csv(pairs, "04_126_pairwise_codetection_by_cluster_GSM_v01.csv")

  qrows <- list()
  for (k in EXPECTED_CLUSTERS) {
    for (g in GSM_IDS) {
      ix <- ci[[k]][mem$gsm[ci[[k]]] == g]
      if (!length(ix))
        next
      cc <- as.numeric(mem$contamination_percent[ix])
      qrows[[length(qrows) + 1L]] <- data.frame(cluster = k, gsm = g, dpi = as.integer(unique(mem$dpi[ix])),
        n_cells = length(ix), median_nCount_RNA = median_safe(as.numeric(md$nCount_RNA[ix])),
        median_nFeature_RNA = median_safe(as.numeric(md$nFeature_RNA[ix])), median_percent_mt = median_safe(as.numeric(md$percent.mt[ix])),
        median_DecontX_contamination_percent = median_safe(cc), n_contamination_ge20 = sum(cc >=
          20), pct_contamination_ge20 = prc(cc >= 20), n_contamination_ge50 = sum(cc >= 50),
        pct_contamination_ge50 = prc(cc >= 50), n_contamination_ge80 = sum(cc >= 80), pct_contamination_ge80 = prc(cc >=
          80), stringsAsFactors = FALSE)
    }
  }
  qc_gsm <- do.call(rbind, qrows)
  write_csv(qc_gsm, "05_126_cluster_GSM_QC_v01.csv")

  nesting <- do.call(rbind, lapply(DISPLAY_CLUSTERS, function(k) {
    ix <- ci[[k]]
    do.call(rbind, lapply(c("cluster_harmony_r06", "cluster_harmony_r08", "cluster_pca_r04"), function(subcol) {
      z <- as.data.frame(table(as.character(mem[[subcol]][ix])), stringsAsFactors = FALSE)
      names(z) <- c("subcluster", "n_cells")
      z$cluster <- k
      z$subcluster_source <- subcol
      z$percent_of_parent_cluster <- 100 * z$n_cells/length(ix)
      z[, c("cluster", "subcluster_source", "subcluster", "n_cells", "percent_of_parent_cluster"),
        drop = FALSE]
    }))
  }))
  write_csv(nesting, "06_126_focus_cross_resolution_partition_v01.csv")

  ev$cluster <- as.character(ev$cluster)
  dec$cluster <- as.character(dec$cluster)
  review_clusters <- EXPECTED_CLUSTERS
  dossier <- do.call(rbind, lapply(review_clusters, function(k) {
    e <- ev[ev$cluster == k, , drop = FALSE]
    d <- dec[dec$cluster == k, , drop = FALSE]
    q <- qc_gsm[qc_gsm$cluster == k, , drop = FALSE]
    ps <- pairs[pairs$cluster == k & pairs$gsm == "ALL", , drop = FALSE]
    getpair <- function(tag) {
      x <- ps$pct_both_detected[ps$pair == tag]
      if (length(x))
        as.numeric(x[[1L]]) else NA_real_
    }
    data.frame(cluster = k, review_priority = if (k %in% FOCUS_CLUSTERS)
      "FOCUS" else "REFERENCE", n_cells = as.integer(e$n_cells), previous_draft = as.character(e$draft_lineage),
      step125_review_status = as.character(d$review_status), step125_manual_lineage = as.character(d$manual_lineage),
      high_contamination_ge50_percent = as.numeric(e$pct_contamination_ge50), max_GSM_contamination_ge50_percent = if (nrow(q))
        max(q$pct_contamination_ge50, na.rm = TRUE) else NA_real_, co_P2ry12_Fcrls_percent = getpair("resident_P2ry12_Fcrls"), co_Fcrls_Clec7a_percent = getpair("resident_Fcrls_Clec7a"),
      co_Ccr2_Ly6c2_percent = getpair("infiltrating_Ccr2_Ly6c2"), co_Pecam1_Cldn5_percent = getpair("endothelial_Pecam1_Cldn5"),
      co_Rgs5_Pdgfrb_percent = getpair("mural_Rgs5_Pdgfrb"), co_Pecam1_Rgs5_percent = getpair("vascular_cross_Pecam1_Rgs5"),
      co_Cldn5_Pdgfrb_percent = getpair("vascular_cross_Cldn5_Pdgfrb"), co_Cd74_H2_Aa_percent = getpair("APC_Cd74_H2_Aa"),
      co_Ly6g_S100a8_percent = getpair("neutrophil_Ly6g_S100a8"), co_Pclaf_Birc5_percent = getpair("cycling_Pclaf_Birc5"),
      focus_question = if (k %in% FOCUS_CLUSTERS)
        unname(REVIEW_QUESTIONS[[k]]) else "Validate existing lineage candidate against independent markers", stringsAsFactors = FALSE)
  }))
  dossier <- dossier[order(as.integer(dossier$cluster)), , drop = FALSE]
  write_csv(dossier, "07_126_cluster_review_dossier_v01.csv")

  worksheet <- dossier[, c("cluster", "review_priority", "n_cells", "previous_draft", "focus_question",
    "step125_review_status")]
  worksheet$manual_lineage <- ""
  worksheet$review_status <- "PENDING"
  worksheet$myeloid_gate <- ""
  worksheet$validation_notes <- ""
  work_path <- file.path(OUT_DIR, "08_126_manual_review_assist_v01.csv")
  if (!file.exists(work_path)) {
    write_csv(worksheet, basename(work_path))
  } else {
    logline("Preserved existing (possibly edited) review worksheet: ", work_path)
  }

  stage_audit <- aggregate(n_cells ~ cluster + dpi, data = qc_gsm, FUN = sum)
  stage_audit$dpi <- as.integer(stage_audit$dpi)
  denoms <- aggregate(n_cells ~ dpi, data = qc_gsm, FUN = sum)
  stage_audit$pct_of_stage_singlets <- 100 * stage_audit$n_cells/denoms$n_cells[match(stage_audit$dpi,
    denoms$dpi)]
  write_csv(stage_audit, "09_126_stage_cluster_descriptive_v01.csv")

  show <- cluster_markers[cluster_markers$cluster %in% DISPLAY_CLUSTERS & cluster_markers$gene_symbol %in%
    c("P2ry12", "Fcrls", "Tmem119", "Ccr2", "Ly6c2", "Clec7a", "Lpl", "Gpnmb", "Pecam1", "Cldn5",
      "Rgs5", "Pdgfrb", "Cd74", "H2-Aa", "Ly6g", "Flt3", "Pclaf", "Birc5"), , drop = FALSE]
  show$cluster <- factor(show$cluster, levels = rev(DISPLAY_CLUSTERS))
  p1 <- ggplot2::ggplot(show, ggplot2::aes(x = gene_symbol, y = cluster, size = detection_percent,
    colour = mean_log_normalized)) + ggplot2::geom_point(alpha = 0.86) + ggplot2::scale_size(range = c(0.2,
    5)) + ggplot2::theme_bw(base_size = 11) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55,
    hjust = 1)) + ggplot2::labs(title = "Focused canonical marker detection by whole-atlas cluster",
    subtitle = "DecontX-corrected LogNormalized expression; not confirmed lineage", x = "Marker",
    y = "Harmony cluster (resolution 0.4)", size = "Detected (%)", colour = "Mean log expression")
  save_plot(p1, "Fig_126_focus_marker_dotplot_v01.png", 13, 8)

  pairs_fig <- pairs[pairs$gsm == "ALL" & pairs$pair %in% c("resident_P2ry12_Fcrls", "resident_Fcrls_Clec7a",
    "infiltrating_Ccr2_Ly6c2", "endothelial_Pecam1_Cldn5", "mural_Rgs5_Pdgfrb", "vascular_cross_Pecam1_Rgs5",
    "APC_Cd74_H2_Aa", "neutrophil_Ly6g_S100a8", "cycling_Pclaf_Birc5"), , drop = FALSE]
  pairs_fig$cluster <- factor(pairs_fig$cluster, levels = rev(DISPLAY_CLUSTERS))
  p2 <- ggplot2::ggplot(pairs_fig, ggplot2::aes(x = pair, y = cluster, fill = pct_both_detected)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 50,
    hjust = 1)) + ggplot2::labs(title = "Cell-level canonical pair co-detection", subtitle = "The two genes detected within the same corrected-expression cell; descriptive",
    x = "Marker pair", y = "Cluster", fill = "Both (%)")
  save_plot(p2, "Fig_126_focus_pairwise_codetection_v01.png", 15, 8)

  qq <- qc_gsm[qc_gsm$cluster %in% DISPLAY_CLUSTERS, , drop = FALSE]
  qq$cluster <- factor(qq$cluster, levels = DISPLAY_CLUSTERS)
  p3 <- ggplot2::ggplot(qq, ggplot2::aes(x = gsm, y = pct_contamination_ge50, fill = factor(dpi))) +
    ggplot2::geom_col() + ggplot2::facet_wrap(~cluster, ncol = 4, scales = "free_y") + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 70, hjust = 1)) + ggplot2::labs(title = "Estimated high contamination by cluster and GSM",
    subtitle = "Model-estimated >=50% contamination; not an automatic exclusion rule", x = "GSM",
    y = "Cells with >=50% DecontX (%)", fill = "Dpi")
  save_plot(p3, "Fig_126_focus_GSM_contamination_v01.png", 16, 11)

  emb_name <- if ("umap_harmony" %in% names(seu@reductions))
    "umap_harmony" else if ("umap.harmony" %in% names(seu@reductions))
    "umap.harmony" else if ("harmony.umap" %in% names(seu@reductions))
    "harmony.umap" else NA_character_
  if (is.na(emb_name)) {

    er <- names(seu@reductions)
    cand <- er[grepl("umap", er, ignore.case = TRUE) & grepl("harm|batch", er, ignore.case = TRUE)]
    if (length(cand))
      emb_name <- cand[[1L]]
  }
  if (!is.na(emb_name)) {
    xy <- SeuratObject::Embeddings(seu[[emb_name]])
    assert(setequal(rownames(xy), colnames(seu)) && ncol(xy) >= 2, "UMAP embedding missing/misaligned cells")
    xy <- xy[colnames(seu), 1:2, drop = FALSE]
    emb <- data.frame(UMAP1 = as.numeric(xy[, 1]), UMAP2 = as.numeric(xy[, 2]), cluster = cl, focus = ifelse(cl %in%
      FOCUS_CLUSTERS, paste0("Cluster ", cl), "Other clusters"), stringsAsFactors = FALSE)
    emb$focus <- factor(emb$focus, levels = c("Other clusters", paste0("Cluster ", FOCUS_CLUSTERS)))
    p4 <- ggplot2::ggplot(emb, ggplot2::aes(x = UMAP1, y = UMAP2)) + ggplot2::geom_point(colour = "grey80",
      size = 0.1, alpha = 0.27) + ggplot2::geom_point(data = emb[emb$cluster %in% FOCUS_CLUSTERS,
      ], ggplot2::aes(colour = focus), size = 0.17, alpha = 0.65) + ggplot2::coord_equal() + ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = "Focus clusters on the existing Harmony UMAP", subtitle = "Exploratory visualization only; 24-cluster identities still unconfirmed",
        colour = "Focus cluster")
    save_plot(p4, "Fig_126_focus_clusters_Harmony_UMAP_v01.png", 12, 8)
  } else {
    logline("Harmony UMAP reduction name not found; omitted UMAP figure.")
  }

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT_DIR, "10_126_REPORT_v01.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "11_126_sessionInfo_v01.txt"), useBytes = TRUE)
  write_csv(data.frame(path = c(IN_RDS, IN_MEM, IN_MARKER, IN_FM, IN_EVID, IN_DEC), source = c("Step123",
    "Step123", "Step123", "Step121", "Step125", "Step125"), used_as = c("corrected expression", "cell identities",
    "marker mapping reference", "full mapping", "annotation evidence", "read-only decisions")), "12_126_provenance_v01.csv")
  logline("SUCCESS | Step126 v01 | targeted evidence generated; no labels changed")
  invisible(TRUE)
}

tryCatch(run_step126(), error = function(e) {
  logline("ERROR | Step126 v01 | ", conditionMessage(e))
  stop(e)
})

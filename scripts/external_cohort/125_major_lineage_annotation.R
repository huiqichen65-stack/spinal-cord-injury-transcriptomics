options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- "D:/zz"
IN_DIR <- file.path(ROOT, "results", "123_GSE281205_major_lineage_clustering")
OUT_DIR <- file.path(ROOT, "results", "125_GSE281205_major_lineage_annotation")
FIG_DIR <- file.path(OUT_DIR, "figures")
RDS_DIR <- file.path(OUT_DIR, "RDS")
LOG_DIR <- file.path(OUT_DIR, "logs")
VERSION <- "v01"
MAIN_CLUSTER <- "cluster_harmony_r04"
EXPECTED_GSM <- paste0("GSM", 8614662:8614667)
EXPECTED_N_CELLS <- 42361L

SAVE_FULL_ANNOTATED_SEURAT <- FALSE
RDS_COMPRESS <- FALSE

IN_MEM <- file.path(IN_DIR, "03_123_per_cell_cluster_membership_v03.csv")
IN_BY_GSM <- file.path(IN_DIR, "04_123_sample_by_cluster_counts_v03.csv")
IN_QC <- file.path(IN_DIR, "06_123_cluster_QC_contamination_v03.csv")
IN_MARKER <- file.path(IN_DIR, "08_123_main_cluster_marker_evidence_v03.csv")
IN_SCORES <- file.path(IN_DIR, "09_123_provisional_major_lineage_scores_v03.csv")
IN_OLD_REVIEW <- file.path(IN_DIR, "10_123_manual_lineage_review_template_v03.csv")
IN_TOP <- file.path(IN_DIR, "13a_123_top25_exploratory_cluster_markers_v03.csv")
IN_RDS <- file.path(IN_DIR, "RDS", "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
IN_124 <- file.path(ROOT, "results", "124_GSE281205_lineage_review_myeloid_extraction", "02_124_cluster_lineage_review_v01.csv")
DECISIONS_PATH <- file.path(OUT_DIR, "04_125_manual_decisions_v01.csv")
LOG_FILE <- file.path(LOG_DIR, "00_125_run_log_v01.log")

fail_unless <- function(ok, message) {
  if (!isTRUE(ok))
    stop(message, call. = FALSE)
}
read_csv <- function(p) {
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}
for (d in c(OUT_DIR, FIG_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
for (p in c(IN_MEM, IN_BY_GSM, IN_QC, IN_MARKER, IN_SCORES, IN_OLD_REVIEW, IN_TOP, IN_RDS)) {
  fail_unless(file.exists(p), paste0("Missing Step123 file: ", p))
}
fail_unless(requireNamespace("ggplot2", quietly = TRUE), "Install ggplot2 first: install.packages('ggplot2')")
cat("", file = LOG_FILE)
log_line <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n", sep = "")
  cat(s, "\n", sep = "", file = LOG_FILE, append = TRUE)
  invisible(s)
}
write_csv <- function(x, basename) {
  fail_unless(is.data.frame(x), paste("Not a data.frame:", basename))
  utils::write.csv(x, file.path(OUT_DIR, basename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
write_plot <- function(p, basename, w = 12, h = 7) {
  ggplot2::ggsave(file.path(FIG_DIR, basename), plot = p, width = w, height = h, dpi = 160, limitsize = FALSE)
}
first_n <- function(x, n = 12L) {
  x <- as.character(x)
  x <- unique(x[!is.na(x) & nzchar(x)])
  paste(utils::head(x, n), collapse = "; ")
}

DRAFT <- data.frame(cluster = as.character(0:23), draft_lineage = c("Myeloid_ambiguous", "Microglia",
  "Microglia", "Myeloid_ambiguous", "Ependymal", "Endothelial", "Monocyte_macrophage", "T_cell", "Unresolved_low_quality",
  "OPC", "Pericyte_vSMC", "Neuron", "NK_cytotoxic", "Vascular_ambiguous", "Fibroblast_meningeal", "Neuron",
  "Ependymal", "Astrocyte", "APC_dendritic", "Unresolved_low_quality", "T_cell", "Neuron", "Microglia",
  "Oligodendrocyte"), review_focus = c("Resident microglial vs infiltrating macrophage ancestry; Clec7a is not sufficient",
  "P2ry12/Tmem119/Sall1/Hexb/Fcrls, verify within-GSM distribution", "IFN-responsive features; inspect microglial identity genes and contrasts",
  "Lpl/Clec7a/Olr1/Acod1; resident vs infiltrating myeloid unclear", "Ciliary/ependymal identity; confirm Foxj1 family",
  "Cldn5/Pecam1/Kdr/Tie1 endothelial evidence", "Mrc1/Pf4/F13a1/H2-Aa macrophage/APC characteristics",
  "Cd3d/Cd3e/Trac T-cell evidence", "Very high contamination and low UMIs; inspect whether a real lineage persists",
  "Pdgfra/Cspg4/Olig1/Olig2 OPC identity", "Vascular mural markers plus high estimated contamination; recheck",
  "Neuronal synaptic genes, rare lineage assessment", "NK cytotoxic genes; distinguish from T cells",
  "Endothelial and mural signatures overlap; do not force one lineage", "Collagen/Dcn/Lum stromal signals; GSM-specific pattern",
  "Neuronal/neuropeptide markers; confirm", "Ependymal-like subgroup; verify marker coherence", "Aqp4/Slc1a3/Aldh1l1 astrocyte identity",
  "Cd209a/Ccl17/Flt3/Cd74 APC/DC-like; not automatically neutrophils", "Low depth/contamination; unresolved until manually adjudicated",
  "Trdc/Trgc family T-cell evidence", "Pkd2l1 neuronal or CSF-contacting subgroup", "Pclaf/Birc5/Ube2c cycling cells; microglial origin requires verification",
  "Mog/Plp1/Mbp/Cnp oligodendrocyte identity"), stringsAsFactors = FALSE)

ALLOWED_LINEAGES <- c("Microglia", "Monocyte_macrophage", "APC_dendritic", "Neutrophil", "T_cell", "NK_cytotoxic",
  "B_cell", "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial", "Pericyte_vSMC",
  "Fibroblast_meningeal", "Myeloid_ambiguous", "Vascular_ambiguous", "Unresolved_low_quality", "Unresolved")

run <- function() {
  mem <- read_csv(IN_MEM)
  qc <- read_csv(IN_QC)
  markers <- read_csv(IN_MARKER)
  scores <- read_csv(IN_SCORES)
  old_review <- read_csv(IN_OLD_REVIEW)
  top <- read_csv(IN_TOP)
  by_gsm <- read_csv(IN_BY_GSM)
  fail_unless(all(c("cell_id", "gsm", "dpi", "replicate", "contamination_percent", MAIN_CLUSTER, "cluster_harmony_r06",
    "cluster_harmony_r08", "cluster_pca_r04") %in% names(mem)), "Step123 per-cell membership columns incorrect")
  fail_unless(all(c("cluster", "n_cells", "pct_contamination_ge50", "median_nCount_RNA", "median_nFeature_RNA") %in%
    names(qc)), "Step123 cluster QC columns incorrect")
  fail_unless(all(c("cluster", "gene_symbol", "mean_log_normalized_expr", "pct_cells_detected", "role") %in%
    names(markers)), "Step123 marker evidence columns incorrect")
  fail_unless(all(c("cluster", "lineage", "mean_z_across_markers") %in% names(scores)), "Step123 score table columns incorrect")
  fail_unless(all(c("cluster", "suggested_lineage") %in% names(old_review)), "Step123 provisional labels columns incorrect")
  fail_unless(all(c("cluster", "gene_symbol", "avg_log2FC") %in% names(top)), "Step123 exploratory cluster marker columns incorrect")
  fail_unless(all(c("dpi", "gsm", "cluster", "n_cells") %in% names(by_gsm)), "Step123 per-sample cluster table columns incorrect")
  fail_unless(nrow(mem) == EXPECTED_N_CELLS, "Step123 per-cell count differs from 42,361: STOP instead of guessing")
  fail_unless(!anyDuplicated(mem$cell_id) && !anyNA(mem$cell_id), "Duplicate/missing Step123 cell IDs")
  fail_unless(setequal(as.character(mem$gsm), EXPECTED_GSM), "Expected six GSM libraries not present")
  fail_unless(all(as.integer(mem$dpi) %in% c(14L, 28L)), "Unexpected stages in Step123 cells")
  for (g in EXPECTED_GSM) {
    ds <- unique(as.integer(mem$dpi[mem$gsm == g]))
    exp_dpi <- if (g %in% EXPECTED_GSM[1:3])
      14L else 28L
    fail_unless(length(ds) == 1L && identical(ds, exp_dpi), paste("GSM/dpi mismatch:", g))
  }
  cl <- as.character(mem[[MAIN_CLUSTER]])
  clusters <- as.character(sort(unique(as.integer(cl))))
  fail_unless(identical(clusters, DRAFT$cluster), "Cluster IDs differ from 0-23: do not reuse the Step123/124 draft")
  fail_unless(!anyDuplicated(qc$cluster) && setequal(as.character(qc$cluster), clusters), "QC cluster mismatch")
  fail_unless(setequal(as.character(old_review$cluster), clusters), "Step123 draft cluster set mismatch")
  fail_unless(sum(qc$n_cells) == EXPECTED_N_CELLS, "QC total does not agree with 42,361 cells")
  n_mem <- table(cl)
  fail_unless(all(as.integer(n_mem[as.character(qc$cluster)]) == as.integer(qc$n_cells)), "Step123 QC cluster counts do not match per-cell labels")
  fail_unless(!anyNA(mem$contamination_percent), "Missing contamination estimate: cannot perform contamination audit")
  fail_unless(sum(by_gsm$n_cells) == EXPECTED_N_CELLS, "GSM cluster counts do not sum to 42,361")
  cross <- table(mem$gsm, cl)
  for (i in seq_len(nrow(by_gsm))) {
    g <- as.character(by_gsm$gsm[i])
    k <- as.character(by_gsm$cluster[i])
    fail_unless(as.integer(cross[g, k]) == as.integer(by_gsm$n_cells[i]), paste("GSM/cluster mismatch:",
      g, k))
  }

  cluster_evidence <- do.call(rbind, lapply(clusters, function(k) {
    q <- qc[as.character(qc$cluster) == k, , drop = FALSE]
    a <- DRAFT[DRAFT$cluster == k, , drop = FALSE]
    old <- old_review[as.character(old_review$cluster) == k, , drop = FALSE]
    s <- scores[as.character(scores$cluster) == k, , drop = FALSE]
    s <- s[order(-s$mean_z_across_markers), , drop = FALSE]
    t <- top[as.character(top$cluster) == k, , drop = FALSE]
    t <- t[order(-t$avg_log2FC), , drop = FALSE]
    m <- markers[as.character(markers$cluster) == k, , drop = FALSE]
    ident_m <- m[m$role == "broad_lineage", , drop = FALSE]
    ident_m <- ident_m[order(-ident_m$pct_cells_detected, -ident_m$mean_log_normalized_expr), , drop = FALSE]
    g <- by_gsm[as.character(by_gsm$cluster) == k, , drop = FALSE]
    g <- g[order(-g$n_cells), , drop = FALSE]
    stage_n <- tapply(mem[[MAIN_CLUSTER]] == k, as.integer(mem$dpi), sum)
    flag <- c()
    if (k %in% c("0", "3"))
      flag <- c(flag, "resident_vs_infiltrating_myeloid")
    if (k %in% c("8", "10", "19"))
      flag <- c(flag, "ambient_RNA_low_depth_review")
    if (k == "13")
      flag <- c(flag, "endothelial_vs_mural_mixed")
    if (k == "18")
      flag <- c(flag, "APC_DC_vs_neutrophil_recheck")
    if (k == "22")
      flag <- c(flag, "cycling_phenotype_identity")
    if (q$pct_contamination_ge50 >= 20)
      flag <- c(flag, "20pct_or_more_cells_ge50pct_DecontX")
    if (q$largest_gsm_fraction >= 50)
      flag <- c(flag, "one_GSM_dominant")
    data.frame(cluster = k, n_cells = as.integer(q$n_cells), draft_lineage = a$draft_lineage, step123_automatic_lineage = old$suggested_lineage,
      annotation_priority = if (length(flag))
        "HIGH" else "STANDARD", review_flags = paste(flag, collapse = "; "), rationale_for_manual_review = a$review_focus,
      top12_enriched_genes = first_n(t$gene_symbol, 12L), top12_detected_canonical = first_n(ident_m$gene_symbol,
        12L), top_lineage_panel = as.character(s$lineage[1]), top_lineage_score = as.numeric(s$mean_z_across_markers[1]),
      runner_up_lineage_panel = if (nrow(s) > 1L)
        as.character(s$lineage[2]) else "", runner_up_score = if (nrow(s) > 1L)
        as.numeric(s$mean_z_across_markers[2]) else NA_real_, median_nFeature_RNA = as.numeric(q$median_nFeature_RNA), median_nCount_RNA = as.numeric(q$median_nCount_RNA),
      median_contamination_percent = as.numeric(q$median_contamination_percent), pct_contamination_ge50 = as.numeric(q$pct_contamination_ge50),
      gsm_with_most_cells = as.character(g$gsm[1]), pct_from_largest_gsm = as.numeric(100 * g$n_cells[1]/q$n_cells),
      n_14dpi = as.integer(stage_n["14"]), n_28dpi = as.integer(stage_n["28"]), stringsAsFactors = FALSE)
  }))
  cluster_evidence <- cluster_evidence[order(as.integer(cluster_evidence$cluster)), , drop = FALSE]
  write_csv(cluster_evidence, "02_125_24cluster_evidence_summary_v01.csv")

  write_csv(markers, "03_125_cluster_canonical_marker_matrix_v01.csv")
  write_csv(by_gsm, "05_125_GSM_by_cluster_audit_v01.csv")
  high <- cluster_evidence[cluster_evidence$annotation_priority == "HIGH", , drop = FALSE]
  write_csv(high, "06_125_priority_manual_review_clusters_v01.csv")
  file_audit <- data.frame(filepath = c(IN_MEM, IN_BY_GSM, IN_QC, IN_MARKER, IN_SCORES, IN_OLD_REVIEW,
    IN_TOP, IN_RDS, IN_124), exists = file.exists(c(IN_MEM, IN_BY_GSM, IN_QC, IN_MARKER, IN_SCORES,
    IN_OLD_REVIEW, IN_TOP, IN_RDS, IN_124)), purpose = c("cell membership", "GSM cluster counts",
    "cluster QC", "canonical marker expression", "lineage panel scores", "automatic unconfirmed suggestions",
    "cluster-enriched top markers", "Seurat object for reviewed annotation UMAP", "optional provisional 124 reference"))
  write_csv(file_audit, "01_125_input_provenance_v01.csv")

  mk <- markers[markers$role == "broad_lineage", , drop = FALSE]
  show_genes <- c("P2ry12", "Tmem119", "Sall1", "Ccr2", "Ly6c2", "Ms4a7", "Flt3", "H2-Aa", "S100a8",
    "S100a9", "Cd3e", "Nkg7", "Aqp4", "Mbp", "Pdgfra", "Snap25", "Foxj1", "Pecam1", "Cldn5", "Rgs5",
    "Pdgfrb", "Col1a1")
  mk <- mk[mk$gene_symbol %in% show_genes, , drop = FALSE]
  mk$cluster <- factor(as.character(mk$cluster), levels = rev(clusters))
  mk$gene_symbol <- factor(mk$gene_symbol, levels = show_genes[show_genes %in% mk$gene_symbol])
  p1 <- ggplot2::ggplot(mk, ggplot2::aes(x = gene_symbol, y = cluster, size = pct_cells_detected, colour = mean_log_normalized_expr)) +
    ggplot2::geom_point(alpha = 0.85) + ggplot2::scale_size(range = c(0.25, 5)) + ggplot2::labs(title = "Whole-atlas canonical marker evidence (unreviewed)",
    subtitle = "Step123 DecontX-normalized expression; cell-level descriptive evidence", x = "Canonical marker",
    y = "Harmony cluster (r=0.4)", size = "% detected", colour = "Mean log expression") + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  write_plot(p1, "Fig_125_canonical_marker_dotplot_v01.png", 15, 9)

  s <- scores
  s$cluster <- factor(as.character(s$cluster), levels = rev(clusters))
  s$lineage <- factor(s$lineage, levels = unique(scores$lineage))
  p2 <- ggplot2::ggplot(s, ggplot2::aes(x = lineage, y = cluster, fill = mean_z_across_markers)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) + ggplot2::labs(title = "Independent broad-lineage marker scores",
    subtitle = "Relative z-scores are NOT probabilities or confirmed lineages", x = "Lineage candidate",
    y = "Harmony cluster (r=0.4)", fill = "Panel z-score") + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 50, hjust = 1))
  write_plot(p2, "Fig_125_lineage_score_heatmap_v01.png", 13, 9)

  if (!file.exists(DECISIONS_PATH)) {
    template <- cluster_evidence[, c("cluster", "n_cells", "draft_lineage", "annotation_priority",
      "review_flags", "top12_enriched_genes", "top12_detected_canonical", "rationale_for_manual_review"),
      drop = FALSE]
    template$manual_lineage <- ""
    template$review_status <- "PENDING"
    template$myeloid_gate <- ""
    template$validation_notes <- ""
    write_csv(template, "04_125_manual_decisions_v01.csv")
  } else {
  }

  dec <- read_csv(DECISIONS_PATH)
  must <- c("cluster", "manual_lineage", "review_status", "myeloid_gate", "validation_notes")
  fail_unless(all(must %in% names(dec)), "Decision file missing required columns; preserve the template layout")
  fail_unless(!anyDuplicated(dec$cluster) && setequal(as.character(dec$cluster), clusters), "Decision file must contain exactly one row for each of 24 clusters")
  dec <- dec[match(clusters, as.character(dec$cluster)), , drop = FALSE]
  dec$review_status <- trimws(as.character(dec$review_status))
  dec$manual_lineage <- trimws(as.character(dec$manual_lineage))
  dec$myeloid_gate <- trimws(as.character(dec$myeloid_gate))
  dec$validation_notes <- trimws(as.character(dec$validation_notes))
  pending <- is.na(dec$review_status) | dec$review_status == "" | dec$review_status == "PENDING"

  if (any(pending)) {
    preview <- data.frame(cluster = clusters, n_cells = cluster_evidence$n_cells, draft_lineage = cluster_evidence$draft_lineage,
      manual_lineage = dec$manual_lineage, review_status = dec$review_status, myeloid_gate = dec$myeloid_gate)
    write_csv(preview, "07_125_review_progress_v01.csv")
    txt <- c("STEP125 v01 | GSE281205 | MANUAL REVIEW REQUIRED", "", paste("Cells inspected in Step123:",
      nrow(mem)), paste("Clusters:", length(clusters)), paste("Reviewed clusters with explicit decision:",
      sum(!pending)), paste("Pending clusters:", sum(pending)), "", "This is a scientific review gate, not a runtime failure.",
      "Edit 04_125_manual_decisions_v01.csv IN PLACE:", " - Enter manual_lineage (one ALLOWED_LINEAGES value; listed in script).",
      " - Enter review_status: CONFIRMED or UNRESOLVED.", " - Enter myeloid_gate: INCLUDE, EXCLUDE or HOLD.",
      " - Enter validation_notes citing marker evidence (especially ambiguous clusters).", " - For UNRESOLVED use manual_lineage=Unresolved, myeloid_gate=HOLD.",
      " - Do not use the 1320-gene candidate signature to choose the labels.", " - Rerun the same Step125 script after editing; it will then apply labels.",
      " - None of the Step123/124 automatic suggestions count as confirmation.", "WARNING: 28dpi GEO title/treatment age field mismatch remains.",
      "WARNING: DecontX high-contamination groups need manual inspection.", "WARNING: Harmony GSM nested within dpi may attenuate stage effects.",
      "", "ALLOWED_LINEAGES:", paste(ALLOWED_LINEAGES, collapse = ", "), "", "No myeloid subclustering or pseudobulk is performed.")
    writeLines(txt, file.path(OUT_DIR, "09_125_REPORT_v01.txt"), useBytes = TRUE)
    writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "10_125_sessionInfo_v01.txt"), useBytes = TRUE)
    return(invisible(TRUE))
  }

  fail_unless(all(dec$review_status %in% c("CONFIRMED", "UNRESOLVED")), "Only CONFIRMED / UNRESOLVED accepted (or PENDING during preparation)")
  fail_unless(!anyNA(dec$manual_lineage) && all(dec$manual_lineage %in% ALLOWED_LINEAGES), paste("Invalid manual_lineage. Allowed:",
    paste(ALLOWED_LINEAGES, collapse = ", ")))
  fail_unless(all(dec$myeloid_gate %in% c("INCLUDE", "EXCLUDE", "HOLD")), "myeloid_gate must be INCLUDE, EXCLUDE or HOLD in each row")
  fail_unless(all(nzchar(dec$validation_notes) & !is.na(dec$validation_notes)), "Every manual decision needs nonempty validation_notes")
  fail_unless(all(dec$manual_lineage[dec$review_status == "UNRESOLVED"] == "Unresolved") && all(dec$myeloid_gate[dec$review_status ==
    "UNRESOLVED"] == "HOLD"), "UNRESOLVED rows require lineage=Unresolved and gate=HOLD")
  fail_unless(all(dec$manual_lineage[dec$review_status == "CONFIRMED"] != "Unresolved"), "CONFIRMED rows cannot use Unresolved lineage")
  fail_unless(all(dec$myeloid_gate[dec$manual_lineage %in% c("T_cell", "NK_cytotoxic", "B_cell", "Astrocyte",
    "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial", "Pericyte_vSMC", "Fibroblast_meningeal",
    "Vascular_ambiguous")] != "INCLUDE"), "Cannot include known non-myeloid lineage in myeloid gate")

  dec$validation_confirmed <- dec$review_status == "CONFIRMED"
  write_csv(dec, "08_125_reviewed_cluster_annotation_v01.csv")
  idx <- match(cl, as.character(dec$cluster))
  fail_unless(!anyNA(idx), "Cell cluster not covered by reviewed decisions")
  annotated <- data.frame(cell_id = as.character(mem$cell_id), gsm = as.character(mem$gsm), dpi = as.integer(mem$dpi),
    replicate = as.integer(mem$replicate), cluster_harmony_r04 = cl, cluster_harmony_r06 = as.character(mem$cluster_harmony_r06),
    cluster_harmony_r08 = as.character(mem$cluster_harmony_r08), cluster_pca_r04 = as.character(mem$cluster_pca_r04),
    contamination_percent = as.numeric(mem$contamination_percent), reviewed_major_lineage_125 = dec$manual_lineage[idx],
    review_status_125 = dec$review_status[idx], myeloid_gate_125 = dec$myeloid_gate[idx], stringsAsFactors = FALSE)
  fail_unless(nrow(annotated) == EXPECTED_N_CELLS && !anyDuplicated(annotated$cell_id), "Reviewed annotations do not cover 42,361 distinct cells")
  write_csv(annotated, "11_125_per_cell_reviewed_lineage_v01.csv")
  cnt <- as.data.frame(table(annotated$gsm, annotated$reviewed_major_lineage_125), stringsAsFactors = FALSE)
  names(cnt) <- c("gsm", "reviewed_major_lineage", "n_cells")
  cnt <- cnt[cnt$n_cells > 0, , drop = FALSE]
  cnt$dpi <- ifelse(cnt$gsm %in% EXPECTED_GSM[1:3], 14L, 28L)
  write_csv(cnt[, c("dpi", "gsm", "reviewed_major_lineage", "n_cells")], "12_125_GSM_reviewed_lineage_counts_v01.csv")
  stage <- as.data.frame(table(dpi = factor(annotated$dpi, levels = c(14, 28)), lineage = annotated$reviewed_major_lineage_125),
    stringsAsFactors = FALSE)
  names(stage) <- c("dpi", "reviewed_major_lineage", "n_cells")
  stage <- stage[stage$n_cells > 0, , drop = FALSE]
  stage$pct_within_stage <- 100 * stage$n_cells/as.numeric(table(annotated$dpi)[as.character(stage$dpi)])
  write_csv(stage, "13_125_stage_reviewed_lineage_composition_v01.csv")
  included <- annotated[annotated$myeloid_gate_125 == "INCLUDE", , drop = FALSE]
  held <- annotated[annotated$myeloid_gate_125 == "HOLD", , drop = FALSE]
  write_csv(included[, c("cell_id", "gsm", "dpi", "cluster_harmony_r04", "reviewed_major_lineage_125")],
    "14_125_myeloid_included_cell_ids_v01.csv")
  write_csv(held[, c("cell_id", "gsm", "dpi", "cluster_harmony_r04", "reviewed_major_lineage_125")],
    "15_125_held_uncertain_cell_ids_v01.csv")
  saveRDS(annotated, file.path(RDS_DIR, "125_GSE281205_reviewed_lineage_metadata_v01.rds"), compress = "gzip")

  fail_unless(requireNamespace("SeuratObject", quietly = TRUE), "Install SeuratObject first: install.packages('SeuratObject')")
  seu <- readRDS(IN_RDS)
  fail_unless(inherits(seu, "Seurat"), "Input Step123 RDS is not Seurat")
  md <- seu[[]]
  fail_unless(setequal(colnames(seu), annotated$cell_id), "Step123 Seurat cell IDs do not match 125 reviewed cell IDs")
  id <- match(colnames(seu), annotated$cell_id)
  fail_unless(!anyNA(id), "Missing Seurat cell annotation after ID mapping")
  fail_unless(all(as.character(md[[MAIN_CLUSTER]]) == annotated$cluster_harmony_r04[id]), "Step123 Seurat cluster identity differs from CSV evidence")
  seu$reviewed_major_lineage_125 <- annotated$reviewed_major_lineage_125[id]
  seu$review_status_125 <- annotated$review_status_125[id]
  seu$myeloid_gate_125 <- annotated$myeloid_gate_125[id]
  seu@misc$step125 <- list(origin = "manually adjudicated Step123 major lineages", decisions_file = DECISIONS_PATH,
    date_applied = as.character(Sys.Date()), resolution = 0.4, signature_1320_used = FALSE, pseudobulk_performed = FALSE,
    incomplete_broad_lineages = sum(dec$review_status == "UNRESOLVED"), age_metadata_warning = "GEO 28dpi young/old annotation conflict",
    original_integer_UMI = "Step119 per-GSM RDS; not this DecontX assay")
  if ("umap_harmony" %in% SeuratObject::Reductions(seu)) {
    e <- SeuratObject::Embeddings(seu[["umap_harmony"]])
    e <- e[colnames(seu), seq_len(2), drop = FALSE]
    dat <- data.frame(UMAP_1 = e[, 1], UMAP_2 = e[, 2], label = seu$reviewed_major_lineage_125, review_status = seu$review_status_125,
      gsm = as.character(md$gsm))
    p3 <- ggplot2::ggplot(dat, ggplot2::aes(x = UMAP_1, y = UMAP_2, colour = label)) + ggplot2::geom_point(size = 0.17,
      alpha = 0.45) + ggplot2::theme_bw(base_size = 11) + ggplot2::labs(title = "Reviewed major lineages on Step123 Harmony UMAP",
      subtitle = "All 42,361 singlets; labels explicitly reviewed; not state-level annotation",
      x = "Harmony UMAP 1", y = "Harmony UMAP 2", colour = "Lineage")
    write_plot(p3, "Fig_125_reviewed_major_lineage_UMAP_v01.png", 13, 8)
  }
  if (SAVE_FULL_ANNOTATED_SEURAT) {
    outfile <- file.path(RDS_DIR, "125_GSE281205_reviewed_major_lineage_Seurat_v01.rds")
    fail_unless(!file.exists(outfile), paste("Will not overwrite existing full RDS:", outfile))
    partial <- paste0(outfile, ".partial")
    if (file.exists(partial))
      unlink(partial)
    saveRDS(seu, partial, compress = RDS_COMPRESS)
    fail_unless(file.exists(partial) && file.info(partial)$size > 1e+06, "Seurat RDS save incomplete")
    fail_unless(isTRUE(file.rename(partial, outfile)), "Could not finalize optional Seurat RDS")
  }
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, "09_125_REPORT_v01.txt"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "10_125_sessionInfo_v01.txt"), useBytes = TRUE)
  invisible(TRUE)
}

tryCatch(run(), error = function(e) {
  log_line("ERROR | Step125 v01 | ", conditionMessage(e))
  stop(e)
})

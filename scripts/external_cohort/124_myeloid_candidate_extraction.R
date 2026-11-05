options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- "D:/zz"
OUT <- file.path(ROOT, "results", "124_GSE281205_lineage_review_myeloid_extraction")
SRC123 <- file.path(ROOT, "results", "123_GSE281205_major_lineage_clustering")
SRC121 <- file.path(ROOT, "results", "121_GSE281205_normalize_PCA_batch_audit")
SRC119 <- file.path(ROOT, "results", "119_GSE281205_scDblFinder")
IN_RDS <- file.path(SRC123, "RDS", "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
IN_MEM <- file.path(SRC123, "03_123_per_cell_cluster_membership_v03.csv")
IN_QC <- file.path(SRC123, "06_123_cluster_QC_contamination_v03.csv")
IN_MARKER <- file.path(SRC123, "08_123_main_cluster_marker_evidence_v03.csv")
IN_REVIEW <- file.path(SRC123, "10_123_manual_lineage_review_template_v03.csv")
IN_TOP <- file.path(SRC123, "13a_123_top25_exploratory_cluster_markers_v03.csv")
IN_FEATURE <- file.path(SRC121, "02_121_gene_feature_mapping_v02.csv")
FIG_DIR <- file.path(OUT, "figures")
RDS_DIR <- file.path(OUT, "RDS")
LOG_DIR <- file.path(OUT, "logs")
VERSION <- "v01"
MAIN_CLUSTER <- "cluster_harmony_r04"

SAVE_EXTENDED_CANDIDATE_SEURAT <- TRUE
RDS_COMPRESS <- FALSE

MANUAL_DECISIONS <- file.path(OUT, "03_124_manual_review_decisions_v01.csv")
EXPECTED_GSM <- paste0("GSM", 8614662:8614667)

stopif <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
for (f in c(IN_RDS, IN_MEM, IN_QC, IN_MARKER, IN_REVIEW, IN_TOP, IN_FEATURE)) {
  stopif(file.exists(f), paste("Missing required direct input:", f))
}
for (pkg in c("Seurat", "SeuratObject", "Matrix", "ggplot2")) {
    stopif(requireNamespace(pkg, quietly = TRUE), paste("Missing R package:", pkg, "(install before running Step124)"))
}
for (d in c(OUT, FIG_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOGFILE <- file.path(OUT, "00_124_run_log_v01.log")
cat("", file = LOGFILE)
log_line <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOGFILE, append = TRUE)
  invisible(msg)
}
write_csv <- function(df, fname) {
    stopif(is.data.frame(df), paste("Not a data.frame:", fname))
    utils::write.csv(df, file.path(OUT, fname), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
read_tab <- function(path) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
plot_save <- function(p, filename, width = 12, height = 7) {
    ggplot2::ggsave(file.path(FIG_DIR, filename), p, width = width, height = height, dpi = 160, limitsize = FALSE)
}
median0 <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}

assignment <- data.frame(
  cluster = as.character(0:23),
  proposed_lineage = c(
    "Myeloid_ambiguous",
    "Microglia",
    "Microglia",
    "Myeloid_ambiguous",
    "Ependymal",
    "Endothelial",
    "Monocyte_macrophage",
    "T_cell",
    "Unresolved_low_quality",
    "OPC",
    "Pericyte_vSMC",
    "Neuron",
    "NK_cytotoxic",
    "Pericyte_vSMC",
    "Fibroblast_meningeal",
    "Neuron",
    "Ependymal",
    "Astrocyte",
    "APC_dendritic",
    "Unresolved_low_quality",
    "T_cell",
    "Neuron",
    "Microglia",
    "Oligodendrocyte"
  ),
  inclusion_tier = c(
    "extended", "core", "core", "extended", "exclude", "exclude", "core",
    "exclude", "exclude", "exclude", "exclude", "exclude", "exclude",
    "exclude", "exclude", "exclude", "exclude", "exclude", "extended",
    "exclude", "exclude", "exclude", "core", "exclude"
  ),
  review_rationale = c(
    "Clec7a/Apoe/Lpl/Gpnmb enriched; partial resident markers; microglia vs macrophage unresolved",
    "P2ry12/Tmem119/Sall1 strongly support microglial identity",
    "Microglial markers with interferon/chemokine-responsive genes",
    "Lpl/Clec7a/Acod1/Olr1; activated myeloid, resident vs recruited uncertain",
    "Ependymal/cilia genes such as Dnah11 and Spag17",
    "Endothelial markers including Cldn5 and Tie1",
    "F13a1/Mrc1/Pf4/Cd74, macrophage-like identity",
    "Cd3d/Cd3e/Trac, T-cell lineage",
    "Low-count and high-background cluster; lineage unresolved (71.6 percent contamination>=50)",
    "Pdgfra/Olig1/Olig2/Gpr17, OPC lineage",
    "Myh11/Myl9/Rgs5; vascular mural lineage, high estimated contamination",
    "Neuron-related synaptic and neuropeptide markers",
    "Ncr1/Gzmc/Klrb1c; NK cytotoxic lineage",
    "Rgs5/Notch3/Abcc9; vascular mural markers; endothelial overlap merits review",
    "Col1a1/Dpt/Mfap5; fibroblast/meningeal lineage",
    "Neurod2/Gal/other neuronal markers",
    "Ependymal-like genes; requires separate confirmation",
    "Gfap/Aqp4/Slc7a10 astrocyte lineage",
    "Cd209a/Ccl17/Clec10a suggest dendritic/APC lineage; previously miscalled neutrophil",
    "Unresolved/low-depth with substantial estimated ambient contamination",
    "Trdc/Trgc1/Il17a gamma-delta-like T lineage",
    "Pkd2l1/neuronal markers; likely CSF-contacting neuronal subgroup",
    "Pclaf/Birc5/Ube2c cell cycle, with microglial identity markers",
    "Mog/Plp1 oligodendrocyte lineage"
  ), stringsAsFactors = FALSE
)
assignment$review_status <- "PROVISIONAL_NOT_MANUALLY_CONFIRMED"
assignment$manual_notes <- ""
stopif(identical(assignment$cluster, as.character(0:23)), "Internal cluster decision table must cover clusters 0..23")
stopif(all(assignment$inclusion_tier %in% c("core", "extended", "exclude")), "Invalid inclusion tier in curated decisions")

run <- function() {

  membership <- read_tab(IN_MEM)
  qc <- read_tab(IN_QC)
  marker <- read_tab(IN_MARKER)
  prevreview <- read_tab(IN_REVIEW)
  top <- read_tab(IN_TOP)
  feature <- read_tab(IN_FEATURE)
  needmem <- c("cell_id", "gsm", "dpi", "replicate", "contamination_percent",
               "cluster_harmony_r04", "cluster_harmony_r06", "cluster_harmony_r08", "cluster_pca_r04")
  stopif(all(needmem %in% names(membership)), "Step123 membership columns missing")
  stopif(all(c("cluster", "n_cells", "pct_contamination_ge50") %in% names(qc)),
         "Step123 QC summary columns missing")
  stopif(all(c("cluster", "gene_symbol", "pct_cells_detected") %in% names(marker)),
         "Step123 marker columns missing")
  stopif(all(c("cluster", "suggested_lineage") %in% names(prevreview)),
         "Step123 review columns missing")
  stopif(all(c("cluster", "gene_symbol") %in% names(top)),
         "Step123 top-marker columns missing")
  stopif(all(c("gene_id", "gene_symbol") %in% names(feature)),
         "Step121 feature mapping columns missing")
  stopif(nrow(membership) == 42361L && !anyDuplicated(membership$cell_id),
         "Step123 cell ID table differs from 42,361 unique singlets")
  stopif(setequal(as.character(membership[[MAIN_CLUSTER]]), assignment$cluster),
         "Harmony r0.4 cluster IDs differ from 0..23; do not reuse assignments")
  stopif(setequal(as.character(membership$gsm), EXPECTED_GSM),
         "Six expected GSM identifiers not present")
  stopif(all(as.integer(membership$dpi) %in% c(14L, 28L)) &&
           all(as.integer(membership$replicate) %in% 1:3),
         "Unexpected dpi or replicate labels")
  for (g in EXPECTED_GSM) {
    s <- unique(membership$dpi[membership$gsm == g])
    expected_stage <- if (g %in% EXPECTED_GSM[1:3]) 14L else 28L
    stopif(length(s) == 1L && as.integer(s[[1]]) == expected_stage,
           paste("Wrong GSM stage metadata:", g))
  }
  stopif(!anyDuplicated(feature$gene_id), "Feature mapping contains duplicate gene IDs")

  if (file.exists(MANUAL_DECISIONS)) {
    decisions <- read_tab(MANUAL_DECISIONS)
    must <- c("cluster", "manual_confirmed", "manual_lineage", "manual_inclusion_tier", "validation_notes")
    stopif(all(must %in% names(decisions)) && !anyDuplicated(decisions$cluster),
           "Manual decisions require columns: cluster,manual_confirmed,manual_lineage,manual_inclusion_tier,validation_notes")
    stopif(setequal(as.character(decisions$cluster), assignment$cluster),
           "Manual decisions must have exactly all 24 clusters")
    decisions <- decisions[match(assignment$cluster, as.character(decisions$cluster)), , drop = FALSE]
    confirmed <- toupper(trimws(as.character(decisions$manual_confirmed))) %in%
      c("TRUE", "T", "YES", "Y", "1")
    goodlineage <- c("Microglia", "Monocyte_macrophage", "APC_dendritic",
                     "Myeloid_ambiguous", "Unresolved_low_quality", "Ependymal",
                     "Endothelial", "OPC", "Pericyte_vSMC", "Neuron", "T_cell",
                     "NK_cytotoxic", "Fibroblast_meningeal", "Astrocyte", "Oligodendrocyte")
    if (any(confirmed)) {
      stopif(all(as.character(decisions$manual_lineage[confirmed]) %in% goodlineage),
             "Confirmed manual_lineage value is not in supported categories")
      stopif(all(as.character(decisions$manual_inclusion_tier[confirmed]) %in%
                   c("core", "extended", "exclude")),
             "Confirmed manual_inclusion_tier must be core/extended/exclude")
      assignment$proposed_lineage[confirmed] <- as.character(decisions$manual_lineage[confirmed])
      assignment$inclusion_tier[confirmed] <- as.character(decisions$manual_inclusion_tier[confirmed])
      assignment$review_status[confirmed] <- "MANUALLY_CONFIRMED_FROM_USER_FILE"
      assignment$manual_notes[confirmed] <- as.character(decisions$validation_notes[confirmed])
    } else {
    }
  } else {
  }

  top$cluster <- as.character(top$cluster)
  top12 <- vapply(assignment$cluster, function(k) {
    v <- as.character(top$gene_symbol[top$cluster == k])
    v <- v[!is.na(v) & nzchar(v)]
    paste(utils::head(v, 12L), collapse = "; ")
  }, character(1))
  assignment$top12_enriched_symbols <- top12
  qc$cluster <- as.character(qc$cluster)
  qc <- qc[match(assignment$cluster, qc$cluster), , drop = FALSE]
  stopif(!anyNA(qc$cluster), "Missing Step123 QC cluster row")
  assignment$n_cells <- as.integer(qc$n_cells)
  assignment$pct_contamination_ge50 <- as.numeric(qc$pct_contamination_ge50)
  assignment$median_nFeature_RNA <- as.numeric(qc$median_nFeature_RNA)
  assignment$median_nCount_RNA <- as.numeric(qc$median_nCount_RNA)
  prevreview$cluster <- as.character(prevreview$cluster)
  assignment$step123_auto_suggestion <- as.character(prevreview$suggested_lineage[
    match(assignment$cluster, prevreview$cluster)])
  stopif(!anyNA(assignment$step123_auto_suggestion),
         "Missing Step123 provisional suggestion")
  interested <- c("P2ry12", "Tmem119", "Sall1", "Clec7a", "Apoe", "Lpl",
                  "Gpnmb", "Ccr2", "Ly6c2", "Flt3", "Itgax", "S100a8", "S100a9",
                  "Cldn5", "Pecam1", "Rgs5", "Pdgfrb")
  marker$cluster <- as.character(marker$cluster)
  evidence <- marker[marker$gene_symbol %in% interested,
                     c("cluster", "gene_symbol", "pct_cells_detected",
                       "mean_log_normalized_expr", "role"), drop = FALSE]
  evidence$review_comment <- "Cluster-level detection is descriptive; not cell-identity proof"
  write_csv(evidence, "01_124_canonical_marker_review_evidence_v01.csv")
  write_csv(assignment, "02_124_cluster_lineage_review_v01.csv")

  template_name <- "03_124_manual_review_template_v01.csv"
  template_file <- file.path(OUT, template_name)
  if (!file.exists(template_file)) {
    temp <- assignment[, c("cluster", "proposed_lineage", "inclusion_tier",
                           "n_cells", "pct_contamination_ge50",
                           "top12_enriched_symbols", "review_rationale"), drop = FALSE]
    names(temp)[1:3] <- c("cluster", "suggested_lineage", "suggested_tier")
    temp$manual_confirmed <- FALSE
    temp$manual_lineage <- ""
    temp$manual_inclusion_tier <- ""
    temp$validation_notes <- ""
    write_csv(temp, template_name)
  } else {
  }

  seu <- readRDS(IN_RDS)
  stopif(inherits(seu, "Seurat"), "Step123 input RDS is not Seurat")
  stopif(ncol(seu) == nrow(membership), "Step123 Seurat ncol != cell manifest")
  stopif(identical(as.character(colnames(seu)), as.character(membership$cell_id)),
         "Step123 cell/order mismatch: cannot safely propagate labels")
  md <- seu[[]]
  stopif(all(c("gsm", "dpi", "replicate", "decontX_contamination_percent",
               MAIN_CLUSTER) %in% colnames(md)),
         "Step123 Seurat lacks necessary provenance columns")
  stopif(identical(as.character(md$gsm), as.character(membership$gsm)) &&
           identical(as.integer(md$dpi), as.integer(membership$dpi)),
         "GSM or stage values disagree between Seurat and CSV")
  stopif(identical(as.character(md[[MAIN_CLUSTER]]),
                   as.character(membership[[MAIN_CLUSTER]])),
         "Cluster assignments differ between RDS and verified CSV")
  stopif(isTRUE(all.equal(as.numeric(md$decontX_contamination_percent),
                            as.numeric(membership$contamination_percent),
                            tolerance = 1e-6, check.attributes = FALSE)),
         "DecontX contamination mismatch between RDS and manifest")

  cl <- as.character(md[[MAIN_CLUSTER]])
  ii <- match(cl, assignment$cluster)
  stopif(!anyNA(ii), "Some cells could not be assigned a reviewed cluster")
  md$lineage_124 <- assignment$proposed_lineage[ii]
  md$myeloid_inclusion_tier_124 <- assignment$inclusion_tier[ii]
  md$lineage_review_status_124 <- assignment$review_status[ii]
  md$is_myeloid_core_124 <- md$myeloid_inclusion_tier_124 == "core"
  md$is_myeloid_extended_124 <- md$myeloid_inclusion_tier_124 %in% c("core", "extended")
  md$high_contam_ge50_124 <- as.numeric(md$decontX_contamination_percent) >= 50
  md$stage_label_124 <- paste0(as.integer(md$dpi), " dpi")
  stopif(!anyNA(md$lineage_124) && !anyNA(md$myeloid_inclusion_tier_124),
         "Missing propagated per-cell labels")

  write_csv(data.frame(cell_id = rownames(md),
      gsm = as.character(md$gsm), dpi = as.integer(md$dpi),
      replicate = as.integer(md$replicate),
      cluster_harmony_r04 = cl,
      lineage_124 = md$lineage_124,
      inclusion_tier_124 = md$myeloid_inclusion_tier_124,
      review_status_124 = md$lineage_review_status_124,
      contamination_percent = as.numeric(md$decontX_contamination_percent),
      high_contam_ge50 = md$high_contam_ge50_124),
    "04_124_per_cell_lineage_manifest_v01.csv")

  core_cells <- rownames(md)[md$is_myeloid_core_124]
  extended_cells <- rownames(md)[md$is_myeloid_extended_124]
  uncertain_high <- rownames(md)[cl %in% c("8", "19") &
                                     !md$is_myeloid_extended_124]
  stopif(length(core_cells) > 0L && length(extended_cells) >= length(core_cells),
         "Empty or inconsistent myeloid cell sets")
  stopif(!anyDuplicated(extended_cells) &&
           all(core_cells %in% extended_cells),
         "Myeloid inclusion sets not nested as expected")
  stopif(length(extended_cells) + sum(!md$is_myeloid_extended_124) == ncol(seu),
         "Myeloid + non-myeloid cell counts fail sum check")
  write_csv(data.frame(cell_id = core_cells),
            "05_124_myeloid_core_cell_ids_v01.csv")
  write_csv(data.frame(cell_id = extended_cells),
            "06_124_myeloid_extended_cell_ids_v01.csv")
  write_csv(data.frame(cell_id = uncertain_high,
      cluster = cl[match(uncertain_high, rownames(md))],
      contamination_percent = as.numeric(md$decontX_contamination_percent)[
        match(uncertain_high, rownames(md))]),
    "07_124_unresolved_high_contamination_review_cells_v01.csv")

  counts <- as.data.frame(table(cluster = cl,
                                inclusion_tier = md$myeloid_inclusion_tier_124),
                          stringsAsFactors = FALSE)
  counts <- counts[counts$Freq > 0, , drop = FALSE]
  names(counts)[3] <- "n_cells"
  write_csv(counts, "08_124_cluster_inclusion_counts_v01.csv")

  gsm_table <- do.call(rbind, lapply(EXPECTED_GSM, function(g) {
    m <- md[as.character(md$gsm) == g, , drop = FALSE]
    data.frame(gsm = g, dpi = as.integer(unique(m$dpi)),
       replicate = as.integer(unique(m$replicate)),
       n_all_singlets = nrow(m),
       n_myeloid_core = sum(m$is_myeloid_core_124),
       n_myeloid_extended = sum(m$is_myeloid_extended_124),
       n_myeloid_extended_only = sum(m$myeloid_inclusion_tier_124 == "extended"),
       n_excluded_from_myeloid = sum(m$myeloid_inclusion_tier_124 == "exclude"),
       n_extended_with_contam_ge50 = sum(m$is_myeloid_extended_124 & m$high_contam_ge50_124),
       pct_myeloid_extended_of_all_singlets = 100 * mean(m$is_myeloid_extended_124),
       original_integer_UMI_exists = file.exists(file.path(SRC119, "RDS",
                                 paste0("119_", g, "_singlet_raw_counts_v02.rds"))),
       age_metadata_note = if (as.integer(unique(m$dpi)) == 28L)
         "GEO title=Young; treatment=Old; BioSample supports Young; verify primary source"
         else "GEO title/treatment consistent")
  }))
  stopif(nrow(gsm_table) == 6L && all(gsm_table$n_myeloid_core <=
                                       gsm_table$n_myeloid_extended),
         "Per-sample count audit failed")
  if (!all(gsm_table$original_integer_UMI_exists)) {
    log_line("WARNING: one or more Step119 raw-UMI files were not found; future pseudobulk BLOCKED")
  }
  write_csv(gsm_table, "09_124_per_sample_myeloid_inclusion_v01.csv")
  stage_table <- do.call(rbind, lapply(c(14L, 28L), function(d) {
    x <- gsm_table[gsm_table$dpi == d, , drop = FALSE]
    data.frame(dpi = d, n_libraries = nrow(x),
      total_singlets = sum(x$n_all_singlets),
      myeloid_core = sum(x$n_myeloid_core),
      myeloid_extended = sum(x$n_myeloid_extended),
      myeloid_extended_only = sum(x$n_myeloid_extended_only),
      pct_extended_of_all = 100 * sum(x$n_myeloid_extended)/sum(x$n_all_singlets),
      median_extended_per_library = median0(x$n_myeloid_extended),
      all_six_have_extended = all(x$n_myeloid_extended > 0L))
  }))
  write_csv(stage_table, "10_124_stage_myeloid_inclusion_v01.csv")
  membership_by_gsm <- as.data.frame(table(
      gsm = as.character(md$gsm),
      cluster = factor(cl, levels = as.character(0:23))),
    stringsAsFactors = FALSE)
  names(membership_by_gsm)[3] <- "n_cells"
  membership_by_gsm <- membership_by_gsm[membership_by_gsm$n_cells > 0, , drop = FALSE]
  membership_by_gsm$inclusion_tier <- assignment$inclusion_tier[
    match(as.character(membership_by_gsm$cluster), assignment$cluster)]
  write_csv(membership_by_gsm, "11_124_GSM_by_cluster_review_v01.csv")

  write_csv(feature, "12_124_gene_feature_mapping_v01.csv")
  write_csv(data.frame(
      file = c(IN_RDS, IN_MEM, IN_QC, IN_MARKER, IN_REVIEW, IN_TOP, IN_FEATURE,
        file.path(SRC119,"RDS",paste0("119_",EXPECTED_GSM,"_singlet_raw_counts_v02.rds"))),
      exists = file.exists(c(IN_RDS, IN_MEM, IN_QC, IN_MARKER, IN_REVIEW, IN_TOP, IN_FEATURE,
        file.path(SRC119,"RDS",paste0("119_",EXPECTED_GSM,"_singlet_raw_counts_v02.rds"))))),
    "13_124_input_provenance_v01.csv")

  dr <- SeuratObject::Reductions(seu)
  if ("umap_harmony" %in% dr) {
    em <- SeuratObject::Embeddings(seu[["umap_harmony"]])
    stopif(all(rownames(md) %in% rownames(em)), "Harmony UMAP cell ID mismatch")
    em <- em[rownames(md), 1:2, drop = FALSE]
    dat <- data.frame(x = em[,1], y = em[,2],
                      lineage = md$lineage_124,
                      tier = md$myeloid_inclusion_tier_124,
                      gsm = md$gsm,
                      contaminated = md$high_contam_ge50_124)
    p1 <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, color = lineage)) +
      ggplot2::geom_point(size = .14, alpha = .45) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = "Provisional whole-atlas lineage review",
        subtitle = "Step123 Harmony UMAP; all labels pending independent manual confirmation",
        x = "Harmony UMAP 1", y = "Harmony UMAP 2")
    plot_save(p1, "Fig_124_provisional_lineage_UMAP_v01.png", 12, 8)
    p2 <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, color = tier)) +
      ggplot2::geom_point(size = .14, alpha = .45) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = "Myeloid extraction tiers (provisional)",
        subtitle = "Core vs expanded vs excluded from proposed myeloid extraction",
        x = "Harmony UMAP 1", y = "Harmony UMAP 2")
    plot_save(p2, "Fig_124_myeloid_inclusion_UMAP_v01.png", 12, 8)
  } else {
    log_line("WARNING: umap_harmony reduction missing; UMAP figures omitted")
  }
  bar <- gsm_table
  bar$gsm <- factor(bar$gsm, levels = EXPECTED_GSM)
  p3 <- ggplot2::ggplot(bar, ggplot2::aes(x = gsm, y = n_myeloid_extended,
                                         fill = factor(dpi))) +
    ggplot2::geom_col() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Candidate myeloid cell counts by library",
      subtitle = "Provisional inclusion, not target-state counts or independent biological validation",
      x = "GSM library", y = "Candidate myeloid singlets", fill = "dpi") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  plot_save(p3, "Fig_124_candidate_myeloid_counts_by_GSM_v01.png", 11, 7)

  saveRDS(md, file.path(RDS_DIR, "124_GSE281205_all_cells_lineage_metadata_v01.rds"),
          compress = "gzip")

  if (SAVE_EXTENDED_CANDIDATE_SEURAT) {
    stopif("DecontX" %in% SeuratObject::Assays(seu),
           "DecontX assay not present in Step123 RDS")
    corrected <- SeuratObject::LayerData(seu[["DecontX"]], layer = "counts")
    stopif(inherits(corrected, "sparseMatrix"),
           "Expected sparse corrected expression matrix; refusing to densify")
    stopif(all(extended_cells %in% colnames(corrected)),
           "Cannot locate all candidate cells in corrected matrix")
    corrected <- corrected[, extended_cells, drop = FALSE]
    stopif(identical(colnames(corrected), extended_cells),
           "Corrected candidate expression cell order mismatch")
    myeloid_md <- md[extended_cells, , drop = FALSE]

    cand <- SeuratObject::CreateSeuratObject(counts = corrected,
      assay = "DecontX", project = "GSE281205_Provisional_Myeloid",
      meta.data = myeloid_md, min.cells = 0L, min.features = 0L)
    stopif(ncol(cand) == length(extended_cells) &&
             identical(colnames(cand), extended_cells),
           "Candidate Seurat creation changed cells or order")
    cand@misc$step124 <- list(
      lineage_status = if (all(assignment$review_status ==
                                 "MANUALLY_CONFIRMED_FROM_USER_FILE"))
                         "all_clusters_manually_confirmed" else
                         "provisional_review_incomplete",
      main_clusters_from_step123 = MAIN_CLUSTER,
      cluster_review_file = file.path(OUT, "02_124_cluster_lineage_review_v01.csv"),
      manual_review_file = MANUAL_DECISIONS,
      extraction_tiers = "core + extended; ambiguity retained; no definitive target microglia",
      input = IN_RDS, corrected_counts_are_fractional_estimates = TRUE,
      raw_integer_UMI_sources = file.path(SRC119, "RDS",
        paste0("119_", EXPECTED_GSM, "_singlet_raw_counts_v02.rds")),
      signature_1320_loaded = FALSE,
      differential_expression_performed = FALSE,
      note = "Recompute cell clustering using corrected expression; use Step119 ORIGINAL INTEGER UMIs for future pseudobulk"
    )
    outfile <- file.path(RDS_DIR,
      "124_GSE281205_extended_myeloid_DecontX_counts_v01.rds")
    tmp <- paste0(outfile, ".partial")
    if (file.exists(tmp)) unlink(tmp)
    saveRDS(cand, tmp, compress = RDS_COMPRESS)
    stopif(file.exists(tmp) && file.info(tmp)$size > 100000,
           "Candidate Seurat output incomplete")
    if (file.exists(outfile)) unlink(outfile)
    stopif(isTRUE(file.rename(tmp, outfile)),
           "Could not finalize candidate RDS (partial retained)")
  }

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "14_124_REPORT_v01.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()),
             file.path(OUT, "15_124_sessionInfo_v01.txt"), useBytes = TRUE)
  index <- data.frame(file = basename(list.files(OUT, recursive = TRUE,
                                               full.names = TRUE)),
                      stringsAsFactors = FALSE)
  write_csv(index, "16_124_OUTPUT_INDEX_v01.csv")
  invisible(list(myeloid_extended = length(extended_cells),
                 myeloid_core = length(core_cells)))
}
tryCatch(run(), error = function(e) {
  log_line("ERROR | Step124 v01 | ", conditionMessage(e))
  stop(e)
})

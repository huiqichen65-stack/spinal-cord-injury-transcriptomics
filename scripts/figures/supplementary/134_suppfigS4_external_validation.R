source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
ROOT <- PROJECT_ROOT
STEM <- "GSE281205"
D123 <- file.path(ROOT, "results", paste0("123_", STEM, "_major_lineage_clustering"))
D128 <- file.path(ROOT, "results", paste0("128_", STEM, "_reviewed_lineage_finalization"))
D130 <- file.path(ROOT, "results", paste0("130_", STEM, "_myeloid_reclustering"))
D131 <- file.path(ROOT, "results", paste0("131_", STEM, "_myeloid_subtype_annotation_and_microglia_groups"))
D132 <- file.path(ROOT, "results", paste0("132_", STEM, "_original_UMI_pseudobulk_signature_validation"))
OUT <- file.path(ROOT, "results", paste0("134_", STEM, "_supplementary_figure_S4"))
FIG_DIR <- file.path(OUT, "figures")
RDS_DIR <- file.path(OUT, "RDS")
LOG_DIR <- file.path(OUT, "logs")
for (d in c(OUT, FIG_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOG <- file.path(LOG_DIR, "00_134_run_log_v03.log")

IN_ATLAS_RDS <- file.path(D123, "RDS", "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
IN_ATLAS_METADATA <- file.path(D128, "03_128_per_cell_curated_lineage_v02.csv")
IN_MYELOID_RDS <- file.path(D130, "RDS", "130_GSE281205_myeloid_LogNorm_PCA_Harmony_clusters_v02.rds")
IN_MYELOID_MEMBERSHIP <- file.path(D130, "05_130_per_cell_myeloid_cluster_membership_v02.csv")
IN_MARKERS <- file.path(D130, "10_130_canonical_marker_evidence_v02.csv")
IN_COMPOSITION <- file.path(D131, "05_131_GSM_primary_group_composition_v01.csv")
IN_PAIRED_DELTA <- file.path(D132, "06_132_per_GSM_frozen_signature_paired_deltas_v02.csv")
IN_STAGE_SUMMARY <- file.path(D132, "07_132_stage_signature_validation_summary_v02.csv")

CACHE_ATLAS <- file.path(RDS_DIR, "134_cached_full_atlas_UMAP_v01.rds")
CACHE_MYELOID <- file.path(RDS_DIR, "134_cached_myeloid_UMAP_v01.rds")
EXPECTED_GSM <- paste0("GSM", 8614662:8614667)
EXPECTED_TOTAL <- 42361L
EXPECTED_MYELOID <- 28716L
FIG_WIDTH <- 18
FIG_HEIGHT <- 12
FIG_DPI <- 600L

PDF_CAIRO <- isTRUE(capabilities("cairo"))
FONT_FAMILY <- if (PDF_CAIRO && .Platform$OS.type == "windows") "Arial" else "sans"

NEEDED <- c("ggplot2", "patchwork", "SeuratObject")
missing_packages <- NEEDED[!vapply(NEEDED, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing R packages: ", paste(missing_packages, collapse = ", "), " (e.g., install.packages(c('ggplot2','patchwork')); SeuratObject via installed Seurat).",
    call. = FALSE)
}

log_it <- function(...) {
  msg <- paste(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG, append = TRUE)
  invisible(msg)
}
assert <- function(condition, message) {
  if (!isTRUE(condition))
    stop(message, call. = FALSE)
}
load_csv <- function(path) {
  assert(file.exists(path), paste("MISSING INPUT:", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
need_columns <- function(df, cols, label) {
  absent <- setdiff(cols, names(df))
  assert(length(absent) == 0L, paste(label, "missing columns:", paste(absent, collapse = ", ")))
}
save_csv <- function(df, filename) {
  path <- file.path(OUT, filename)
  utils::write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  log_it("Saved data: ", path, " [", nrow(df), " rows]")
  invisible(path)
}

STATE_LABELS <- c(`0` = "P2ry12-high microglia", `1` = "Clec7a-high lipid-associated microglia", `2` = "Monocyte/macrophage-like",
  `3` = "APC/DC-like", `4` = "Proliferating microglia")
STATE_COLORS <- c(`P2ry12-high microglia` = "#418A9B", `Clec7a-high lipid-associated microglia` = "#D98946",
  `Monocyte/macrophage-like` = "#8166A2", `APC/DC-like` = "#4B957B", `Proliferating microglia` = "#C76170",
  `QC-excluded / ambiguous` = "#B9BEC5")
ATLAS_LABELS <- c(Microglia = "Microglia", Monocyte_macrophage = "Monocyte/macrophage", APC_dendritic = "APC/DC-like",
  Ependymal = "Ependymal", Endothelial = "Endothelial", T_cell = "T cells", Unresolved = "Unresolved",
  OPC = "OPCs", Pericyte_vSMC = "Pericyte/VSMC", Neuron = "Neurons", NK_cytotoxic = "NK/cytotoxic",
  Vascular_ambiguous = "Vascular ambiguous", Fibroblast_meningeal = "Fibroblast/meningeal", Astrocyte = "Astrocytes",
  Oligodendrocyte = "Oligodendrocytes")
ATLAS_COLORS <- c(Microglia = "#418A9B", `Monocyte/macrophage` = "#8166A2", `APC/DC-like` = "#4B957B",
  Ependymal = "#82A7BD", Endothelial = "#3976AD", `T cells` = "#D56B80", Unresolved = "#C1C4C9", OPCs = "#B59CCB",
  `Pericyte/VSMC` = "#607E91", Neurons = "#6B6F9F", `NK/cytotoxic` = "#E09B9C", `Vascular ambiguous` = "#989B9F",
  `Fibroblast/meningeal` = "#B59272", Astrocytes = "#79A474", Oligodendrocytes = "#9DA75C")
GSM_LABELS <- stats::setNames(paste0("466", 2:7), EXPECTED_GSM)
STAGE_COLORS <- c(`14 dpi` = "#626971", `28 dpi` = "#1E2228")

base_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size, base_family = FONT_FAMILY) + ggplot2::theme(plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = base_size +
      0.2), axis.text = ggplot2::element_text(size = base_size - 0.1), plot.margin = ggplot2::margin(12,
      8, 8, 8), legend.title = ggplot2::element_text(size = base_size - 0.3), legend.text = ggplot2::element_text(size = base_size -
      1.2), panel.grid = ggplot2::element_blank())
}

extract_umap <- function(rds_path, cache_path, reduction, reference_ids, expected_n, label) {
  if (file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (is.data.frame(cached) && all(c("cell_id", "UMAP1", "UMAP2") %in% names(cached)) && nrow(cached) ==
      expected_n && !anyDuplicated(cached$cell_id) && setequal(as.character(cached$cell_id), as.character(reference_ids)) &&
      all(is.finite(cached$UMAP1)) && all(is.finite(cached$UMAP2))) {
      log_it("Reusing cached ", label, " UMAP; large Seurat RDS not reloaded")
      return(cached)
    }
    log_it("Cached UMAP missing/inconsistent; extracting from Step", label, " RDS")
  }
  assert(file.exists(rds_path), paste("MISSING RDS for", label, ":", rds_path))
  log_it("Loading ", label, " RDS (may consume considerable RAM): ", rds_path)
  object <- readRDS(rds_path)
  assert(inherits(object, "Seurat"), paste(label, "input is not Seurat"))
  assert(reduction %in% SeuratObject::Reductions(object), paste(label, "UMAP reduction not found:",
    reduction))
  emb <- SeuratObject::Embeddings(object[[reduction]])
  assert(ncol(emb) >= 2L, paste(label, "UMAP has less than 2 dimensions"))
  coords <- data.frame(cell_id = rownames(emb), UMAP1 = as.numeric(emb[, 1]), UMAP2 = as.numeric(emb[,
    2]), stringsAsFactors = FALSE)
  rm(emb, object)
  invisible(gc())
  assert(nrow(coords) == expected_n, paste(label, "UMAP cell count mismatch"))
  assert(!anyDuplicated(coords$cell_id) && setequal(as.character(coords$cell_id), as.character(reference_ids)),
    paste(label, "UMAP cell IDs differ from validated metadata"))
  assert(all(is.finite(coords$UMAP1)) && all(is.finite(coords$UMAP2)), paste(label, "nonfinite UMAP coordinates"))
  saveRDS(coords, cache_path, compress = "gzip")
  log_it("Saved ", label, " coordinate cache: ", cache_path)
  coords
}

save_panel <- function(p, tag) {
  dest <- file.path(FIG_DIR, paste0("Fig_134_S4_panel_", tag, "_v03.png"))
  ggplot2::ggsave(dest, plot = p, device = "png", width = 8.5, height = 6, units = "in", dpi = FIG_DPI,
    bg = "white", limitsize = FALSE)
  log_it("Saved panel ", tag, ": ", dest)
}

run_134 <- function() {
  log_it("============================================================")
  log_it("Step134 v03 | Supplementary Figure S4 | ", format(Sys.time()))
  log_it("Output: ", OUT)
  log_it("No ZIP; direct local inputs; NO new DE, CAMERA, clustering or UMAP")
  log_it("============================================================")
  paths <- c(IN_ATLAS_METADATA, IN_MYELOID_MEMBERSHIP, IN_MARKERS, IN_COMPOSITION, IN_PAIRED_DELTA,
    IN_STAGE_SUMMARY)
  missing <- paths[!file.exists(paths)]
  assert(!length(missing), paste("Missing required input(s):", paste(missing, collapse = "; ")))
  save_csv(data.frame(purpose = c("123_whole_atlas_UMAP_RDS", "128_full_atlas_lineage_CSV", "130_myeloid_UMAP_RDS",
    "130_myeloid_membership_CSV", "130_canonical_markers", "131_GSM_composition", "132_v02_paired_deltas",
    "132_v02_stage_summary"), path = c(IN_ATLAS_RDS, IN_ATLAS_METADATA, IN_MYELOID_RDS, IN_MYELOID_MEMBERSHIP,
    IN_MARKERS, IN_COMPOSITION, IN_PAIRED_DELTA, IN_STAGE_SUMMARY), stringsAsFactors = FALSE), "01_134_input_provenance_v03.csv")

  atlas <- load_csv(IN_ATLAS_METADATA)
  need_columns(atlas, c("cell_id", "gsm", "dpi", "curated_major_lineage_128"), "Step128 full atlas")
  assert(nrow(atlas) == EXPECTED_TOTAL && !anyDuplicated(atlas$cell_id), "Step128 atlas must contain exactly 42,361 unique cell IDs")
  assert(setequal(as.character(unique(atlas$gsm)), EXPECTED_GSM), "Unexpected GSM IDs in Step128 atlas")
  assert(all(as.character(atlas$curated_major_lineage_128) %in% names(ATLAS_LABELS)), "Step128 contains a major lineage not listed in ATLAS_LABELS")
  atlas_xy <- extract_umap(IN_ATLAS_RDS, CACHE_ATLAS, "umap_harmony", atlas$cell_id, EXPECTED_TOTAL,
    "123/full-atlas")
  ai <- match(atlas_xy$cell_id, atlas$cell_id)
  assert(!anyNA(ai), "Full-atlas UMAP could not match Step128 annotations")
  atlas_xy$lineage <- factor(unname(ATLAS_LABELS[atlas$curated_major_lineage_128[ai]]), levels = names(ATLAS_COLORS))
  assert(!anyNA(atlas_xy$lineage), "Unmapped atlas annotations")
  save_csv(as.data.frame(table(lineage = atlas_xy$lineage)), "02_134_panelA_atlas_lineage_counts_v03.csv")
  pA <- ggplot2::ggplot(atlas_xy, ggplot2::aes(x = UMAP1, y = UMAP2, color = lineage)) + ggplot2::geom_point(size = 0.16,
    alpha = 0.94, shape = 16, stroke = 0) + ggplot2::scale_color_manual(values = ATLAS_COLORS, drop = FALSE) +
    ggplot2::coord_equal() + ggplot2::labs(x = "UMAP 1", y = "UMAP 2", color = NULL) + ggplot2::guides(color = ggplot2::guide_legend(ncol = 3,
    override.aes = list(size = 2, alpha = 1))) + base_theme(10.3) + ggplot2::theme(legend.position = "bottom",
    legend.key.height = grid::unit(0.27, "cm"), legend.spacing.y = grid::unit(0.06, "cm"))
  rm(atlas, atlas_xy)
  invisible(gc())
  log_it("Panel A ready")

  membership <- load_csv(IN_MYELOID_MEMBERSHIP)
  need_columns(membership, c("cell_id", "gsm", "dpi", "myeloid_harmony_r010"), "Step130 membership")
  assert(nrow(membership) == EXPECTED_MYELOID && !anyDuplicated(membership$cell_id), "Step130 must contain exactly 28,716 unique candidate myeloid cells")
  membership$cluster <- as.character(membership$myeloid_harmony_r010)
  assert(all(membership$cluster %in% names(STATE_LABELS)), "Unexpected Step130 r0.10 myeloid cluster")
  my_xy <- extract_umap(IN_MYELOID_RDS, CACHE_MYELOID, "umap_myeloid_harmony", membership$cell_id,
    EXPECTED_MYELOID, "130/myeloid")
  bi <- match(my_xy$cell_id, membership$cell_id)
  assert(!anyNA(bi), "Myeloid UMAP could not match Step130 membership")
  my_xy$state <- factor(unname(STATE_LABELS[membership$cluster[bi]]), levels = unname(STATE_LABELS))
  save_csv(as.data.frame(table(state = my_xy$state)), "03_134_panelB_myeloid_state_counts_v03.csv")
  pB <- ggplot2::ggplot(my_xy, ggplot2::aes(x = UMAP1, y = UMAP2, color = state)) + ggplot2::geom_point(size = 0.24,
    alpha = 0.91, stroke = 0) + ggplot2::scale_color_manual(values = STATE_COLORS[unname(STATE_LABELS)],
    drop = FALSE) + ggplot2::coord_equal() + ggplot2::labs(x = "UMAP 1", y = "UMAP 2", color = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 2.2, alpha = 1))) +
    base_theme(10.3) + ggplot2::theme(legend.position = "bottom", legend.key.height = grid::unit(0.27,
    "cm"))
  rm(my_xy, membership)
  invisible(gc())
  log_it("Panel B ready")

  mk <- load_csv(IN_MARKERS)
  need_columns(mk, c("cluster", "gene_symbol", "pct_cells_corrected_counts_ge1", "mean_log_normalized_expr"),
    "Step130 marker evidence")
  MARKERS <- c("P2ry12", "Tmem119", "Fcrls", "Hexb", "Cx3cr1", "Clec7a", "Lpl", "Spp1", "Gpnmb", "Cst7",
    "Apoe", "Lyz2", "Ms4a7", "Mrc1", "Cd74", "H2-Aa", "Flt3", "Cd209a", "Mki67", "Pclaf", "Birc5")
  missing_mk <- setdiff(MARKERS, unique(as.character(mk$gene_symbol)))
  assert(!length(missing_mk), paste("Step130 canonical marker evidence missing:", paste(missing_mk,
    collapse = ", ")))
  mk <- mk[mk$gene_symbol %in% MARKERS & as.character(mk$cluster) %in% names(STATE_LABELS), , drop = FALSE]
  assert(nrow(mk) == 5L * length(MARKERS) && !anyDuplicated(paste(mk$cluster, mk$gene_symbol)), "Marker evidence must have exactly one row per gene x 5 myeloid clusters")
  assert(all(is.finite(mk$mean_log_normalized_expr)) && all(is.finite(mk$pct_cells_corrected_counts_ge1)) &&
    all(mk$pct_cells_corrected_counts_ge1 >= 0 & mk$pct_cells_corrected_counts_ge1 <= 100), "Step130 marker evidence has non-finite or out-of-range values")

  mk$scaled_mean <- ave(as.numeric(mk$mean_log_normalized_expr), mk$gene_symbol, FUN = function(x) {
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0)
      return(rep(0, length(x)))
    pmax(-2, pmin(2, (x - mean(x))/s))
  })
  mk$gene_symbol <- factor(mk$gene_symbol, levels = MARKERS)
  mk$state <- factor(unname(STATE_LABELS[as.character(mk$cluster)]), levels = rev(unname(STATE_LABELS)))
  save_csv(mk, "04_134_panelC_marker_dotplot_data_v03.csv")
  pC <- ggplot2::ggplot(mk, ggplot2::aes(x = gene_symbol, y = state)) + ggplot2::geom_point(ggplot2::aes(size = pct_cells_corrected_counts_ge1,
    color = scaled_mean)) + ggplot2::scale_size_area(max_size = 4, limits = c(0, 100), breaks = c(0,
    25, 50, 75, 100), name = "Detected (%)") + ggplot2::scale_color_gradient2(low = "#3979A4", mid = "#F4F4F4",
    high = "#BF643B", midpoint = 0, limits = c(-2, 2), name = "Scaled mean") + ggplot2::labs(x = NULL,
    y = NULL) + base_theme(10.3) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
    hjust = 1, vjust = 1, size = 9.4), axis.text.y = ggplot2::element_text(size = 9), legend.position = "bottom",
    legend.box = "horizontal", panel.grid.major = ggplot2::element_line(color = "#ECEFF2", linewidth = 0.25)) +
    ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = grid::unit(2.3, "cm"), barheight = grid::unit(0.28,
      "cm")), size = ggplot2::guide_legend(nrow = 1))
  log_it("Panel C ready")

  comp <- load_csv(IN_COMPOSITION)
  need_columns(comp, c("gsm", "group_131", "n_cells", "dpi"), "Step131 specimen composition")
  group_map <- c(TARGET_CLEC7A_HIGH_MICROGLIA = STATE_LABELS[["1"]], REFERENCE_P2RY12_HIGH_MICROGLIA = STATE_LABELS[["0"]],
    OTHER_MONOCYTE_MACROPHAGE = STATE_LABELS[["2"]], OTHER_APC_DC_LIKE = STATE_LABELS[["3"]], PROLIFERATING_MICROGLIA = STATE_LABELS[["4"]],
    EXCLUDED_AMBIGUOUS_OR_QC = "QC-excluded / ambiguous")
  assert(setequal(unique(as.character(comp$group_131)), names(group_map)), "Unexpected Step131 composition group labels")
  assert(setequal(unique(as.character(comp$gsm)), EXPECTED_GSM), "Step131 composition must represent all six GSMs")
  assert(nrow(comp) == 36L && !anyDuplicated(paste(comp$gsm, comp$group_131)), "Step131 should have 6 GSM x 6 state/QC groups")
  comp$gsm <- as.character(comp$gsm)
  comp$dpi <- as.integer(as.character(comp$dpi))
  assert(all(comp$dpi == ifelse(comp$gsm %in% EXPECTED_GSM[1:3], 14L, 28L)), "Step131 GSM vs dpi mismatch")
  comp$n_cells <- as.integer(comp$n_cells)
  assert(all(is.finite(comp$n_cells)) && all(comp$n_cells >= 0), "Invalid state count in Step131 composition")
  assert(sum(comp$n_cells) == EXPECTED_MYELOID, "Step131 specimen composition does not sum to 28,716")
  comp$state <- factor(unname(group_map[comp$group_131]), levels = c(unname(STATE_LABELS), "QC-excluded / ambiguous"))
  comp$sample <- factor(comp$gsm, levels = EXPECTED_GSM, labels = paste0(unname(GSM_LABELS[EXPECTED_GSM]),
    "\n", rep(c("14 dpi", "28 dpi"), each = 3)))
  comp$percent <- 100 * comp$n_cells/as.numeric(ave(comp$n_cells, comp$gsm, FUN = sum))
  save_csv(comp, "05_134_panelD_sample_state_composition_v03.csv")
  pD <- ggplot2::ggplot(comp, ggplot2::aes(x = sample, y = percent, fill = state)) + ggplot2::geom_col(width = 0.78,
    color = "white", linewidth = 0.18) + ggplot2::scale_fill_manual(values = STATE_COLORS, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = ggplot2::expansion(mult = c(0,
      0.02))) + ggplot2::labs(x = "GSM (short ID) and injury stage", y = "Candidate myeloid cells (%)",
    fill = NULL) + base_theme(10.6) + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 9),
    legend.key.height = grid::unit(0.3, "cm")) + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3))
  log_it("Panel D ready")

  delta <- load_csv(IN_PAIRED_DELTA)
  summ <- load_csv(IN_STAGE_SUMMARY)
  A <- "contrast_A_target_vs_other_microglia_131"
  B <- "contrast_B_target_vs_all_other_highconf_myeloid_131"
  desired <- c(A, B)
  need_columns(delta, c("comparison", "dpi", "gsm", "mean_signature_target_minus_other_logCPM", "n_signature_genes_mapped"),
    "Step132 v02 paired deltas")
  need_columns(summ, c("comparison", "dpi", "n_GSM_pairs", "mean_of_three_GSM_signature_deltas", "n_GSM_with_positive_signature_delta",
    "n_signature_mapped"), "Step132 v02 stage summary")
  delta <- delta[delta$comparison %in% desired, , drop = FALSE]
  summ <- summ[summ$comparison %in% desired, , drop = FALSE]
  delta$dpi <- as.integer(as.character(delta$dpi))
  summ$dpi <- as.integer(as.character(summ$dpi))
  assert(nrow(delta) == 12L && !anyDuplicated(paste(delta$comparison, delta$dpi, delta$gsm)), "Step132 must give 2 comparisons x 2 dpi x 3 GSM = 12 specimen deltas")
  assert(nrow(summ) == 4L && !anyDuplicated(paste(summ$comparison, summ$dpi)), "Step132 must give 2 comparisons x 2 dpi = 4 stage summaries")
  assert(all(delta$gsm %in% EXPECTED_GSM) && all(delta$dpi == ifelse(delta$gsm %in% EXPECTED_GSM[1:3],
    14L, 28L)), "Step132 delta specimen/stage mismatch")
  assert(all(summ$n_GSM_pairs == 3L) && all(summ$n_signature_mapped == 1301L) && all(delta$n_signature_genes_mapped ==
    1301L), "Step132 v02 replicate or frozen-signature mapping audit failed")
  test_means <- stats::aggregate(mean_signature_target_minus_other_logCPM ~ comparison + dpi, delta,
    mean)
  idx <- match(paste(summ$comparison, summ$dpi), paste(test_means$comparison, test_means$dpi))
  assert(!anyNA(idx) && max(abs(summ$mean_of_three_GSM_signature_deltas - test_means$mean_signature_target_minus_other_logCPM[idx])) <
    1e-06, "Step132 summary mean and specimen deltas do not agree")
  nice <- c(contrast_A_target_vs_other_microglia_131 = "Target vs other microglia", contrast_B_target_vs_all_other_highconf_myeloid_131 = "Target vs other high-confidence myeloid")
  delta$comparison_label <- factor(unname(nice[delta$comparison]), levels = unname(nice[desired]))
  summ$comparison_label <- factor(unname(nice[summ$comparison]), levels = unname(nice[desired]))
  delta$stage <- factor(paste0(delta$dpi, " dpi"), levels = names(STAGE_COLORS))
  summ$stage <- factor(paste0(summ$dpi, " dpi"), levels = names(STAGE_COLORS))
  save_csv(delta, "06_134_panelE_specimen_signature_deltas_v03.csv")

  drop_ci <- intersect(c("exploratory_95CI_lower_three_GSM", "exploratory_95CI_upper_three_GSM"), names(summ))
  if (length(drop_ci))
    summ <- summ[, !(names(summ) %in% drop_ci), drop = FALSE]
  save_csv(summ, "07_134_panelE_stage_effect_and_consistency_v03.csv")

  set.seed(13401L)
  pE <- ggplot2::ggplot() + ggplot2::geom_hline(yintercept = 0, linewidth = 0.38, linetype = "dashed",
    color = "#8D939A") + ggplot2::geom_point(data = delta, ggplot2::aes(x = stage, y = mean_signature_target_minus_other_logCPM,
    color = stage), position = ggplot2::position_jitter(width = 0.105, height = 0, seed = 13401L),
    size = 3.1, alpha = 0.95) + ggplot2::geom_point(data = summ, ggplot2::aes(x = stage, y = mean_of_three_GSM_signature_deltas),
    shape = 23, fill = "white", color = "#24292E", stroke = 1.2, size = 4.5) + ggplot2::facet_wrap(~comparison_label,
    ncol = 2) + ggplot2::scale_color_manual(values = STAGE_COLORS, guide = "none") + ggplot2::labs(x = NULL,
    y = "Target - comparator: mean signature delta (logCPM)") + base_theme(10.8) + ggplot2::theme(strip.background = ggplot2::element_rect(fill = "#F2F4F5",
    color = NA), strip.text = ggplot2::element_text(face = "bold", size = 10.1), panel.spacing.x = grid::unit(0.35,
    "cm"), axis.text.x = ggplot2::element_text(size = 10))

  pE <- pE + ggplot2::geom_text(data = summ, ggplot2::aes(x = stage, y = mean_of_three_GSM_signature_deltas,
    label = sprintf("%+.3f", mean_of_three_GSM_signature_deltas)), nudge_y = 0.012, size = 3.5, color = "#222222",
    fontface = "bold", family = FONT_FAMILY, show.legend = FALSE) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12,
    0.23)))
  log_it("Panel E ready")

  design <- paste(c("AAABBBCCCC", "AAABBBCCCC", "DDDDDEEEEE", "DDDDDEEEEE"), collapse = "\n")

  no_panel_heading <- ggplot2::theme(plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
    plot.tag.position = c(0.015, 0.995), plot.tag = ggplot2::element_text(face = "bold", size = 21,
      hjust = 0, vjust = 1, family = FONT_FAMILY))
  pA <- pA + no_panel_heading
  pB <- pB + no_panel_heading
  pC <- pC + no_panel_heading
  pD <- pD + no_panel_heading
  pE <- pE + no_panel_heading
  fig <- patchwork::wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE, design = design) + patchwork::plot_annotation(tag_levels = "A",
    theme = ggplot2::theme(plot.tag.position = c(0.015, 0.995), plot.tag = ggplot2::element_text(face = "bold",
      size = 21, hjust = 0, vjust = 1, family = FONT_FAMILY)))
  out_base <- file.path(FIG_DIR, "Supplementary_Figure_S4_GSE281205_v03")
  if (PDF_CAIRO) {
    ggplot2::ggsave(paste0(out_base, ".pdf"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT,
      units = "in", bg = "white", device = grDevices::cairo_pdf, limitsize = FALSE)
  } else {
    ggplot2::ggsave(paste0(out_base, ".pdf"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT,
      units = "in", bg = "white", device = grDevices::pdf, useDingbats = FALSE, limitsize = FALSE)
  }
  log_it("Saved vector PDF")
  ggplot2::ggsave(paste0(out_base, ".png"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in",
    dpi = FIG_DPI, bg = "white", device = "png", limitsize = FALSE)
  log_it("Saved ", FIG_DPI, " dpi PNG")
  ggplot2::ggsave(paste0(out_base, ".tiff"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in",
    dpi = FIG_DPI, bg = "white", device = "tiff", compression = "lzw", limitsize = FALSE)
  log_it("Saved ", FIG_DPI, " dpi LZW TIFF")
  save_panel(pA, "A")
  save_panel(pB, "B")
  save_panel(pC, "C")
  save_panel(pD, "D")
  save_panel(pE, "E")

  legend <- c("Supplementary Figure S4. Independent single-cell characterization and specimen-level validation of Clec7a-high lipid-associated microglia in GSE281205.",
    "(A) Harmony UMAP of 42,361 retained cells from six independent specimens, colored by evidence-based working major-lineage annotation (Step128). Unresolved/vascular-ambiguous groups are retained as such.",
    "(B) Harmony UMAP of 28,716 candidate myeloid cells at resolution 0.10, colored by five working transcriptomic states (Steps130-131); state identity is not a lineage-tracing result.",
    "(C) Canonical-marker DotPlot across the five candidate myeloid states. Dot size indicates the percentage of cells with DecontX-estimated counts >=1; color indicates per-gene standardized mean log-normalized expression. These descriptive markers are not specimen-level differential expression statistics.",
    "(D) Per-specimen composition of Step131 high-confidence myeloid states and QC-excluded/ambiguous candidate cells. Denominator: all candidate myeloid cells in each GSM; display is descriptive and does not establish a significant 14-to-28 dpi change.",
    "(E) Stage-matched within-specimen effects of the frozen 1,320-gene signature, using original-integer-UMI specimen-level pseudobulk results (Step132 v02). Each dot is a GSM-specific target-minus-comparator mean logCPM difference over 1,301 mapped signature genes; white diamonds and adjacent numeric labels mark the three-GSM stage mean. Confidence intervals are not plotted; all three independent specimen effects at each stage are shown. Left: Clec7a-high vs other P2ry12-high microglia; right: Clec7a-high vs all other high-confidence myeloid cells. Positive means are not evidence that all signature genes replicate.",
    "All 14 dpi and 28 dpi comparisons were performed within stage, with three independent GSMs at each stage. GSM is nested within injury stage in the Harmony embedding; stage interpretations should therefore remain cautious. The GEO 28-dpi age-metadata discrepancy is documented in analysis provenance.")
  writeLines(legend, file.path(OUT, "08_134_Supplementary_Figure_S4_legend_v03.txt"), useBytes = TRUE)
  log_it("Saved figure legend")
  audit <- data.frame(check = c("atlas_n_cells", "candidate_myeloid_n_cells", "n_specimens", "stages",
    "Step131_composition_n_cells", "figure_panels", "Step132_primary_comparisons", "Step132_stage_summaries",
    "mapped_signature_genes", "camera_shown", "workflow_shown", "ci_shown", "dpi"), value = c(EXPECTED_TOTAL,
    EXPECTED_MYELOID, length(EXPECTED_GSM), "14dpi;28dpi", sum(comp$n_cells), "A;B;C;D;E", 2, 4,
    unique(summ$n_signature_mapped), "NO", "NO", "NO", FIG_DPI))
  save_csv(audit, "09_134_figure_validation_audit_v03.csv")
  capture.output(sessionInfo(), file = file.path(OUT, "10_134_sessionInfo_v03.txt"))
  log_it("SUCCESS | Step134 v03 | Supplementary Figure S4 completed")
  invisible(fig)
}

tryCatch(run_134(), error = function(e) {
  log_it("FAILED | Step134 v03 | ", conditionMessage(e))
  stop(e)
})

source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
ROOT <- PROJECT_ROOT
STEM <- "Reviewer1_Major4_GSE281205"
D123 <- file.path(ROOT, "results", paste0("123_", STEM, "_major_lineage_clustering"))
D128 <- file.path(ROOT, "results", paste0("128_", STEM, "_reviewed_lineage_finalization"))
D130 <- file.path(ROOT, "results", paste0("130_", STEM, "_myeloid_reclustering"))
D131 <- file.path(ROOT, "results", paste0("131_", STEM, "_myeloid_subtype_annotation_and_microglia_groups"))
D132 <- file.path(ROOT, "results", paste0("132_", STEM, "_original_UMI_pseudobulk_signature_validation"))
D145 <- file.path(ROOT, "results", "145_Reviewer1_RemainingMajor1_GSE281205_naive_extension")
OUT <- file.path(ROOT, "results", paste0("134_", STEM, "_supplementary_figure_S4"))
FIG_DIR <- file.path(OUT, "figures")
RDS_DIR <- file.path(OUT, "RDS")
LOG_DIR <- file.path(OUT, "logs")
for (d in c(OUT, FIG_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOG <- file.path(LOG_DIR, "00_134_run_log_v05.log")

IN_ATLAS_RDS <- file.path(D123, "RDS", "123_GSE281205_Harmony_major_lineage_clustering_v03.rds")
IN_ATLAS_METADATA <- file.path(D128, "03_128_per_cell_curated_lineage_v02.csv")
IN_MYELOID_RDS <- file.path(D130, "RDS", "130_GSE281205_myeloid_LogNorm_PCA_Harmony_clusters_v02.rds")
IN_MYELOID_MEMBERSHIP <- file.path(D130, "05_130_per_cell_myeloid_cluster_membership_v02.csv")
IN_MARKERS <- file.path(D130, "10_130_canonical_marker_evidence_v02.csv")
IN_COMPOSITION <- file.path(D131, "05_131_GSM_primary_group_composition_v01.csv")
IN_PAIRED_DELTA <- file.path(D132, "06_132_per_GSM_frozen_signature_paired_deltas_v02.csv")
IN_STAGE_SUMMARY <- file.path(D132, "07_132_stage_signature_validation_summary_v02.csv")
IN145_SPECIMEN <- file.path(D145, "09_145_specimen_Clec7a_high_proportions_v01.csv")
IN145_STAGE <- file.path(D145, "10_145_stage_descriptive_summary_v01.csv")
IN145_PAIRWISE <- file.path(D145, "11_145_pairwise_effect_size_exact_permutation_v01.csv")
IN145_OMNIBUS <- file.path(D145, "12_145_three_group_exact_permutation_v01.csv")

CACHE_ATLAS <- file.path(RDS_DIR, "134_cached_full_atlas_UMAP_v01.rds")
CACHE_MYELOID <- file.path(RDS_DIR, "134_cached_myeloid_UMAP_v01.rds")
EXPECTED_GSM <- paste0("GSM", 8614662:8614667)
EXPECTED_TOTAL <- 42361L
EXPECTED_MYELOID <- 28716L
FIG_WIDTH <- 18.0
FIG_HEIGHT <- 13.4
FIG_DPI <- 600L

PDF_CAIRO <- isTRUE(capabilities("cairo"))
FONT_FAMILY <- if (PDF_CAIRO && .Platform$OS.type == "windows") "Arial" else "sans"

NEEDED <- c("ggplot2", "patchwork", "SeuratObject")
missing_packages <- NEEDED[!vapply(NEEDED, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing R packages: ", paste(missing_packages, collapse = ", "), " (e.g., install.packages(c('ggplot2','patchwork')); SeuratObject via installed Seurat).", call. = FALSE)
}

log_it <- function(...) {
  msg <- paste(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG, append = TRUE)
  invisible(msg)
}
assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
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

STATE_LABELS <- c( "0" = "P2ry12-high microglia", "1" = "Clec7a-high candidate myeloid state", "2" = "Monocyte/macrophage-like", "3" = "APC/DC-like", "4" = "Proliferating microglia" )
STATE_COLORS <- c( "P2ry12-high microglia" = "#418A9B", "Clec7a-high candidate myeloid state" = "#D98946", "Monocyte/macrophage-like" = "#8166A2", "APC/DC-like" = "#4B957B",
  "Proliferating microglia" = "#C76170", "QC-excluded / ambiguous" = "#B9BEC5" )
ATLAS_LABELS <- c( "Microglia" = "Microglia", "Monocyte_macrophage" = "Monocyte/macrophage", "APC_dendritic" = "APC/DC-like", "Ependymal" = "Ependymal", "Endothelial" = "Endothelial",
  "T_cell" = "T cells", "Unresolved" = "Unresolved", "OPC" = "OPCs", "Pericyte_vSMC" = "Pericyte/VSMC", "Neuron" = "Neurons", "NK_cytotoxic" = "NK/cytotoxic",
  "Vascular_ambiguous" = "Vascular ambiguous", "Fibroblast_meningeal" = "Fibroblast/meningeal", "Astrocyte" = "Astrocytes", "Oligodendrocyte" = "Oligodendrocytes" )
ATLAS_COLORS <- c( "Microglia" = "#418A9B", "Monocyte/macrophage" = "#8166A2", "APC/DC-like" = "#4B957B", "Ependymal" = "#82A7BD", "Endothelial" = "#3976AD", "T cells" = "#D56B80",
  "Unresolved" = "#C1C4C9", "OPCs" = "#B59CCB", "Pericyte/VSMC" = "#607E91", "Neurons" = "#6B6F9F", "NK/cytotoxic" = "#E09B9C", "Vascular ambiguous" = "#989B9F",
  "Fibroblast/meningeal" = "#B59272", "Astrocytes" = "#79A474", "Oligodendrocytes" = "#9DA75C" )
GSM_LABELS <- stats::setNames(paste0("466", 2:7), EXPECTED_GSM)
STAGE_COLORS <- c("14 dpi" = "#626971", "28 dpi" = "#1E2228")
EXT_STAGE_COLORS <- c("Naive" = "#B9BEC5", "14 dpi" = "#626971", "28 dpi" = "#1E2228")

base_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    ggplot2::theme( plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = base_size + 0.2),
      axis.text = ggplot2::element_text(size = base_size - 0.1), plot.margin = ggplot2::margin(12, 8, 8, 8), legend.title = ggplot2::element_text(size = base_size - 0.3),
      legend.text = ggplot2::element_text(size = base_size - 1.2), panel.grid = ggplot2::element_blank() )
}

extract_umap <- function(rds_path, cache_path, reduction, reference_ids,
                         expected_n, label) {
  if (file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (is.data.frame(cached) &&
        all(c("cell_id", "UMAP1", "UMAP2") %in% names(cached)) &&
        nrow(cached) == expected_n &&
        !anyDuplicated(cached$cell_id) &&
        setequal(as.character(cached$cell_id), as.character(reference_ids)) &&
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
  assert(reduction %in% SeuratObject::Reductions(object), paste(label, "UMAP reduction not found:", reduction))
  emb <- SeuratObject::Embeddings(object[[reduction]])
  assert(ncol(emb) >= 2L, paste(label, "UMAP has less than 2 dimensions"))
  coords <- data.frame(cell_id = rownames(emb), UMAP1 = as.numeric(emb[, 1]), UMAP2 = as.numeric(emb[, 2]), stringsAsFactors = FALSE)
  rm(emb, object)
  invisible(gc())
  assert(nrow(coords) == expected_n, paste(label, "UMAP cell count mismatch"))
  assert(!anyDuplicated(coords$cell_id) && setequal(as.character(coords$cell_id), as.character(reference_ids)), paste(label, "UMAP cell IDs differ from validated metadata"))
  assert(all(is.finite(coords$UMAP1)) && all(is.finite(coords$UMAP2)), paste(label, "nonfinite UMAP coordinates"))
  saveRDS(coords, cache_path, compress = "gzip")
  log_it("Saved ", label, " coordinate cache: ", cache_path)
  coords
}

save_panel <- function(p, tag) {
  dest <- file.path(FIG_DIR, paste0("Fig_134_S4_panel_", tag, "_v05.png"))
  ggplot2::ggsave(dest, plot = p, device = "png", width = 8.5, height = 6.0, units = "in", dpi = FIG_DPI, bg = "white", limitsize = FALSE)
  log_it("Saved panel ", tag, ": ", dest)
}

run_134 <- function() {
  log_it("============================================================")
  log_it("Step134 v05 | Supplementary Figure S4 | ", format(Sys.time()))
  log_it("Output: ", OUT)
  log_it("No ZIP; direct local inputs; A-E unchanged; Panel F reads frozen Step145 outputs")
  log_it("============================================================")
  assert(identical(unname(STATE_LABELS[["1"]]), "Clec7a-high candidate myeloid state"), "Reviewer 2 focal-state display terminology is inconsistent")
  assert("Clec7a-high candidate myeloid state" %in% names(STATE_COLORS), "Reviewer 2 focal-state color key is missing")
  paths <- c(IN_ATLAS_METADATA, IN_MYELOID_MEMBERSHIP, IN_MARKERS, IN_COMPOSITION, IN_PAIRED_DELTA, IN_STAGE_SUMMARY, IN145_SPECIMEN, IN145_STAGE, IN145_PAIRWISE, IN145_OMNIBUS)
  missing <- paths[!file.exists(paths)]
  assert(!length(missing), paste("Missing required input(s):", paste(missing, collapse = "; ")))
  save_csv(data.frame( purpose = c("123_whole_atlas_UMAP_RDS", "128_full_atlas_lineage_CSV", "130_myeloid_UMAP_RDS", "130_myeloid_membership_CSV",
    "130_canonical_markers", "131_GSM_composition", "132_v02_paired_deltas", "132_v02_stage_summary", "145_v01_nine_specimen_proportions", "145_v01_stage_summary",
    "145_v01_pairwise_effects", "145_v01_three_group_omnibus"), path = c(IN_ATLAS_RDS, IN_ATLAS_METADATA, IN_MYELOID_RDS, IN_MYELOID_MEMBERSHIP, IN_MARKERS, IN_COMPOSITION,
    IN_PAIRED_DELTA, IN_STAGE_SUMMARY, IN145_SPECIMEN, IN145_STAGE, IN145_PAIRWISE, IN145_OMNIBUS), stringsAsFactors = FALSE ), "01_134_input_provenance_v05.csv")

  atlas <- load_csv(IN_ATLAS_METADATA)
  need_columns(atlas, c("cell_id", "gsm", "dpi", "curated_major_lineage_128"), "Step128 full atlas")
  assert(nrow(atlas) == EXPECTED_TOTAL && !anyDuplicated(atlas$cell_id), "Step128 atlas must contain exactly 42,361 unique cell IDs")
  assert(setequal(as.character(unique(atlas$gsm)), EXPECTED_GSM), "Unexpected GSM IDs in Step128 atlas")
  assert(all(as.character(atlas$curated_major_lineage_128) %in% names(ATLAS_LABELS)), "Step128 contains a major lineage not listed in ATLAS_LABELS")
  atlas_xy <- extract_umap(IN_ATLAS_RDS, CACHE_ATLAS, "umap_harmony", atlas$cell_id, EXPECTED_TOTAL, "123/full-atlas")
  ai <- match(atlas_xy$cell_id, atlas$cell_id)
  assert(!anyNA(ai), "Full-atlas UMAP could not match Step128 annotations")
  atlas_xy$lineage <- factor(unname(ATLAS_LABELS[ atlas$curated_major_lineage_128[ai]]), levels = names(ATLAS_COLORS))
  assert(!anyNA(atlas_xy$lineage), "Unmapped atlas annotations")
  save_csv(as.data.frame(table(lineage = atlas_xy$lineage)), "02_134_panelA_atlas_lineage_counts_v05.csv")
  pA <- ggplot2::ggplot(atlas_xy, ggplot2::aes(x = UMAP1, y = UMAP2, color = lineage)) +
    ggplot2::geom_point(size = 0.16, alpha = 0.94, shape = 16, stroke = 0) +
    ggplot2::scale_color_manual(values = ATLAS_COLORS, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2", color = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 3, override.aes = list(size = 2.0, alpha = 1))) +
    base_theme(10.3) +
    ggplot2::theme(legend.position = "bottom", legend.key.height = grid::unit(0.27, "cm"), legend.spacing.y = grid::unit(0.06, "cm"))
  rm(atlas, atlas_xy)
  invisible(gc())
  log_it("Panel A ready")

  membership <- load_csv(IN_MYELOID_MEMBERSHIP)
  need_columns(membership, c("cell_id", "gsm", "dpi", "myeloid_harmony_r010"), "Step130 membership")
  assert(nrow(membership) == EXPECTED_MYELOID && !anyDuplicated(membership$cell_id), "Step130 must contain exactly 28,716 unique candidate myeloid cells")
  membership$cluster <- as.character(membership$myeloid_harmony_r010)
  assert(all(membership$cluster %in% names(STATE_LABELS)), "Unexpected Step130 r0.10 myeloid cluster")
  my_xy <- extract_umap(IN_MYELOID_RDS, CACHE_MYELOID, "umap_myeloid_harmony", membership$cell_id, EXPECTED_MYELOID, "130/myeloid")
  bi <- match(my_xy$cell_id, membership$cell_id)
  assert(!anyNA(bi), "Myeloid UMAP could not match Step130 membership")
  my_xy$state <- factor(unname(STATE_LABELS[membership$cluster[bi]]), levels = unname(STATE_LABELS))
  save_csv(as.data.frame(table(state = my_xy$state)), "03_134_panelB_myeloid_state_counts_v05.csv")
  pB <- ggplot2::ggplot(my_xy, ggplot2::aes(x = UMAP1, y = UMAP2, color = state)) +
    ggplot2::geom_point(size = 0.24, alpha = 0.91, stroke = 0) +
    ggplot2::scale_color_manual(values = STATE_COLORS[unname(STATE_LABELS)], drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2", color = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 2.2, alpha = 1))) +
    base_theme(10.3) +
    ggplot2::theme(legend.position = "bottom", legend.key.height = grid::unit(0.27, "cm"))
  rm(my_xy, membership)
  invisible(gc())
  log_it("Panel B ready")

  mk <- load_csv(IN_MARKERS)
  need_columns(mk, c("cluster", "gene_symbol", "pct_cells_corrected_counts_ge1", "mean_log_normalized_expr"), "Step130 marker evidence")
  MARKERS <- c("P2ry12", "Tmem119", "Fcrls", "Hexb", "Cx3cr1", "Clec7a", "Lpl", "Spp1", "Gpnmb", "Cst7", "Apoe", "Lyz2", "Ms4a7", "Mrc1", "Cd74", "H2-Aa", "Flt3", "Cd209a",
    "Mki67", "Pclaf", "Birc5")
  missing_mk <- setdiff(MARKERS, unique(as.character(mk$gene_symbol)))
  assert(!length(missing_mk), paste("Step130 canonical marker evidence missing:", paste(missing_mk, collapse = ", ")))
  mk <- mk[mk$gene_symbol %in% MARKERS & as.character(mk$cluster) %in% names(STATE_LABELS), , drop = FALSE]
  assert(nrow(mk) == 5L * length(MARKERS) && !anyDuplicated(paste(mk$cluster, mk$gene_symbol)), "Marker evidence must have exactly one row per gene x 5 myeloid clusters")
  assert(all(is.finite(mk$mean_log_normalized_expr)) && all(is.finite(mk$pct_cells_corrected_counts_ge1)) && all(mk$pct_cells_corrected_counts_ge1 >= 0 &
    mk$pct_cells_corrected_counts_ge1 <= 100), "Step130 marker evidence has non-finite or out-of-range values")

  mk$scaled_mean <- ave(as.numeric(mk$mean_log_normalized_expr),
    mk$gene_symbol, FUN = function(x) {
      s <- stats::sd(x)
      if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
      pmax(-2, pmin(2, (x - mean(x)) / s))
    })
  mk$gene_symbol <- factor(mk$gene_symbol, levels = MARKERS)
  mk$state <- factor(unname(STATE_LABELS[as.character(mk$cluster)]), levels = rev(unname(STATE_LABELS)))
  save_csv(mk, "04_134_panelC_marker_dotplot_data_v05.csv")
  pC <- ggplot2::ggplot(mk, ggplot2::aes(x = gene_symbol, y = state)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_cells_corrected_counts_ge1, color = scaled_mean)) +
    ggplot2::scale_size_area(max_size = 4.0, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = "Detected (%)") +
    ggplot2::scale_color_gradient2(low = "#3979A4", mid = "#F4F4F4", high = "#BF643B", midpoint = 0, limits = c(-2, 2), name = "Scaled mean") +
    ggplot2::labs(x = NULL, y = NULL) +
    base_theme(10.3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 9.4), axis.text.y = ggplot2::element_text(size = 9.0), legend.position = "bottom",
      legend.box = "horizontal", panel.grid.major = ggplot2::element_line(color = "#ECEFF2", linewidth = 0.25)) +
    ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = grid::unit(2.3, "cm"), barheight = grid::unit(0.28, "cm")), size = ggplot2::guide_legend(nrow = 1))
  log_it("Panel C ready")

  comp <- load_csv(IN_COMPOSITION)
  need_columns(comp, c("gsm", "group_131", "n_cells", "dpi"), "Step131 specimen composition")
  group_map <- c( "TARGET_CLEC7A_HIGH_MICROGLIA" = STATE_LABELS[["1"]], "REFERENCE_P2RY12_HIGH_MICROGLIA" = STATE_LABELS[["0"]], "OTHER_MONOCYTE_MACROPHAGE" = STATE_LABELS[["2"]],
    "OTHER_APC_DC_LIKE" = STATE_LABELS[["3"]], "PROLIFERATING_MICROGLIA" = STATE_LABELS[["4"]], "EXCLUDED_AMBIGUOUS_OR_QC" = "QC-excluded / ambiguous" )
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
  comp$sample <- factor(comp$gsm, levels = EXPECTED_GSM, labels = paste0(unname(GSM_LABELS[EXPECTED_GSM]), "\n", rep(c("14 dpi", "28 dpi"), each = 3)))
  comp$percent <- 100 * comp$n_cells /
    as.numeric(ave(comp$n_cells, comp$gsm, FUN = sum))
  save_csv(comp, "05_134_panelD_sample_state_composition_v05.csv")
  pD <- ggplot2::ggplot(comp, ggplot2::aes(x = sample, y = percent, fill = state)) +
    ggplot2::geom_col(width = 0.78, color = "white", linewidth = 0.18) +
    ggplot2::scale_fill_manual(values = STATE_COLORS, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "GSM (short ID) and injury stage", y = "Candidate myeloid cells (%)", fill = NULL) +
    base_theme(10.6) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 9.0), legend.key.height = grid::unit(0.3, "cm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3))
  log_it("Panel D ready")

  delta <- load_csv(IN_PAIRED_DELTA)
  summ <- load_csv(IN_STAGE_SUMMARY)
  A <- "contrast_A_target_vs_other_microglia_131"
  B <- "contrast_B_target_vs_all_other_highconf_myeloid_131"
  desired <- c(A, B)
  need_columns(delta, c("comparison", "dpi", "gsm", "mean_signature_target_minus_other_logCPM", "n_signature_genes_mapped"), "Step132 v02 paired deltas")
  need_columns(summ, c("comparison", "dpi", "n_GSM_pairs", "mean_of_three_GSM_signature_deltas", "n_GSM_with_positive_signature_delta", "n_signature_mapped"), "Step132 v02 stage summary")
  delta <- delta[delta$comparison %in% desired, , drop = FALSE]
  summ <- summ[summ$comparison %in% desired, , drop = FALSE]
  delta$dpi <- as.integer(as.character(delta$dpi))
  summ$dpi <- as.integer(as.character(summ$dpi))
  assert(nrow(delta) == 12L && !anyDuplicated(paste(delta$comparison, delta$dpi, delta$gsm)), "Step132 must give 2 comparisons x 2 dpi x 3 GSM = 12 specimen deltas")
  assert(nrow(summ) == 4L && !anyDuplicated(paste(summ$comparison, summ$dpi)), "Step132 must give 2 comparisons x 2 dpi = 4 stage summaries")
  assert(all(delta$gsm %in% EXPECTED_GSM) && all(delta$dpi == ifelse(delta$gsm %in% EXPECTED_GSM[1:3], 14L, 28L)), "Step132 delta specimen/stage mismatch")
  assert(all(summ$n_GSM_pairs == 3L) && all(summ$n_signature_mapped == 1301L) && all(delta$n_signature_genes_mapped == 1301L),
    "Step132 v02 replicate or frozen-signature mapping audit failed")
  test_means <- stats::aggregate(mean_signature_target_minus_other_logCPM ~ comparison + dpi, delta, mean)
  idx <- match(paste(summ$comparison, summ$dpi), paste(test_means$comparison, test_means$dpi))
  assert(!anyNA(idx) && max(abs(summ$mean_of_three_GSM_signature_deltas - test_means$mean_signature_target_minus_other_logCPM[idx])) < 1e-6,
    "Step132 summary mean and specimen deltas do not agree")
  nice <- c( "contrast_A_target_vs_other_microglia_131" = "Candidate state vs P2ry12-high microglia", "contrast_B_target_vs_all_other_highconf_myeloid_131" =
    "Candidate state vs other high-confidence myeloid" )
  delta$comparison_label <- factor(unname(nice[delta$comparison]), levels = unname(nice[desired]))
  summ$comparison_label <- factor(unname(nice[summ$comparison]), levels = unname(nice[desired]))
  delta$stage <- factor(paste0(delta$dpi, " dpi"), levels = names(STAGE_COLORS))
  summ$stage <- factor(paste0(summ$dpi, " dpi"), levels = names(STAGE_COLORS))
  save_csv(delta, "06_134_panelE_specimen_signature_deltas_v05.csv")

  drop_ci <- intersect(c("exploratory_95CI_lower_three_GSM", "exploratory_95CI_upper_three_GSM"), names(summ))
  if (length(drop_ci)) summ <- summ[, !(names(summ) %in% drop_ci), drop = FALSE]
  save_csv(summ, "07_134_panelE_stage_effect_and_consistency_v05.csv")

  set.seed(13401L)
  pE <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.38, linetype = "dashed", color = "#8D939A") +
    ggplot2::geom_point(data = delta, ggplot2::aes(x = stage, y = mean_signature_target_minus_other_logCPM, color = stage),
      position = ggplot2::position_jitter(width = 0.105, height = 0, seed = 13401L), size = 3.1, alpha = 0.95) +
    ggplot2::geom_point(data = summ, ggplot2::aes(x = stage, y = mean_of_three_GSM_signature_deltas), shape = 23, fill = "white", color = "#24292E", stroke = 1.2, size = 4.5) +
    ggplot2::facet_wrap(~comparison_label, ncol = 2) +
    ggplot2::scale_color_manual(values = STAGE_COLORS, guide = "none") +
    ggplot2::labs(x = NULL, y = "Candidate state - comparator: mean signature delta (logCPM)") +
    base_theme(10.8) +
    ggplot2::theme( strip.background = ggplot2::element_rect(fill = "#F2F4F5", color = NA), strip.text = ggplot2::element_text(face = "bold", size = 10.1),
      panel.spacing.x = grid::unit(0.35, "cm"), axis.text.x = ggplot2::element_text(size = 10.0))

  pE <- pE +
    ggplot2::geom_text(data = summ, ggplot2::aes(x = stage, y = mean_of_three_GSM_signature_deltas, label = sprintf("%+.3f", mean_of_three_GSM_signature_deltas)),
      nudge_y = 0.012, size = 3.5, color = "#222222", fontface = "bold", family = FONT_FAMILY, show.legend = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.23)))
  log_it("Panel E ready")

  spec145 <- load_csv(IN145_SPECIMEN)
  stage145 <- load_csv(IN145_STAGE)
  pair145 <- load_csv(IN145_PAIRWISE)
  omni145 <- load_csv(IN145_OMNIBUS)

  need_columns(spec145, c("gsm", "stage", "n_candidate_myeloid", "n_primary_target", "percent", "source"), "Step145 nine-specimen proportions")
  need_columns(stage145, c("stage", "n_specimens", "mean_percent", "median_percent", "min_percent", "max_percent"), "Step145 stage descriptive summary")
  need_columns(pair145, c("group_A", "group_B", "n_A", "n_B", "absolute_mean_difference_B_minus_A_percentage_points",
    "welch_95CI_lower_percentage_points", "welch_95CI_upper_percentage_points", "exact_two_sided_permutation_P"), "Step145 pairwise effect-size table")
  need_columns(omni145, c("statistic", "observed_statistic", "exact_two_sided_or_omnibus_P", "n_unique_group_allocations"), "Step145 three-group exact permutation table")

  EXPECTED_ALL9 <- paste0("GSM", 8614659:8614667)
  EXPECTED_NAIVE <- paste0("GSM", 8614659:8614661)
  assert(nrow(spec145) == 9L && !anyDuplicated(spec145$gsm) && setequal(as.character(spec145$gsm), EXPECTED_ALL9),
    "Step145 specimen table must contain exactly GSM8614659-GSM8614667 once each")
  assert(setequal(as.character(spec145$stage), c("Naive", "14dpi", "28dpi")), "Unexpected Step145 stage labels")
  assert(all(table(spec145$stage) == 3L), "Step145 must contain three independent GSMs in each stage")
  assert(nrow(stage145) == 3L && all(stage145$n_specimens == 3L) && setequal(as.character(stage145$stage), c("Naive", "14dpi", "28dpi")),
    "Step145 stage summary must contain Naive/14dpi/28dpi with n=3 each")
  assert(nrow(pair145) == 3L && all(pair145$n_A == 3L) && all(pair145$n_B == 3L), "Step145 pairwise table must contain three 3-vs-3 comparisons")
  assert(nrow(omni145) == 1L && is.finite(omni145$exact_two_sided_or_omnibus_P[1]), "Step145 omnibus result is missing or invalid")
  assert(all(is.finite(spec145$percent)) && all(spec145$percent >= 0 & spec145$percent <= 100), "Step145 specimen percentages are invalid")

  old_target <- comp[as.character(comp$group_131) == "TARGET_CLEC7A_HIGH_MICROGLIA", c("gsm", "n_cells"), drop = FALSE]
  old_den <- stats::aggregate(n_cells ~ gsm, comp, sum)
  old_target$denominator <- old_den$n_cells[match(old_target$gsm, old_den$gsm)]
  old_target$percent_131 <- 100 * old_target$n_cells / old_target$denominator
  si <- match(as.character(old_target$gsm), as.character(spec145$gsm))
  assert(!anyNA(si) && max(abs(old_target$percent_131 - spec145$percent[si])) < 1e-6, "Step145 14/28-dpi percentages differ from frozen Step131 values")
  assert(all(grepl("Frozen Step131", spec145$source[spec145$gsm %in% EXPECTED_GSM], fixed = TRUE)), "Step145 14/28 source labels do not identify frozen Step131 results")
  assert(all(grepl("Step145 naive frozen-reference transfer", spec145$source[spec145$gsm %in% EXPECTED_NAIVE], fixed = TRUE)),
    "Step145 naive source labels do not identify frozen-reference transfer")

  stage_map145 <- c("Naive" = "Naive", "14dpi" = "14 dpi", "28dpi" = "28 dpi")
  spec145$stage_plot <- factor(unname(stage_map145[as.character(spec145$stage)]), levels = names(EXT_STAGE_COLORS))
  stage145$stage_plot <- factor(unname(stage_map145[as.character(stage145$stage)]), levels = names(EXT_STAGE_COLORS))
  assert(!anyNA(spec145$stage_plot) && !anyNA(stage145$stage_plot), "Step145 stage plotting labels could not be mapped")

  calc145 <- stats::aggregate(percent ~ stage_plot, spec145, mean)
  mi <- match(as.character(stage145$stage_plot), as.character(calc145$stage_plot))
  assert(!anyNA(mi) && max(abs(stage145$mean_percent - calc145$percent[mi])) < 1e-6, "Step145 stage means do not match the nine specimen values")

  save_csv(spec145, "08_134_panelF_nine_specimen_proportions_v05.csv")
  save_csv(stage145, "09_134_panelF_stage_descriptive_summary_v05.csv")
  save_csv(pair145, "10_134_panelF_pairwise_effects_uncertainty_exactP_v05.csv")
  save_csv(omni145, "11_134_panelF_three_group_exact_permutation_v05.csv")

  set.seed(13404L)
  ymaxF <- max(60, ceiling((max(spec145$percent, stage145$mean_percent, na.rm = TRUE) + 5) / 10) * 10)
  pF <- ggplot2::ggplot() +
    ggplot2::geom_point(data = spec145, ggplot2::aes(x = stage_plot, y = percent, color = stage_plot), position = ggplot2::position_jitter(width = 0.09, height = 0, seed = 13404L),
      size = 3.2, alpha = 0.96) +
    ggplot2::geom_point(data = stage145, ggplot2::aes(x = stage_plot, y = mean_percent), shape = 23, fill = "white", color = "#24292E", stroke = 1.2, size = 4.7) +
    ggplot2::geom_text(data = stage145, ggplot2::aes(x = stage_plot, y = mean_percent, label = sprintf("%.2f%%", mean_percent)),
      nudge_y = 3.2, size = 3.5, color = "#222222", fontface = "bold", family = FONT_FAMILY, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = EXT_STAGE_COLORS, guide = "none", drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, ymaxF), breaks = pretty(c(0, ymaxF), n = 6), expand = ggplot2::expansion(mult = c(0, 0.03))) +
    ggplot2::labs(x = NULL, y = "Clec7a-high candidate myeloid state (%)") +
    base_theme(10.8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10.0))
  log_it("Panel F ready: nine independent specimens (Naive/14 dpi/28 dpi; n=3 each)")

  design <- paste(c("AAABBBCCCCCC", "AAABBBCCCCCC", "DDDDEEEEFFFF", "DDDDEEEEFFFF"), collapse = "\n")

  no_panel_heading <- ggplot2::theme( plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(), plot.tag.position = c(0.015, 0.995),
    plot.tag = ggplot2::element_text(face = "bold", size = 21, hjust = 0, vjust = 1, family = FONT_FAMILY) )
  pA <- pA + no_panel_heading
  pB <- pB + no_panel_heading
  pC <- pC + no_panel_heading
  pD <- pD + no_panel_heading
  pE <- pE + no_panel_heading
  pF <- pF + no_panel_heading
  fig <- patchwork::wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE, F = pF, design = design) +
    patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme( plot.tag.position = c(0.015, 0.995), plot.tag = ggplot2::element_text(face = "bold", size = 21,
      hjust = 0, vjust = 1, family = FONT_FAMILY)))
  out_base <- file.path(FIG_DIR, "Supplementary_Figure_S4_GSE281205_v05")
  if (PDF_CAIRO) {
    ggplot2::ggsave(paste0(out_base, ".pdf"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", bg = "white", device = grDevices::cairo_pdf, limitsize = FALSE)
  } else {
    ggplot2::ggsave(paste0(out_base, ".pdf"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", bg = "white", device = grDevices::pdf, useDingbats = FALSE, limitsize = FALSE)
  }
  log_it("Saved vector PDF")
  ggplot2::ggsave(paste0(out_base, ".png"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", dpi = FIG_DPI, bg = "white", device = "png", limitsize = FALSE)
  log_it("Saved ", FIG_DPI, " dpi PNG")
  ggplot2::ggsave(paste0(out_base, ".tiff"), plot = fig, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", dpi = FIG_DPI,
    bg = "white", device = "tiff", compression = "lzw", limitsize = FALSE)
  log_it("Saved ", FIG_DPI, " dpi LZW TIFF")
  save_panel(pA, "A")
  save_panel(pB, "B")
  save_panel(pC, "C")
  save_panel(pD, "D")
  save_panel(pE, "E")
  save_panel(pF, "F")

  legend <- c( "Supplementary Figure S4. External single-cell characterization and specimen-level evaluation of the Clec7a-high candidate myeloid state in GSE281205.",
    "(A) Harmony UMAP of 42,361 retained cells from the original six injured specimens (14 and 28 dpi), colored by evidence-based working major-lineage annotation (Step128). Unresolved/vascular-ambiguous groups are retained as such.",
    "(B) Harmony UMAP of 28,716 candidate myeloid cells from the original six injured specimens at resolution 0.10, colored by five working transcriptomic states (Steps130-131); state identity is an operational transcriptional annotation rather than lineage-tracing evidence.",
    "(C) Canonical-marker DotPlot across the five candidate myeloid states. Dot size indicates the percentage of cells with DecontX-estimated counts >=1; color indicates per-gene standardized mean log-normalized expression. These descriptive markers are not specimen-level differential-expression statistics.",
    "(D) Per-specimen composition of Step131 high-confidence myeloid states and QC-excluded/ambiguous candidate cells in the original six injured GSMs. Denominator: all candidate myeloid cells in each GSM; display is descriptive and does not by itself establish a statistically significant 14-to-28-dpi change.",
    "(E) Stage-matched within-specimen effects of the frozen 1,320-gene signature using original-integer-UMI specimen-level pseudobulk results (Step132 v02). Each dot is a GSM-specific target-minus-comparator mean logCPM difference over 1,301 mapped signature genes; white diamonds and adjacent numeric labels mark the three-GSM stage mean. Confidence intervals are not plotted; all independent specimen effects are shown. Left: Clec7a-high candidate myeloid state vs P2ry12-high microglia; right: Clec7a-high candidate myeloid state vs all other high-confidence myeloid cells. Positive means are not evidence that all signature genes replicate.",
    "(F) Clec7a-high candidate myeloid-state proportions after addition of three previously unused Young-naive biological replicates from the same GSE281205 cohort (Step145 v01). Each point is one independent GSM (Naive n=3, 14 dpi n=3, 28 dpi n=3); white diamonds and labels indicate stage means. The Young-naive specimens were assigned using the frozen major-lineage/myeloid-state references and fixed Step131 rules without naive-informed reclustering, marker selection, or threshold optimization. The original six 14/28-dpi specimen values were retained unchanged. Pairwise effect sizes, 95% confidence intervals, and exact permutation tests are reported in Supplementary Table S24; the three-group exact permutation test is an omnibus between-stage test, not a trajectory model.",
    "Panels A-E therefore retain the original six-specimen external analysis, whereas Panel F specifically evaluates enrichment relative to uninjured tissue using all nine specimens. GSE281205 still lacks acute early post-injury specimens, so this figure does not independently establish emergence across the complete early-to-late trajectory. The GEO 28-dpi age-metadata discrepancy remains documented in analysis provenance."
    )
  writeLines(legend, file.path(OUT, "12_134_Supplementary_Figure_S4_legend_v05.txt"), useBytes = TRUE)
  log_it("Saved figure legend")
  audit <- data.frame( check = c("atlas_n_cells_A", "candidate_myeloid_n_cells_BD", "n_injured_specimens_AE", "injured_stages_AE", "Step131_composition_n_cells", "n_total_specimens_F",
    "stages_F", "n_specimens_per_stage_F", "figure_panels", "Step132_primary_comparisons", "Step132_stage_summaries", "mapped_signature_genes", "Step145_pairwise_comparisons",
    "Step145_omnibus_exact_P", "Step145_unique_group_allocations", "naive_rule_refit", "reviewer2_focal_state_label", "camera_shown", "workflow_shown", "ci_bars_shown", "dpi"),
    value = c(EXPECTED_TOTAL, EXPECTED_MYELOID, length(EXPECTED_GSM), "14dpi;28dpi", sum(comp$n_cells), nrow(spec145), "Naive;14dpi;28dpi", 3, "A;B;C;D;E;F", 2, 4,
    unique(summ$n_signature_mapped), nrow(pair145), omni145$exact_two_sided_or_omnibus_P[1], omni145$n_unique_group_allocations[1], "NO - frozen reference/rules from Step145",
    STATE_LABELS[["1"]], "NO", "NO", "NO", FIG_DPI), stringsAsFactors = FALSE)
  save_csv(audit, "13_134_figure_validation_audit_v05.csv")
  capture.output(sessionInfo(), file = file.path(OUT, "14_134_sessionInfo_v05.txt"))
  log_it("SUCCESS | Step134 v05 | Supplementary Figure S4 completed")
  invisible(fig)
}

tryCatch(run_134(), error = function(e) {
  log_it("FAILED | Step134 v05 | ", conditionMessage(e))
  stop(e)
})

source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
set.seed(1234)

ROOT <- PROJECT_ROOT

STEP62_DIR <- file.path(ROOT, "results", "62_GSE195783_joint_embedding_batch_assessment_clustering")

STEP64_DIR <- file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation")

INPUT_OBJECT <- file.path(STEP62_DIR, "13_62_GSE195783_joint_unintegrated_clustered_v01.rds")

MARKER_FILE <- file.path(STEP64_DIR, "02_64_top_markers_sample_consistency_v01.csv")

AUTO_ANNOTATION_FILE <- file.path(STEP64_DIR, "05_64_provisional_transcriptional_region_annotation_v01.csv")

VERSION <- "v03"
CLUSTER_COL <- "prelim_cluster_r0.4"
FINAL_LABEL_COL <- "spatial_region_v03"

required_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "dplyr", "patchwork")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

write_csv_utf8 <- write_csv_utf8

save_pdf_png <- function(p, stem, width, height) {
    ggsave(filename = file.path(STEP64_DIR, paste0(stem, "_", VERSION, ".pdf")), plot = p, width = width, 
        height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
    ggsave(filename = file.path(STEP64_DIR, paste0(stem, "_", VERSION, ".png")), plot = p, width = width, 
        height = height, units = "in", dpi = 300, bg = "white")
}

collapse_markers <- function(df, n = 12L) {
    if (nrow(df) == 0L) {
        return("")
    }
    paste(head(df$gene, n), collapse = "; ")
}

required_files <- c(INPUT_OBJECT, MARKER_FILE, AUTO_ANNOTATION_FILE)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop("Required file(s) missing:\n", paste(missing_files, collapse = "\n"))
}

obj <- readRDS(INPUT_OBJECT)

markers <- read.csv(MARKER_FILE, stringsAsFactors = FALSE, check.names = FALSE)

auto_ann <- read.csv(AUTO_ANNOTATION_FILE, stringsAsFactors = FALSE, check.names = FALSE)

if (!inherits(obj, "Seurat")) {
  stop("Step 62 input is not a Seurat object.")
}

required_md <- c(CLUSTER_COL, "sample_label", "timepoint", "array_x", "array_y", "histology_registered")

missing_md <- setdiff(required_md, colnames(obj[[]]))

if (length(missing_md) > 0L) {
    stop("Step 62 object is missing metadata column(s): ", paste(missing_md, collapse = ", "))
}

if (ncol(obj) != 15204L) {
    stop("Expected 15,204 Step 62 spots, found ", ncol(obj), ".")
}

if (any(obj$histology_registered)) {
  stop("H&E registration must remain FALSE.")
}

required_marker_cols <- c("cluster", "gene", "avg_log2FC_std", "n_samples_evaluable", "fraction_samples_positive", 
    "robust_directional_support")

missing_marker_cols <- setdiff(required_marker_cols, colnames(markers))

if (length(missing_marker_cols) > 0L) {
    stop("Marker file is missing column(s): ", paste(missing_marker_cols, collapse = ", "))
}

required_auto_cols <- c("cluster", "provisional_program", "top_program_score", "second_program", "second_program_score", 
    "score_margin", "annotation_status", "dominant_timepoint", "dominant_timepoint_fraction")

missing_auto_cols <- setdiff(required_auto_cols, colnames(auto_ann))

if (length(missing_auto_cols) > 0L) {
    stop("Automatic annotation file is missing column(s): ", paste(missing_auto_cols, collapse = ", "))
}

curated_map <- data.frame(cluster = as.character(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)), spatial_region_v03 = c("Oligodendrocyte / white-matter-like", 
    "Oligodendrocyte / white-matter-like", "Myeloid/immune–lipid-handling enriched", "Cholinergic/motor-neuron-like neuronal", 
    "Inhibitory/interneuron-like neuronal", "Acute proliferative injury-associated myeloid", "Peripheral-myelin / Schwann-like region", 
    "Neuronal subtype region, 7 dpi-enriched", "Fibrotic/ECM–myeloid mixed injury niche", "IFN-responsive inflammatory/astroglial region", 
    "Perivascular / vSMC–pericyte region"), short_label = c("Oligo/WM-like", "Oligo/WM-like", "Myeloid + lipid handling", 
    "Cholinergic/motor neuron-like", "Interneuron-like", "Acute proliferative injury myeloid", "Peripheral myelin/Schwann-like", 
    "7 dpi neuronal subtype", "Fibrotic/ECM–myeloid niche", "IFN-inflammatory/astroglial", "Perivascular vSMC/pericyte"), 
    confidence = c("moderate", "moderate", "moderate", "high", "high", "high", "high", "moderate", "moderate", 
        "moderate", "high"), rationale = c(paste0("Automatic oligodendrocyte/white-matter program retained; ", 
        "robust myelin-associated markers include Mog and Plp1."), paste0("Automatic oligodendrocyte/white-matter program retained; ", 
        "broad oligodendroglial/myelin transcriptional pattern."), paste0("Myeloid/immune and lipid-handling scores are both high; robust markers ", 
        "include Cst7, Apoc1, Lpl, Apoe, C1qa, Trem2 and Tyrobp. ", "Kept at region/program level rather than a pure microglial subtype."), 
        paste0("Strong neuronal program with cholinergic/motor-neuron-associated markers ", "including Mnx1, Chat, Slc5a7 and Slc18a3."), 
        paste0("Strong neuronal program with inhibitory/interneuron-associated markers ", "including Slc32a1, Slc6a5, Lamp5, Pax2 and Kcnc1."), 
        paste0("Strong injury-response program with Arg1, Hmox1, Ccl2 and Spp1 plus ", "proliferation genes Mki67, Ube2c, Top2a and Ccna2; predominantly 3 dpi."), 
        paste0("Mpz, Prx, Ncmap, Cldn19, Pmp2 and Dhh support a peripheral-myelin / ", "Schwann-like pattern more strongly than the automatic fibrotic/ECM label."), 
        paste0("Strong neuronal subtype markers including Tac2, Vip, Tac1, Sst and Penk; ", "strongly restricted to 7 dpi, so the label remains descriptive."), 
        paste0("Strongest markers are dominated by ECM/fibrotic genes (Col10a1, Mmp13, ", "Lrrc15, Postn, Col1a1, Col3a1, Fap) with concurrent myeloid markers ", 
            "(Lyz2, Ms4a7, Clec12a, Mmp12), supporting a mixed injury niche."), paste0("Automatic astroglial score retained only as one component; strongest ", 
            "markers are IFN/inflammatory genes including Cxcl10, Iigp1, Rsad2, Oas3, ", "Ifi44, Ifit3, Isg15 and Irf7."), 
        paste0("Strong vascular smooth-muscle/perivascular markers including Myh11, ", "Acta2, Lmod1, Myl9, Gja4 and Cnn1.")), 
    stringsAsFactors = FALSE)

markers$cluster <- as.character(markers$cluster)
markers$robust_directional_support <- as.logical(markers$robust_directional_support)

marker_evidence <- markers %>% group_by(cluster) %>% arrange(desc(robust_directional_support), desc(avg_log2FC_std), 
    .by_group = TRUE) %>% summarise(n_top_markers_tested = n(), n_robust_directional = sum(robust_directional_support, 
    na.rm = TRUE), robust_marker_fraction = mean(robust_directional_support, na.rm = TRUE), robust_top_markers_consistent = collapse_markers(data.frame(gene = gene[robust_directional_support %in% 
    TRUE], stringsAsFactors = FALSE), n = 12L), .groups = "drop")

auto_ann$cluster <- as.character(auto_ann$cluster)

if ("robust_top_markers" %in% colnames(auto_ann)) {
    auto_ann <- auto_ann %>% rename(robust_top_markers_v01 = robust_top_markers)
}

final_table <- curated_map %>% left_join(auto_ann, by = "cluster") %>% left_join(marker_evidence, by = "cluster") %>% 
    select(cluster, spatial_region_v03, short_label, confidence, rationale, provisional_program, top_program_score, 
        second_program, second_program_score, score_margin, annotation_status, dominant_timepoint, dominant_timepoint_fraction, 
        n_top_markers_tested, n_robust_directional, robust_marker_fraction, any_of("robust_top_markers_v01"), 
        robust_top_markers_consistent) %>% arrange(as.numeric(cluster))

required_final_cols <- c("cluster", FINAL_LABEL_COL, "confidence", "dominant_timepoint", "dominant_timepoint_fraction", 
    "robust_top_markers_consistent")

missing_final_cols <- setdiff(required_final_cols, colnames(final_table))

if (length(missing_final_cols) > 0L) {
    stop("Final annotation table is missing column(s): ", paste(missing_final_cols, collapse = ", "), 
        "\nColumns present: ", paste(colnames(final_table), collapse = ", "))
}

if (nrow(final_table) != 11L) {
    stop("Curated annotation table does not contain exactly 11 clusters.")
}

if (any(is.na(final_table$spatial_region_v03))) {
    stop("At least one cluster is missing a curated v03 label.")
}

FINAL_TABLE_FILE <- file.path(STEP64_DIR, paste0("11_64_curated_transcriptional_region_annotation_", 
    VERSION, ".csv"))

write_csv_utf8(final_table, FINAL_TABLE_FILE)

cluster_vec <- as.character(obj[[CLUSTER_COL]][, 1])

label_lookup <- setNames(curated_map$spatial_region_v03, curated_map$cluster)

short_lookup <- setNames(curated_map$short_label, curated_map$cluster)

confidence_lookup <- setNames(curated_map$confidence, curated_map$cluster)

obj[[FINAL_LABEL_COL]] <- unname(label_lookup[cluster_vec])

obj$spatial_region_short_v03 <- unname(short_lookup[cluster_vec])

obj$spatial_region_confidence_v03 <- unname(confidence_lookup[cluster_vec])

if (any(is.na(obj[[FINAL_LABEL_COL]][, 1]))) {
    stop("Failed to map at least one spot to a curated v03 label.")
}

obj$spatial_cluster_step62 <- cluster_vec

obj@misc$GSE195783_step64_v03 <- list(step = 64L, version = VERSION, purpose = "Curated refinement of v01 automatic transcriptional-region labels", 
    marker_testing_rerun = FALSE, program_scoring_rerun = FALSE, clustering_rerun = FALSE, primary_representation = "Step62 unintegrated", 
    harmony_primary = FALSE, histology_registered = FALSE, label_column = FINAL_LABEL_COL, interpretation = "transcriptionally defined spatial regions/programs, not H&E-defined compartments")

OBJECT_FILE <- file.path(STEP64_DIR, paste0("12_64_GSE195783_spatial_region_annotated_", VERSION, ".rds"))

saveRDS(obj, file = OBJECT_FILE, compress = TRUE)

region_levels <- unique(curated_map$spatial_region_v03)

region_palette <- c(`Oligodendrocyte / white-matter-like` = "#4C78A8", `Myeloid/immune–lipid-handling enriched` = "#B279A2", 
    `Cholinergic/motor-neuron-like neuronal` = "#F58518", `Inhibitory/interneuron-like neuronal` = "#E45756", 
    `Acute proliferative injury-associated myeloid` = "#D62728", `Peripheral-myelin / Schwann-like region` = "#72B7B2", 
    `Neuronal subtype region, 7 dpi-enriched` = "#FF9DA6", `Fibrotic/ECM–myeloid mixed injury niche` = "#8C564B", 
    `IFN-responsive inflammatory/astroglial region` = "#54A24B", `Perivascular / vSMC–pericyte region` = "#9D755D")

map_df <- obj[[]]

map_df[[FINAL_LABEL_COL]] <- factor(map_df[[FINAL_LABEL_COL]], levels = region_levels)

sample_order <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1", 
    "28dpi_rep2")

sample_order <- sample_order[sample_order %in% unique(as.character(map_df$sample_label))]

map_df$sample_label <- factor(as.character(map_df$sample_label), levels = sample_order)

p_spatial <- ggplot(map_df, aes(x = array_x, y = array_y, color = .data[[FINAL_LABEL_COL]])) + geom_point(size = 0.48, 
    alpha = 0.95) + coord_equal() + scale_y_reverse() + scale_color_manual(values = region_palette, drop = FALSE) + 
    facet_wrap(~sample_label, ncol = 4) + labs(title = "Curated transcriptional regions across GSE195783", 
    subtitle = "Visium V1 array coordinates only; no H&E/pixel registration", x = "Visium array X", y = "Visium array Y", 
    color = NULL) + theme_classic(base_family = "Arial", base_size = 9) + theme(plot.title = element_text(size = 12), 
    plot.subtitle = element_text(size = 8.5), strip.text = element_text(size = 8), legend.text = element_text(size = 7.5), 
    legend.position = "right")

save_pdf_png(p = p_spatial, stem = "13_64_curated_transcriptional_region_spatial_map", width = 12.5, 
    height = 7.5)

if (!"umap" %in% names(obj@reductions)) {
    stop("Step 62 UMAP is missing from the object.")
}

p_umap <- DimPlot(object = obj, reduction = "umap", group.by = FINAL_LABEL_COL, cols = region_palette, 
    pt.size = 0.28, raster = TRUE, label = FALSE) + labs(title = "Curated transcriptional-region labels on Step 62 UMAP", 
    color = NULL) + theme_classic(base_family = "Arial", base_size = 9) + theme(plot.title = element_text(size = 12), 
    legend.text = element_text(size = 7.5), legend.position = "right")

save_pdf_png(p = p_umap, stem = "14_64_curated_transcriptional_region_UMAP", width = 9, height = 6.2)

cluster_count_before <- length(unique(cluster_vec))

cluster_count_after <- length(unique(obj$spatial_cluster_step62))

final_checks <- data.frame(check = c("Step62 unintegrated object was used as primary input", "Input contains 15204 spots", 
    "Original Step62 cluster IDs were preserved", "No marker testing was rerun", "No program scoring was rerun", 
    "No clustering was rerun", "All 11 clusters received curated v03 labels", "Curated label metadata were added to all spots", 
    "Spatial map was generated", "UMAP label map was generated", "No H&E registration is claimed"), pass = c(grepl("62_GSE195783", 
    INPUT_OBJECT), ncol(obj) == 15204L, cluster_count_before == cluster_count_after, TRUE, TRUE, TRUE, 
    nrow(final_table) == 11L, all(!is.na(obj[[FINAL_LABEL_COL]][, 1])), file.exists(file.path(STEP64_DIR, 
        paste0("13_64_curated_transcriptional_region_spatial_map_", VERSION, ".pdf"))), file.exists(file.path(STEP64_DIR, 
        paste0("14_64_curated_transcriptional_region_UMAP_", VERSION, ".pdf"))), all(obj$histology_registered == 
        FALSE)), stringsAsFactors = FALSE)

QC_FILE <- file.path(STEP64_DIR, paste0("15_64_GSE195783_annotation_refinement_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_checks, QC_FILE)

annotation_lines <- apply(final_table, 1, function(z) {
    paste0("  Cluster ", z[["cluster"]], ": ", z[["spatial_region_v03"]], " [confidence=", z[["confidence"]], 
        "]", "; dominant_timepoint=", z[["dominant_timepoint"]], " (", sprintf("%.1f%%", 100 * as.numeric(z[["dominant_timepoint_fraction"]])), 
        ")", "; robust markers=", z[["robust_top_markers_consistent"]])
})

REPORT_FILE <- file.path(STEP64_DIR, paste0("16_64_GSE195783_annotation_refinement_REPORT_", VERSION, 
    ".txt"))

report_con <- file(REPORT_FILE, open = "wt", encoding = "UTF-8")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = report_con)

close(report_con)

print(final_table[, c("cluster", "spatial_region_v03", "confidence", "dominant_timepoint", "dominant_timepoint_fraction")])
print(final_checks)

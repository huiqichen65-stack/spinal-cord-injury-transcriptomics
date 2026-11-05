source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234

FINAL_CLUSTER_COL <- "myeloid_res.0.1"

MIN_PCT <- 0.1
LOGFC_THRESHOLD <- 0.25
PADJ_THRESHOLD <- 0.05

MAX_PLOT_CELLS <- 50000
TOP_MARKERS_PER_CLUSTER_HEATMAP <- 5
TOP_MARKERS_PER_CLUSTER_SUMMARY <- 10

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "dplyr", "tidyr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\nPlease install them before rerunning.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "16_myeloid_cluster_resolution_assessment", "16_myeloid_cluster_resolution_assessment_3scRNA_v03.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^16_myeloid_cluster_resolution_assessment_3scRNA_v03\\.rds$", 
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Input file not found:\n", "16_myeloid_cluster_resolution_assessment_3scRNA_v03.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "17_myeloid_marker_annotation_validation")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "17_myeloid_marker_annotation_validation_v02.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

md <- obj[[]]

if (!FINAL_CLUSTER_COL %in% colnames(md)) {
  stop("The selected clustering column is absent: ", FINAL_CLUSTER_COL)
}

required_meta <- c("source_dataset", "parent_major_celltype")

missing_meta <- setdiff(required_meta, colnames(md))

if (length(missing_meta) > 0) {
  stop("Missing required metadata field(s): ", paste(missing_meta, collapse = ", "))
}

obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) {
  log_msg("JoinLayers skipped/failed safely: ", conditionMessage(e))
  obj
})

DefaultAssay(obj) <- "RNA"

rna_layers <- tryCatch(Layers(obj[["RNA"]]), error = function(e) character())

if (!"data" %in% rna_layers) {
    obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, 
        verbose = TRUE)
}

cluster_raw <- as.character(obj[[FINAL_CLUSTER_COL]][, 1])

cluster_levels <- as.character(sort(unique(as.integer(cluster_raw))))

obj$myeloid_cluster_raw <- factor(cluster_raw, levels = cluster_levels)

Idents(obj) <- "myeloid_cluster_raw"

capture.output(print(table(obj$myeloid_cluster_raw)), file = log_file, append = TRUE)

markers <- FindAllMarkers(object = obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox", min.pct = MIN_PCT,
  logfc.threshold = LOGFC_THRESHOLD, return.thresh = PADJ_THRESHOLD, densify = FALSE, verbose = TRUE)

if (nrow(markers) == 0) {
  stop("FindAllMarkers returned no significant positive markers.")
}

markers$cluster <- as.character(markers$cluster)

fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")

fc_col <- fc_candidates[fc_candidates %in% colnames(markers)][1]

if (is.na(fc_col)) {
  stop("Could not identify the fold-change column in FindAllMarkers output.")
}

if (!"gene" %in% colnames(markers)) {
  markers$gene <- rownames(markers)
}

markers <- markers %>% arrange(as.integer(cluster), p_val_adj, desc(.data[[fc_col]]))

markers_file <- file.path(out_dir, "17_myeloid_FindAllMarkers_res0.10_v02.csv")

write.csv(markers, markers_file, row.names = FALSE, fileEncoding = "UTF-8")

technical_pattern <- paste0("^mt-|", "^Rpl|^Rps|", "^Gm[0-9]+$|", "^Hba-|^Hbb-|", "^Malat1$")

markers_annotation <- markers %>% filter(!grepl(technical_pattern, gene, ignore.case = FALSE))

top_summary <- markers_annotation %>% group_by(cluster) %>% arrange(p_val_adj, desc(.data[[fc_col]]), 
    .by_group = TRUE) %>% slice_head(n = TOP_MARKERS_PER_CLUSTER_SUMMARY) %>% summarise(top_markers = paste(gene, 
    collapse = "; "), .groups = "drop")

top_heat <- markers_annotation %>% group_by(cluster) %>% arrange(p_val_adj, desc(.data[[fc_col]]), .by_group = TRUE) %>% 
    slice_head(n = TOP_MARKERS_PER_CLUSTER_HEATMAP) %>% ungroup()

heat_genes <- unique(top_heat$gene)

heat_genes <- heat_genes[heat_genes %in% rownames(obj)]

md <- obj[[]] %>% mutate(cluster = as.character(myeloid_cluster_raw))

cluster_size <- md %>% count(cluster, name = "cluster_n")

parent_detail <- md %>% count(cluster, parent_major_celltype, name = "n_cells") %>% group_by(cluster) %>% 
    mutate(fraction = n_cells/sum(n_cells)) %>% ungroup()

parent_dom <- parent_detail %>% group_by(cluster) %>% slice_max(order_by = fraction, n = 1, with_ties = FALSE) %>% 
    ungroup() %>% transmute(cluster, dominant_parent = as.character(parent_major_celltype), dominant_parent_fraction = fraction)

dataset_detail <- md %>% count(cluster, source_dataset, name = "n_cells") %>% group_by(cluster) %>% mutate(fraction = n_cells/sum(n_cells)) %>% 
    ungroup()

dataset_dom <- dataset_detail %>% group_by(cluster) %>% slice_max(order_by = fraction, n = 1, with_ties = FALSE) %>% 
    ungroup() %>% transmute(cluster, dominant_dataset = as.character(source_dataset), dominant_dataset_fraction = fraction)

dataset_support <- dataset_detail %>% group_by(cluster) %>% summarise(n_datasets_present = sum(n_cells > 
    0), n_datasets_ge1pct = sum(fraction >= 0.01), .groups = "drop")

has_timepoint <- ("timepoint_main" %in% colnames(md) && any(!is.na(md$timepoint_main)))

if (has_timepoint) {
    time_detail <- md %>% filter(!is.na(timepoint_main)) %>% count(cluster, timepoint_main, name = "n_cells") %>% 
        group_by(cluster) %>% mutate(fraction = n_cells/sum(n_cells)) %>% ungroup()
    time_dom <- time_detail %>% group_by(cluster) %>% slice_max(order_by = fraction, n = 1, with_ties = FALSE) %>% 
        ungroup() %>% transmute(cluster, dominant_timepoint = as.character(timepoint_main), dominant_timepoint_fraction = fraction)
    time_support <- time_detail %>% group_by(cluster) %>% summarise(n_timepoints_present = sum(n_cells > 
        0), n_timepoints_ge1pct = sum(fraction >= 0.01), .groups = "drop")
} else {
    time_detail <- data.frame()
    time_dom <- data.frame(cluster = cluster_levels, dominant_timepoint = NA_character_, dominant_timepoint_fraction = NA_real_, 
        stringsAsFactors = FALSE)
    time_support <- data.frame(cluster = cluster_levels, n_timepoints_present = NA_integer_, n_timepoints_ge1pct = NA_integer_, 
        stringsAsFactors = FALSE)
}

validation_summary <- cluster_size %>% left_join(parent_dom, by = "cluster") %>% left_join(dataset_dom, 
    by = "cluster") %>% left_join(dataset_support, by = "cluster") %>% left_join(time_dom, by = "cluster") %>% 
    left_join(time_support, by = "cluster") %>% left_join(top_summary, by = "cluster") %>% mutate(flag_small_lt200 = cluster_n < 
    200, flag_small_lt500 = cluster_n < 500, flag_dataset_ge90 = dominant_dataset_fraction >= 0.9, flag_dataset_ge99 = dominant_dataset_fraction >= 
    0.99, flag_timepoint_ge90 = ifelse(is.na(dominant_timepoint_fraction), FALSE, dominant_timepoint_fraction >= 
    0.9), flag_parent_mixed = dominant_parent_fraction < 0.8, review_priority = ifelse(flag_small_lt500 | 
    flag_dataset_ge90 | flag_timepoint_ge90 | flag_parent_mixed, "Review", "Standard")) %>% arrange(as.integer(cluster))

validation_file <- file.path(out_dir, "17_myeloid_cluster_validation_summary_v02.csv")

write.csv(validation_summary, validation_file, row.names = FALSE, fileEncoding = "UTF-8")

capture.output(print(validation_summary %>% filter(review_priority == "Review") %>% select(cluster, cluster_n, 
    dominant_parent, dominant_parent_fraction, dominant_dataset, dominant_dataset_fraction, dominant_timepoint, 
    dominant_timepoint_fraction, top_markers)), file = log_file, append = TRUE)

marker_sets <- list(`Microglia identity` = c("P2ry12", "Tmem119", "Sall1", "Hexb", "Cx3cr1"), `Monocyte/Macrophage` = c("Ccr2",
  "Ly6c2", "Ms4a7", "Lgals3", "Fn1"), `APC/DC-like` = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Itgax"),
  `Inflammatory/acute` = c("Il1b", "Tnf", "Cxcl2", "Nfkbia", "S100a4"), `Interferon response` = c("Isg15",
    "Ifit1", "Ifit3", "Irf7", "Stat1"), Proliferating = c("Mki67", "Top2a", "Ccnb2", "Pclaf", "Ube2c"),
  `Repair/phagolysosomal` = c("Gpnmb", "Arg1", "Mrc1", "C1qa", "C1qc"), `Lipid-handling reference` = c("Lpl",
    "Apoc4", "Apoe", "Abca1", "Trem2"), `Contamination check` = c("S100a8", "S100a9", "Ly6g", "Retnlg",
    "Cd3d", "Nkg7", "Ms4a1", "Cd79a", "Mbp", "Pecam1", "Hba-a1", "Hbb-bs"))

marker_sets_present <- lapply(marker_sets, function(x) {
  x[x %in% rownames(obj)]
})

marker_sets_present <- marker_sets_present[lengths(marker_sets_present) > 0]

marker_order <- unique(unlist(marker_sets_present, use.names = FALSE))

marker_group_lengths <- vapply(marker_sets_present, length, FUN.VALUE = integer(1))

marker_group_ends <- cumsum(marker_group_lengths)

separator_x <- marker_group_ends[-length(marker_group_ends)] + 0.5

if (length(heat_genes) == 0) {
  stop("No top-marker genes were available for the heatmap.")
}

data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

avg_mat <- matrix(NA_real_, nrow = length(heat_genes), ncol = length(cluster_levels), dimnames = list(heat_genes,
  cluster_levels))

for (cl in cluster_levels) {
    cl_cells <- WhichCells(obj, idents = cl)
    avg_mat[, cl] <- Matrix::rowMeans(data_mat[heat_genes, cl_cells, drop = FALSE])
}

z_mat <- t(scale(t(avg_mat)))

z_mat[!is.finite(z_mat)] <- 0

z_mat <- pmax(pmin(z_mat, 2), -2)

heat_df <- as.data.frame(as.table(z_mat), stringsAsFactors = FALSE)

colnames(heat_df) <- c("gene", "cluster", "z")

heat_df$cluster <- factor(heat_df$cluster, levels = cluster_levels)

heat_df$gene <- factor(heat_df$gene, levels = rev(heat_genes))

rm(data_mat, avg_mat, z_mat)

gc(verbose = FALSE)

set.seed(SEED)

plot_cells <- Cells(obj)

if (length(plot_cells) > MAX_PLOT_CELLS) {
  plot_cells <- sample(plot_cells, MAX_PLOT_CELLS)
}

obj_plot <- subset(obj, cells = plot_cells)

Idents(obj_plot) <- "myeloid_cluster_raw"

diagnostic_pdf <- file.path(out_dir, "17_myeloid_marker_annotation_diagnostics_v02.pdf")

pdf(diagnostic_pdf, width = 12, height = 8.5, onefile = TRUE)

p_cluster <- DimPlot(obj_plot, reduction = "umap_myeloid_harmony30", group.by = "myeloid_cluster_raw",
  label = TRUE, repel = TRUE, pt.size = 0.12, raster = TRUE) + ggtitle("Candidate raw myeloid clustering | resolution = 0.10") +
  theme_classic(base_size = 11) + theme(legend.position = "none")

p_parent <- DimPlot(obj_plot, reduction = "umap_myeloid_harmony30", group.by = "parent_major_celltype",
  pt.size = 0.12, raster = TRUE) + ggtitle("Parent major-cell-type annotation") + theme_classic(base_size = 11)

print(p_cluster + p_parent + plot_layout(widths = c(1, 1)))

p_dot <- DotPlot(obj, features = marker_order, assay = "RNA", group.by = "myeloid_cluster_raw", dot.scale = 5.2,
  scale = TRUE) + RotatedAxis() + labs(title = "Canonical and state-reference markers across resolution-0.10 clusters",
  x = NULL, y = "Raw cluster", color = "Average\nexpression", size = "Percent\nexpressed") + theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(size = 7.3), axis.text.y = element_text(size = 8))

for (xx in separator_x) {
  p_dot <- p_dot + geom_vline(xintercept = xx, linetype = 3, linewidth = 0.3, color = "grey55")
}

print(p_dot)

p_heat <- ggplot(heat_df, aes(x = cluster, y = gene, fill = z)) + geom_tile() + scale_fill_gradient2(low = "#3B6FB6",
  mid = "white", high = "#D94F4F", midpoint = 0, limits = c(-2, 2), name = "Row Z") + labs(title = paste0("Top ",
  TOP_MARKERS_PER_CLUSTER_HEATMAP, " annotation-oriented markers per raw cluster"), x = "Raw cluster",
  y = NULL) + theme_classic(base_size = 8) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5.6),
  axis.ticks = element_blank())

print(p_heat)

dataset_plot_df <- dataset_detail %>% mutate(cluster = factor(cluster, levels = cluster_levels))

p_dataset <- ggplot(dataset_plot_df, aes(x = cluster, y = fraction, fill = source_dataset)) + geom_col(width = 0.82) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) {
    paste0(round(x * 100), "%")
  }, expand = c(0, 0)) + labs(title = "Source-dataset support of resolution-0.10 clusters", x = "Raw cluster",
  y = "Within-cluster fraction", fill = "Dataset") + theme_classic(base_size = 9)

if (has_timepoint) {
    time_plot_df <- time_detail %>% mutate(cluster = factor(cluster, levels = cluster_levels))
    p_time <- ggplot(time_plot_df, aes(x = cluster, y = fraction, fill = timepoint_main)) + geom_col(width = 0.82) + 
        scale_y_continuous(limits = c(0, 1), labels = function(x) {
            paste0(round(x * 100), "%")
        }, expand = c(0, 0)) + labs(title = "Time-point distribution of resolution-0.10 clusters", x = "Raw cluster", 
        y = "Within-cluster fraction", fill = "Time point") + theme_classic(base_size = 9)
    print(p_dataset/p_time + plot_layout(heights = c(1, 1)))
} else {
    print(p_dataset)
}

parent_plot_df <- parent_detail %>% mutate(cluster = factor(cluster, levels = cluster_levels))

p_parent_comp <- ggplot(parent_plot_df, aes(x = cluster, y = fraction, fill = parent_major_celltype)) +
  geom_col(width = 0.82) + scale_y_continuous(limits = c(0, 1), labels = function(x) {
  paste0(round(x * 100), "%")
}, expand = c(0, 0)) + labs(title = "Parent-annotation support of resolution-0.10 clusters", x = "Raw cluster",
  y = "Within-cluster fraction", fill = "Parent annotation") + theme_classic(base_size = 10)

print(p_parent_comp)

feature_candidates <- c("P2ry12", "Ccr2", "Cd74", "Il1b", "Isg15", "Mki67", "Gpnmb", "Lpl", "Apoc4")

feature_present <- feature_candidates[feature_candidates %in% rownames(obj_plot)]

if (length(feature_present) > 0) {
    p_feature <- FeaturePlot(obj_plot, features = feature_present, reduction = "umap_myeloid_harmony30", 
        ncol = 3, order = TRUE, pt.size = 0.12, raster = TRUE) & theme_classic(base_size = 8)
    print(p_feature)
}

dev.off()

rm(obj_plot, plot_cells)

gc(verbose = FALSE)

out_rds <- file.path(out_dir, "17_myeloid_marker_annotation_validation_3scRNA_v02.rds")

saveRDS(obj, out_rds, compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)

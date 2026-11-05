source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
RAW_CLUSTER_COL <- "myeloid_cluster_raw"
UMAP_REDUCTION <- "umap_myeloid_harmony30"

MIN_PCT <- 0.10
LOGFC_THRESHOLD <- 0.25
PADJ_THRESHOLD <- 0.05

MAX_PLOT_CELLS <- 50000

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

expected_input <- file.path(project_dir, "results", "17_myeloid_marker_annotation_validation", "17_myeloid_marker_annotation_validation_3scRNA_v02.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^17_myeloid_marker_annotation_validation_3scRNA_v02\\.rds$", 
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Input file not found:\n", "17_myeloid_marker_annotation_validation_3scRNA_v02.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "18_myeloid_annotation_cleanup")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "18_myeloid_annotation_cleanup_v01.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

md <- obj[[]]

if (!RAW_CLUSTER_COL %in% colnames(md)) {
    stop("Missing raw cluster metadata column: ", RAW_CLUSTER_COL)
}

if (!UMAP_REDUCTION %in% names(obj@reductions)) {
    stop("Missing UMAP reduction: ", UMAP_REDUCTION)
}

required_meta <- c("source_dataset", "parent_major_celltype")

missing_meta <- setdiff(required_meta, colnames(md))

if (length(missing_meta) > 0) {
    stop("Missing required metadata field(s): ", paste(missing_meta, collapse = ", "))
}

raw_cluster <- as.character(obj[[RAW_CLUSTER_COL]][, 1])

raw_levels <- as.character(sort(unique(as.integer(raw_cluster))))

expected_raw <- as.character(0:14)

if (!setequal(raw_levels, expected_raw)) {
    stop("Unexpected raw-cluster set.\n", "Expected 0:14 from resolution 0.10, but found: ", paste(raw_levels, 
        collapse = ", "))
}

capture.output(print(table(obj[[RAW_CLUSTER_COL]][, 1])), file = log_file, append = TRUE)

mapping_df <- data.frame(
  raw_cluster = as.character(0:14),
  action = c(
    "keep_merge",
    "keep",
    "keep",
    "keep_merge",
    "keep_merge",
    "keep",
    "keep_merge",
    "keep",
    "keep",
    "keep_merge",
    "remove",
    "keep_merge",
    "remove",
    "keep_merge",
    "keep_merge"
  ),
  myeloid_state = c(
    "Clec7a+ lipid-associated Microglia",
    "MHC-II+ antigen-presenting Macrophage",
    "Gpnmb+ tissue-repair Macrophage",
    "Activated inflammatory Microglia",
    "Activated inflammatory Microglia",
    "Ccl12+ acute Microglia",
    "Homeostatic Microglia",
    "Proliferating Microglia",
    "Tnf+ inflammatory Microglia",
    "Acute-response Macrophage",
    NA,
    "Homeostatic Microglia",
    NA,
    "Clec7a+ lipid-associated Microglia",
    "Acute-response Macrophage"
  ),
  rationale = c(
    "Clec7a/Lpl/Apoc4/Apoe-associated activated microglial program",
    "H2-Aa/H2-Ab1/H2-Eb1/Cd74/Mrc1 antigen-presentation program",
    "Gpnmb/Atp6v0d2/Fabp4 tissue-repair/phagolysosomal program",
    "Activated inflammatory microglial branch; merged with raw cluster 4",
    "Activated inflammatory microglial branch; merged with raw cluster 3",
    "Ccl12-enriched acute microglial state",
    "P2ry12/Tmem119-enriched homeostatic microglia; merged with raw cluster 11",
    "Mki67/Pclaf/Top2a proliferating microglia",
    "Tnf-enriched inflammatory microglial state",
    "S100a8/S100a9 acute-response macrophage; merged with raw cluster 14",
    "B-lineage contamination (Ms4a1/Cd79a/Cd19/Pou2af1)",
    "P2ry12/Tmem119/Hexb homeostatic microglia; merged with raw cluster 6",
    "Erythroid/vascular contamination (Hb/Alas2/Spta1/Cxcl12-related)",
    "Small Apoe+ activated microglial branch; merged with raw cluster 0",
    "Nos2/Vcan acute-response macrophage; merged with raw cluster 9"
  ),
  stringsAsFactors = FALSE
)

write.csv(mapping_df, file.path(out_dir, "18_myeloid_annotation_mapping_v01.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

capture.output(print(mapping_df), file = log_file, append = TRUE)

remove_clusters <- mapping_df$raw_cluster[mapping_df$action == "remove"]

remove_cells <- Cells(obj)[as.character(obj[[RAW_CLUSTER_COL]][, 1]) %in% remove_clusters]

removed_md <- obj[[]][remove_cells, , drop = FALSE]

removed_df <- data.frame(cell = rownames(removed_md), raw_cluster = as.character(removed_md[[RAW_CLUSTER_COL]]), 
    parent_major_celltype = as.character(removed_md$parent_major_celltype), source_dataset = as.character(removed_md$source_dataset), 
    stringsAsFactors = FALSE)

if ("timepoint_main" %in% colnames(removed_md)) {
    removed_df$timepoint_main <- as.character(removed_md$timepoint_main)
}

write.csv(removed_df, file.path(out_dir, "18_myeloid_removed_cells_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

keep_cells <- setdiff(Cells(obj), remove_cells)

obj_clean <- subset(obj, cells = keep_cells)

rm(obj, md, removed_md, keep_cells)

gc(verbose = FALSE)

capture.output(print(table(removed_df$raw_cluster)), file = log_file, append = TRUE)

raw_clean <- as.character(obj_clean[[RAW_CLUSTER_COL]][, 1])

state_lookup <- setNames(mapping_df$myeloid_state, mapping_df$raw_cluster)

state_vec <- unname(state_lookup[raw_clean])

if (anyNA(state_vec)) {
    bad <- unique(raw_clean[is.na(state_vec)])
    stop("Some retained raw clusters have no myeloid-state mapping: ", paste(bad, collapse = ", "))
}

state_levels <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

obj_clean$myeloid_state <- factor(state_vec, levels = state_levels)

Idents(obj_clean) <- "myeloid_state"

capture.output(print(table(obj_clean$myeloid_state)), file = log_file, append = TRUE)

md_clean <- obj_clean[[]] %>% mutate(myeloid_state = as.character(myeloid_state), raw_cluster = as.character(.data[[RAW_CLUSTER_COL]]))

state_count <- md_clean %>% count(myeloid_state, name = "n_cells") %>% mutate(fraction_of_clean_myeloid = n_cells/sum(n_cells))

raw_mapping_count <- md_clean %>% count(myeloid_state, raw_cluster, name = "n_cells") %>% group_by(myeloid_state) %>% 
    mutate(fraction_within_state = n_cells/sum(n_cells)) %>% ungroup()

dataset_summary <- md_clean %>% count(myeloid_state, source_dataset, name = "n_cells") %>% group_by(myeloid_state) %>% 
    mutate(fraction_within_state = n_cells/sum(n_cells)) %>% ungroup()

has_timepoint <- ("timepoint_main" %in% colnames(md_clean) && any(!is.na(md_clean$timepoint_main)))

if (has_timepoint) {
    time_summary <- md_clean %>% filter(!is.na(timepoint_main)) %>% count(myeloid_state, timepoint_main, 
        name = "n_cells") %>% group_by(myeloid_state) %>% mutate(fraction_within_state = n_cells/sum(n_cells)) %>% 
        ungroup()
} else {
    time_summary <- data.frame()
}

dominant_dataset <- dataset_summary %>% group_by(myeloid_state) %>% slice_max(fraction_within_state, 
    n = 1, with_ties = FALSE) %>% ungroup() %>% transmute(myeloid_state, dominant_dataset = as.character(source_dataset), 
    dominant_dataset_fraction = fraction_within_state)

if (has_timepoint) {
    dominant_time <- time_summary %>% group_by(myeloid_state) %>% slice_max(fraction_within_state, n = 1, 
        with_ties = FALSE) %>% ungroup() %>% transmute(myeloid_state, dominant_timepoint = as.character(timepoint_main), 
        dominant_timepoint_fraction = fraction_within_state)
} else {
    dominant_time <- data.frame(myeloid_state = state_count$myeloid_state, dominant_timepoint = NA_character_, 
        dominant_timepoint_fraction = NA_real_, stringsAsFactors = FALSE)
}

annotation_summary <- state_count %>% left_join(dominant_dataset, by = "myeloid_state") %>% left_join(dominant_time, 
    by = "myeloid_state") %>% mutate(myeloid_state = factor(myeloid_state, levels = state_levels)) %>% 
    arrange(myeloid_state)

write.csv(annotation_summary, file.path(out_dir, "18_myeloid_annotation_summary_v01.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

support_long <- bind_rows(raw_mapping_count %>% transmute(myeloid_state, variable = "raw_cluster", category = raw_cluster, 
    n_cells, fraction = fraction_within_state), dataset_summary %>% transmute(myeloid_state, variable = "source_dataset", 
    category = as.character(source_dataset), n_cells, fraction = fraction_within_state), if (has_timepoint) {
    time_summary %>% transmute(myeloid_state, variable = "timepoint_main", category = as.character(timepoint_main), 
        n_cells, fraction = fraction_within_state)
} else {
    NULL
})

write.csv(support_long, file.path(out_dir, "18_myeloid_annotation_support_v01.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

obj_clean <- tryCatch(JoinLayers(obj_clean, assay = "RNA"), error = function(e) {
    log_msg("JoinLayers skipped/failed safely: ", conditionMessage(e))
    obj_clean
})

DefaultAssay(obj_clean) <- "RNA"

rna_layers <- tryCatch(Layers(obj_clean[["RNA"]]), error = function(e) character())

if (!"data" %in% rna_layers) {
    obj_clean <- NormalizeData(obj_clean, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, 
        verbose = TRUE)
}

Idents(obj_clean) <- "myeloid_state"

state_markers <- FindAllMarkers(object = obj_clean, assay = "RNA", only.pos = TRUE, test.use = "wilcox", 
    min.pct = MIN_PCT, logfc.threshold = LOGFC_THRESHOLD, return.thresh = PADJ_THRESHOLD, densify = FALSE, 
    verbose = TRUE)

if (nrow(state_markers) == 0) {
    stop("FindAllMarkers returned no significant state-level markers.")
}

if (!"gene" %in% colnames(state_markers)) {
  state_markers$gene <- rownames(state_markers)
}

state_markers$cluster <- as.character(state_markers$cluster)

fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")

fc_col <- fc_candidates[fc_candidates %in% colnames(state_markers)][1]

if (is.na(fc_col)) {
    stop("Could not identify fold-change column in state-marker output.")
}

state_markers <- state_markers %>% mutate(state_order = match(cluster, state_levels)) %>% arrange(state_order, 
    p_val_adj, desc(.data[[fc_col]])) %>% select(-state_order)

write.csv(state_markers, file.path(out_dir, "18_myeloid_state_markers_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

marker_sets <- list(`Homeostatic Microglia` = c("P2ry12", "Tmem119", "Sall1", "Hexb"), `Lipid-associated Microglia` = c("Clec7a", 
    "Lpl", "Apoc4", "Apoe", "Trem2"), `Activated inflammatory Microglia` = c("Il1b", "Nfkbia", "Tnfaip3", 
    "Cxcl10"), `Ccl12+ acute Microglia` = c("Ccl12", "Il1rn", "Lgals3"), `Tnf+ inflammatory Microglia` = c("Tnf", 
    "Slc7a5", "F3"), `Proliferating Microglia` = c("Mki67", "Pclaf", "Top2a", "Birc5"), `MHC-II+ antigen-presenting Macrophage` = c("H2-Aa", 
    "H2-Ab1", "H2-Eb1", "Cd74", "Mrc1"), `Gpnmb+ tissue-repair Macrophage` = c("Gpnmb", "Atp6v0d2", "Fabp4", 
    "Arg1"), `Acute-response Macrophage` = c("S100a8", "S100a9", "Nos2", "Vcan", "Cxcl3"))

marker_sets_present <- lapply(marker_sets, function(x) {
    x[x %in% rownames(obj_clean)]
})

marker_sets_present <- marker_sets_present[lengths(marker_sets_present) > 0]

marker_order <- unique(unlist(marker_sets_present, use.names = FALSE))

marker_group_lengths <- vapply(marker_sets_present, length, FUN.VALUE = integer(1))

marker_group_ends <- cumsum(marker_group_lengths)

separator_x <- marker_group_ends[-length(marker_group_ends)] + 0.5

set.seed(SEED)

plot_cells <- Cells(obj_clean)

if (length(plot_cells) > MAX_PLOT_CELLS) {
    plot_cells <- sample(plot_cells, MAX_PLOT_CELLS)
}

obj_plot <- subset(obj_clean, cells = plot_cells)

Idents(obj_plot) <- "myeloid_state"

state_colors <- c(`Homeostatic Microglia` = "#3B6FB6", `Clec7a+ lipid-associated Microglia` = "#E08B2C", 
    `Activated inflammatory Microglia` = "#C94C4C", `Ccl12+ acute Microglia` = "#D66BA0", `Tnf+ inflammatory Microglia` = "#8E63B6", 
    `Proliferating Microglia` = "#2F9E8F", `MHC-II+ antigen-presenting Macrophage` = "#4F8A3B", `Gpnmb+ tissue-repair Macrophage` = "#8A6D5A", 
    `Acute-response Macrophage` = "#D9A514")

diagnostic_pdf <- file.path(out_dir, "18_myeloid_annotation_diagnostics_v01.pdf")

pdf(diagnostic_pdf, width = 12, height = 8.5, onefile = TRUE)

p_raw <- DimPlot(obj_plot, reduction = UMAP_REDUCTION, group.by = RAW_CLUSTER_COL, label = TRUE, repel = TRUE, 
    pt.size = 0.12, raster = TRUE) + ggtitle("Resolution-0.10 raw clusters after contaminant removal") + 
    theme_classic(base_size = 10) + theme(legend.position = "none")

p_state <- DimPlot(obj_plot, reduction = UMAP_REDUCTION, group.by = "myeloid_state", label = FALSE, pt.size = 0.12, 
    raster = TRUE, cols = state_colors) + ggtitle("Consolidated myeloid states") + theme_classic(base_size = 10)

print(p_raw + p_state + plot_layout(widths = c(1, 1.18)))

mapping_plot_df <- md_clean %>% count(raw_cluster, myeloid_state, name = "n_cells") %>% mutate(raw_cluster = factor(raw_cluster, 
    levels = as.character(sort(unique(as.integer(raw_cluster))))), myeloid_state = factor(myeloid_state, 
    levels = state_levels))

p_mapping <- ggplot(mapping_plot_df, aes(x = raw_cluster, y = n_cells, fill = myeloid_state)) + geom_col(width = 0.8) + 
    scale_fill_manual(values = state_colors, drop = FALSE) + labs(title = "Raw-cluster to consolidated-state mapping", 
    x = "Resolution-0.10 raw cluster", y = "Cell count", fill = "Myeloid state") + theme_classic(base_size = 10)

print(p_mapping)

p_dot <- DotPlot(obj_clean, features = marker_order, assay = "RNA", group.by = "myeloid_state", dot.scale = 5.5, 
    scale = TRUE) + RotatedAxis() + labs(title = "Marker validation of consolidated myeloid states", 
    x = NULL, y = NULL, color = "Average\nexpression", size = "Percent\nexpressed") + theme_classic(base_size = 8.5) + 
    theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7.5))

for (xx in separator_x) {
    p_dot <- p_dot + geom_vline(xintercept = xx, linetype = 3, linewidth = 0.3, color = "grey55")
}

print(p_dot)

dataset_plot_df <- dataset_summary %>% mutate(myeloid_state = factor(myeloid_state, levels = state_levels))

p_dataset <- ggplot(dataset_plot_df, aes(x = myeloid_state, y = fraction_within_state, fill = source_dataset)) + 
    geom_col(width = 0.82) + scale_y_continuous(limits = c(0, 1), labels = function(x) {
    paste0(round(x * 100), "%")
}, expand = c(0, 0)) + labs(title = "Source-dataset support of consolidated myeloid states", x = NULL, 
    y = "Within-state fraction", fill = "Dataset") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35, 
    hjust = 1, vjust = 1))

if (has_timepoint) {
    time_plot_df <- time_summary %>% mutate(myeloid_state = factor(myeloid_state, levels = state_levels))
    p_time <- ggplot(time_plot_df, aes(x = myeloid_state, y = fraction_within_state, fill = timepoint_main)) + 
        geom_col(width = 0.82) + scale_y_continuous(limits = c(0, 1), labels = function(x) {
        paste0(round(x * 100), "%")
    }, expand = c(0, 0)) + labs(title = "Time-point distribution of consolidated myeloid states", x = NULL, 
        y = "Within-state fraction", fill = "Time point") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35, 
        hjust = 1, vjust = 1))
    print(p_dataset/p_time + plot_layout(heights = c(1, 1)))
} else {
    print(p_dataset)
}

feature_candidates <- c("P2ry12", "Clec7a", "Lpl", "Ccl12", "Tnf", "Mki67", "Cd74", "Gpnmb", "S100a9")

feature_present <- feature_candidates[feature_candidates %in% rownames(obj_plot)]

if (length(feature_present) > 0) {
    p_feature <- FeaturePlot(obj_plot, features = feature_present, reduction = UMAP_REDUCTION, ncol = 3, 
        order = TRUE, pt.size = 0.12, raster = TRUE) & theme_classic(base_size = 8)
    print(p_feature)
}

dev.off()

rm(obj_plot, plot_cells, remove_cells)

gc(verbose = FALSE)

out_rds <- file.path(out_dir, "18_myeloid_annotation_cleanup_3scRNA_v01.rds")

saveRDS(obj_clean, out_rds, compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)

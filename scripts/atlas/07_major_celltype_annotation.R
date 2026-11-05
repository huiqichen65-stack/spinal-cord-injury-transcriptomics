source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
CLUSTER_COL <- "harmony30_res.0.4"
MAX_FEATUREPLOT_CELLS <- 50000

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "patchwork", "data.table")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("缺少R包：", paste(missing_pkgs, collapse = ", "), "\n请先安装后再运行。")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

in_file <- file.path(project_dir, "results", "06_Harmony_cluster_resolution", "06_Harmony_cluster_resolution_3scRNA_v01.rds")

out_dir <- file.path(project_dir, "results", "07_major_celltype_annotation")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "07_major_celltype_annotation_3scRNA_v02.log")

if (file.exists(log_file)) file.remove(log_file)

log_msg <- write_log

if (!file.exists(in_file)) {
  stop("未找到输入文件：\n", in_file)
}

obj <- readRDS(in_file)
DefaultAssay(obj) <- "RNA"

if (!CLUSTER_COL %in% colnames(obj[[]])) {
  stop("对象中没有cluster字段：", CLUSTER_COL)
}

raw_cluster <- as.character(obj[[CLUSTER_COL]][, 1])

cluster_levels_num <- sort(unique(as.numeric(raw_cluster)))

safe_levels <- paste0("c", cluster_levels_num)

obj$cluster07 <- factor(paste0("c", raw_cluster), levels = safe_levels)

PLOT_CLUSTER_COL <- "cluster07"
Idents(obj) <- PLOT_CLUSTER_COL

log_msg("Clusters = ", length(levels(Idents(obj))))

cluster_tb <- sort(table(Idents(obj)), decreasing = TRUE)

cluster_size_df <- data.frame(cluster = names(cluster_tb), n_cells = as.integer(cluster_tb), percent = as.integer(cluster_tb)/ncol(obj) *
  100, stringsAsFactors = FALSE)

write.csv(cluster_size_df, file.path(out_dir, "07_major_celltype_annotation_3scRNA_cluster_size_v02.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

marker_sets <- list(Microglia = c("P2ry12", "Tmem119", "Sall1", "Hexb", "Cx3cr1", "C1qa", "C1qb"), Macrophage_Monocyte = c("Lyz2",
  "Ccr2", "Ly6c2", "Ms4a7", "Fcgr3", "Cd68", "Adgre1"), Neutrophil = c("S100a8", "S100a9", "Ly6g",
  "Retnlg", "Camp"), T_NK = c("Cd3d", "Cd3e", "Trbc1", "Cd4", "Cd8a", "Nkg7", "Klrd1"), B_cell = c("Cd79a",
  "Ms4a1", "Cd74", "Cd37", "H2-Aa"), Astrocyte = c("Aldoc", "Aqp4", "Slc1a3", "Gfap", "S100b", "Sparcl1"),
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp"), OPC = c("Pdgfra", "Cspg4", "Olig1",
    "Olig2", "Sox10"), Neuron = c("Snap25", "Syt1", "Rbfox3", "Tubb3", "Gria2"), Endothelial = c("Pecam1",
    "Cldn5", "Kdr", "Emcn", "Esam", "Vwf"), Pericyte_VSMC = c("Pdgfrb", "Rgs5", "Cspg4", "Des", "Acta2",
    "Tagln"), Fibroblast_Meningeal = c("Col1a1", "Col1a2", "Dcn", "Lum", "Col3a1", "Pdgfra"), Ependymal = c("Foxj1",
    "Pifo", "Tppp3", "Cfap126"))

all_markers <- unique(unlist(marker_sets, use.names = FALSE))
present_markers <- all_markers[all_markers %in% rownames(obj)]
missing_markers <- setdiff(all_markers, present_markers)

log_msg("Markers missing = ", length(missing_markers))

if (length(missing_markers) > 0) {
  log_msg("Missing markers: ", paste(missing_markers, collapse = ", "))
}

marker_sets_present <- lapply(marker_sets, function(x) x[x %in% rownames(obj)])

marker_sets_present <- marker_sets_present[lengths(marker_sets_present) > 0]

flat_markers <- unlist(marker_sets_present, use.names = FALSE)
dup_markers <- unique(flat_markers[duplicated(flat_markers)])

if (length(dup_markers) > 0) {
}

seen_markers <- character(0)

marker_sets_overview <- lapply(marker_sets_present, function(x) {
  keep <- x[!x %in% seen_markers]
  seen_markers <<- c(seen_markers, keep)
  keep
})

marker_sets_overview <- marker_sets_overview[lengths(marker_sets_overview) > 0]

avg_list <- AggregateExpression(obj, assays = "RNA", features = present_markers, group.by = PLOT_CLUSTER_COL,
  slot = "data", return.seurat = FALSE, verbose = FALSE)

avg_mat <- avg_list$RNA

avg_df <- as.data.frame(avg_mat)
avg_df$gene <- rownames(avg_df)

avg_df <- avg_df[, c("gene", setdiff(colnames(avg_df), "gene"))]

write.csv(avg_df, file.path(out_dir, "07_major_celltype_annotation_3scRNA_marker_matrix_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

plot_cells <- Cells(obj)

if (length(plot_cells) > MAX_FEATUREPLOT_CELLS) {
  set.seed(SEED)
  plot_cells <- sample(plot_cells, MAX_FEATUREPLOT_CELLS)
}

obj_plot <- subset(obj, cells = plot_cells)

Idents(obj_plot) <- PLOT_CLUSTER_COL

pdf_file <- file.path(out_dir, "07_major_celltype_annotation_3scRNA_diagnostics_v02.pdf")

pdf(pdf_file, width = 14, height = 9, onefile = TRUE)

p_cluster <- DimPlot(obj_plot, reduction = "umap_harmony30", group.by = PLOT_CLUSTER_COL, label = TRUE,
  repel = TRUE, pt.size = 0.12, raster = TRUE) + ggtitle("Main atlas | resolution = 0.4") + theme_classic(base_size = 11) +
  theme(legend.position = "none")

print(p_cluster)

p_dot_all <- DotPlot(obj, features = marker_sets_overview, group.by = PLOT_CLUSTER_COL, dot.scale = 5) +
  RotatedAxis() + labs(title = "Canonical marker screening across 34 clusters", x = NULL, y = "Cluster") +
  theme_classic(base_size = 9)

print(p_dot_all)

neural_markers <- unique(c(marker_sets_present$Astrocyte, marker_sets_present$Oligodendrocyte, marker_sets_present$OPC,
  marker_sets_present$Neuron, marker_sets_present$Ependymal))

neural_markers <- neural_markers[neural_markers %in% rownames(obj)]

if (length(neural_markers) > 0) {
    p_dot_neural <- DotPlot(obj, features = neural_markers, group.by = PLOT_CLUSTER_COL, dot.scale = 6) + 
        RotatedAxis() + labs(title = "Neural / glial canonical markers", x = NULL, y = "Cluster") + theme_classic(base_size = 10)
    print(p_dot_neural)
}

immune_markers <- unique(c(marker_sets_present$Microglia, marker_sets_present$Macrophage_Monocyte, marker_sets_present$Neutrophil,
  marker_sets_present$T_NK, marker_sets_present$B_cell))

immune_markers <- immune_markers[immune_markers %in% rownames(obj)]

if (length(immune_markers) > 0) {
    p_dot_immune <- DotPlot(obj, features = immune_markers, group.by = PLOT_CLUSTER_COL, dot.scale = 6) + 
        RotatedAxis() + labs(title = "Immune canonical markers", x = NULL, y = "Cluster") + theme_classic(base_size = 10)
    print(p_dot_immune)
}

vascular_markers <- unique(c(marker_sets_present$Endothelial, marker_sets_present$Pericyte_VSMC, marker_sets_present$Fibroblast_Meningeal))

vascular_markers <- vascular_markers[vascular_markers %in% rownames(obj)]

if (length(vascular_markers) > 0) {
    p_dot_vascular <- DotPlot(obj, features = vascular_markers, group.by = PLOT_CLUSTER_COL, dot.scale = 6) + 
        RotatedAxis() + labs(title = "Vascular / stromal canonical markers", x = NULL, y = "Cluster") + 
        theme_classic(base_size = 10)
    print(p_dot_vascular)
}

feature_groups <- list(Major_lineages_1 = c("P2ry12", "Lyz2", "S100a8", "Cd3d", "Aldoc", "Mbp", "Pdgfra",
  "Snap25"), Major_lineages_2 = c("Pecam1", "Pdgfrb", "Col1a1", "Foxj1", "Tmem119", "Ccr2", "Nkg7",
  "Cd79a"), Myeloid_focus = c("P2ry12", "Tmem119", "Cx3cr1", "C1qa", "Lyz2", "Ccr2", "Ly6c2", "Ms4a7"))

for (nm in names(feature_groups)) {
    feats <- feature_groups[[nm]]
    feats <- feats[feats %in% rownames(obj_plot)]
    if (length(feats) == 0) 
        next
    p <- FeaturePlot(obj_plot, features = feats, reduction = "umap_harmony30", ncol = 4, pt.size = 0.1, 
        raster = TRUE, order = TRUE) + plot_annotation(title = nm)
    print(p)
}

dev.off()

log_msg("Clusters = ", length(levels(Idents(obj))))
capture.output(sessionInfo(), file = log_file, append = TRUE)

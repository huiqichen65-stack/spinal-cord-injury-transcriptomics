source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
CLUSTER_COL <- "harmony30_res.0.4"
TARGET_CLUSTERS <- c("c8", "c23", "c27", "c31", "c32", "c33")
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

out_dir <- file.path(project_dir, "results", "08_ambiguous_cluster_validation")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "08_ambiguous_cluster_validation_3scRNA_v01.log")

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

obj$cluster08 <- factor(paste0("c", raw_cluster), levels = safe_levels)

Idents(obj) <- "cluster08"

missing_targets <- setdiff(TARGET_CLUSTERS, levels(Idents(obj)))

if (length(missing_targets) > 0) {
  stop("以下目标cluster不存在：", paste(missing_targets, collapse = ", "))
}

log_msg("Total clusters = ", length(levels(Idents(obj))))

target_counts <- table(Idents(obj))[TARGET_CLUSTERS]
capture.output(print(target_counts), file = log_file, append = TRUE)

markers_c8 <- c("Itgax", "Flt3", "Clec10a", "Xcr1", "Ccr7", "H2-Ab1", "H2-Aa", "H2-Eb1", "Cd74", "Cd79a",
  "Cd79b", "Ms4a1", "Cd19", "Cd37", "Csf1r", "Lyz2", "Fcgr3", "Ms4a7", "Cd68")

markers_c23 <- c("Cd3d", "Cd3e", "Trbc1", "Trbc2", "Cd4", "Cd8a", "Nkg7", "Klrd1", "Klrk1", "Prf1", "Gzmb",
  "Cd79a", "Cd79b", "Ms4a1", "Cd19", "Jchain", "Mzb1", "Sdc1")

markers_c27 <- c("Pdgfra", "Cspg4", "Olig1", "Olig2", "Sox10", "Bcan", "Vcan", "Snap25", "Syt1", "Rbfox3",
  "Tubb3", "Gria2", "Slc17a7", "Gad1", "Gad2")

markers_c31 <- c("S100a8", "S100a9", "Ly6g", "Retnlg", "Camp", "Ngp", "Mmp8", "Mmp9", "Lyz2", "Ccr2",
  "Ly6c2", "Fcgr3", "Ms4a7", "Cd68", "Adgre1", "Csf1r")

markers_c32c33 <- c("P2ry12", "Tmem119", "Sall1", "Hexb", "Cx3cr1", "C1qa", "C1qb", "Lyz2", "Ccr2", "Ly6c2",
  "Ms4a7", "Adgre1", "Cd68", "Mrc1", "Lyve1", "Pf4", "Apoe", "Lpl", "Mki67", "Top2a", "Pclaf", "Ccnb2",
  "Cdk1", "Stmn1", "Isg15", "Ifit1", "Ifit3", "Irf7", "Stat1")

panel_list <- list(c8_APC_DC_B_Myeloid = markers_c8, c23_T_NK_B_Plasma = markers_c23, c27_OPC_Neuron = markers_c27,
  c31_Neutrophil_Monocyte = markers_c31, c32_c33_Myeloid_State = markers_c32c33)

panel_list_present <- lapply(panel_list, function(x) unique(x[x %in% rownames(obj)]))

all_requested <- unique(unlist(panel_list, use.names = FALSE))
all_present <- unique(unlist(panel_list_present, use.names = FALSE))
all_missing <- setdiff(all_requested, all_present)

log_msg("Missing markers = ", length(all_missing))

if (length(all_missing) > 0) {
  log_msg("Missing markers: ", paste(all_missing, collapse = ", "))
}

obj_target <- subset(obj, idents = TARGET_CLUSTERS)

obj_target$cluster08 <- factor(as.character(obj_target$cluster08), levels = TARGET_CLUSTERS)

Idents(obj_target) <- "cluster08"

all_present_unique <- unique(all_present)

dot_all <- DotPlot(obj_target, features = all_present_unique, group.by = "cluster08")

marker_data <- dot_all$data

gene_to_panel <- do.call(rbind, lapply(names(panel_list_present), function(pn) {
  data.frame(features.plot = panel_list_present[[pn]], marker_panel = pn, stringsAsFactors = FALSE)
}))

gene_to_panel <- gene_to_panel[!duplicated(gene_to_panel$features.plot), , drop = FALSE]

marker_data <- merge(marker_data, gene_to_panel, by = "features.plot", all.x = TRUE, sort = FALSE)

marker_data <- marker_data[order(match(marker_data$id, TARGET_CLUSTERS), match(marker_data$features.plot,
  all_present_unique)), ]

write.csv(marker_data, file.path(out_dir, "08_ambiguous_cluster_validation_3scRNA_marker_data_v01.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

set.seed(SEED)

plot_cells <- Cells(obj)

if (length(plot_cells) > MAX_FEATUREPLOT_CELLS) {
  plot_cells <- sample(plot_cells, MAX_FEATUREPLOT_CELLS)
}

obj_plot <- subset(obj, cells = plot_cells)

obj_plot$validation_group <- ifelse(as.character(obj_plot$cluster08) %in% TARGET_CLUSTERS, as.character(obj_plot$cluster08),
  "Other")

obj_plot$validation_group <- factor(obj_plot$validation_group, levels = c("Other", TARGET_CLUSTERS))

pdf_file <- file.path(out_dir, "08_ambiguous_cluster_validation_3scRNA_diagnostics_v01.pdf")

pdf(pdf_file, width = 14, height = 9, onefile = TRUE)

p_loc <- DimPlot(obj_plot, reduction = "umap_harmony30", group.by = "validation_group", pt.size = 0.12,
  raster = TRUE) + ggtitle("Ambiguous clusters highlighted in the global atlas") + theme_classic(base_size = 11)

print(p_loc)

p_all <- DotPlot(obj_target, features = all_present_unique, group.by = "cluster08", dot.scale = 6) +
  RotatedAxis() + labs(title = "Targeted marker screening of ambiguous clusters", x = NULL, y = "Cluster") +
  theme_classic(base_size = 9)

print(p_all)

for (pn in names(panel_list_present)) {

  feats <- panel_list_present[[pn]]

  if (length(feats) == 0)
    next

  clusters_this <- switch(pn, c8_APC_DC_B_Myeloid = "c8", c23_T_NK_B_Plasma = "c23", c27_OPC_Neuron = "c27",
    c31_Neutrophil_Monocyte = "c31", c32_c33_Myeloid_State = c("c32", "c33"), TARGET_CLUSTERS)

  tmp <- subset(obj_target, idents = clusters_this)

  p <- DotPlot(tmp, features = feats, group.by = "cluster08", dot.scale = 7) + RotatedAxis() + labs(title = pn,
    x = NULL, y = "Cluster") + theme_classic(base_size = 10)

  print(p)
}

key_features <- c("H2-Ab1", "Cd74", "Itgax", "Cd3d", "Nkg7", "Cd79a", "Pdgfra", "Snap25", "Ly6g", "Ccr2",
  "P2ry12", "Mki67")

key_features <- key_features[key_features %in% rownames(obj_plot)]

if (length(key_features) > 0) {
    p_feat <- FeaturePlot(obj_plot, features = key_features, reduction = "umap_harmony30", ncol = 4, 
        pt.size = 0.1, raster = TRUE, order = TRUE) + plot_annotation(title = "Representative validation markers")
    print(p_feat)
}

dev.off()

capture.output(sessionInfo(), file = log_file, append = TRUE)

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggrastr)
})

FONT_FAMILY <- "sans"
USE_SHOWTEXT <- FALSE

arial_regular <- "C:/Windows/Fonts/arial.ttf"
arial_bold    <- "C:/Windows/Fonts/arialbd.ttf"

if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE) && file.exists(arial_regular)) {
    if (file.exists(arial_bold)) {
        sysfonts::font_add(family = "ArialLocal", regular = arial_regular, bold = arial_bold)
    }
    else {
        sysfonts::font_add(family = "ArialLocal", regular = arial_regular)
    }
    showtext::showtext_auto(enable = TRUE)
    FONT_FAMILY <- "ArialLocal"
    USE_SHOWTEXT <- TRUE
} else {
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

in_file <- file.path(project_dir, "results", "06_Harmony_cluster_resolution", "06_Harmony_cluster_resolution_3scRNA_v01.rds")

out_dir <- file.path(project_dir, "results", "09_Figure1A_major_celltype_UMAP")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(in_file)

cluster_col   <- "harmony30_res.0.4"
reduction_use <- "umap_harmony30"

if (!cluster_col %in% colnames(obj[[]])) {
  stop("缺少cluster字段: ", cluster_col)
}

if (!reduction_use %in% names(obj@reductions)) {
  stop("缺少 reduction: ", reduction_use)
}

cluster_raw  <- as.character(obj[[cluster_col]][, 1])
cluster_safe <- paste0("c", cluster_raw)

major_map <- c(c0 = "Microglia", c1 = "Macrophage/Monocyte", c2 = "Endothelial", c3 = "Macrophage/Monocyte", 
    c4 = "Ependymal", c5 = "Microglia", c6 = "Microglia", c7 = "Microglia", c8 = "APC/DC-like myeloid", 
    c9 = "Microglia", c10 = "OPC", c11 = "Endothelial", c12 = "Endothelial", c13 = "Neutrophil", c14 = "Pericyte/VSMC", 
    c15 = "Microglia", c16 = "Astrocyte", c17 = "Neuron", c18 = "Endothelial", c19 = "OPC", c20 = "Ependymal", 
    c21 = "Fibroblast/Meningeal", c22 = "Ependymal", c23 = "T/NK cell", c24 = "Oligodendrocyte", c25 = "Ependymal", 
    c26 = "Microglia", c27 = "OPC", c28 = "Microglia", c29 = "Oligodendrocyte", c30 = "Macrophage/Monocyte", 
    c31 = "Neutrophil", c32 = "Macrophage/Monocyte", c33 = "Microglia")

missing_cluster <- setdiff(unique(cluster_safe), names(major_map))

if (length(missing_cluster) > 0) {
    stop("以下cluster没有major cell type映射: ", paste(missing_cluster, collapse = ", "))
}

obj$major_celltype <- unname(major_map[cluster_safe])

major_levels <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid", "Neutrophil", "T/NK cell", 
    "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

obj$major_celltype <- factor(obj$major_celltype, levels = major_levels)

celltype_cols <- c(Microglia = "#3B6FB6", `Macrophage/Monocyte` = "#D94F4F", `APC/DC-like myeloid` = "#8E63B6", 
    Neutrophil = "#D9A514", `T/NK cell` = "#3E9B50", Astrocyte = "#44A6A6", Oligodendrocyte = "#8A6D5A", 
    OPC = "#7E8A97", Neuron = "#E47D8E", Ependymal = "#C8A600", Endothelial = "#169EB3", `Pericyte/VSMC` = "#8B4A43", 
    `Fibroblast/Meningeal` = "#626262")

umap_df <- as.data.frame(Embeddings(obj, reduction_use)[, 1:2, drop = FALSE])

colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$major_celltype <- obj$major_celltype

label_df <- umap_df %>% group_by(major_celltype) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), 
    .groups = "drop")

label_map <- c(Microglia = "Microglia", `Macrophage/Monocyte` = "Macrophage/Monocyte", `APC/DC-like myeloid` = "APC/DC-like", 
    Neutrophil = "Neutrophil", `T/NK cell` = "T/NK cell", Astrocyte = "Astrocyte", Oligodendrocyte = "Oligodendrocyte", 
    OPC = "OPC", Neuron = "Neuron", Ependymal = "Ependymal", Endothelial = "Endothelial", `Pericyte/VSMC` = "Pericyte/VSMC", 
    `Fibroblast/Meningeal` = "Fibroblast")

label_df$label <- unname(label_map[as.character(label_df$major_celltype)])

label_df <- label_df %>% mutate(UMAP_1 = case_when(major_celltype == "Macrophage/Monocyte" ~ UMAP_1 + 
    1.35, major_celltype == "APC/DC-like myeloid" ~ UMAP_1 + 0.5, major_celltype == "Fibroblast/Meningeal" ~ 
    UMAP_1 + 0.2, major_celltype == "Pericyte/VSMC" ~ UMAP_1 + 0.65, major_celltype == "Endothelial" ~ 
    UMAP_1 - 0.45, TRUE ~ UMAP_1), UMAP_2 = case_when(major_celltype == "Macrophage/Monocyte" ~ UMAP_2 + 
    0.15, major_celltype == "APC/DC-like myeloid" ~ UMAP_2 + 0.2, major_celltype == "Fibroblast/Meningeal" ~ 
    UMAP_2 + 0.85, major_celltype == "Pericyte/VSMC" ~ UMAP_2 + 0.55, major_celltype == "Endothelial" ~ 
    UMAP_2 - 0.55, TRUE ~ UMAP_2))

p <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = major_celltype
  )
) +
  ggrastr::geom_point_rast(
    size = 0.30,
    alpha = 0.90,
    raster.dpi = 600
  ) +
  scale_color_manual(
    values = celltype_cols,
    drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = label
    ),
    inherit.aes = FALSE,
    family = FONT_FAMILY,
    fontface = "bold",
    size = 4.8,
    color = "black",
    box.padding = 0.15,
    point.padding = 0.05,
    min.segment.length = Inf,
    max.overlaps = Inf,
    seed = 1234
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme_classic(
    base_family = FONT_FAMILY,
    base_size = 12
  ) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(
      family = FONT_FAMILY,
      size = 14,
      face = "plain",
      color = "black"
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    legend.position = "right",
    legend.text = element_text(
      family = FONT_FAMILY,
      size = 11.5,
      color = "black"
    ),
    legend.key.height = grid::unit(0.36, "cm"),
    legend.key.width = grid::unit(0.36, "cm"),
    legend.spacing.y = grid::unit(0.02, "cm"),
    panel.border = element_blank(),
    plot.margin = margin(5, 7, 5, 5)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 2.8, alpha = 1),
      ncol = 1,
      byrow = TRUE
    )
  )

ggsave(filename = file.path(out_dir, "09_Figure1A_major_celltype_UMAP_v05.pdf"), plot = p, width = 9, 
    height = 7, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = file.path(out_dir, "09_Figure1A_major_celltype_UMAP_v05.svg"), width = 9, 
        height = 7)
    print(p)
    dev.off()
} else {
    ggsave(filename = file.path(out_dir, "09_Figure1A_major_celltype_UMAP_v05.svg"), plot = p, width = 9, 
        height = 7, units = "in")
}

if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = file.path(out_dir, "09_Figure1A_major_celltype_UMAP_v05.png"), width = 9, 
        height = 7, units = "in", res = 600, background = "white")
    print(p)
    dev.off()
} else {
    ggsave(filename = file.path(out_dir, "09_Figure1A_major_celltype_UMAP_v05.png"), plot = p, width = 9, 
        height = 7, units = "in", dpi = 600, bg = "white")
}

frozen_object_file <- file.path(out_dir, "09_Figure1A_major_celltype_UMAP_object_v04.rds")
if (!file.exists(frozen_object_file)) {
  saveRDS(obj, frozen_object_file, compress = FALSE)
} else {
}

print(table(obj$major_celltype))

if (USE_SHOWTEXT) {
  showtext::showtext_auto(enable = FALSE)
}

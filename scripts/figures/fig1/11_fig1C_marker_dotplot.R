rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- "D:/zz"

IN_FILE <- file.path(ROOT, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

OUT_DIR <- file.path(ROOT, "results", "12_Figure1C_major_celltype_marker_DotPlot_v03")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "12_Figure1C_major_celltype_marker_DotPlot_v03"

FIG_W <- 13.5
FIG_H <- 6.8
DPI <- 600

FONT_GENE <- 17.0
FONT_CELLTYPE <- 16.0
FONT_LEGEND_TITLE <- 15.0
FONT_LEGEND_TEXT <- 15.0
GENE_ANGLE <- 75

COLOR_MIN <- -1
COLOR_MAX <- 1

required_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "svglite", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(svglite)
  library(scales)
})

if (!file.exists(IN_FILE)) {
 stop("Figure 1A object not found:\n", IN_FILE)
}

obj <- readRDS(IN_FILE)

DefaultAssay(obj) <- "RNA"

obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) {
    message("JoinLayers skipped: ", conditionMessage(e))
    obj
})

if (!"major_celltype" %in% colnames(obj[[]])) {
  stop("Metadata column 'major_celltype' was not found.")
}

rna_layers <- SeuratObject::Layers(obj[["RNA"]])

if (!"data" %in% rna_layers) {
    stop("RNA normalized 'data' layer was not found.\n", "This plotting script will NOT rerun normalization.")
}

major_levels <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid", "Neutrophil", "T/NK cell", 
    "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

marker_sets <- list(Microglia = c("P2ry12", "Tmem119", "Sall1"), `Macrophage/Monocyte` = c("Ccr2", "Ly6c2", 
    "Ms4a7"), `APC/DC-like myeloid` = c("Cd74", "H2-Aa", "Itgax"), Neutrophil = c("S100a8", "Ly6g", "Retnlg"), 
    `T/NK cell` = c("Cd3d", "Nkg7", "Klrd1"), Astrocyte = c("Aqp4", "Slc1a3", "Gfap"), Oligodendrocyte = c("Mbp", 
        "Mog", "Plp1"), OPC = c("Pdgfra", "Olig1", "Olig2"), Neuron = c("Snap25", "Rbfox3", "Tubb3"), 
    Ependymal = c("Foxj1", "Pifo", "Tppp3"), Endothelial = c("Pecam1", "Cldn5", "Kdr"), `Pericyte/VSMC` = c("Rgs5", 
        "Acta2", "Tagln"), `Fibroblast/Meningeal` = c("Col1a1", "Dcn", "Lum"))

obj$major_celltype <- factor(as.character(obj$major_celltype), levels = major_levels)

all_requested <- unlist(marker_sets, use.names = FALSE)

missing_genes <- setdiff(all_requested, rownames(obj))

if (length(missing_genes) > 0L) {
    stop("The following canonical marker gene(s) are missing from the object:\n", paste(missing_genes, 
        collapse = ", "), "\nThe script stops rather than silently changing the marker panel.")
}

gene_order <- all_requested

dp0 <- Seurat::DotPlot(object = obj, features = gene_order, group.by = "major_celltype", assay = "RNA", 
    scale = TRUE, scale.by = "radius", dot.scale = 6)

dot_df <- dp0$data

required_cols <- c("features.plot", "id", "pct.exp", "avg.exp", "avg.exp.scaled")

missing_cols <- setdiff(required_cols, colnames(dot_df))

if (length(missing_cols) > 0L) {
    stop("Unexpected Seurat DotPlot output; missing column(s): ", paste(missing_cols, collapse = ", "))
}

dot_df$gene <- factor(dot_df$features.plot, levels = gene_order)

dot_df$major_celltype_display <- factor(as.character(dot_df$id), levels = rev(major_levels))

dot_df$avg_exp_plot <- pmax(COLOR_MIN, pmin(COLOR_MAX, dot_df$avg.exp.scaled))

group_ends <- cumsum(lengths(marker_sets))

separator_x <- group_ends[-length(group_ends)] + 0.5

p <- ggplot(
  dot_df,
  aes(
    x = gene,
    y = major_celltype_display
  )
) +

  geom_vline(
    xintercept = separator_x,
    color = "#EFEFEF",
    linewidth = 0.30
  ) +

  geom_point(
    aes(
      size = pct.exp,
      color = avg_exp_plot
    ),
    alpha = 1
  ) +

  scale_size_continuous(
    name = "Percent Expressed",
    range = c(0.45, 6.2),
    limits = c(0, 100),
    breaks = c(25, 50, 75)
  ) +

  scale_color_gradient2(
    name = "Average Expression",
    low = "#3B6FB6",
    mid = "#FFF3B0",
    high = "#D73027",
    midpoint = 0,
    limits = c(COLOR_MIN, COLOR_MAX),
    breaks = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    oob = scales::squish
  ) +

  scale_x_discrete(
    expand = expansion(add = 0.45)
  ) +

  scale_y_discrete(
    expand = expansion(add = 0.35)
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  guides(
    size = guide_legend(
      order = 1,
      title.position = "left",
      title.hjust = 0,
      direction = "horizontal",
      override.aes = list(
        color = "black"
      )
    ),

    color = guide_colorbar(
      order = 2,
      title.position = "left",
      title.hjust = 0,
      direction = "horizontal",
      barwidth = grid::unit(1.80, "cm"),
      barheight = grid::unit(0.30, "cm"),
      ticks = TRUE
    )
  ) +

  theme_bw(
    base_family = "Arial",
    base_size = 15
  ) +

  theme(
    text = element_text(
      family = "Arial",
      color = "black"
    ),

    axis.text.x = element_text(
      family = "Arial",
      face = "italic",
      size = FONT_GENE,
      angle = GENE_ANGLE,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),

    axis.text.y = element_text(
      family = "Arial",
      size = FONT_CELLTYPE,
      color = "black"
    ),

    axis.ticks = element_blank(),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      color = "#666666",
      fill = NA,
      linewidth = 0.55
    ),

    legend.position = "top",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.spacing.x = grid::unit(0.15, "cm"),

    legend.title = element_text(
      family = "Arial",
      size = FONT_LEGEND_TITLE,
      color = "black"
    ),

    legend.text = element_text(
      family = "Arial",
      size = FONT_LEGEND_TEXT,
      color = "black"
    ),

    plot.margin = margin(
      t = 4,
      r = 7,
      b = 4,
      l = 4,
      unit = "pt"
    )
  )

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))

PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

DATA_FILE <- file.path(OUT_DIR, paste0(PREFIX, "_data.csv"))

MARKER_MAP_FILE <- file.path(OUT_DIR, paste0(PREFIX, "_marker_map.csv"))

ggsave(filename = SVG_FILE, plot = p, device = svglite::svglite, width = FIG_W, height = FIG_H, units = "in", 
    bg = "white")

ggsave(filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = FIG_W, height = FIG_H, units = "in", 
    family = "Arial", bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = FIG_W, height = FIG_H, units = "in", dpi = DPI, bg = "white")

write.csv(dot_df[, c("major_celltype_display", "features.plot", "pct.exp", "avg.exp", "avg.exp.scaled", 
    "avg_exp_plot")], DATA_FILE, row.names = FALSE, fileEncoding = "UTF-8")

marker_map <- do.call(rbind, lapply(names(marker_sets), function(ct) {
    data.frame(major_celltype = ct, gene = marker_sets[[ct]], stringsAsFactors = FALSE)
}))

write.csv(marker_map, MARKER_MAP_FILE, row.names = FALSE, fileEncoding = "UTF-8")

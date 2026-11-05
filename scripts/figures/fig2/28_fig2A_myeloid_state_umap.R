options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
UMAP_REDUCTION <- "umap_myeloid_harmony30"

POINT_SIZE_RASTER <- 0.30
POINT_ALPHA_RASTER <- 1.00
RASTER_DPI <- 600

FIG_W <- 9.1
FIG_H <- 6.6
DPI <- 600

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "dplyr", "svglite", "ggrastr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nPlease install them before running Figure 2A.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(ggrastr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^20_AUCell_lipid_metabolism_scoring_3scRNA_v03\\.rds$", 
        recursive = TRUE, full.names = TRUE)
    if (length(hits) == 0) {
        stop("Cannot find Step 20 v03 RDS:\n", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "28_Figure2A_myeloid_state_UMAP")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(in_file)

if (!"myeloid_state" %in% colnames(obj[[]])) {
  stop("Metadata column 'myeloid_state' was not found.")
}

if (!UMAP_REDUCTION %in% names(obj@reductions)) {
    stop("Reduction '", UMAP_REDUCTION, "' was not found.\nAvailable reductions: ", paste(names(obj@reductions), 
        collapse = ", "))
}

state_levels <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

display_labels <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb⁺ macrophages", 
    `Acute-response Macrophage` = "Acute-response Macrophage")

label_short <- c(`Homeostatic Microglia` = "Homeostatic", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling", 
    `Activated inflammatory Microglia` = "Activated inflammatory", `Ccl12+ acute Microglia` = "Ccl12+ acute", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory", `Proliferating Microglia` = "Proliferating", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb⁺ macrophages", 
    `Acute-response Macrophage` = "Acute-response")

state_colors <- c(`Homeostatic Microglia` = "#4C78A8", `Clec7a+ lipid-associated Microglia` = "#D64F4F", 
    `Activated inflammatory Microglia` = "#F2B134", `Ccl12+ acute Microglia` = "#8E6BBE", `Tnf+ inflammatory Microglia` = "#E77AA8", 
    `Proliferating Microglia` = "#4FA6A2", `MHC-II+ antigen-presenting Macrophage` = "#59A14F", `Gpnmb+ tissue-repair Macrophage` = "#D5A832", 
    `Acute-response Macrophage` = "#9C755F")

observed_states <- unique(as.character(obj$myeloid_state))

missing_states <- setdiff(state_levels, observed_states)

extra_states <- setdiff(observed_states, state_levels)

if (length(missing_states) > 0) {
    warning("Expected state(s) not found: ", paste(missing_states, collapse = ", "))
}

if (length(extra_states) > 0) {
    warning("Unexpected state(s) found: ", paste(extra_states, collapse = ", "))
}

umap <- Embeddings(obj, reduction = UMAP_REDUCTION)

plot_df <- data.frame(cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], myeloid_state = as.character(obj$myeloid_state[match(rownames(umap), 
    colnames(obj))]), stringsAsFactors = FALSE)

plot_df <- plot_df[!is.na(plot_df$myeloid_state), , drop = FALSE]

plot_df$myeloid_state <- factor(plot_df$myeloid_state, levels = state_levels)

state_n <- plot_df %>% dplyr::count(myeloid_state, name = "n")

plot_df <- plot_df %>% dplyr::left_join(state_n, by = "myeloid_state") %>% dplyr::arrange(dplyr::desc(n))

label_df <- plot_df %>% dplyr::group_by(myeloid_state) %>% dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), 
    UMAP_2 = median(UMAP_2, na.rm = TRUE), .groups = "drop")

xr <- range(plot_df$UMAP_1, finite = TRUE)

yr <- range(plot_df$UMAP_2, finite = TRUE)

xspan <- diff(xr)
yspan <- diff(yr)

label_offsets <- data.frame(
  myeloid_state = state_levels,
  dx = c(
     0.000,
     0.000,
     0.000,
    -0.010,
    -0.005,
     0.005,
     0.000,
     0.000,
     0.008
  ),
  dy = c(
     0.006,
     0.000,
     0.003,
    -0.003,
    -0.006,
     0.006,
     0.003,
     0.006,
    -0.005
  ),
  stringsAsFactors = FALSE
)

label_df <- label_df %>% dplyr::left_join(label_offsets, by = "myeloid_state") %>% dplyr::mutate(UMAP_1 = UMAP_1 + 
    dx * xspan, UMAP_2 = UMAP_2 + dy * yspan, label = unname(label_short[as.character(myeloid_state)]))

axis_x0 <- xr[1] + 0.115 * xspan
axis_y0 <- yr[1] + 0.145 * yspan

axis_x1 <- axis_x0 + 0.105 * xspan
axis_y1 <- axis_y0 + 0.105 * yspan

p <- ggplot(
  plot_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = myeloid_state
  )
) +
  ggrastr::geom_point_rast(
    size = POINT_SIZE_RASTER,
    alpha = POINT_ALPHA_RASTER,
    stroke = 0,
    raster.dpi = RASTER_DPI
  ) +

  geom_text(
    data = label_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = label
    ),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "plain",
    size = 4.0,
    color = "black",
    lineheight = 0.92,
    check_overlap = FALSE
  ) +

  scale_color_manual(
    values = state_colors,
    breaks = state_levels,
    labels = unname(
      display_labels[
        state_levels
      ]
    ),
    drop = FALSE
  ) +

  annotate(
    "segment",
    x = axis_x0,
    xend = axis_x1,
    y = axis_y0,
    yend = axis_y0,
    linewidth = 0.40,
    color = "black",
    arrow = grid::arrow(
      length = grid::unit(
        0.070,
        "inches"
      ),
      type = "closed"
    )
  ) +
  annotate(
    "segment",
    x = axis_x0,
    xend = axis_x0,
    y = axis_y0,
    yend = axis_y1,
    linewidth = 0.40,
    color = "black",
    arrow = grid::arrow(
      length = grid::unit(
        0.070,
        "inches"
      ),
      type = "closed"
    )
  ) +
  annotate(
    "text",
    x = axis_x1,
    y = axis_y0 - 0.032 * yspan,
    label = "UMAP1",
    family = "Arial",
    size = 3.65,
    hjust = 1
  ) +
  annotate(
    "text",
    x = axis_x0 - 0.024 * xspan,
    y = axis_y1,
    label = "UMAP2",
    family = "Arial",
    size = 3.65,
    angle = 90,
    hjust = 1
  ) +

  coord_equal(
    clip = "off"
  ) +

  guides(
    color = guide_legend(
      override.aes = list(
        size = 3.2,
        alpha = 1
      ),
      ncol = 1,
      byrow = TRUE
    )
  ) +

  labs(
    x = NULL,
    y = NULL,
    color = NULL
  ) +

  theme_void(
    base_family = "Arial"
  ) +

  theme(
    text = element_text(
      family = "Arial",
      color = "black"
    ),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(
      family = "Arial",
      size = 14.2,
      color = "black",
      lineheight = 0.95
    ),
    legend.key.height = grid::unit(
      0.43,
      "cm"
    ),
    legend.key.width = grid::unit(
      0.35,
      "cm"
    ),
    legend.spacing.y = grid::unit(
      0.03,
      "cm"
    ),
    plot.margin = margin(
      4,
      8,
      4,
      4,
      unit = "pt"
    )
  )

svg_file <- file.path(out_dir, "28_Figure2A_myeloid_state_UMAP_v09.svg")

pdf_file <- file.path(out_dir, "28_Figure2A_myeloid_state_UMAP_v09.pdf")

png_file <- file.path(out_dir, "28_Figure2A_myeloid_state_UMAP_v09.png")

tiff_file <- file.path(out_dir, "28_Figure2A_myeloid_state_UMAP_v09.tif")

ggsave(filename = svg_file, plot = p, device = svglite::svglite, width = FIG_W, height = FIG_H, units = "in", 
    bg = "white")

ggsave(filename = pdf_file, plot = p, device = grDevices::cairo_pdf, width = FIG_W, height = FIG_H, units = "in", 
    bg = "white")

ggsave(filename = png_file, plot = p, width = FIG_W, height = FIG_H, units = "in", dpi = DPI, bg = "white")

ggsave(filename = tiff_file, plot = p, device = "tiff", width = FIG_W, height = FIG_H, units = "in", 
    dpi = DPI, compression = "lzw", bg = "white")

count_table <- plot_df %>% dplyr::count(myeloid_state, name = "n_cells") %>% dplyr::mutate(display_label = unname(display_labels[as.character(myeloid_state)]), 
    pct = 100 * n_cells/sum(n_cells))

write.csv(count_table, file.path(out_dir, "28_Figure2A_myeloid_state_counts_v09.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

print(count_table)

source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(ggrastr)
})

FONT_FAMILY <- "sans"
USE_SHOWTEXT <- FALSE

arial_regular <- FONT_REGULAR
arial_bold <- FONT_BOLD

if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE) && file.exists(arial_regular)) {
  if (file.exists(arial_bold)) {
    sysfonts::font_add(family = "ArialLocal", regular = arial_regular, bold = arial_bold)
  } else {
    sysfonts::font_add(family = "ArialLocal", regular = arial_regular)
  }

  showtext::showtext_auto(enable = TRUE)
  FONT_FAMILY <- "ArialLocal"
  USE_SHOWTEXT <- TRUE
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

in_file <- file.path(project_dir, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

out_dir <- file.path(project_dir, "results", "13_Figure1E_major_celltype_timepoint_UMAP")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop("未找到Figure 1A对象：\n", in_file)
}

obj <- readRDS(in_file)

reduction_use <- "umap_harmony30"

if (!reduction_use %in% names(obj@reductions)) {
  stop("对象中缺少 reduction: ", reduction_use)
}

if (!"major_celltype" %in% colnames(obj[[]])) {
  stop("对象中缺少 major_celltype 字段。")
}

major_levels <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid", "Neutrophil", "T/NK cell",
  "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

obj$major_celltype <- factor(as.character(obj$major_celltype), levels = major_levels)

celltype_cols <- c(Microglia = "#3B6FB6", `Macrophage/Monocyte` = "#D94F4F", `APC/DC-like myeloid` = "#8E63B6",
  Neutrophil = "#D9A514", `T/NK cell` = "#3E9B50", Astrocyte = "#44A6A6", Oligodendrocyte = "#8A6D5A",
  OPC = "#7E8A97", Neuron = "#E47D8E", Ependymal = "#C8A600", Endothelial = "#169EB3", `Pericyte/VSMC` = "#8B4A43",
  `Fibroblast/Meningeal` = "#626262")

md <- obj[[]]

candidate_cols <- c("timepoint_main", "timepoint", "time", "Time", "injury_time", "day", "dpi", "group",
  "Group", "condition", "Condition", "sample", "sample_id", "sampleID", "library_id", "orig.ident")

candidate_cols <- candidate_cols[candidate_cols %in% colnames(md)]

if (length(candidate_cols) == 0) {
  stop("metadata中未找到可用于识别时间点的字段。")
}

standardize_timepoint <- function(x) {

  s <- tolower(trimws(as.character(x)))

  if (is.na(s) || s == "")
    return(NA_character_)

  if (str_detect(s, "(^|[^0-9])38\\s*(d|dpi|day|days)\\b|\\bd38\\b"))
    return("38 dpi")

  if (str_detect(s, "(^|[^0-9])28\\s*(d|dpi|day|days)\\b|\\bd28\\b"))
    return("28 dpi")

  if (str_detect(s, "(^|[^0-9])14\\s*(d|dpi|day|days)\\b|\\bd14\\b"))
    return("14 dpi")

  if (str_detect(s, "(^|[^0-9])0?7\\s*(d|dpi|day|days)\\b|\\bd0?7\\b"))
    return("7 dpi")

  if (str_detect(s, "(^|[^0-9])0?3\\s*(d|dpi|day|days)\\b|\\bd0?3\\b"))
    return("3 dpi")

  if (str_detect(s, "(^|[^0-9])0?1\\s*(d|dpi|day|days)\\b|\\bd0?1\\b"))
    return("1 dpi")

  if (str_detect(s, "(^|[^0-9])4\\s*(h|hr|hrs|hour|hours)\\b"))
    return("4 h")

  if (str_detect(s, paste0("\\buninjured\\b|\\buninj\\b|\\bnaive\\b|", "\\bsham\\b|\\bcontrol\\b|\\bctrl\\b|",
    "(^|[^0-9])0{1,2}\\s*(d|dpi|day|days)\\b")))
    return("Control")

  NA_character_
}

md$timepoint_main <- NA_character_

for (cc in candidate_cols) {
    vals <- vapply(md[[cc]], standardize_timepoint, FUN.VALUE = character(1))
    idx <- is.na(md$timepoint_main) & !is.na(vals)
    md$timepoint_main[idx] <- vals[idx]
}

time_text <- apply(md[, candidate_cols, drop = FALSE], 1, function(x) {
  paste(unique(as.character(x)), collapse = " | ")
})

fallback_vals <- vapply(time_text, standardize_timepoint, FUN.VALUE = character(1))

idx_fb <- is.na(md$timepoint_main) & !is.na(fallback_vals)

md$timepoint_main[idx_fb] <- fallback_vals[idx_fb]

if (any(is.na(md$timepoint_main))) {
    bad_n <- sum(is.na(md$timepoint_main))
    bad_show <- unique(time_text[is.na(md$timepoint_main)])
    bad_show <- head(bad_show, 20)
    stop("仍有 ", bad_n, " 个细胞无法识别timepoint。\n", "metadata示例：\n", paste(bad_show, 
        collapse = "\n"))
}

obj$timepoint_main <- md$timepoint_main

key_levels <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

keep_cells <- Cells(obj)[obj$timepoint_main %in% key_levels]

obj_key <- subset(obj, cells = keep_cells)

obj_key$timepoint_main <- factor(obj_key$timepoint_main, levels = key_levels)

umap_df <- as.data.frame(Embeddings(obj_key, reduction_use)[, 1:2, drop = FALSE])

colnames(umap_df) <- c("UMAP_1", "UMAP_2")

umap_df$cell <- rownames(umap_df)

umap_df$timepoint_main <- obj_key$timepoint_main

umap_df$major_celltype <- factor(as.character(obj_key$major_celltype), levels = major_levels)

full_umap <- Embeddings(obj, reduction_use)

x_range <- range(full_umap[, 1], finite = TRUE)

y_range <- range(full_umap[, 2], finite = TRUE)

x_pad <- diff(x_range) * 0.025
y_pad <- diff(y_range) * 0.025

x_lim <- c(x_range[1] - x_pad, x_range[2] + x_pad)

y_lim <- c(y_range[1] - y_pad, y_range[2] + y_pad)

write.csv(umap_df, file.path(out_dir, "13_Figure1E_major_celltype_timepoint_UMAP_data_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

count_df <- umap_df %>% count(timepoint_main, major_celltype, name = "n") %>% group_by(timepoint_main) %>% 
    mutate(total_cells = sum(n), proportion = n/total_cells) %>% ungroup()

write.csv(count_df, file.path(out_dir, "13_Figure1E_major_celltype_timepoint_cellcount_v02.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = major_celltype)) + ggrastr::geom_point_rast(size = 0.25,
  alpha = 0.88, raster.dpi = 600) + scale_color_manual(values = celltype_cols, drop = FALSE) + facet_wrap(~timepoint_main,
  ncol = 3, drop = FALSE) + coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE, clip = "on") +
  labs(x = "UMAP 1", y = "UMAP 2", color = NULL) + theme_classic(base_family = FONT_FAMILY, base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(family = FONT_FAMILY, size = 14,
    face = "bold", color = "black", margin = margin(b = 4)), axis.title = element_text(family = FONT_FAMILY,
    size = 13, color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_line(color = "black",
    linewidth = 0.45), panel.spacing = grid::unit(0.5, "cm"), legend.position = "right", legend.text = element_text(family = FONT_FAMILY,
    size = 12, color = "black"), legend.key.height = grid::unit(0.36, "cm"), legend.key.width = grid::unit(0.36,
    "cm"), legend.spacing.y = grid::unit(0.015, "cm"), plot.margin = margin(5, 7, 5, 5)) + guides(color = guide_legend(override.aes = list(size = 2.6,
  alpha = 1), ncol = 1, byrow = TRUE))

ggsave(filename = file.path(out_dir, "13_Figure1E_major_celltype_timepoint_UMAP_v02.pdf"), plot = p,
  width = 11.2, height = 7.1, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = file.path(out_dir, "13_Figure1E_major_celltype_timepoint_UMAP_v02.svg"), 
        width = 11.2, height = 7.1)
    print(p)
    dev.off()
} else {
    ggsave(filename = file.path(out_dir, "13_Figure1E_major_celltype_timepoint_UMAP_v02.svg"), plot = p, 
        width = 11.2, height = 7.1, units = "in")
}

if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = file.path(out_dir, "13_Figure1E_major_celltype_timepoint_UMAP_v02.png"), 
        width = 11.2, height = 7.1, units = "in", res = 600, background = "white")
    print(p)
    dev.off()
} else {
    ggsave(filename = file.path(out_dir, "13_Figure1E_major_celltype_timepoint_UMAP_v02.png"), plot = p, 
        width = 11.2, height = 7.1, units = "in", dpi = 600, bg = "white")
}

print(umap_df %>% count(timepoint_main, name = "n_cells"))

if (USE_SHOWTEXT) {
  showtext::showtext_auto(enable = FALSE)
}

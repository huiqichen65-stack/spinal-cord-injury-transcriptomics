options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

FONT_FAMILY <- "sans"
USE_SHOWTEXT <- FALSE

arial_regular <- "C:/Windows/Fonts/arial.ttf"
arial_bold <- "C:/Windows/Fonts/arialbd.ttf"

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

in_file1 <- file.path(project_dir, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

in_file2 <- file.path(project_dir, "results", "06_Harmony_cluster_resolution", "06_Harmony_cluster_resolution_3scRNA_v01.rds")

out_dir <- file.path(project_dir, "results", "10_Figure1B_major_celltype_composition")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(in_file1)) {
  obj <- readRDS(in_file1)
} else if (file.exists(in_file2)) {
  obj <- readRDS(in_file2)
} else {
  stop("未找到输入对象。")
}

cluster_col <- "harmony30_res.0.4"

if (!"major_celltype" %in% colnames(obj[[]])) {
    if (!cluster_col %in% colnames(obj[[]])) {
        stop("对象中既没有 major_celltype，也没有 ", cluster_col)
    }
    cluster_raw <- as.character(obj[[cluster_col]][, 1])
    cluster_safe <- paste0("c", cluster_raw)
    major_map <- c(c0 = "Microglia", c1 = "Macrophage/Monocyte", c2 = "Endothelial", c3 = "Macrophage/Monocyte", 
        c4 = "Ependymal", c5 = "Microglia", c6 = "Microglia", c7 = "Microglia", c8 = "APC/DC-like myeloid", 
        c9 = "Microglia", c10 = "OPC", c11 = "Endothelial", c12 = "Endothelial", c13 = "Neutrophil", 
        c14 = "Pericyte/VSMC", c15 = "Microglia", c16 = "Astrocyte", c17 = "Neuron", c18 = "Endothelial", 
        c19 = "OPC", c20 = "Ependymal", c21 = "Fibroblast/Meningeal", c22 = "Ependymal", c23 = "T/NK cell", 
        c24 = "Oligodendrocyte", c25 = "Ependymal", c26 = "Microglia", c27 = "OPC", c28 = "Microglia", 
        c29 = "Oligodendrocyte", c30 = "Macrophage/Monocyte", c31 = "Neutrophil", c32 = "Macrophage/Monocyte", 
        c33 = "Microglia")
    obj$major_celltype <- unname(major_map[cluster_safe])
}

major_levels <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid", "Neutrophil", "T/NK cell",
  "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

obj$major_celltype <- factor(obj$major_celltype, levels = major_levels)

celltype_cols <- c(Microglia = "#3B6FB6", `Macrophage/Monocyte` = "#D94F4F", `APC/DC-like myeloid` = "#8E63B6",
  Neutrophil = "#D9A514", `T/NK cell` = "#3E9B50", Astrocyte = "#44A6A6", Oligodendrocyte = "#8A6D5A",
  OPC = "#7E8A97", Neuron = "#E47D8E", Ependymal = "#C8A600", Endothelial = "#169EB3", `Pericyte/VSMC` = "#8B4A43",
  `Fibroblast/Meningeal` = "#626262")

md <- obj[[]]

candidate_cols <- c("timepoint_main", "timepoint", "time", "Time", "injury_time", "day", "dpi", "group",
  "Group", "condition", "Condition", "sample", "sample_id", "sampleID", "library_id", "orig.ident")

candidate_cols <- candidate_cols[candidate_cols %in% colnames(md)]

if (length(candidate_cols) == 0) {
  stop("未找到可用于判定时间点的metadata字段。")
}

time_text <- apply(md[, candidate_cols, drop = FALSE], 1, function(x) {
  paste(unique(as.character(x)), collapse = " | ")
})

standardize_timepoint <- function(x) {
  s <- tolower(trimws(as.character(x)))

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

  if (str_detect(s, "\\buninjured\\b|\\buninj\\b|\\bnaive\\b|\\bsham\\b|\\bcontrol\\b|\\bctrl\\b|(^|[^0-9])0{1,2}\\s*(d|dpi|day|days)\\b"))
    return("Control")

  return(NA_character_)
}

md$timepoint_main <- NA_character_

for (cc in candidate_cols) {
    vals <- vapply(md[[cc]], standardize_timepoint, FUN.VALUE = character(1))
    idx <- is.na(md$timepoint_main) & !is.na(vals)
    md$timepoint_main[idx] <- vals[idx]
}

fallback_vals <- vapply(time_text, standardize_timepoint, FUN.VALUE = character(1))

idx_fb <- is.na(md$timepoint_main) & !is.na(fallback_vals)
md$timepoint_main[idx_fb] <- fallback_vals[idx_fb]

if ("source_dataset" %in% colnames(md)) {
    is_196928 <- md$source_dataset == "GSE196928"
    source_txt <- tolower(time_text)
    md$timepoint_main[is_196928 & is.na(md$timepoint_main) & str_detect(source_txt, "\\bb1\\b")] <- "Control"
    md$timepoint_main[is_196928 & is.na(md$timepoint_main) & str_detect(source_txt, "\\bb3\\b")] <- "7 dpi"
    md$timepoint_main[is_196928 & is.na(md$timepoint_main) & str_detect(source_txt, "\\bb7\\b")] <- "28 dpi"
}

if (any(is.na(md$timepoint_main))) {
  bad_n <- sum(is.na(md$timepoint_main))
  bad_show <- unique(time_text[is.na(md$timepoint_main)])
  bad_show <- head(bad_show, 20)
  stop("仍有 ", bad_n, " 个细胞无法识别timepoint。\n", "请检查以下metadata示例：\n",
    paste(bad_show, collapse = "\n"))
}

obj$timepoint_main <- md$timepoint_main

full_levels <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")
key_levels <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

comp_full <- obj[[]] %>% transmute(timepoint_main = factor(timepoint_main, levels = full_levels), major_celltype = factor(major_celltype, 
    levels = major_levels)) %>% filter(!is.na(timepoint_main), !is.na(major_celltype)) %>% count(timepoint_main, 
    major_celltype, .drop = FALSE) %>% group_by(timepoint_main) %>% mutate(total_cells = sum(n), proportion = ifelse(total_cells > 
    0, n/total_cells, 0), percent = proportion * 100) %>% ungroup()

write.csv(comp_full, file.path(out_dir, "10_Figure1B_major_celltype_composition_full_timepoints_v04.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

comp_key <- comp_full %>% filter(as.character(timepoint_main) %in% key_levels) %>% mutate(timepoint_main = factor(timepoint_main, 
    levels = key_levels))

write.csv(comp_key, file.path(out_dir, "10_Figure1B_major_celltype_composition_data_v04.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

p <- ggplot(comp_key, aes(x = timepoint_main, y = proportion, fill = major_celltype)) + geom_col(color = "white",
  linewidth = 0.25, width = 0.8) + scale_fill_manual(values = celltype_cols, drop = FALSE) + scale_y_continuous(labels = function(x) paste0(round(x *
  100), "%"), expand = c(0, 0)) + labs(x = NULL, y = "Cell proportion", fill = NULL) + theme_classic(base_family = FONT_FAMILY,
  base_size = 12) + theme(axis.title.y = element_text(size = 14, color = "black"), axis.text.x = element_text(size = 13,
  color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.line = element_line(color = "black",
  linewidth = 0.6), axis.ticks = element_line(color = "black", linewidth = 0.5), legend.position = "right",
  legend.text = element_text(size = 12.5, color = "black"), legend.key.height = grid::unit(0.37, "cm"),
  legend.key.width = grid::unit(0.37, "cm"), legend.spacing.y = grid::unit(0.02, "cm"), panel.grid = element_blank(),
  plot.margin = margin(5, 7, 5, 5)) + guides(fill = guide_legend(reverse = TRUE, ncol = 1, byrow = TRUE))

ggsave(filename = file.path(out_dir, "10_Figure1B_major_celltype_composition_v04.pdf"), plot = p, width = 8.6,
  height = 5.8, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file = file.path(out_dir, "10_Figure1B_major_celltype_composition_v04.svg"), width = 8.6,
    height = 5.8)
  print(p)
  dev.off()
} else {
  ggsave(filename = file.path(out_dir, "10_Figure1B_major_celltype_composition_v04.svg"), plot = p,
    width = 8.6, height = 5.8, units = "in")
}

if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(filename = file.path(out_dir, "10_Figure1B_major_celltype_composition_v04.png"), width = 8.6,
    height = 5.8, units = "in", res = 600, background = "white")
  print(p)
  dev.off()
} else {
  ggsave(filename = file.path(out_dir, "10_Figure1B_major_celltype_composition_v04.png"), plot = p,
    width = 8.6, height = 5.8, units = "in", dpi = 600, bg = "white")
}

print(comp_key %>% distinct(timepoint_main, total_cells))

if (USE_SHOWTEXT) {
  showtext::showtext_auto(enable = FALSE)
}

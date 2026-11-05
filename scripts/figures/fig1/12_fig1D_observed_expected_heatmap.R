source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

FONT_FAMILY <- "sans"
USE_SHOWTEXT <- FALSE

arial_regular <- FONT_REGULAR
arial_bold    <- FONT_BOLD

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
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

in_candidates <- file.path(project_dir, "results", "10_Figure1B_major_celltype_composition", c("10_Figure1B_major_celltype_composition_data_v04.csv", 
    "10_Figure1B_major_celltype_composition_data_v03.csv"))
in_existing <- in_candidates[file.exists(in_candidates)]
if (length(in_existing) == 0L) {
  stop("Figure 1B数据不存在。请先运行Figure 1B v03或v04。\n", paste(in_candidates, collapse = "\n"))
}
in_file <- in_existing[1L]

out_dir <- file.path(project_dir, "results", "11_Figure1D_major_celltype_Roe_heatmap_v02")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
    stop("未找到Figure 1B数据文件：\n", in_file, "\n请先确认10_Figure1B v03已成功运行。")
}

dat <- read.csv(in_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("timepoint_main", "major_celltype", "n", "total_cells")

missing_cols <- setdiff(required_cols, colnames(dat))

if (length(missing_cols) > 0) {
    stop("输入CSV缺少字段：", paste(missing_cols, collapse = ", "))
}

key_levels <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

major_levels <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid", "Neutrophil", "T/NK cell", 
    "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

dat <- dat %>% filter(timepoint_main %in% key_levels, major_celltype %in% major_levels) %>% mutate(timepoint_main = factor(timepoint_main, 
    levels = key_levels), major_celltype = factor(major_celltype, levels = major_levels))

tab <- dat %>% select(timepoint_main, major_celltype, n) %>% tidyr::complete(timepoint_main = factor(key_levels, 
    levels = key_levels), major_celltype = factor(major_levels, levels = major_levels), fill = list(n = 0))

grand_total <- sum(tab$n)

if (grand_total <= 0) {
  stop("总细胞数为0，无法计算R/oe。")
}

time_totals <- tab %>% group_by(timepoint_main) %>% summarise(time_total = sum(n), .groups = "drop")

celltype_totals <- tab %>% group_by(major_celltype) %>% summarise(celltype_total = sum(n), .groups = "drop")

roe_df <- tab %>% left_join(time_totals, by = "timepoint_main") %>% left_join(celltype_totals, by = "major_celltype") %>% 
    mutate(expected = time_total * celltype_total/grand_total, Roe = ifelse(expected > 0, n/expected, 
        NA_real_), log2_Roe = ifelse(Roe > 0, log2(Roe), NA_real_), log2_Roe_plot = pmax(pmin(log2_Roe, 
        2), -2), Roe_label = ifelse(is.na(Roe), "", sprintf("%.2f", Roe)))

write.csv(roe_df, file.path(out_dir, "11_Figure1D_major_celltype_Roe_data_v02.csv"), row.names = FALSE, 
    fileEncoding = "UTF-8")

roe_df$major_celltype_plot <- factor(as.character(roe_df$major_celltype), levels = rev(major_levels))

p <- ggplot(
  roe_df,
  aes(
    x = timepoint_main,
    y = major_celltype_plot,
    fill = log2_Roe_plot
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.7,
    width = 0.98,
    height = 0.98
  ) +
  geom_text(
    aes(label = Roe_label),
    family = FONT_FAMILY,
    size = 3.9,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#D94F4F",
    midpoint = 0,
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    name = expression(log[2]~"(R/oe)"),
    na.value = "#D9D9D9"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_fixed(ratio = 0.75) +
  theme_minimal(
    base_family = FONT_FAMILY,
    base_size = 11
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 12.5,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 11.5,
      color = "black"
    ),
    legend.text = element_text(
      size = 10.5,
      color = "black"
    ),
    plot.margin = margin(5, 7, 5, 5)
  )

ggsave(filename = file.path(out_dir, "11_Figure1D_major_celltype_Roe_heatmap_v02.pdf"), plot = p, width = 8.4, 
    height = 6.4, units = "in", device = cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = file.path(out_dir, "11_Figure1D_major_celltype_Roe_heatmap_v02.svg"), width = 8.4, 
        height = 6.4)
    print(p)
    dev.off()
} else {
    ggsave(filename = file.path(out_dir, "11_Figure1D_major_celltype_Roe_heatmap_v02.svg"), plot = p, 
        width = 8.4, height = 6.4, units = "in")
}

if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = file.path(out_dir, "11_Figure1D_major_celltype_Roe_heatmap_v02.png"), width = 8.4, 
        height = 6.4, units = "in", res = 600, background = "white")
    print(p)
    dev.off()
} else {
    ggsave(filename = file.path(out_dir, "11_Figure1D_major_celltype_Roe_heatmap_v02.png"), plot = p, 
        width = 8.4, height = 6.4, units = "in", dpi = 600, bg = "white")
}

if (USE_SHOWTEXT) {
  showtext::showtext_auto(enable = FALSE)
}

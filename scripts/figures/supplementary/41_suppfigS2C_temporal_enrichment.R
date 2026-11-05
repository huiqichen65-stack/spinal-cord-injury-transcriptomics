options(stringsAsFactors = FALSE)

required_pkgs <- c("ggplot2", "dplyr", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

ROOT <- "D:/zz"

IN_DIR <- file.path(ROOT, "results", "19_myeloid_state_robustness_temporal_dynamics")

IN_FILE <- file.path(IN_DIR, "19_myeloid_Roe_main_v02.csv")

OUT_DIR <- file.path(ROOT, "results", "41_SuppFigure2C_Roe_temporal_heatmap")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "41_SuppFigure2C_Roe_temporal_heatmap_v02")

if (!file.exists(IN_FILE)) {
    stop("Cannot find fixed Step19 R/oe file:\n", IN_FILE)
}

TIME_LEVELS <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

DISPLAY_LABELS <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ tissue-repair Macrophage", 
    `Acute-response Macrophage` = "Acute-response Macrophage")

roe <- read.csv(IN_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("timepoint_main", "myeloid_state", "observed", "expected", "Roe", "log2_Roe")

missing_cols <- setdiff(required_cols, colnames(roe))

if (length(missing_cols) > 0) {
    stop("R/oe file is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

roe_plot <- roe %>%
  filter(
    timepoint_main %in% TIME_LEVELS,
    myeloid_state %in% STATE_LEVELS
  ) %>%
  mutate(
    timepoint_main = factor(
      timepoint_main,
      levels = TIME_LEVELS
    ),

    myeloid_state = factor(
      myeloid_state,
      levels = rev(STATE_LEVELS)
    ),

    plot_log2_Roe = ifelse(
      is.finite(log2_Roe),
      pmax(
        pmin(log2_Roe, 4),
        -4
      ),
      NA_real_
    )
  )

if (nrow(roe_plot) == 0) {
 stop("No valid rows remained for the main timeline.")
}

write.csv(roe_plot, file.path(OUT_DIR, "41_SuppFigure2C_Roe_temporal_heatmap_plotdata_v02.csv"), row.names = FALSE)

p <- ggplot(roe_plot, aes(x = timepoint_main, y = myeloid_state, fill = plot_log2_Roe)) + geom_tile(color = "white", 
    linewidth = 0.7, width = 0.98, height = 0.98) + scale_fill_gradient2(low = "#2B5DA8", mid = "#F7F7F7", 
    high = "#C43C39", midpoint = 0, limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4), na.value = "#BDBDBD", 
    name = expression(log[2](R/oe))) + scale_y_discrete(labels = function(x) {
    unname(DISPLAY_LABELS[as.character(x)])
}, drop = FALSE) + scale_x_discrete(drop = FALSE) + coord_fixed(ratio = 0.72) + labs(x = NULL, y = NULL) + 
    theme_minimal(base_family = "Arial", base_size = 9) + theme(text = element_text(family = "Arial", 
    color = "black"), panel.grid = element_blank(), axis.text.x = element_text(family = "Arial", size = 8.2, 
    color = "black", angle = 0, hjust = 0.5), axis.text.y = element_text(family = "Arial", size = 7.6, 
    color = "black"), axis.ticks = element_blank(), legend.position = "right", legend.title = element_text(family = "Arial", 
    size = 8.3), legend.text = element_text(family = "Arial", size = 7.5), legend.key.height = grid::unit(6, 
    "mm"), plot.margin = margin(6, 6, 6, 6))

W_IN <- 7.3
H_IN <- 5.4

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

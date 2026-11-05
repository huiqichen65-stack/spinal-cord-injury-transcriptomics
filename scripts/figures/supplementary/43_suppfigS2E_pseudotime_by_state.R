options(stringsAsFactors = FALSE)

required_pkgs <- c("monocle3", "SingleCellExperiment", "ggplot2", "dplyr", "ggridges", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(monocle3)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(ggridges)
})

ROOT <- "D:/zz"

STEP25_DIR <- file.path(ROOT, "results", "25_Monocle3_all_myeloid_trajectory")

CDS_FILE <- file.path(STEP25_DIR, "25_Monocle3_all_myeloid_trajectory_3scRNA_v03.rds")

OUT_DIR <- file.path(ROOT, "results", "43_SuppFigure2E_pseudotime_ridge_by_state")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "43_SuppFigure2E_pseudotime_ridge_by_state_v01")

if (!file.exists(CDS_FILE)) {
    stop("Cannot find Step25 CDS:\n", CDS_FILE)
}

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

DISPLAY_LABELS <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ tissue-repair Macrophage", 
    `Acute-response Macrophage` = "Acute-response Macrophage")

STATE_COLORS <- c(`Homeostatic Microglia` = "#4C78A8", `Clec7a+ lipid-associated Microglia` = "#D64F4F", 
    `Activated inflammatory Microglia` = "#F2B134", `Ccl12+ acute Microglia` = "#8E6BBE", `Tnf+ inflammatory Microglia` = "#E77AA8", 
    `Proliferating Microglia` = "#4FA6A2", `MHC-II+ antigen-presenting Macrophage` = "#59A14F", `Gpnmb+ tissue-repair Macrophage` = "#D5A832", 
    `Acute-response Macrophage` = "#9C755F")

cds <- readRDS(CDS_FILE)

if (!inherits(cds, "cell_data_set")) {
    stop("Input RDS is not a Monocle3 cell_data_set.")
}

meta <- as.data.frame(colData(cds))

required_cols <- c("myeloid_state")
missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
    stop("CDS is missing required colData column(s): ", paste(missing_cols, collapse = ", "))
}

pt <- monocle3::pseudotime(cds)

if (is.null(names(pt))) {
  names(pt) <- colnames(cds)
}

meta$cell <- rownames(meta)

plot_df <- data.frame(cell = colnames(cds), pseudotime = as.numeric(pt[colnames(cds)]), stringsAsFactors = FALSE) %>% 
    left_join(meta %>% select(cell, myeloid_state), by = "cell") %>% filter(is.finite(pseudotime), myeloid_state %in% 
    STATE_LEVELS) %>% mutate(myeloid_state = factor(myeloid_state, levels = rev(STATE_LEVELS)))

if (nrow(plot_df) == 0) {
    stop("No cells with finite pseudotime and valid myeloid state were found.")
}

summary_df <- plot_df %>% mutate(myeloid_state = factor(as.character(myeloid_state), levels = STATE_LEVELS)) %>% 
    group_by(myeloid_state) %>% summarise(n_cells = n(), median_pseudotime = median(pseudotime, na.rm = TRUE), 
    mean_pseudotime = mean(pseudotime, na.rm = TRUE), q25_pseudotime = quantile(pseudotime, 0.25, na.rm = TRUE), 
    q75_pseudotime = quantile(pseudotime, 0.75, na.rm = TRUE), .groups = "drop") %>% mutate(display_label = unname(DISPLAY_LABELS[as.character(myeloid_state)]))

write.csv(summary_df, file.path(OUT_DIR, "43_SuppFigure2E_pseudotime_summary_by_state_v01.csv"), row.names = FALSE)

print(summary_df)

p <- ggplot(plot_df, aes(x = pseudotime, y = myeloid_state, fill = myeloid_state)) + ggridges::geom_density_ridges(scale = 1.2, 
    rel_min_height = 0.035, linewidth = 0.45, color = "#444444", alpha = 0.8) + geom_segment(data = summary_df %>% 
    mutate(myeloid_state = factor(as.character(myeloid_state), levels = rev(STATE_LEVELS)), y_num = as.numeric(myeloid_state)), 
    aes(x = median_pseudotime, xend = median_pseudotime, y = y_num - 0.04, yend = y_num + 0.34), inherit.aes = FALSE, 
    linewidth = 0.55, color = "black") + scale_fill_manual(values = STATE_COLORS, guide = "none") + scale_y_discrete(labels = function(x) {
    unname(DISPLAY_LABELS[as.character(x)])
}) + scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) + labs(x = "Monocle3 pseudotime", y = NULL) + 
    theme_classic(base_family = "Arial", base_size = 9) + theme(text = element_text(family = "Arial", 
    color = "black"), axis.title.x = element_text(family = "Arial", size = 9.2, margin = margin(t = 5)), 
    axis.text.x = element_text(family = "Arial", size = 7.8, color = "black"), axis.text.y = element_text(family = "Arial", 
        size = 7.6, color = "black"), axis.line = element_line(linewidth = 0.48, color = "black"), axis.ticks = element_line(linewidth = 0.42, 
        color = "black"), panel.grid = element_blank(), plot.margin = margin(6, 10, 6, 6))

W_IN <- 7.6
H_IN <- 6.2

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

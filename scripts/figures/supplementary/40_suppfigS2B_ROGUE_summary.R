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

IN_FILE <- file.path(IN_DIR, "19_myeloid_ROGUE_v02.csv")

OUT_DIR <- file.path(ROOT, "results", "40_SuppFigure2B_ROGUE_summary")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "40_SuppFigure2B_ROGUE_summary_v01")

if (!file.exists(IN_FILE)) {
    stop("Cannot find fixed Step19 ROGUE file:\n", IN_FILE)
}

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

DISPLAY_LABELS <- c(`Homeostatic Microglia` = "Homeostatic", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling", 
    `Activated inflammatory Microglia` = "Activated inflammatory", `Ccl12+ acute Microglia` = "Ccl12+ acute", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory", `Proliferating Microglia` = "Proliferating", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ tissue-repair", 
    `Acute-response Macrophage` = "Acute-response")

STATE_COLORS <- c(`Homeostatic Microglia` = "#4C78A8", `Clec7a+ lipid-associated Microglia` = "#D64F4F", 
    `Activated inflammatory Microglia` = "#F2B134", `Ccl12+ acute Microglia` = "#8E6BBE", `Tnf+ inflammatory Microglia` = "#E77AA8", 
    `Proliferating Microglia` = "#4FA6A2", `MHC-II+ antigen-presenting Macrophage` = "#59A14F", `Gpnmb+ tissue-repair Macrophage` = "#D5A832", 
    `Acute-response Macrophage` = "#9C755F")

rogue <- read.csv(IN_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("source_dataset", "myeloid_state", "ROGUE")

missing_cols <- setdiff(required_cols, colnames(rogue))

if (length(missing_cols) > 0) {
    stop("ROGUE file missing required column(s): ", paste(missing_cols, collapse = ", "))
}

rogue <- rogue %>% filter(myeloid_state %in% STATE_LEVELS, is.finite(ROGUE)) %>% mutate(myeloid_state = factor(myeloid_state, 
    levels = STATE_LEVELS))

if (nrow(rogue) == 0) {
  stop("No valid ROGUE estimates were found.")
}

summary_df <- rogue %>% group_by(myeloid_state) %>% summarise(n_dataset_estimates = n(), median_ROGUE = median(ROGUE, 
    na.rm = TRUE), mean_ROGUE = mean(ROGUE, na.rm = TRUE), min_ROGUE = min(ROGUE, na.rm = TRUE), max_ROGUE = max(ROGUE, 
    na.rm = TRUE), .groups = "drop") %>% mutate(myeloid_state = factor(myeloid_state, levels = STATE_LEVELS), 
    display_label = unname(DISPLAY_LABELS[as.character(myeloid_state)]))

write.csv(summary_df, paste0(OUT_BASE, ".csv"), row.names = FALSE)

print(summary_df)

p <- ggplot(summary_df, aes(x = myeloid_state, y = median_ROGUE, fill = myeloid_state)) + geom_col(width = 0.68, 
    color = "black", linewidth = 0.38) + geom_errorbar(aes(ymin = min_ROGUE, ymax = max_ROGUE), width = 0.18, 
    linewidth = 0.55, color = "black") + geom_point(size = 1.7, shape = 21, fill = "white", color = "black", 
    stroke = 0.45) + scale_fill_manual(values = STATE_COLORS, guide = "none", drop = FALSE) + scale_x_discrete(labels = unname(DISPLAY_LABELS[STATE_LEVELS]), 
    drop = FALSE) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 
    0.02))) + labs(x = NULL, y = "Median ROGUE score") + theme_classic(base_family = "Arial", base_size = 9) + 
    theme(text = element_text(family = "Arial", color = "black"), axis.title.y = element_text(size = 9.2, 
        margin = margin(r = 5)), axis.text.y = element_text(size = 7.8, color = "black"), axis.text.x = element_text(size = 7.4, 
        color = "black", angle = 32, hjust = 1, vjust = 1), axis.line = element_line(linewidth = 0.48, 
        color = "black"), axis.ticks = element_line(linewidth = 0.42, color = "black"), plot.margin = margin(6, 
        6, 6, 6))

W_IN <- 6.9
H_IN <- 4.5

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

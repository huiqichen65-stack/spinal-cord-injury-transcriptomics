rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

UP_FILE <- file.path(ROOT, "results", "48_Bulk_functional_enrichment", "02_48_Bulk_functional_enrichment_v04_ORA_UP_GO_BP.csv")

DOWN_FILE <- file.path(ROOT, "results", "48_Bulk_functional_enrichment", "04_48_Bulk_functional_enrichment_v04_ORA_DOWN_GO_BP.csv")

OUT_DIR <- file.path(ROOT, "results", "55_Figure3E_Bulk_GO_BP_barplot")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "55_Figure3E_Bulk_GO_BP_barplot_v04"

PLOT_TITLE <- "Established late vs Early GO biological process enrichment"
BAR_BORDER <- NA

EXPECTED_UP_GO_SIG <- 1362L
EXPECTED_DOWN_GO_SIG <- 26L

FIG_WIDTH <- 5.8
FIG_HEIGHT <- 7.3

MANUAL_DOWN <- c("cholesterol biosynthetic process via lathosterol", "cholesterol biosynthetic process via desmosterol", 
    "zymosterol biosynthetic process", "cholesterol biosynthetic process", "sterol biosynthetic process")

MANUAL_UP <- c("response to lipoprotein particle", "cellular response to lipoprotein particle stimulus", 
    "lipid localization", "lipid transport", "positive regulation of endocytosis", "lysosome organization", 
    "receptor-mediated endocytosis", "phagocytosis", "microglial cell activation", "myeloid leukocyte activation")

install_if_missing <- require_package

for (p in c("ggplot2", "dplyr", "svglite", "scales")) {
    install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(scales)
})

write_csv <- write_result_csv

wrap_label <- function(x, width = 36L) {
    vapply(x, function(z) {
        paste(strwrap(z, width = width), collapse = "\n")
    }, character(1))
}

if (!file.exists(UP_FILE)) {
  stop("UP GO BP file not found:\n", UP_FILE)
}

if (!file.exists(DOWN_FILE)) {
  stop("DOWN GO BP file not found:\n", DOWN_FILE)
}

up_go <- read.csv(UP_FILE, check.names = FALSE, stringsAsFactors = FALSE)

down_go <- read.csv(DOWN_FILE, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("ID", "Description", "pvalue", "p.adjust", "Count")

for (nm in c("up_go", "down_go")) {
    tab <- get(nm)
    miss <- setdiff(required_cols, names(tab))
    if (length(miss) > 0L) {
        stop(nm, " lacks required column(s): ", paste(miss, collapse = ", "))
    }
}

n_up_sig <- sum(!is.na(up_go$p.adjust) & up_go$p.adjust < 0.05)
n_down_sig <- sum(!is.na(down_go$p.adjust) & down_go$p.adjust < 0.05)

if (n_up_sig != EXPECTED_UP_GO_SIG) {
    stop("UP GO BP significant count mismatch: observed ", n_up_sig, ", expected ", EXPECTED_UP_GO_SIG)
}

if (n_down_sig != EXPECTED_DOWN_GO_SIG) {
    stop("DOWN GO BP significant count mismatch: observed ", n_down_sig, ", expected ", EXPECTED_DOWN_GO_SIG)
}

select_manual <- function(df, pathway_vector, direction_label) {
    idx <- match(pathway_vector, df$Description)
    if (any(is.na(idx))) {
        missing_names <- pathway_vector[is.na(idx)]
        stop("The following manually selected pathways were not found in ", direction_label, " table:\n", 
            paste(missing_names, collapse = "\n"))
    }
    z <- df[idx, , drop = FALSE]
    z$direction <- direction_label
    z$minus_log10_pvalue <- -log10(pmax(z$pvalue, .Machine$double.xmin))
    z$minus_log10_padj <- -log10(pmax(z$p.adjust, .Machine$double.xmin))
    z$Description_plot <- wrap_label(z$Description, width = 36L)
    z
}

plot_down <- select_manual(down_go, MANUAL_DOWN, "Down")
plot_up <- select_manual(up_go, MANUAL_UP, "Up")

plot_df <- dplyr::bind_rows(plot_down, plot_up)

plot_df <- plot_df %>% dplyr::group_by(direction) %>% dplyr::arrange(dplyr::desc(minus_log10_pvalue), 
    p.adjust, pvalue, dplyr::desc(Count), Description, .by_group = TRUE) %>% dplyr::mutate(sort_rank = dplyr::row_number()) %>% 
    dplyr::ungroup()

plot_df$facet_label <- factor(plot_df$direction, levels = c("Down", "Up"))

plot_df$facet_pathway <- paste0(plot_df$direction, "___", plot_df$sort_rank, "___", plot_df$Description_plot)

down_levels <- plot_df$facet_pathway[plot_df$direction == "Down"]
up_levels <- plot_df$facet_pathway[plot_df$direction == "Up"]

plot_df$facet_pathway <- factor(plot_df$facet_pathway, levels = c(rev(down_levels), rev(up_levels)))

max_x <- max(plot_df$minus_log10_pvalue, na.rm = TRUE)

base_label_x <- 0.12

plot_df$label_x <- pmin(pmax(base_label_x, plot_df$minus_log10_pvalue * 0.06), plot_df$minus_log10_pvalue * 
    0.18)

plot_df$label_x <- pmax(plot_df$label_x, 0.08)

plot_df$label_color <- "black"

label_map <- setNames(plot_df$Description_plot, plot_df$facet_pathway)

p <- ggplot(plot_df, aes(x = minus_log10_pvalue, y = facet_pathway, fill = minus_log10_padj)) + geom_col(width = 0.78, 
    color = BAR_BORDER) + geom_text(aes(x = label_x, label = Description_plot), hjust = 0, lineheight = 0.92, 
    family = "Arial", size = 3.1, color = "black", show.legend = FALSE) + facet_grid(facet_label ~ ., 
    scales = "free_y", space = "free_y", switch = "y") + scale_fill_gradient(low = "#DDD3CE", high = "#D86445", 
    name = expression(-Log[10] * "(P adjust)")) + scale_x_continuous(expand = expansion(mult = c(0, 0.04)), 
    breaks = scales::pretty_breaks(n = 5)) + scale_y_discrete(labels = function(x) rep("", length(x))) + 
    labs(title = PLOT_TITLE, x = expression(-Log[10] * "(pvalue)"), y = NULL) + coord_cartesian(clip = "off") + 
    theme_minimal(base_family = "Arial", base_size = 10) + theme(plot.title = if (is.null(PLOT_TITLE)) {
    element_blank()
} else {
    element_text(family = "Arial", size = 11.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 8))
}, axis.title.x = element_text(family = "Arial", size = 10, face = "bold", color = "black"), axis.text.x = element_text(family = "Arial", 
    size = 8.5, color = "black"), axis.text.y = element_blank(), axis.title.y = element_blank(), panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank(), panel.grid.major.x = element_line(linewidth = 0.35, color = "#D8D8D8"), 
    panel.border = element_rect(color = "#3A3A3A", fill = NA, linewidth = 0.6), strip.background = element_blank(), 
    strip.text.y.left = element_text(angle = 90, family = "Arial", size = 12, face = "bold", color = "black"), 
    strip.placement = "outside", legend.position = "right", legend.title = element_text(family = "Arial", 
        size = 10, face = "bold"), legend.text = element_text(family = "Arial", size = 8.5), plot.margin = margin(t = 8, 
        r = 18, b = 8, l = 8, unit = "pt")) + guides(fill = guide_colorbar(barheight = grid::unit(32, 
    "mm"), frame.colour = NA, ticks.colour = "black"))

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))
PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))
PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave(filename = SVG_FILE, plot = p, device = svglite::svglite, width = FIG_WIDTH, height = FIG_HEIGHT, 
    units = "in", bg = "white")

ggsave(filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = FIG_WIDTH, height = FIG_HEIGHT, 
    units = "in", family = "Arial", bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", dpi = 600, 
    bg = "white")

qc <- data.frame(metric = c("UP_source_file", "DOWN_source_file", "database", "panel_layout", "x_axis", 
    "fill_metric", "labels_inside_bar", "bar_order_within_panel", "significant_UP_GO_BP_n", "significant_DOWN_GO_BP_n", 
    "displayed_UP_n", "displayed_DOWN_n", "analysis_rerun", "plot_title_inside"), value = c(UP_FILE, 
    DOWN_FILE, "GO Biological Process only", "Top = Down; Bottom = Up", "-log10(pvalue)", "-log10(p.adjust)", 
    "YES", "Descending -log10(pvalue) from top to bottom within each panel", n_up_sig, n_down_sig, nrow(plot_up), 
    nrow(plot_down), "NO", ifelse(is.null(PLOT_TITLE), "NO", "YES")), stringsAsFactors = FALSE)

write_csv(qc, paste0("01_", PREFIX, "_plot_QC.csv"))

selected_export <- plot_df[, intersect(c("direction", "sort_rank", "ID", "Description", "Count", "GeneRatio", 
    "pvalue", "p.adjust", "minus_log10_pvalue", "minus_log10_padj"), names(plot_df)), drop = FALSE]

write_csv(selected_export, paste0("02_", PREFIX, "_selected_pathways.csv"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, paste0("03_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

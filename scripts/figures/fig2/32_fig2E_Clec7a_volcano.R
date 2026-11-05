source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)

SEED <- 1234
set.seed(SEED)

LOG2FC_CUTOFF <- 0.25
PADJ_CUTOFF <- 0.05

required_pkgs <- c("ggplot2", "dplyr", "ggrastr", "ggrepel", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrastr)
  library(ggrepel)
})

ROOT <- PROJECT_ROOT
RES_DIR <- file.path(ROOT, "results")

STEP21_DIR <- file.path(RES_DIR, "21_Clec7a_DEG_volcano_myeloid")

DEG_FILE <- file.path(STEP21_DIR, "21_Clec7a_vs_all_myeloid_DEG_all_v01.csv")

OUT_DIR <- file.path(RES_DIR, "32_Figure2E_Clec7a_DEG_volcano")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "32_Figure2E_Clec7a_DEG_volcano_v09")

if (!file.exists(DEG_FILE)) {
    stop("Step 21 DEG table was not found:\n", DEG_FILE, "\nDo NOT rerun DEG analysis. Check the Step 21 folder/file name.")
}

deg <- read.csv(DEG_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("gene", "avg_log2FC", "p_val_adj")

missing_cols <- setdiff(required_cols, colnames(deg))

if (length(missing_cols) > 0) {
    stop("Required DEG column(s) missing: ", paste(missing_cols, collapse = ", "))
}

pct_candidates_1 <- c("pct.1", "pct1", "pct_1", "pct_in_group1")
pct_candidates_2 <- c("pct.2", "pct2", "pct_2", "pct_in_group2")

pct1_col <- pct_candidates_1[pct_candidates_1 %in% colnames(deg)][1]
pct2_col <- pct_candidates_2[pct_candidates_2 %in% colnames(deg)][1]

if (is.na(pct1_col) || is.na(pct2_col) || is.null(pct1_col) || is.null(pct2_col)) {
    stop("The DEG table does not contain usable percentage columns ", "(expected pct.1 and pct.2 or equivalents). ", 
        "v03 requires Delta percentage difference for the x-axis.")
}

deg$gene <- as.character(deg$gene)
deg$avg_log2FC <- as.numeric(deg$avg_log2FC)
deg$p_val_adj <- as.numeric(deg$p_val_adj)
deg$pct1_raw <- as.numeric(deg[[pct1_col]])
deg$pct2_raw <- as.numeric(deg[[pct2_col]])

scale_pct <- function(x) {
  mx <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(mx)) return(x)
  if (mx <= 1.01) {
    return(100 * x)
  } else {
    return(x)
  }
}

deg$pct1 <- scale_pct(deg$pct1_raw)
deg$pct2 <- scale_pct(deg$pct2_raw)
deg$delta_pct <- deg$pct1 - deg$pct2

deg <- deg %>% mutate(significance = case_when(!is.na(p_val_adj) & p_val_adj < PADJ_CUTOFF & avg_log2FC >= 
    LOG2FC_CUTOFF ~ "Up", !is.na(p_val_adj) & p_val_adj < PADJ_CUTOFF & avg_log2FC <= -LOG2FC_CUTOFF ~ 
    "Down", TRUE ~ "NS"))

deg$significance <- factor(deg$significance, levels = c("Down", "NS", "Up"))

n_up <- sum(deg$significance == "Up", na.rm = TRUE)
n_down <- sum(deg$significance == "Down", na.rm = TRUE)

if (n_up != 1320 || n_down != 2497) {
    warning("Counts differ from the Step 21 result ", "(expected Up=1320, Down=2497). ", "Check that the correct Step 21 DEG table is being used.")
}

core_up <- c("Clec7a", "Apoc4", "Lpl", "Apoe", "Trem2", "Abca1")

core_down <- c("Ccr2", "S100a4", "Il1b", "Tnf", "Cd36", "Plin2")

volcano_colors <- c(Down = "#4C78A8", NS = "#C9C9C9", Up = "#D64F4F")

plot_df <- bind_rows(deg %>% filter(significance == "NS"), deg %>% filter(significance == "Down"), deg %>% 
    filter(significance == "Up"))

x_q <- stats::quantile(plot_df$delta_pct, probs = c(0.005, 0.995), na.rm = TRUE)
x_display <- c(floor((x_q[1] - 4)/5) * 5, ceiling((x_q[2] + 4)/5) * 5)
x_display[1] <- min(x_display[1], -5)
x_display[2] <- max(x_display[2], 5)

y_q <- stats::quantile(plot_df$avg_log2FC, probs = c(0.005, 0.995), na.rm = TRUE)
y_display <- c(floor((y_q[1] - 0.4) * 2)/2, ceiling((y_q[2] + 0.4) * 2)/2)
y_display[1] <- min(y_display[1], -0.5)
y_display[2] <- max(y_display[2], 0.5)

clec7a_row_for_window <- plot_df %>% filter(gene == "Clec7a")
if (nrow(clec7a_row_for_window) > 0) {
    x_display[2] <- max(x_display[2], ceiling((clec7a_row_for_window$delta_pct[1] + 4)/5) * 5)
    y_display[2] <- max(y_display[2], ceiling((clec7a_row_for_window$avg_log2FC[1] + 0.4) * 2)/2)
}

core_label_df <- plot_df %>% filter(gene %in% unique(c(core_up, core_down)), significance != "NS", delta_pct >= 
    x_display[1], delta_pct <= x_display[2], avg_log2FC >= y_display[1], avg_log2FC <= y_display[2])

extra_down <- plot_df %>% filter(significance == "Down", !gene %in% core_down, delta_pct >= x_display[1], 
    delta_pct <= x_display[2], avg_log2FC >= y_display[1], avg_log2FC <= y_display[2], is.finite(avg_log2FC), 
    is.finite(p_val_adj)) %>% arrange(avg_log2FC, p_val_adj) %>% slice_head(n = 8) %>% pull(gene)

extra_up <- plot_df %>% filter(significance == "Up", !gene %in% core_up, delta_pct >= x_display[1], delta_pct <= 
    x_display[2], avg_log2FC >= y_display[1], avg_log2FC <= y_display[2], is.finite(avg_log2FC), is.finite(p_val_adj)) %>% 
    arrange(desc(avg_log2FC), p_val_adj) %>% slice_head(n = 2) %>% pull(gene)

label_genes_requested <- unique(c(core_up, extra_up, core_down, extra_down))

label_df_visible <- plot_df %>% filter(gene %in% label_genes_requested, significance != "NS", delta_pct >= 
    x_display[1], delta_pct <= x_display[2], avg_log2FC >= y_display[1], avg_log2FC <= y_display[2])

clec7a_force <- plot_df %>% filter(gene == "Clec7a", significance != "NS")

label_df <- bind_rows(label_df_visible, clec7a_force) %>% distinct(gene, .keep_all = TRUE) %>% mutate(label = gene, 
    label_group = case_when(gene %in% core_up ~ "Core up", gene %in% core_down ~ "Core down", significance == 
        "Up" ~ "Extra up", significance == "Down" ~ "Extra down", TRUE ~ "Other"))

message("  Core up visible:   ", paste(intersect(core_up, label_df$gene), collapse = ", "))
message("  Core down visible: ", paste(intersect(core_down, label_df$gene), collapse = ", "))

write.csv(label_df %>% select(gene, avg_log2FC, p_val_adj, pct1, pct2, delta_pct, significance, label_group), 
    paste0(OUT_BASE, "_labeled_genes.csv"), row.names = FALSE)

p <- ggplot(plot_df, aes(x = delta_pct, y = avg_log2FC)) + ggrastr::rasterise(geom_point(aes(color = significance), 
    size = 0.95, alpha = 0.74, stroke = 0), dpi = 600) + geom_vline(xintercept = 0, linetype = "solid", 
    linewidth = 0.45, color = "black") + geom_hline(yintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", 
    linewidth = 0.4, color = "black") + ggrepel::geom_label_repel(data = label_df, aes(x = delta_pct, 
    y = avg_log2FC, label = label), inherit.aes = FALSE, family = "Arial", size = 3.65, color = "black", 
    fill = "white", label.padding = grid::unit(0.07, "lines"), label.r = grid::unit(0.06, "lines"), label.size = 0.2, 
    box.padding = 0.3, point.padding = 0.18, min.segment.length = 0, segment.size = 0.24, segment.color = "black", 
    force = 2, force_pull = 0.55, max.overlaps = Inf, seed = SEED, show.legend = FALSE) + scale_color_manual(values = volcano_colors, 
    breaks = c("Up", "Down", "NS"), labels = c(paste0("Up (", n_up, ")"), paste0("Down (", n_down, ")"), 
        "NS"), name = NULL, drop = FALSE) + scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) + coord_cartesian(xlim = x_display, 
    ylim = y_display, clip = "off") + labs(x = expression(Delta ~ "Percentage difference"), y = "Average log2 fold change") + 
    theme_classic(base_family = "Arial", base_size = 13) + theme(text = element_text(family = "Arial", 
    color = "black"), axis.title = element_text(size = 13.4), axis.text = element_text(size = 12.1, color = "black"), 
    axis.line = element_line(linewidth = 0.55, color = "black"), axis.ticks = element_line(linewidth = 0.45, 
        color = "black"), legend.position = "right", legend.justification = "center", legend.background = element_blank(), 
    legend.key = element_blank(), legend.text = element_text(size = 12.4), plot.margin = margin(14, 26, 
        13, 13)) + guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

W_IN <- 6.7
H_IN <- 5.6

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

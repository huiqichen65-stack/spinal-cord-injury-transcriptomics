source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(20260930)

ROOT <- "D:/zz"

STEP46_RDS <- file.path(ROOT, "results", "46_Bulk_limma_temporal_DEG_analysis", "30_46_Bulk_limma_temporal_DEG_analysis_v06_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "52_Figure3B_Bulk_multiVolcano")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "52_Figure3B_Bulk_multiVolcano_v02"

ADJ_P_CUTOFF <- 0.05
ABS_LOGFC_CUTOFF <- 0.5
LABEL_N_EACH_DIRECTION <- 2L
JITTER_WIDTH <- 0.32

CONTRASTS <- c("15min_vs_Uninjured", "1dpi_vs_Uninjured", "3dpi_vs_Uninjured", "7dpi_vs_Uninjured", "14dpi_vs_Uninjured", 
    "28dpi_vs_Uninjured", "42dpi_vs_Uninjured")

DISPLAY_LABELS <- c("15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi")

names(DISPLAY_LABELS) <- CONTRASTS

TILE_COLORS <- c(`15 min` = "#B8E0D2", `1 dpi` = "#D6EADF", `3 dpi` = "#F6E7A7", `7 dpi` = "#F4C38A", 
    `14 dpi` = "#EFA1A1", `28 dpi` = "#D9B7DB", `42 dpi` = "#CDB7A2")

POINT_COLORS <- c(Up = "#D84A4A", Down = "#3F7FBF", `Not significant` = "#BDBDBD")

install_if_missing <- require_package

for (p in c("ggplot2", "ggrepel", "svglite")) {
    install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(svglite)
})

write_csv <- write_result_csv

select_labels_one_group <- function(df) {
    out <- list()
    for (dir_i in c("Up", "Down")) {
        z <- df[df$status == dir_i, , drop = FALSE]
        if (nrow(z) == 0L) {
            next
        }
        z <- z[order(z$adj.P.Val, -abs(z$logFC), z$gene_symbol_restored), , drop = FALSE]
        out[[dir_i]] <- head(z, LABEL_N_EACH_DIRECTION)
    }
    if (length(out) == 0L) {
        return(df[FALSE, , drop = FALSE])
    }
    do.call(rbind, out)
}

if (!file.exists(STEP46_RDS)) {
    stop("Step46 object not found:\n", STEP46_RDS)
}

obj46 <- readRDS(STEP46_RDS)

if (!"full_time" %in% names(obj46) || !"contrast_tables" %in% names(obj46$full_time)) {
    stop("Step46 object lacks full_time$contrast_tables.")
}

time_tables <- obj46$full_time$contrast_tables

missing_contrasts <- setdiff(CONTRASTS, names(time_tables))

if (length(missing_contrasts) > 0L) {
    stop("Step46 object is missing time contrast(s): ", paste(missing_contrasts, collapse = ", "))
}

required_cols <- c("gene_id", "gene_symbol_restored", "logFC", "P.Value", "adj.P.Val")

all_list <- vector("list", length(CONTRASTS))

names(all_list) <- CONTRASTS

for (i in seq_along(CONTRASTS)) {
    cn <- CONTRASTS[i]
    tt <- time_tables[[cn]]
    missing_cols <- setdiff(required_cols, names(tt))
    if (length(missing_cols) > 0L) {
        stop(cn, " table is missing column(s): ", paste(missing_cols, collapse = ", "))
    }
    if (nrow(tt) != 12977L) {
        stop(cn, " contains ", nrow(tt), " tested genes; expected 12977.")
    }
    if (anyDuplicated(tt$gene_id)) {
        stop("Duplicated gene_id values in ", cn)
    }
    x <- tt[, required_cols, drop = FALSE]
    x$contrast <- cn
    x$timepoint <- DISPLAY_LABELS[[cn]]
    x$group_index <- i
    x$status <- "Not significant"
    x$status[!is.na(x$adj.P.Val) & x$adj.P.Val < ADJ_P_CUTOFF & x$logFC >= ABS_LOGFC_CUTOFF] <- "Up"
    x$status[!is.na(x$adj.P.Val) & x$adj.P.Val < ADJ_P_CUTOFF & x$logFC <= -ABS_LOGFC_CUTOFF] <- "Down"
    x$status <- factor(x$status, levels = c("Up", "Down", "Not significant"))
    all_list[[i]] <- x
}

all_df <- do.call(rbind, all_list)

rownames(all_df) <- NULL

plot_df <- all_df[is.finite(all_df$logFC) & abs(all_df$logFC) >= ABS_LOGFC_CUTOFF, , drop = FALSE]

if (nrow(plot_df) == 0L) {
    stop("No genes remain after |logFC| >= ", ABS_LOGFC_CUTOFF)
}

plot_df$jitter_x <- plot_df$group_index + stats::runif(nrow(plot_df), min = -JITTER_WIDTH, max = JITTER_WIDTH)

count_rows <- lapply(seq_along(CONTRASTS), function(i) {
    cn <- CONTRASTS[i]
    x_all <- all_df[all_df$contrast == cn, , drop = FALSE]
    x_plot <- plot_df[plot_df$contrast == cn, , drop = FALSE]
    data.frame(contrast = cn, display_label = DISPLAY_LABELS[[cn]], tested_genes_n = nrow(x_all), plotted_abs_logFC_ge_0.5_n = nrow(x_plot), 
        Up_effect_significant_n = sum(x_all$status == "Up", na.rm = TRUE), Down_effect_significant_n = sum(x_all$status == 
            "Down", na.rm = TRUE), FDR_lt_0.05_any_effect_n = sum(!is.na(x_all$adj.P.Val) & x_all$adj.P.Val < 
            ADJ_P_CUTOFF, na.rm = TRUE), role = if (cn %in% c("15min_vs_Uninjured", "1dpi_vs_Uninjured")) {
            "Exploratory due to operational-batch imbalance"
        }
        else {
            "Temporal landscape"
        }, stringsAsFactors = FALSE)
})

count_summary <- do.call(rbind, count_rows)

label_list <- lapply(CONTRASTS, function(cn) {
    z <- plot_df[plot_df$contrast == cn, , drop = FALSE]
    select_labels_one_group(z)
})

label_df <- do.call(rbind, label_list)

rownames(label_df) <- NULL

group_ranges <- lapply(seq_along(CONTRASTS), function(i) {
    cn <- CONTRASTS[i]
    z <- plot_df[plot_df$contrast == cn, , drop = FALSE]
    data.frame(contrast = cn, timepoint = DISPLAY_LABELS[[cn]], group_index = i, xmin = i - 0.43, xmax = i + 
        0.43, ymin = min(z$logFC, na.rm = TRUE), ymax = max(z$logFC, na.rm = TRUE), background_fill = ifelse(i%%2L == 
        1L, "#F3F3F3", "#FAFAFA"), stringsAsFactors = FALSE)
})

group_ranges <- do.call(rbind, group_ranges)

tile_df <- data.frame(timepoint = DISPLAY_LABELS, group_index = seq_along(DISPLAY_LABELS), xmin = seq_along(DISPLAY_LABELS) - 
    0.43, xmax = seq_along(DISPLAY_LABELS) + 0.43, ymin = -0.3, ymax = 0.3, stringsAsFactors = FALSE)

y_min_data <- min(plot_df$logFC, na.rm = TRUE)

y_max_data <- max(plot_df$logFC, na.rm = TRUE)

if (!is.finite(y_min_data) || !is.finite(y_max_data)) {
    stop("Invalid y-axis range.")
}

lower_pad <- max(0.35, 0.08 * abs(y_min_data))

upper_pad <- max(0.45, 0.08 * abs(y_max_data))

y_lower <- floor((y_min_data - lower_pad) * 2)/2

y_upper <- ceiling((y_max_data + upper_pad) * 2)/2

if (!is.finite(y_lower) || !is.finite(y_upper) || y_lower >= 0 || y_upper <= 0 || y_lower >= y_upper) {
    stop("Invalid asymmetric Y-axis limits: ", y_lower, " to ", y_upper)
}

p <- ggplot() +

  geom_rect(
    data = group_ranges,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = group_ranges$background_fill,
    color = NA
  ) +

  geom_point(
    data = plot_df[
      plot_df$status == "Not significant",
      ,
      drop = FALSE
    ],
    aes(
      x = jitter_x,
      y = logFC,
      color = status
    ),
    size = 0.62,
    alpha = 0.22,
    stroke = 0
  ) +

  geom_point(
    data = plot_df[
      plot_df$status != "Not significant",
      ,
      drop = FALSE
    ],
    aes(
      x = jitter_x,
      y = logFC,
      color = status
    ),
    size = 0.88,
    alpha = 0.88,
    stroke = 0
  ) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    color = "black"
  ) +

  geom_rect(
    data = tile_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = timepoint
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +

  geom_text(
    data = tile_df,
    aes(
      x = group_index,
      y = 0,
      label = timepoint
    ),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3.0,
    color = "black"
  ) +

  ggrepel::geom_text_repel(
    data = label_df,
    aes(
      x = jitter_x,
      y = logFC,
      label = gene_symbol_restored
    ),
    family = "Arial",
    size = 2.45,
    box.padding = 0.30,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.color = "#555555",
    max.overlaps = Inf,
    force = 0.8,
    seed = 20260930,
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = POINT_COLORS,
    breaks = c(
      "Up",
      "Down",
      "Not significant"
    ),
    name = NULL,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = TILE_COLORS,
    guide = "none"
  ) +

  scale_x_continuous(
    breaks = seq_along(
      DISPLAY_LABELS
    ),
    labels = rep(
      "",
      length(DISPLAY_LABELS)
    ),
    limits = c(
      0.50,
      length(DISPLAY_LABELS) + 0.50
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_y_continuous(
    breaks = scales::pretty_breaks(
      n = 6
    ),
    expand = expansion(
      mult = c(
        0.04,
        0.08
      )
    )
  ) +

  coord_cartesian(
    ylim = c(
      y_lower,
      y_upper
    ),
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = expression(log[2] * "FC")
  ) +

  theme_classic(
    base_family = "Arial",
    base_size = 10
  ) +

  theme(
    plot.title = element_blank(),
    plot.tag = element_blank(),

    axis.title.x = element_blank(),

    axis.title.y = element_text(
      family = "Arial",
      size = 10,
      face = "plain",
      color = "black",
      margin = margin(
        r = 5
      )
    ),

    axis.text.y = element_text(
      family = "Arial",
      size = 8.5,
      color = "black"
    ),

    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),

    axis.line.y = element_line(
      linewidth = 0.55,
      color = "black"
    ),

    axis.ticks.y = element_line(
      linewidth = 0.45,
      color = "black"
    ),

    panel.grid = element_blank(),

    legend.position = "right",
    legend.title = element_blank(),

    legend.text = element_text(
      family = "Arial",
      size = 8.5
    ),

    legend.key.height = grid::unit(
      4,
      "mm"
    ),

    legend.key.width = grid::unit(
      4,
      "mm"
    ),

    plot.margin = margin(
      t = 10,
      r = 8,
      b = 6,
      l = 6,
      unit = "pt"
    )
  ) +

  guides(
    color = guide_legend(
      override.aes = list(
        size = 2.5,
        alpha = 1
      )
    )
  )

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))

PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave(filename = SVG_FILE, plot = p, device = svglite::svglite, width = 9, height = 5.8, units = "in", 
    bg = "white")

ggsave(filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = 9, height = 5.8, units = "in", 
    family = "Arial", bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = 9, height = 5.8, units = "in", dpi = 600, bg = "white")

write_csv(count_summary, paste0("01_", PREFIX, "_DEG_count_summary.csv"))

write_csv(plot_df, paste0("02_", PREFIX, "_plot_data_absLogFC_ge_0.5.csv"))

write_csv(label_df, paste0("03_", PREFIX, "_labeled_genes.csv"))

plot_qc <- data.frame(metric = c("Step46_source", "contrasts_n", "tested_genes_per_contrast", "FDR_cutoff", 
    "absolute_logFC_cutoff", "displayed_gene_rule", "Up_rule", "Down_rule", "Not_significant_rule", "label_rule", 
    "y_axis_rule", "y_axis_lower_limit", "y_axis_upper_limit", "15min_vs_Uninjured_role", "1dpi_vs_Uninjured_role", 
    "panel_letter_inside_plot", "internal_plot_title"), value = c(STEP46_RDS, length(CONTRASTS), 12977, 
    ADJ_P_CUTOFF, ABS_LOGFC_CUTOFF, "|logFC| >= 0.5", "adj.P.Val < 0.05 & logFC >= 0.5", "adj.P.Val < 0.05 & logFC <= -0.5", 
    "|logFC| >= 0.5 but not effect-significant", paste0("Top ", LABEL_N_EACH_DIRECTION, " Up + Top ", 
        LABEL_N_EACH_DIRECTION, " Down per comparison by smallest adjusted P; tie-break larger |logFC|"), 
    "Asymmetric data-driven range based on displayed logFC values", y_lower, y_upper, "Exploratory because operational batch is incompletely balanced", 
    "Exploratory because operational batch is incompletely balanced", "NO", "NO"), stringsAsFactors = FALSE)

write_csv(plot_qc, paste0("04_", PREFIX, "_plot_QC.csv"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = file.path(OUT_DIR, paste0("05_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

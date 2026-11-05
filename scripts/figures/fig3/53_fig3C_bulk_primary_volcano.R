rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

INPUT_CSV <- file.path(ROOT, "results", "46_Bulk_limma_temporal_DEG_analysis", "04_46_Bulk_limma_temporal_DEG_analysis_v06_PRIMARY_EstablishedLate_vs_Early.csv")

OUT_DIR <- file.path(ROOT, "results", "53_Figure3C_Bulk_PRIMARY_volcano")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "53_Figure3C_Bulk_PRIMARY_volcano_v01"

EXPECTED_TESTED <- 12977L
EXPECTED_UP <- 1801L
EXPECTED_DOWN <- 249L

ADJ_P_CUTOFF <- 0.05
ABS_LOGFC_CUTOFF <- 0.5

CORE_LIPID_GENES <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

STEROL_GENES <- c("Cyp51", "Sqle", "Hmgcr", "Dhcr7", "Dhcr24", "Sc5d")

POINT_COLORS <- c(Up = "#D84A4A", Down = "#3F7FBF", `Not significant` = "#C7C7C7")

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

if (!file.exists(INPUT_CSV)) {
    stop("Step46 Primary table not found:\n", INPUT_CSV)
}

df <- read.csv(INPUT_CSV, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("gene_id", "gene_symbol_restored", "logFC", "P.Value", "adj.P.Val")

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0L) {
    stop("Primary table lacks required column(s): ", paste(missing_cols, collapse = ", "))
}

if (nrow(df) != EXPECTED_TESTED) {
    stop("Unexpected tested-gene count: ", nrow(df), "; expected ", EXPECTED_TESTED)
}

if (anyDuplicated(df$gene_id)) {
    stop("Duplicated gene_id values found in Primary table.")
}

df$status <- "Not significant"

df$status[!is.na(df$adj.P.Val) & df$adj.P.Val < ADJ_P_CUTOFF & df$logFC >= ABS_LOGFC_CUTOFF] <- "Up"

df$status[!is.na(df$adj.P.Val) & df$adj.P.Val < ADJ_P_CUTOFF & df$logFC <= -ABS_LOGFC_CUTOFF] <- "Down"

df$status <- factor(df$status, levels = c("Up", "Down", "Not significant"))

n_up <- sum(df$status == "Up", na.rm = TRUE)

n_down <- sum(df$status == "Down", na.rm = TRUE)

if (n_up != EXPECTED_UP || n_down != EXPECTED_DOWN) {
    stop("Primary DEG count validation failed.\n", "Observed Up = ", n_up, ", expected = ", EXPECTED_UP, 
        "\nObserved Down = ", n_down, ", expected = ", EXPECTED_DOWN)
}

adj_for_plot <- pmax(df$adj.P.Val, .Machine$double.xmin)

df$minus_log10_FDR <- -log10(adj_for_plot)

if (any(!is.finite(df$minus_log10_FDR))) {
    stop("Non-finite -log10(FDR) values remain after numerical protection.")
}

label_candidates <- data.frame(gene_symbol_restored = c(CORE_LIPID_GENES, STEROL_GENES), biological_group = c(rep("Lipid handling", 
    length(CORE_LIPID_GENES)), rep("Sterol biosynthesis", length(STEROL_GENES))), stringsAsFactors = FALSE)

idx <- match(label_candidates$gene_symbol_restored, df$gene_symbol_restored)

label_audit <- label_candidates

label_audit$found <- !is.na(idx)

label_audit$logFC <- NA_real_
label_audit$adj.P.Val <- NA_real_
label_audit$status <- NA_character_

found_idx <- which(label_audit$found)

if (length(found_idx) > 0L) {
    label_audit$logFC[found_idx] <- df$logFC[idx[found_idx]]
    label_audit$adj.P.Val[found_idx] <- df$adj.P.Val[idx[found_idx]]
    label_audit$status[found_idx] <- as.character(df$status[idx[found_idx]])
}

label_audit$effect_significant <- label_audit$status %in% c("Up", "Down")

label_df <- df[df$gene_symbol_restored %in% label_audit$gene_symbol_restored[label_audit$effect_significant], 
    , drop = FALSE]

label_df$biological_group <- ifelse(label_df$gene_symbol_restored %in% CORE_LIPID_GENES, "Lipid handling", 
    "Sterol biosynthesis")

core_audit <- label_audit[label_audit$gene_symbol_restored %in% CORE_LIPID_GENES, , drop = FALSE]

if (nrow(core_audit) != length(CORE_LIPID_GENES) || any(!core_audit$found) || any(core_audit$status != 
    "Up")) {
    stop("six-core-gene validation failed. ", "All Clec7a/Trem2/Apoe/Lpl/Apoc4/Abca1 must be significantly UP.")
}

x_abs <- max(abs(df$logFC), na.rm = TRUE)

x_lim <- ceiling(x_abs * 1.08 * 2)/2

y_max <- max(df$minus_log10_FDR, na.rm = TRUE)

y_upper <- ceiling(y_max * 1.08)

if (!is.finite(x_lim) || !is.finite(y_upper) || x_lim <= 0 || y_upper <= 0) {
    stop("Invalid volcano axis limits.")
}

p <- ggplot(
  df,
  aes(
    x = logFC,
    y = minus_log10_FDR
  )
) +

  geom_point(
    data = df[
      df$status == "Not significant",
      ,
      drop = FALSE
    ],
    aes(
      color = status
    ),
    size = 0.85,
    alpha = 0.28,
    stroke = 0
  ) +

  geom_point(
    data = df[
      df$status != "Not significant",
      ,
      drop = FALSE
    ],
    aes(
      color = status
    ),
    size = 1.05,
    alpha = 0.82,
    stroke = 0
  ) +

  geom_vline(
    xintercept = c(
      -ABS_LOGFC_CUTOFF,
      ABS_LOGFC_CUTOFF
    ),
    linetype = "dashed",
    linewidth = 0.40,
    color = "#6E6E6E"
  ) +

  geom_hline(
    yintercept = -log10(
      ADJ_P_CUTOFF
    ),
    linetype = "dashed",
    linewidth = 0.40,
    color = "#6E6E6E"
  ) +

  ggrepel::geom_text_repel(
    data = label_df,
    aes(
      label = gene_symbol_restored
    ),
    family = "Arial",
    size = 3.0,
    color = "black",
    box.padding = 0.40,
    point.padding = 0.18,
    min.segment.length = 0,
    segment.size = 0.30,
    segment.color = "#555555",
    max.overlaps = Inf,
    force = 0.9,
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

  scale_x_continuous(
    limits = c(
      -x_lim,
      x_lim
    ),
    breaks = scales::pretty_breaks(
      n = 6
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.02
      )
    )
  ) +

  scale_y_continuous(
    limits = c(
      0,
      y_upper
    ),
    breaks = scales::pretty_breaks(
      n = 6
    ),
    expand = expansion(
      mult = c(
        0,
        0.04
      )
    )
  ) +

  labs(
    x = expression(log[2] * "FC"),
    y = expression(-log[10] * "(FDR)")
  ) +

  theme_classic(
    base_family = "Arial",
    base_size = 10
  ) +

  theme(
    plot.title = element_blank(),
    plot.tag = element_blank(),

    axis.title = element_text(
      family = "Arial",
      size = 10,
      face = "plain",
      color = "black"
    ),

    axis.text = element_text(
      family = "Arial",
      size = 8.5,
      color = "black"
    ),

    axis.line = element_line(
      linewidth = 0.55,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),

    panel.grid = element_blank(),

    legend.position = "right",

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
      t = 8,
      r = 8,
      b = 6,
      l = 6,
      unit = "pt"
    )
  ) +

  guides(
    color = guide_legend(
      override.aes = list(
        size = 2.8,
        alpha = 1
      )
    )
  )

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))

PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave(filename = SVG_FILE, plot = p, device = svglite::svglite, width = 5.7, height = 5, units = "in", 
    bg = "white")

ggsave(filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = 5.7, height = 5, units = "in", 
    family = "Arial", bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = 5.7, height = 5, units = "in", dpi = 600, bg = "white")

count_summary <- data.frame(metric = c("tested_genes_n", "Up_effect_significant_n", "Down_effect_significant_n", 
    "Not_significant_n", "FDR_cutoff", "absolute_logFC_cutoff"), value = c(nrow(df), n_up, n_down, sum(df$status == 
    "Not significant", na.rm = TRUE), ADJ_P_CUTOFF, ABS_LOGFC_CUTOFF), stringsAsFactors = FALSE)

write_csv(count_summary, paste0("01_", PREFIX, "_DEG_count_QC.csv"))

plot_export <- df[, c("gene_id", "gene_symbol_restored", "logFC", "P.Value", "adj.P.Val", "minus_log10_FDR", 
    "status"), drop = FALSE]

write_csv(plot_export, paste0("02_", PREFIX, "_plot_data.csv"))

write_csv(label_audit, paste0("03_", PREFIX, "_highlight_gene_audit.csv"))

write_csv(label_df[, c("gene_id", "gene_symbol_restored", "biological_group", "logFC", "adj.P.Val", "minus_log10_FDR", 
    "status"), drop = FALSE], paste0("04_", PREFIX, "_labeled_genes.csv"))

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT_DIR, paste0("05_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

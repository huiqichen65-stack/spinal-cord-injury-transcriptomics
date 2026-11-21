rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(20260930)

ROOT <- PROJECT_ROOT
INPUT_GZ <- file.path( ROOT, "results", "147_Reviewer1_RemainingMajor4_bulk_scale_revision", paste0( "05_147_Reviewer1_RemainingMajor4_bulk_scale_revision_v01_",
  "TIMEPOINT_all_genes_transformed_scale.csv.gz" ) )

OUT_DIR <- file.path( ROOT, "results", "52_Figure3B_Bulk_multiVolcano" )

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "52_Figure3B_Bulk_multiVolcano_v03"

ADJ_P_CUTOFF <- 0.05
LABEL_N_EACH_DIRECTION <- 2L
JITTER_WIDTH <- 0.32
EXPECTED_TESTED_PER_CONTRAST <- 12977L

CONTRASTS <- c( "15min_vs_Uninjured", "1dpi_vs_Uninjured", "3dpi_vs_Uninjured", "7dpi_vs_Uninjured", "14dpi_vs_Uninjured", "28dpi_vs_Uninjured", "42dpi_vs_Uninjured" )

DISPLAY_LABELS <- c( "15min_vs_Uninjured" = "15 min", "1dpi_vs_Uninjured" = "1 dpi", "3dpi_vs_Uninjured" = "3 dpi", "7dpi_vs_Uninjured" = "7 dpi", "14dpi_vs_Uninjured" = "14 dpi",
  "28dpi_vs_Uninjured" = "28 dpi", "42dpi_vs_Uninjured" = "42 dpi" )

TILE_COLORS <- c( "15 min" = "#B8E0D2", "1 dpi"  = "#D6EADF", "3 dpi"  = "#F6E7A7", "7 dpi"  = "#F4C38A", "14 dpi" = "#EFA1A1", "28 dpi" = "#D9B7DB", "42 dpi" = "#CDB7A2" )

POINT_COLORS <- c( "Higher" = "#D84A4A", "Lower" = "#3F7FBF", "Not significant" = "#BDBDBD" )

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages( pkg, repos = "https://cloud.r-project.org", dependencies = TRUE )
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package unavailable: ", pkg)
  }
}

for (p in c("ggplot2", "ggrepel", "svglite", "scales")) {
  install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(svglite)
  library(scales)
})

msg <- function(...) {
  cat( sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "" )
}

write_csv <- function(x, filename) {
  utils::write.csv( x, file.path(OUT_DIR, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8" )
}

select_labels_one_group <- function(df) {
  out <- list()

  for (dir_i in c("Higher", "Lower")) {
    z <- df[df$status == dir_i, , drop = FALSE]
    if (nrow(z) == 0L) next

    z <- z[ order( z$adj.P.Val, -abs(z$transformed_scale_difference), z$gene_symbol_restored ), , drop = FALSE ]

    out[[dir_i]] <- head(z, LABEL_N_EACH_DIRECTION)
  }

  if (length(out) == 0L) {
    return(df[FALSE, , drop = FALSE])
  }

  do.call(rbind, out)
}

if (!file.exists(INPUT_GZ)) {
  stop( "Script-147 time-point result not found:\n", INPUT_GZ, "\nRun script 147 first." )
}

msg("Reading script-147 time-point transformed-scale statistics...")

con <- gzfile(INPUT_GZ, open = "rt")
all_df <- tryCatch( utils::read.csv( con, check.names = FALSE, stringsAsFactors = FALSE ), finally = close(con) )

required_cols <- c( "contrast", "display_label", "gene_id", "gene_symbol_restored", "transformed_scale_difference", "P.Value", "adj.P.Val", "FDR_lt_0.05", "scale_direction" )

missing_cols <- setdiff(required_cols, names(all_df))
if (length(missing_cols) > 0L) {
  stop( "Script-147 input lacks required column(s): ", paste(missing_cols, collapse = ", ") )
}

all_df <- all_df[all_df$contrast %in% CONTRASTS, , drop = FALSE]

if (nrow(all_df) != EXPECTED_TESTED_PER_CONTRAST * length(CONTRASTS)) {
  stop( "Unexpected number of rows after contrast filtering: ", nrow(all_df), "; expected ", EXPECTED_TESTED_PER_CONTRAST * length(CONTRASTS) )
}

for (cn in CONTRASTS) {
  n_i <- sum(all_df$contrast == cn)
  if (n_i != EXPECTED_TESTED_PER_CONTRAST) {
    stop( cn, " has ", n_i, " genes; expected ", EXPECTED_TESTED_PER_CONTRAST )
  }
}

if (anyDuplicated(all_df[c("contrast", "gene_id")])) {
  stop("Duplicated contrast + gene_id combinations found.")
}

if (any(!is.finite(all_df$transformed_scale_difference))) {
  stop("Non-finite transformed-scale differences detected.")
}

all_df$status <- "Not significant"

all_df$status[ !is.na(all_df$adj.P.Val) & all_df$adj.P.Val < ADJ_P_CUTOFF & all_df$transformed_scale_difference > 0 ] <- "Higher"

all_df$status[ !is.na(all_df$adj.P.Val) & all_df$adj.P.Val < ADJ_P_CUTOFF & all_df$transformed_scale_difference < 0 ] <- "Lower"

all_df$status <- factor( all_df$status, levels = c("Higher", "Lower", "Not significant") )

expected_direction <- ifelse( all_df$adj.P.Val < ADJ_P_CUTOFF & all_df$transformed_scale_difference > 0, "Higher", ifelse( all_df$adj.P.Val < ADJ_P_CUTOFF &
  all_df$transformed_scale_difference < 0, "Lower", "Not significant" ) )

if (any(as.character(all_df$status) != expected_direction)) {
  stop("Internal FDR-only status reconstruction failed.")
}

plot_df <- all_df

group_index_map <- setNames(seq_along(CONTRASTS), CONTRASTS)
plot_df$group_index <- unname(group_index_map[plot_df$contrast])
plot_df$timepoint <- unname(DISPLAY_LABELS[plot_df$contrast])

plot_df$jitter_x <- plot_df$group_index +
  stats::runif( nrow(plot_df), min = -JITTER_WIDTH, max = JITTER_WIDTH )

count_summary <- do.call(
  rbind,
  lapply(
    CONTRASTS,
    function(cn) {
      z <- plot_df[plot_df$contrast == cn, , drop = FALSE]

      data.frame(
        contrast = cn,
        display_label = DISPLAY_LABELS[[cn]],
        tested_genes_n = nrow(z),
        Higher_FDR_lt_0.05_n = sum(z$status == "Higher", na.rm = TRUE),
        Lower_FDR_lt_0.05_n = sum(z$status == "Lower", na.rm = TRUE),
        FDR_lt_0.05_total_n = sum( z$status %in% c("Higher", "Lower"), na.rm = TRUE ),
        absolute_effect_threshold_used = FALSE,
        effect_interpretation =
          "Difference on supplied transformed expression scale",
        role = if (
          cn %in% c("15min_vs_Uninjured", "1dpi_vs_Uninjured")
        ) {
          "Exploratory due to operational-group imbalance"
        } else {
          "Temporal landscape"
        },
        stringsAsFactors = FALSE
      )
    }
  )
)

label_df <- do.call(
  rbind,
  lapply(
    CONTRASTS,
    function(cn) {
      select_labels_one_group( plot_df[plot_df$contrast == cn, , drop = FALSE] )
    }
  )
)

rownames(label_df) <- NULL

group_ranges <- do.call(
  rbind,
  lapply(
    seq_along(CONTRASTS),
    function(i) {
      cn <- CONTRASTS[i]
      z <- plot_df[plot_df$contrast == cn, , drop = FALSE]

      data.frame( contrast = cn, timepoint = DISPLAY_LABELS[[cn]], group_index = i, xmin = i - 0.43, xmax = i + 0.43, ymin = min(z$transformed_scale_difference, na.rm = TRUE),
        ymax = max(z$transformed_scale_difference, na.rm = TRUE), background_fill = ifelse( i %% 2L == 1L, "#F3F3F3", "#FAFAFA" ), stringsAsFactors = FALSE )
    }
  )
)

tile_df <- data.frame( timepoint = unname(DISPLAY_LABELS[CONTRASTS]), group_index = seq_along(CONTRASTS), xmin = seq_along(CONTRASTS) - 0.43, xmax = seq_along(CONTRASTS) + 0.43,
  ymin = -0.30, ymax = 0.30, stringsAsFactors = FALSE )

y_min_data <- min(plot_df$transformed_scale_difference, na.rm = TRUE)
y_max_data <- max(plot_df$transformed_scale_difference, na.rm = TRUE)

lower_pad <- max(0.35, 0.08 * abs(y_min_data))
upper_pad <- max(0.45, 0.08 * abs(y_max_data))

y_lower <- floor((y_min_data - lower_pad) * 2) / 2
y_upper <- ceiling((y_max_data + upper_pad) * 2) / 2

if (
  !is.finite(y_lower) ||
  !is.finite(y_upper) ||
  y_lower >= 0 ||
  y_upper <= 0 ||
  y_lower >= y_upper
) {
  stop("Invalid transformed-scale y-axis limits.")
}

p <- ggplot() +

  geom_rect( data = group_ranges, aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax ), inherit.aes = FALSE, fill = group_ranges$background_fill, color = NA ) +

  geom_point( data = plot_df[ plot_df$status == "Not significant", , drop = FALSE ], aes( x = jitter_x, y = transformed_scale_difference, color = status ), size = 0.50, alpha = 0.13,
    stroke = 0 ) +

  geom_point( data = plot_df[ plot_df$status != "Not significant", , drop = FALSE ], aes( x = jitter_x, y = transformed_scale_difference, color = status ), size = 0.78, alpha = 0.78,
    stroke = 0 ) +

  geom_hline( yintercept = 0, linewidth = 0.35, color = "black" ) +

  geom_rect( data = tile_df, aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = timepoint ), inherit.aes = FALSE, linewidth = 0.35, color = "black" ) +

  geom_text( data = tile_df, aes( x = group_index, y = 0, label = timepoint ), inherit.aes = FALSE, family = "Arial", size = 3.0, color = "black" ) +

  ggrepel::geom_text_repel( data = label_df, aes( x = jitter_x, y = transformed_scale_difference, label = gene_symbol_restored ), family = "Arial", size = 2.35, box.padding = 0.28,
    point.padding = 0.10, min.segment.length = 0, segment.size = 0.23, segment.color = "#555555", max.overlaps = Inf, force = 0.75, seed = 20260930, show.legend = FALSE ) +

  scale_color_manual( values = POINT_COLORS, breaks = c("Higher", "Lower", "Not significant"), name = NULL, drop = FALSE ) +

  scale_fill_manual( values = TILE_COLORS, guide = "none" ) +

  scale_x_continuous( breaks = seq_along(CONTRASTS), labels = rep("", length(CONTRASTS)), limits = c(0.50, length(CONTRASTS) + 0.50), expand = expansion(mult = c(0, 0)) ) +

  scale_y_continuous( breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0.04, 0.08)) ) +

  coord_cartesian( ylim = c(y_lower, y_upper), clip = "off" ) +

  labs( x = NULL, y = "Difference on supplied transformed\nexpression scale" ) +

  theme_classic( base_family = "Arial", base_size = 10 ) +

  theme( plot.title = element_blank(), plot.tag = element_blank(), axis.title.y = element_text( family = "Arial", size = 9.4, face = "plain", color = "black", margin = margin(r = 6) ),
    axis.text.y = element_text( family = "Arial", size = 8.4, color = "black" ), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line( linewidth = 0.55,
    color = "black" ), panel.grid = element_blank(), legend.position = "right", legend.text = element_text( family = "Arial", size = 8.5 ), plot.margin = margin( t = 8, r = 10, b = 6, l = 8,
    unit = "pt" ) ) +

  guides( color = guide_legend( override.aes = list( size = 2.6, alpha = 1 ) ) )

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))
PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))
PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave( filename = SVG_FILE, plot = p, device = svglite::svglite, width = 7.8, height = 5.5, units = "in", bg = "white" )

ggsave( filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = 7.8, height = 5.5, units = "in", family = "Arial", bg = "white" )

ggsave( filename = PNG_FILE, plot = p, width = 7.8, height = 5.5, units = "in", dpi = 600, bg = "white" )

write_csv( count_summary, paste0("01_", PREFIX, "_FDR_count_summary.csv") )

plot_export <- plot_df[ , c( "contrast", "timepoint", "gene_id", "gene_symbol_restored", "transformed_scale_difference", "P.Value", "adj.P.Val", "status" ), drop = FALSE ]

write_csv( plot_export, paste0("02_", PREFIX, "_plot_data_all_genes.csv") )

write_csv( label_df[ , c( "contrast", "timepoint", "gene_id", "gene_symbol_restored", "transformed_scale_difference", "adj.P.Val", "status" ), drop = FALSE ],
  paste0("03_", PREFIX, "_labeled_genes.csv") )

qc <- data.frame( metric = c( "input_file", "statistics_source", "tested_genes_per_contrast", "display_rule", "absolute_effect_threshold_used", "effect_axis", "RDS_read",
  "bulk_expression_matrix_loaded", "statistical_analysis_rerun" ), value = c( INPUT_GZ, "Script 147 transformed-scale outputs derived from frozen Step46 statistics",
  EXPECTED_TESTED_PER_CONTRAST, "All tested genes plotted; significance defined by BH-FDR < 0.05 and coefficient sign", "NO", "Difference on supplied transformed expression scale", "NO",
  "NO", "NO" ), stringsAsFactors = FALSE )

write_csv( qc, paste0("04_", PREFIX, "_plot_QC.csv") )

report_lines <- c( "Figure 3B — Multi-volcano transformed-scale revision (v03)", "===========================================================",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "", "Input:", INPUT_GZ, "", "Reviewer-1 Remaining Major Comment 4 revision:",
  "The upstream normalization/transformation of the deposited bulk matrix remains undocumented.",
  "Therefore coefficients are displayed as differences on the supplied transformed expression scale.", "No log2 fold-change interpretation is made.",
  "No absolute effect-size threshold is applied.", "", "Display rule:", "Higher = BH-FDR < 0.05 and transformed-scale difference > 0.",
  "Lower = BH-FDR < 0.05 and transformed-scale difference < 0.", "Not significant = BH-FDR >= 0.05.", "All tested genes are plotted.", "", "Counts:", paste(
  capture.output(print(count_summary, row.names = FALSE)), collapse = "\n" ), "", "Interpretation boundary:",
  "15 min vs Uninjured and 1 dpi vs Uninjured remain exploratory because of the documented operational-group imbalance.", "",
  "No RDS was read; the bulk expression matrix was not loaded; Step46 was not rerun." )

writeLines( report_lines, con = file.path(OUT_DIR, paste0("05_", PREFIX, "_REPORT.txt")), useBytes = TRUE )

msg("Figure 3B v03 completed successfully.")
msg("Output directory: ", OUT_DIR)
msg("Displayed all tested genes: ", nrow(plot_df))
msg("Labeled genes: ", nrow(label_df))

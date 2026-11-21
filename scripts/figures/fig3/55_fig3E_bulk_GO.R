rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT
REV147_DIR <- file.path( ROOT, "results", "147_Reviewer1_RemainingMajor4_bulk_scale_revision" )

HIGHER_FILE <- file.path( REV147_DIR, paste0( "08_147_Reviewer1_RemainingMajor4_bulk_scale_revision_v01_", "ORA_HIGHER_GO_BP.csv" ) )

LOWER_FILE <- file.path( REV147_DIR, paste0( "09_147_Reviewer1_RemainingMajor4_bulk_scale_revision_v01_", "ORA_LOWER_GO_BP.csv" ) )

OUT_DIR <- file.path( ROOT, "results", "55_Figure3E_Bulk_GO_BP_barplot" )

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "55_Figure3E_Bulk_GO_BP_barplot_v05"

PLOT_TITLE <- "Established late vs Early: GO BP enrichment"
BAR_BORDER <- NA

FIG_WIDTH <- 6.8
FIG_HEIGHT <- 7.5

MANUAL_LOWER <- c( "cholesterol biosynthetic process via lathosterol", "cholesterol biosynthetic process via desmosterol", "zymosterol biosynthetic process",
  "cholesterol biosynthetic process", "sterol biosynthetic process" )

MANUAL_HIGHER <- c( "response to lipoprotein particle", "cellular response to lipoprotein particle stimulus", "lipid localization", "lipid transport", "positive regulation of endocytosis",
  "lysosome organization", "receptor-mediated endocytosis", "phagocytosis", "microglial cell activation", "myeloid leukocyte activation" )

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages( pkg, repos = "https://cloud.r-project.org", dependencies = TRUE )
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package unavailable: ", pkg)
  }
}

for (p in c("ggplot2", "dplyr", "svglite", "scales")) {
  install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(scales)
})

msg <- function(...) {
  cat( sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "" )
}

write_csv <- function(x, filename) {
  utils::write.csv( x, file.path(OUT_DIR, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8" )
}

wrap_label <- function(x, width = 36L) {
  vapply(
    x,
    function(z) {
      paste( strwrap(z, width = width), collapse = "\n" )
    },
    character(1)
  )
}

if (!file.exists(HIGHER_FILE)) {
  stop( "Revised Higher GO BP file not found:\n", HIGHER_FILE, "\nRun script 147 first." )
}

if (!file.exists(LOWER_FILE)) {
  stop( "Revised Lower GO BP file not found:\n", LOWER_FILE, "\nRun script 147 first." )
}

msg("Reading revised script-147 GO BP ORA tables...")

higher_go <- utils::read.csv( HIGHER_FILE, check.names = FALSE, stringsAsFactors = FALSE )

lower_go <- utils::read.csv( LOWER_FILE, check.names = FALSE, stringsAsFactors = FALSE )

required_cols <- c( "ID", "Description", "pvalue", "p.adjust", "Count" )

for (nm in c("higher_go", "lower_go")) {
  tab <- get(nm)
  miss <- setdiff(required_cols, names(tab))
  if (length(miss) > 0L) {
    stop( nm, " lacks required column(s): ", paste(miss, collapse = ", ") )
  }
}

select_manual <- function(df, pathway_vector, direction_label) {
  idx <- match(pathway_vector, df$Description)

  if (any(is.na(idx))) {
    missing_names <- pathway_vector[is.na(idx)]
    stop( "The following manually selected pathways were not found in ", direction_label, " revised ORA table:\n", paste(missing_names, collapse = "\n") )
  }

  z <- df[idx, , drop = FALSE]

  if (any(is.na(z$p.adjust) | z$p.adjust >= 0.05)) {
    bad <- z$Description[is.na(z$p.adjust) | z$p.adjust >= 0.05]
    stop( "Selected pathway(s) are not FDR < 0.05 in the revised ORA:\n", paste(bad, collapse = "\n") )
  }

  z$direction <- direction_label
  z$minus_log10_pvalue <- -log10( pmax(z$pvalue, .Machine$double.xmin) )
  z$minus_log10_padj <- -log10( pmax(z$p.adjust, .Machine$double.xmin) )
  z$Description_plot <- wrap_label( z$Description, width = 36L )

  z
}

plot_lower <- select_manual( lower_go, MANUAL_LOWER, "Lower" )

plot_higher <- select_manual( higher_go, MANUAL_HIGHER, "Higher" )

plot_df <- dplyr::bind_rows( plot_lower, plot_higher )

plot_df <- plot_df %>%
  dplyr::group_by(direction) %>%
  dplyr::arrange( dplyr::desc(minus_log10_pvalue), p.adjust, pvalue, dplyr::desc(Count), Description, .by_group = TRUE ) %>%
  dplyr::mutate( sort_rank = dplyr::row_number() ) %>%
  dplyr::ungroup()

plot_df$facet_label <- factor( plot_df$direction, levels = c("Lower", "Higher") )

plot_df$facet_pathway <- paste0( plot_df$direction, "___", plot_df$sort_rank, "___", plot_df$Description_plot )

lower_levels <- plot_df$facet_pathway[ plot_df$direction == "Lower" ]

higher_levels <- plot_df$facet_pathway[ plot_df$direction == "Higher" ]

plot_df$facet_pathway <- factor( plot_df$facet_pathway, levels = c( rev(lower_levels), rev(higher_levels) ) )

base_label_x <- 0.12

plot_df$label_x <- pmin( pmax( base_label_x, plot_df$minus_log10_pvalue * 0.06 ), plot_df$minus_log10_pvalue * 0.18 )

plot_df$label_x <- pmax( plot_df$label_x, 0.08 )

p <- ggplot( plot_df, aes( x = minus_log10_pvalue, y = facet_pathway, fill = minus_log10_padj ) ) +

  geom_col( width = 0.78, color = BAR_BORDER ) +

  geom_text( aes( x = label_x, label = Description_plot ), hjust = 0, lineheight = 0.92, family = "Arial", size = 3.0, color = "black", show.legend = FALSE ) +

  facet_grid( facet_label ~ ., scales = "free_y", space = "free_y", switch = "y" ) +

  scale_fill_gradient( low = "#DDD3CE", high = "#D86445", name = expression(-Log[10] * "(P adjust)") ) +

  scale_x_continuous( expand = expansion(mult = c(0, 0.05)), breaks = scales::pretty_breaks(n = 5) ) +

  scale_y_discrete( labels = function(x) rep("", length(x)) ) +

  labs( title = PLOT_TITLE, x = expression(-Log[10] * "(pvalue)"), y = NULL ) +

  coord_cartesian(clip = "off") +

  theme_minimal( base_family = "Arial", base_size = 10 ) +

  theme( plot.title = element_text( family = "Arial", size = 10.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 8) ), axis.title.x = element_text( family = "Arial",
    size = 10, face = "bold", color = "black" ), axis.text.x = element_text( family = "Arial", size = 8.5, color = "black" ), axis.text.y = element_blank(), axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line( linewidth = 0.35, color = "#D8D8D8" ), panel.border = element_rect(
    color = "#3A3A3A", fill = NA, linewidth = 0.6 ), strip.background = element_blank(), strip.text.y.left = element_text( angle = 90, family = "Arial", size = 12, face = "bold",
    color = "black" ), strip.placement = "outside", legend.position = "right", legend.title = element_text( family = "Arial", size = 10, face = "bold" ), legend.text = element_text(
    family = "Arial", size = 8.5 ), plot.margin = margin( t = 10, r = 26, b = 8, l = 10, unit = "pt" ) ) +

  guides( fill = guide_colorbar( barheight = grid::unit(32, "mm"), frame.colour = NA, ticks.colour = "black" ) )

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))
PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))
PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave( filename = SVG_FILE, plot = p, device = svglite::svglite, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", bg = "white" )

ggsave( filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", family = "Arial", bg = "white" )

ggsave( filename = PNG_FILE, plot = p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", dpi = 600, bg = "white" )

qc <- data.frame( metric = c( "Higher_source_file", "Lower_source_file", "foreground_definition", "absolute_effect_threshold_used", "panel_layout", "x_axis", "fill_metric",
  "displayed_Higher_n", "displayed_Lower_n", "old_Step48_ORA_used", "enrichGO_rerun_in_this_script", "RDS_read", "bulk_expression_matrix_loaded", "plot_title" ), value = c( HIGHER_FILE,
  LOWER_FILE, "Script 147 FDR-only Higher/Lower gene sets", "NO", "Top = Lower; Bottom = Higher", "-log10(pvalue)", "-log10(p.adjust)", nrow(plot_higher), nrow(plot_lower), "NO", "NO", "NO",
  "NO", PLOT_TITLE ), stringsAsFactors = FALSE )

write_csv( qc, paste0("01_", PREFIX, "_plot_QC.csv") )

selected_export <- plot_df[ , intersect( c( "direction", "sort_rank", "ID", "Description", "Count", "GeneRatio", "pvalue", "p.adjust", "minus_log10_pvalue", "minus_log10_padj" ),
  names(plot_df) ), drop = FALSE ]

write_csv( selected_export, paste0("02_", PREFIX, "_selected_pathways.csv") )

report <- c( "Figure 3E — GO BP ORA barplot after transformed-scale revision (v05)", "====================================================================",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "", "Revised source files:", HIGHER_FILE, LOWER_FILE, "", "Foreground definition:",
  "Script 147 FDR-only gene sets; direction determined by the sign of the transformed-scale coefficient.", "No |logFC| or other absolute effect-size threshold is used.", "", "Display:",
  paste0("Title = ", PLOT_TITLE), "Top panel = Lower", "Bottom panel = Higher", "x = -log10(pvalue)", "fill = -log10(p.adjust)", "pathway labels are placed inside bars",
  "bars are ordered from longest to shortest within each panel", "", paste0("Displayed Higher pathways = ", nrow(plot_higher)), paste0("Displayed Lower pathways = ", nrow(plot_lower)), "",
  "Important:", "This script does not use the old Step48 ORA tables.", "It does not rerun enrichGO; revised ORA was generated by script 147.",
  "No RDS was read and the bulk expression matrix was not loaded." )

writeLines( report, file.path( OUT_DIR, paste0("03_", PREFIX, "_REPORT.txt") ), useBytes = TRUE )

msg("Figure 3E v05 completed successfully.")
msg( "Displayed pathways: Higher=", nrow(plot_higher), "; Lower=", nrow(plot_lower) )
msg("Output directory: ", OUT_DIR)

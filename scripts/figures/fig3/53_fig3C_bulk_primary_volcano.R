rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT
INPUT_CSV <- file.path( ROOT, "results", "147_Reviewer1_RemainingMajor4_bulk_scale_revision", paste0( "03_147_Reviewer1_RemainingMajor4_bulk_scale_revision_v01_",
  "PRIMARY_all_genes_transformed_scale.csv" ) )

OUT_DIR <- file.path( ROOT, "results", "53_Figure3C_Bulk_PRIMARY_volcano" )

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "53_Figure3C_Bulk_PRIMARY_volcano_v02"

EXPECTED_TESTED <- 12977L
EXPECTED_HIGHER <- 3221L
EXPECTED_LOWER <- 2420L

ADJ_P_CUTOFF <- 0.05

CORE_LIPID_GENES <- c( "Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1" )

STEROL_GENES <- c( "Cyp51", "Sqle", "Hmgcr", "Dhcr7", "Dhcr24", "Sc5d" )

POINT_COLORS <- c( "Higher" = "#D84A4A", "Lower" = "#3F7FBF", "Not significant" = "#C7C7C7" )

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

if (!file.exists(INPUT_CSV)) {
  stop( "Script-147 Primary result not found:\n", INPUT_CSV, "\nRun script 147 first." )
}

msg("Reading script-147 Primary transformed-scale result...")

df <- utils::read.csv( INPUT_CSV, check.names = FALSE, stringsAsFactors = FALSE )

required_cols <- c( "gene_id", "gene_symbol_restored", "transformed_scale_difference", "P.Value", "adj.P.Val", "FDR_lt_0.05", "scale_direction" )

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0L) {
  stop( "Script-147 Primary table lacks required column(s): ", paste(missing_cols, collapse = ", ") )
}

if (nrow(df) != EXPECTED_TESTED) {
  stop( "Unexpected tested-gene count: ", nrow(df), "; expected ", EXPECTED_TESTED )
}

if (anyDuplicated(df$gene_id)) {
  stop("Duplicated gene_id values found in Primary table.")
}

if (any(!is.finite(df$transformed_scale_difference))) {
  stop("Non-finite transformed-scale differences detected.")
}

df$status <- "Not significant"

df$status[ !is.na(df$adj.P.Val) & df$adj.P.Val < ADJ_P_CUTOFF & df$transformed_scale_difference > 0 ] <- "Higher"

df$status[ !is.na(df$adj.P.Val) & df$adj.P.Val < ADJ_P_CUTOFF & df$transformed_scale_difference < 0 ] <- "Lower"

df$status <- factor( df$status, levels = c("Higher", "Lower", "Not significant") )

n_higher <- sum(df$status == "Higher", na.rm = TRUE)
n_lower <- sum(df$status == "Lower", na.rm = TRUE)

if (
  n_higher != EXPECTED_HIGHER ||
  n_lower != EXPECTED_LOWER
) {
  stop( "Script-147 Primary FDR count validation failed.\n", "Observed Higher = ", n_higher, ", expected = ", EXPECTED_HIGHER, "\nObserved Lower = ", n_lower, ", expected = ", EXPECTED_LOWER
    )
}

adj_for_plot <- pmax(df$adj.P.Val, .Machine$double.xmin)
df$minus_log10_FDR <- -log10(adj_for_plot)

if (any(!is.finite(df$minus_log10_FDR))) {
  stop("Non-finite -log10(FDR) values remain.")
}

label_candidates <- data.frame( gene_symbol_restored = c( CORE_LIPID_GENES, STEROL_GENES ), biological_group = c( rep("Lipid handling", length(CORE_LIPID_GENES)),
  rep("Sterol biosynthesis", length(STEROL_GENES)) ), stringsAsFactors = FALSE )

idx <- match( label_candidates$gene_symbol_restored, df$gene_symbol_restored )

label_audit <- label_candidates
label_audit$found <- !is.na(idx)
label_audit$transformed_scale_difference <- NA_real_
label_audit$adj.P.Val <- NA_real_
label_audit$status <- NA_character_

found_idx <- which(label_audit$found)

if (length(found_idx) > 0L) {
  label_audit$transformed_scale_difference[found_idx] <-
    df$transformed_scale_difference[idx[found_idx]]

  label_audit$adj.P.Val[found_idx] <-
    df$adj.P.Val[idx[found_idx]]

  label_audit$status[found_idx] <-
    as.character(df$status[idx[found_idx]])
}

label_audit$FDR_significant <-
  label_audit$status %in% c("Higher", "Lower")

label_df <- df[ df$gene_symbol_restored %in% label_audit$gene_symbol_restored[ label_audit$FDR_significant ], , drop = FALSE ]

label_df$biological_group <- ifelse( label_df$gene_symbol_restored %in% CORE_LIPID_GENES, "Lipid handling", "Sterol biosynthesis" )

core_audit <- label_audit[ label_audit$gene_symbol_restored %in% CORE_LIPID_GENES, , drop = FALSE ]

if (
  nrow(core_audit) != length(CORE_LIPID_GENES) ||
  any(!core_audit$found) ||
  any(core_audit$status != "Higher")
) {
  stop( "Six-core-gene validation failed. ", "All Clec7a/Trem2/Apoe/Lpl/Apoc4/Abca1 must be FDR-significant Higher." )
}

x_abs <- max( abs(df$transformed_scale_difference), na.rm = TRUE )

x_lim <- ceiling(x_abs * 1.08 * 2) / 2

y_max <- max(df$minus_log10_FDR, na.rm = TRUE)
y_upper <- ceiling(y_max * 1.08)

if (
  !is.finite(x_lim) ||
  !is.finite(y_upper) ||
  x_lim <= 0 ||
  y_upper <= 0
) {
  stop("Invalid volcano axis limits.")
}

p <- ggplot( df, aes( x = transformed_scale_difference, y = minus_log10_FDR ) ) +

  geom_point( data = df[ df$status == "Not significant", , drop = FALSE ], aes(color = status), size = 0.78, alpha = 0.22, stroke = 0 ) +

  geom_point( data = df[ df$status != "Not significant", , drop = FALSE ], aes(color = status), size = 1.00, alpha = 0.78, stroke = 0 ) +

  geom_hline( yintercept = -log10(ADJ_P_CUTOFF), linetype = "dashed", linewidth = 0.40, color = "#6E6E6E" ) +

  ggrepel::geom_text_repel( data = label_df, aes(label = gene_symbol_restored), family = "Arial", size = 3.0, color = "black", box.padding = 0.40, point.padding = 0.18,
    min.segment.length = 0, segment.size = 0.30, segment.color = "#555555", max.overlaps = Inf, force = 0.9, seed = 20260930, show.legend = FALSE ) +

  scale_color_manual( values = POINT_COLORS, breaks = c("Higher", "Lower", "Not significant"), name = NULL, drop = FALSE ) +

  scale_x_continuous( limits = c(-x_lim, x_lim), breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0.02, 0.02)) ) +

  scale_y_continuous( limits = c(0, y_upper), breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0, 0.04)) ) +

  labs( x = "Difference on supplied transformed expression scale", y = expression(-log[10] * "(FDR)") ) +

  theme_classic( base_family = "Arial", base_size = 10 ) +

  theme( plot.title = element_blank(), plot.tag = element_blank(), axis.title = element_text( family = "Arial", size = 9.5, face = "plain", color = "black" ), axis.text = element_text(
    family = "Arial", size = 8.5, color = "black" ), axis.line = element_line( linewidth = 0.55, color = "black" ), axis.ticks = element_line( linewidth = 0.45, color = "black" ),
    panel.grid = element_blank(), legend.position = "right", legend.text = element_text( family = "Arial", size = 8.5 ), legend.key.height = grid::unit(4, "mm"),
    legend.key.width = grid::unit(4, "mm"), plot.margin = margin( t = 8, r = 8, b = 10, l = 6, unit = "pt" ) ) +

  guides( color = guide_legend( override.aes = list( size = 2.8, alpha = 1 ) ) )

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))
PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))
PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave( filename = SVG_FILE, plot = p, device = svglite::svglite, width = 6.2, height = 5.1, units = "in", bg = "white" )

ggsave( filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = 6.2, height = 5.1, units = "in", family = "Arial", bg = "white" )

ggsave( filename = PNG_FILE, plot = p, width = 6.2, height = 5.1, units = "in", dpi = 600, bg = "white" )

count_summary <- data.frame( metric = c( "tested_genes_n", "Higher_FDR_lt_0.05_n", "Lower_FDR_lt_0.05_n", "Not_significant_n", "FDR_cutoff", "absolute_effect_threshold_used",
  "effect_interpretation", "RDS_read", "bulk_expression_matrix_loaded", "limma_rerun" ), value = c( nrow(df), n_higher, n_lower, sum(df$status == "Not significant", na.rm = TRUE),
  ADJ_P_CUTOFF, "NO", "Difference on supplied transformed expression scale", "NO", "NO", "NO" ), stringsAsFactors = FALSE )

write_csv( count_summary, paste0("01_", PREFIX, "_FDR_count_QC.csv") )

plot_export <- df[ , c( "gene_id", "gene_symbol_restored", "transformed_scale_difference", "P.Value", "adj.P.Val", "minus_log10_FDR", "status" ), drop = FALSE ]

write_csv( plot_export, paste0("02_", PREFIX, "_plot_data.csv") )

write_csv( label_audit, paste0("03_", PREFIX, "_highlight_gene_audit.csv") )

write_csv( label_df[ , c( "gene_id", "gene_symbol_restored", "biological_group", "transformed_scale_difference", "adj.P.Val", "minus_log10_FDR", "status" ), drop = FALSE ],
  paste0("04_", PREFIX, "_labeled_genes.csv") )

report <- c( "Figure 3C — PRIMARY transformed-scale volcano (v02)", "====================================================", paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
  "Input:", INPUT_CSV, "", "Comparison:", "Established late (14 + 28 dpi) vs Early (Uninjured + 15 min + 1 + 3 dpi)", "", "Significance rule:",
  "BH-adjusted P < 0.05; direction determined only by the sign of the transformed-scale coefficient.", "No absolute effect-size threshold is used.", "", paste0("Tested genes = ", nrow(df)),
  paste0("Higher = ", n_higher), paste0("Lower = ", n_lower), "", "Plot:", "X = Difference on supplied transformed expression scale", "Y = -log10(BH-adjusted P value)",
  "Higher = red; Lower = blue; Not significant = gray.", "No +/-0.5 vertical threshold lines.", "Dashed horizontal line = FDR 0.05.", "", "Interpretation boundary:",
  "The coefficient is not described as log2 fold change because the deposited bulk matrix transformation is undocumented.",
  "No RDS was read; the bulk expression matrix was not loaded; limma was not rerun." )

writeLines( report, file.path(OUT_DIR, paste0("05_", PREFIX, "_REPORT.txt")), useBytes = TRUE )

msg("Figure 3C v02 completed successfully.")
msg( "Tested = ", nrow(df), "; Higher = ", n_higher, "; Lower = ", n_lower )
msg("Labels shown = ", nrow(label_df))
msg("Output directory: ", OUT_DIR)

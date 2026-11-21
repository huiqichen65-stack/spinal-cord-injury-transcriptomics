rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT
INPUT_FILE <- file.path( ROOT, "results", "73_GSE195783_Clec7a_lipid_spatial_coupling", "02_73_Clec7a_lipid_spot_scores_section_aware_v03.csv.gz" )

OUT_DIR <- file.path( ROOT, "results", "76_Figure4G_Clec7a_lipid_density_scatter_28dpi" )

dir.create( OUT_DIR, recursive = TRUE, showWarnings = FALSE )

VERSION <- "v04"
TARGET_TIMEPOINT <- "28dpi"

COMPARISON_TIMEPOINTS <- c( "14dpi", "28dpi" )

required_pkgs <- c( "ggplot2", "dplyr", "MASS" )

missing_pkgs <- required_pkgs[ !vapply( required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1) ) ]

if (length(missing_pkgs) > 0L) {
  stop( "Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step76." )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(MASS)
})

read_csv_auto <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

write_csv_utf8 <- function(x, filename) {
  write.csv( x, file = filename, row.names = FALSE, fileEncoding = "UTF-8", na = "" )
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(z))
  out[z %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[z %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found:\n", INPUT_FILE)
}

message("============================================================")
message("Step76 / Figure4G")
message("28 dpi Clec7a-state / lipid-program density scatter")
message("Version: ", VERSION)
message("============================================================")
message("Input file: ", INPUT_FILE)
message("No RCTD or scoring will be rerun.")

raw0 <- read_csv_auto(INPUT_FILE)

required_cols <- c( "barcode", "GSM", "sample_label", "timepoint", "section_uid_internal", "Clec7a_RCTDnorm", "score_Core_lipid_handling" )

missing_cols <- setdiff(required_cols, colnames(raw0))
if (length(missing_cols) > 0L) {
  stop( "Input file is missing required column(s): ", paste(missing_cols, collapse = " | ") )
}

if ("section_eligible" %in% colnames(raw0)) {
  raw0$section_eligible <- as_logical_safe(raw0$section_eligible)
  if (any(is.na(raw0$section_eligible))) {
    stop("Could not safely parse section_eligible.")
  }
} else {
  raw0$section_eligible <- !is.na(raw0$section_uid_internal) &
    nzchar(trimws(as.character(raw0$section_uid_internal)))
}

all_data <- raw0 %>%
  transmute( barcode = as.character(barcode), GSM = as.character(GSM), sample_label = as.character(sample_label), timepoint = as.character(timepoint),
    section_uid_internal = as.character(section_uid_internal), section_eligible = section_eligible, Clec7a_RCTDnorm = as.numeric(Clec7a_RCTDnorm),
    score_Core_lipid_handling = as.numeric(score_Core_lipid_handling) ) %>%
  filter( section_eligible, timepoint %in% COMPARISON_TIMEPOINTS, is.finite(Clec7a_RCTDnorm), is.finite(score_Core_lipid_handling) )

dat <- all_data %>%
  filter( timepoint == TARGET_TIMEPOINT )

if (nrow(dat) == 0L) {
  stop("No eligible spots found for ", TARGET_TIMEPOINT, ".")
}

if (n_distinct(dat$sample_label) != 2L) {
  stop( "Expected 2 capture areas at 28 dpi; found ", n_distinct(dat$sample_label), "." )
}

if (n_distinct(dat$section_uid_internal) != 8L) {
  stop( "Expected 8 major sections at 28 dpi; found ", n_distinct(dat$section_uid_internal), "." )
}

rho_obj <- suppressWarnings( cor.test( dat$Clec7a_RCTDnorm, dat$score_Core_lipid_handling, method = "spearman", exact = FALSE ) )

rho <- unname(rho_obj$estimate)
pval <- rho_obj$p.value
n_spots <- nrow(dat)

fit <- lm( score_Core_lipid_handling ~ Clec7a_RCTDnorm, data = dat )

dens_palette <- colorRampPalette( c( "#2D2B8E", "#5E3C99", "#A12CC8", "#E84A8A", "#F98E52", "#F2E85E" ) )

dat$point_col <- densCols( x = dat$Clec7a_RCTDnorm, y = dat$score_Core_lipid_handling, colramp = dens_palette )

x_rng_shared <- range( all_data$Clec7a_RCTDnorm, na.rm = TRUE )

y_rng_shared <- range( all_data$score_Core_lipid_handling, na.rm = TRUE )

x_span_shared <- diff( x_rng_shared )

y_span_shared <- diff( y_rng_shared )

if (!is.finite(x_span_shared) || x_span_shared <= 0) {
  stop("Invalid shared X-axis range.")
}

if (!is.finite(y_span_shared) || y_span_shared <= 0) {
  stop("Invalid shared Y-axis range.")
}

X_LOWER <- 0
X_UPPER <- x_rng_shared[2] + 0.06 * x_span_shared

Y_LOWER <- y_rng_shared[1] - 0.05 * y_span_shared
Y_UPPER <- y_rng_shared[2] + 0.05 * y_span_shared

label_text <- paste0( "\u03C1 = ", sprintf("%.2f", rho) )

LABEL_X <- X_UPPER - 0.035 * (X_UPPER - X_LOWER)
LABEL_Y <- Y_UPPER - 0.045 * (Y_UPPER - Y_LOWER)

p <- ggplot( dat, aes( x = Clec7a_RCTDnorm, y = score_Core_lipid_handling ) ) +
  geom_point( color = dat$point_col, size = 1.8, alpha = 0.90, stroke = 0 ) +
  geom_smooth( method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.9 ) +
  annotate( "text", x = LABEL_X, y = LABEL_Y, hjust = 1, vjust = 1, size = 3.7, fontface = "bold", label = label_text ) +
  labs( title = "28 dpi", x = "Normalized Clec7a+ state RCTD weight", y = "Lipid-handling program score" ) +
  coord_cartesian( xlim = c(X_LOWER, X_UPPER), ylim = c(Y_LOWER, Y_UPPER), expand = FALSE ) +
  theme_classic( base_family = "Arial", base_size = 12 ) +
  theme( plot.title = element_text(face = "bold", size = 14, hjust = 0.5), axis.title = element_text(face = "bold", size = 11.5), axis.text = element_text(size = 9.5, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"), axis.ticks = element_line(linewidth = 0.6, color = "black"), plot.margin = margin(7, 7, 7, 7, unit = "mm") )

SUMMARY_FILE <- file.path( OUT_DIR, paste0("01_76_Figure4G_spot_summary_", VERSION, ".csv") )

summary_df <- data.frame( timepoint = TARGET_TIMEPOINT, n_spots = n_spots, n_capture_areas = n_distinct(dat$sample_label), n_sections = n_distinct(dat$section_uid_internal),
  spearman_rho = rho, spearman_p = pval, lm_intercept = coef(fit)[1], lm_slope = coef(fit)[2], shared_x_lower = X_LOWER, shared_x_upper = X_UPPER, shared_y_lower = Y_LOWER,
  shared_y_upper = Y_UPPER, stringsAsFactors = FALSE )

write_csv_utf8(summary_df, SUMMARY_FILE)

PDF_FILE <- file.path( OUT_DIR, paste0("02_76_Figure4G_Clec7a_lipid_density_scatter_", VERSION, ".pdf") )

PNG_FILE <- file.path( OUT_DIR, paste0("02_76_Figure4G_Clec7a_lipid_density_scatter_", VERSION, ".png") )

SVG_FILE <- file.path( OUT_DIR, paste0("02_76_Figure4G_Clec7a_lipid_density_scatter_", VERSION, ".svg") )

ggsave( filename = PDF_FILE, plot = p, width = 4.5, height = 4.1, units = "in", device = grDevices::cairo_pdf, bg = "white" )

ggsave( filename = PNG_FILE, plot = p, width = 4.5, height = 4.1, units = "in", dpi = 600, bg = "white" )

ggsave( filename = SVG_FILE, plot = p, width = 4.5, height = 4.1, units = "in", device = "svg", bg = "white" )

final_qc <- data.frame( check = c( "Frozen Step73 v03 spot-score file was used", "Target timepoint is 28dpi", "Exactly 2 capture areas are present", "Exactly 8 major sections are present",
  "All plotted spots are section-eligible", "All X values are finite", "All Y values are finite", "Spearman correlation computed", "LM fit computed for display line only",
  "Shared 14/28 dpi axis limits were used", "Only rho is shown in the main figure", "Spot-level P value is omitted from the main figure", "Spot count is omitted from the main figure",
  "PDF generated", "PNG generated", "SVG generated" ), pass = c( identical(basename(INPUT_FILE), "02_73_Clec7a_lipid_spot_scores_section_aware_v03.csv.gz"),
  identical(TARGET_TIMEPOINT, "28dpi"), n_distinct(dat$sample_label) == 2L, n_distinct(dat$section_uid_internal) == 8L, all(dat$section_eligible), all(is.finite(dat$Clec7a_RCTDnorm)),
  all(is.finite(dat$score_Core_lipid_handling)), is.finite(rho), is.finite(coef(fit)[2]), is.finite(X_LOWER) & is.finite(X_UPPER) & is.finite(Y_LOWER) & is.finite(Y_UPPER) &
  X_LOWER == 0 & X_UPPER > X_LOWER & Y_UPPER > Y_LOWER, TRUE, TRUE, TRUE, file.exists(PDF_FILE), file.exists(PNG_FILE), file.exists(SVG_FILE) ), stringsAsFactors = FALSE )

QC_FILE <- file.path( OUT_DIR, paste0("03_76_Figure4G_FINAL_QC_", VERSION, ".csv") )
write_csv_utf8(final_qc, QC_FILE)

REPORT_FILE <- file.path( OUT_DIR, paste0("04_76_Figure4G_REPORT_", VERSION, ".txt") )

report_lines <- c( "============================================================", "Figure4G / Step76", "28 dpi Clec7a-state / lipid-program density scatter", paste0("Version: ", VERSION),
  "============================================================", "", "INPUT", paste0("  ", INPUT_FILE), "", "DISPLAY", "  Each point = one eligible 28 dpi Visium spot.",
  "  X = normalized Clec7a+ state RCTD weight (Clec7a_RCTDnorm).", "  Y = lipid-handling program score (score_Core_lipid_handling).", "  Point color = local point density via densCols.",
  "  Black line = linear trend for visualization only.", "  Only Spearman rho is shown in the panel, at the top-right.", "  Spot-level P and n are NOT shown in the main figure.",
  "  Figure4F and Figure4G use the SAME X/Y axis limits.", "  RCTD weight is a model-derived spatial estimate, not directly measured cell abundance.",
  "  The plotted spatial signal represents predicted transcriptomic enrichment rather than histologically confirmed anatomical localization.", "", "DESCRIPTIVE SUMMARY",
  paste0("  n spots = ", format(n_spots, big.mark = ",")), paste0("  capture areas = ", n_distinct(dat$sample_label)), paste0("  sections = ", n_distinct(dat$section_uid_internal)),
  paste0("  Spearman rho = ", sprintf("%.4f", rho)), paste0("  descriptive spot-level P = ", format.pval(pval, digits = 3, eps = 1e-300)),
  paste0("  LM slope = ", sprintf("%.6f", coef(fit)[2])), paste0("  shared X limits = ", sprintf("%.4f", X_LOWER), " to ", sprintf("%.4f", X_UPPER)),
  paste0("  shared Y limits = ", sprintf("%.4f", Y_LOWER), " to ", sprintf("%.4f", Y_UPPER)), "", "INTERPRETATION LIMIT", "  This is a descriptive spot-level association plot.",
  "  Visium spots are multicellular mixtures and not independent animals.", "  The displayed rho is descriptive; spot-level P is not used as manuscript inference.",
  "  Therefore the panel supports association, not causal inference.",
  "  It also does not establish anatomically validated localization because matched histological registration is unavailable.", "", "FINAL QC",
  paste0(ifelse(final_qc$pass, "[PASS] ", "[CHECK] "), final_qc$check), "", "============================================================" )

writeLines(report_lines, con = REPORT_FILE, useBytes = TRUE)

message("")
message("Step76 finished.")
print(summary_df)
print(final_qc)
message("")
message("Please send me these 3 files:")
message("  ", basename(PDF_FILE))
message("  ", basename(QC_FILE))
message("  ", basename(REPORT_FILE))

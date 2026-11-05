source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
source("config.R", local = TRUE)


options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(1234)

ROOT <- PROJECT_ROOT

AUDIT_DIR <- file.path(ROOT, "results", "67_GSE195783_spatial_section_identity_audit")

INPUT_SPOT_COMPONENT_FILE <- file.path(AUDIT_DIR, "05_67_GSE195783_spatial_section_identity_audit_spot_component_map_v01.csv.gz")

INPUT_SAMPLE_SUMMARY_FILE <- file.path(AUDIT_DIR, "01_67_GSE195783_spatial_section_identity_audit_sample_summary_v01.csv")

INPUT_STABILITY_FILE <- file.path(AUDIT_DIR, "04_67_GSE195783_spatial_section_identity_audit_stability_summary_v01.csv")

OUT_DIR <- file.path(ROOT, "results", "68_GSE195783_spatial_section_map_finalize")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v01"
STEP <- 68L

EXPECTED_CAPTURE_AREAS <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", 
    "28dpi_rep1", "28dpi_rep2")

EXPECTED_MAJOR_SECTIONS_PER_CAPTURE <- 4L
EXPECTED_TOTAL_MAJOR_SECTIONS <- 32L

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

POSITION_LEVELS <- c("upper_left", "upper_right", "lower_left", "lower_right")

required_pkgs <- c("dplyr", "ggplot2")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

write_csv_utf8 <- write_csv_utf8

read_csv_gz_base <- read_gzip_csv

write_csv_gz_base <- write_gzip_csv

assign_component_positions <- function(
    component_summary_one_sample) {

  z <- component_summary_one_sample

  if (nrow(z) != 4L) {

    stop(
      "Expected exactly 4 major components for sample ",
      unique(z$sample_label),
      ", found ",
      nrow(z),
      "."
    )
  }

  if (anyDuplicated(
    z$component_id_internal
  ) > 0L) {

    stop(
      "Duplicated component_id_internal detected for sample ",
      unique(z$sample_label),
      "."
    )
  }

  if (
    any(
      !is.finite(
        z$centroid_x
      )
    ) ||
    any(
      !is.finite(
        z$centroid_y
      )
    )
  ) {

    stop(
      "Non-finite component centroid detected for sample ",
      unique(z$sample_label),
      "."
    )
  }

  ord_y <- order(
    z$centroid_y,
    z$centroid_x
  )

  upper_idx <- ord_y[1:2]
  lower_idx <- ord_y[3:4]

  upper <- z[
    upper_idx,
    ,
    drop = FALSE
  ]

  lower <- z[
    lower_idx,
    ,
    drop = FALSE
  ]

  upper <- upper[
    order(
      upper$centroid_x
    ),
    ,
    drop = FALSE
  ]

  lower <- lower[
    order(
      lower$centroid_x
    ),
    ,
    drop = FALSE
  ]

  upper$section_position_internal <- c(
    "upper_left",
    "upper_right"
  )

  lower$section_position_internal <- c(
    "lower_left",
    "lower_right"
  )

  out <- bind_rows(
    upper,
    lower
  )

  out$section_position_internal <- factor(
    out$section_position_internal,
    levels = POSITION_LEVELS
  )

  out <- out %>%
    arrange(
      section_position_internal
    )

  out$section_position_internal <- as.character(
    out$section_position_internal
  )

  out
}

required_files <- c(INPUT_SPOT_COMPONENT_FILE, INPUT_SAMPLE_SUMMARY_FILE, INPUT_STABILITY_FILE)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop("Required audit file(s) missing:\n", paste(missing_files, collapse = "\n"), "\nRun the completed section-identity audit first.")
}

spot_map <- read_csv_gz_base(INPUT_SPOT_COMPONENT_FILE)

sample_summary <- read.csv(INPUT_SAMPLE_SUMMARY_FILE, stringsAsFactors = FALSE, check.names = FALSE)

stability_summary <- read.csv(INPUT_STABILITY_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_spot_cols <- c("barcode", "GSM", "timepoint", "replicate", "sample_label", "array_x", "array_y", 
    "component_id_internal", "is_major_component", "major_component_rank_internal")

missing_spot_cols <- setdiff(required_spot_cols, colnames(spot_map))

if (length(missing_spot_cols) > 0L) {
    stop("Spot component map is missing required column(s): ", paste(missing_spot_cols, collapse = ", "))
}

required_sample_cols <- c("sample_label", "timepoint", "replicate", "GSM", "n_spots_total", "n_major_components", 
    "n_tiny_components", "n_spots_in_major_components", "n_spots_in_tiny_components", "fraction_spots_in_major_components", 
    "expected_four_major_components")

missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_summary))

if (length(missing_sample_cols) > 0L) {
    stop("Sample summary is missing required column(s): ", paste(missing_sample_cols, collapse = ", "))
}

required_stability_cols <- c("sample_label", "timepoint", "replicate", "GSM", "n_thresholds_tested", 
    "n_thresholds_with_four_major_components", "all_thresholds_support_four", "major_component_count_range")

missing_stability_cols <- setdiff(required_stability_cols, colnames(stability_summary))

if (length(missing_stability_cols) > 0L) {
    stop("Stability summary is missing required column(s): ", paste(missing_stability_cols, collapse = ", "))
}

observed_samples <- sort(unique(as.character(spot_map$sample_label)))

if (!setequal(observed_samples, EXPECTED_CAPTURE_AREAS)) {
    stop("Capture-area labels do not match the expected 8 samples.\n", "Observed: ", paste(observed_samples, 
        collapse = " | "))
}

if (any(sample_summary$n_major_components != EXPECTED_MAJOR_SECTIONS_PER_CAPTURE)) {
    bad <- sample_summary$sample_label[sample_summary$n_major_components != EXPECTED_MAJOR_SECTIONS_PER_CAPTURE]
    stop("Primary audit did not identify exactly four major components in: ", paste(bad, collapse = " | "))
}

major_spots <- spot_map %>% filter(is_major_component %in% TRUE)

component_summary <- major_spots %>% group_by(sample_label, GSM, timepoint, replicate, component_id_internal) %>% 
    summarise(n_spots = n(), centroid_x = mean(array_x), centroid_y = mean(array_y), min_x = min(array_x), 
        max_x = max(array_x), min_y = min(array_y), max_y = max(array_y), .groups = "drop")

component_counts <- component_summary %>% count(sample_label, name = "n_major_components")

if (any(component_counts$n_major_components != EXPECTED_MAJOR_SECTIONS_PER_CAPTURE)) {
    stop("Component table does not contain exactly four major components ", "for every capture area.")
}

manifest_list <- list()

for (samp in EXPECTED_CAPTURE_AREAS) {
    z <- component_summary %>% filter(sample_label == samp)
    z <- assign_component_positions(z)
    manifest_list[[samp]] <- z
}

section_manifest <- bind_rows(manifest_list)

section_manifest$section_uid_internal <- paste0(section_manifest$sample_label, "__", section_manifest$section_position_internal)

section_manifest$timepoint <- factor(section_manifest$timepoint, levels = TIME_LEVELS)

section_manifest$sample_label <- factor(section_manifest$sample_label, levels = EXPECTED_CAPTURE_AREAS)

section_manifest$section_position_internal <- factor(section_manifest$section_position_internal, levels = POSITION_LEVELS)

section_manifest <- section_manifest %>% arrange(timepoint, sample_label, section_position_internal)

section_manifest$timepoint <- as.character(section_manifest$timepoint)

section_manifest$sample_label <- as.character(section_manifest$sample_label)

section_manifest$section_position_internal <- as.character(section_manifest$section_position_internal)

geometry_qc_rows <- list()

for (samp in EXPECTED_CAPTURE_AREAS) {
    z <- section_manifest %>% filter(sample_label == samp)
    z_ul <- z[z$section_position_internal == "upper_left", , drop = FALSE]
    z_ur <- z[z$section_position_internal == "upper_right", , drop = FALSE]
    z_ll <- z[z$section_position_internal == "lower_left", , drop = FALSE]
    z_lr <- z[z$section_position_internal == "lower_right", , drop = FALSE]
    geometry_qc_rows[[samp]] <- data.frame(sample_label = samp, upper_row_x_order_ok = z_ul$centroid_x < 
        z_ur$centroid_x, lower_row_x_order_ok = z_ll$centroid_x < z_lr$centroid_x, left_column_y_order_ok = z_ul$centroid_y < 
        z_ll$centroid_y, right_column_y_order_ok = z_ur$centroid_y < z_lr$centroid_y, upper_lower_centroid_y_gap = mean(c(z_ll$centroid_y, 
        z_lr$centroid_y)) - mean(c(z_ul$centroid_y, z_ur$centroid_y)), left_right_upper_x_gap = z_ur$centroid_x - 
        z_ul$centroid_x, left_right_lower_x_gap = z_lr$centroid_x - z_ll$centroid_x, stringsAsFactors = FALSE)
}

geometry_qc <- bind_rows(geometry_qc_rows)

geometry_qc$all_position_checks_pass <- with(geometry_qc, upper_row_x_order_ok & lower_row_x_order_ok & 
    left_column_y_order_ok & right_column_y_order_ok & upper_lower_centroid_y_gap > 0 & left_right_upper_x_gap > 
    0 & left_right_lower_x_gap > 0)

if (!all(geometry_qc$all_position_checks_pass)) {
    bad <- geometry_qc$sample_label[!geometry_qc$all_position_checks_pass]
    stop("Deterministic 2x2 positional assignment failed geometry QC for: ", paste(bad, collapse = " | "))
}

section_lookup <- section_manifest %>% select(sample_label, component_id_internal, section_position_internal, 
    section_uid_internal)

spot_section_map <- spot_map %>% left_join(section_lookup, by = c("sample_label", "component_id_internal"))

spot_section_map$section_status <- ifelse(spot_section_map$is_major_component %in% TRUE, "major_coordinate_defined_section", 
    "unassigned_tiny_component")

spot_section_map$section_position_internal[!spot_section_map$is_major_component] <- NA_character_

spot_section_map$section_uid_internal[!spot_section_map$is_major_component] <- NA_character_

spot_section_map$section_join_key <- paste(as.character(spot_section_map$GSM), as.character(spot_section_map$array_x), 
    as.character(spot_section_map$array_y), sep = "||")

if (anyDuplicated(spot_section_map$section_join_key) > 0L) {
    stop("Duplicate GSM + array_x + array_y join keys detected.")
}

if (any(spot_section_map$is_major_component & is.na(spot_section_map$section_uid_internal))) {
    stop("At least one major-component spot failed section-ID assignment.")
}

if (any(!spot_section_map$is_major_component & !is.na(spot_section_map$section_uid_internal))) {
    stop("At least one tiny-component spot was incorrectly assigned to a section.")
}

capture_area_summary <- spot_section_map %>% group_by(sample_label, GSM, timepoint, replicate) %>% summarise(n_spots_total = n(), 
    n_spots_major_sections = sum(is_major_component), n_spots_unassigned_tiny = sum(!is_major_component), 
    fraction_spots_major_sections = mean(is_major_component), n_major_sections = n_distinct(section_uid_internal[!is.na(section_uid_internal)]), 
    .groups = "drop")

capture_area_summary <- capture_area_summary %>% left_join(stability_summary %>% select(sample_label, 
    n_thresholds_tested, n_thresholds_with_four_major_components, all_thresholds_support_four, major_component_count_range), 
    by = "sample_label")

SECTION_MANIFEST_FILE <- file.path(OUT_DIR, paste0("01_68_GSE195783_spatial_section_map_section_manifest_", 
    VERSION, ".csv"))

SPOT_SECTION_MAP_FILE <- file.path(OUT_DIR, paste0("02_68_GSE195783_spatial_section_map_spot_to_section_", 
    VERSION, ".csv.gz"))

CAPTURE_SUMMARY_FILE <- file.path(OUT_DIR, paste0("03_68_GSE195783_spatial_section_map_capture_area_summary_", 
    VERSION, ".csv"))

GEOMETRY_QC_FILE <- file.path(OUT_DIR, paste0("04_68_GSE195783_spatial_section_map_geometry_QC_", VERSION, 
    ".csv"))

write_csv_utf8(section_manifest, SECTION_MANIFEST_FILE)

write_csv_gz_base(spot_section_map, SPOT_SECTION_MAP_FILE)

write_csv_utf8(capture_area_summary, CAPTURE_SUMMARY_FILE)

write_csv_utf8(geometry_qc, GEOMETRY_QC_FILE)

plot_df <- spot_section_map

plot_df$plot_group <- ifelse(plot_df$is_major_component, plot_df$section_position_internal, "unassigned_tiny")

plot_df$sample_label <- factor(plot_df$sample_label, levels = EXPECTED_CAPTURE_AREAS)

diagnostic_palette <- c(upper_left = "#4C78A8", upper_right = "#F58518", lower_left = "#54A24B", lower_right = "#B279A2", 
    unassigned_tiny = "#BDBDBD")

p_diag <- ggplot(plot_df, aes(x = array_x, y = array_y, color = plot_group)) + geom_point(size = 0.55, 
    alpha = 0.95) + coord_equal() + scale_y_reverse() + scale_color_manual(values = diagnostic_palette, 
    drop = FALSE) + facet_wrap(~sample_label, ncol = 4) + labs(title = "Finalized coordinate-defined tissue-section map", 
    subtitle = paste0("Four major disconnected tissue pieces per Visium capture area; ", "internal section labels are not displayed"), 
    x = "Visium array X", y = "Visium array Y") + theme_classic(base_family = "Arial", base_size = 9) + 
    theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5), strip.text = element_text(size = 8), 
        legend.position = "none")

DIAGNOSTIC_PDF_FILE <- file.path(OUT_DIR, paste0("05_68_GSE195783_spatial_section_map_diagnostic_", VERSION, 
    ".pdf"))

DIAGNOSTIC_PNG_FILE <- file.path(OUT_DIR, paste0("05_68_GSE195783_spatial_section_map_diagnostic_", VERSION, 
    ".png"))

ggsave(filename = DIAGNOSTIC_PDF_FILE, plot = p_diag, width = 11.5, height = 6.8, units = "in", device = grDevices::cairo_pdf, 
    bg = "white")

ggsave(filename = DIAGNOSTIC_PNG_FILE, plot = p_diag, width = 11.5, height = 6.8, units = "in", dpi = 300, 
    bg = "white")

n_final_sections <- nrow(section_manifest)

sections_per_capture <- section_manifest %>% count(sample_label, name = "n_sections")

unstable_samples <- stability_summary$sample_label[!stability_summary$all_thresholds_support_four]

final_qc <- data.frame(check = c("Eight expected Visium capture areas are present", "Primary audit contains four major tissue pieces per capture area", 
    "Final section manifest contains 32 major tissue sections", "Each capture area contributes exactly four finalized major sections", 
    "All finalized positions pass deterministic geometry QC", "Every major-component spot has a finalized internal section UID", 
    "Tiny disconnected components remain unassigned", "GSM + array_x + array_y join key is unique", "No animal identity is inferred", 
    "No Step60-66 object was modified"), pass = c(setequal(unique(section_manifest$sample_label), EXPECTED_CAPTURE_AREAS), 
    all(sample_summary$n_major_components == EXPECTED_MAJOR_SECTIONS_PER_CAPTURE), n_final_sections == 
        EXPECTED_TOTAL_MAJOR_SECTIONS, all(sections_per_capture$n_sections == EXPECTED_MAJOR_SECTIONS_PER_CAPTURE), 
    all(geometry_qc$all_position_checks_pass), all(!is.na(spot_section_map$section_uid_internal[spot_section_map$is_major_component])), 
    all(is.na(spot_section_map$section_uid_internal[!spot_section_map$is_major_component])), anyDuplicated(spot_section_map$section_join_key) == 
        0L, TRUE, TRUE), stringsAsFactors = FALSE)

FINAL_QC_FILE <- file.path(OUT_DIR, paste0("06_68_GSE195783_spatial_section_map_FINAL_QC_", VERSION, 
    ".csv"))

write_csv_utf8(final_qc, FINAL_QC_FILE)

REPORT_FILE <- file.path(OUT_DIR, paste0("07_68_GSE195783_spatial_section_map_REPORT_", VERSION, ".txt"))

if (length(unstable_samples) == 0L) {
    stability_line <- "  All capture areas retained four major components across all tested thresholds."
} else {
    stability_line <- paste0("  Threshold sensitivity was not perfectly stable for: ", paste(unstable_samples, 
        collapse = " | "), ". Primary audit still identified four major components in every capture area.")
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(section_manifest[, c("sample_label", "timepoint", "section_position_internal", "section_uid_internal", 
    "n_spots", "centroid_x", "centroid_y")])

print(capture_area_summary)

print(final_qc)

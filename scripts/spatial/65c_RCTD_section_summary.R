rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)


options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(1234)

ROOT <- PROJECT_ROOT

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

STEP68_DIR <- file.path(ROOT, "results", "68_GSE195783_spatial_section_map_finalize")

INPUT_SPOT_WEIGHTS <- file.path(STEP65_DIR, "03_65_RCTD_normalized_spot_weights_v04.csv.gz")

INPUT_SECTION_MANIFEST <- file.path(STEP68_DIR, "01_68_GSE195783_spatial_section_map_section_manifest_v01.csv")

INPUT_SPOT_SECTION_MAP <- file.path(STEP68_DIR, "02_68_GSE195783_spatial_section_map_spot_to_section_v01.csv.gz")

VERSION <- "v06"

TARGET_STATE <- "Clec7a+ lipid-handling Microglia"

EXPECTED_STEP65_SPOTS <- 15204L
EXPECTED_CAPTURE_AREAS <- 8L
EXPECTED_MAJOR_SECTIONS <- 32L
EXPECTED_REFERENCE_TYPES <- 19L

SUM_TOL <- 1e-8

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

SAMPLE_ORDER <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1", 
    "28dpi_rep2")

required_pkgs <- c("dplyr", "tidyr", "ggplot2")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step65 v06.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

write_csv_utf8 <- write_csv_utf8

read_csv_gz_base <- read_gzip_csv

write_csv_gz_base <- write_gzip_csv

safe_q <- function(x, p) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

as_logical_safe <- function(x) {
    if (is.logical(x)) {
        return(x)
    }
    y <- toupper(trimws(as.character(x)))
    out <- rep(NA, length(y))
    out[y %in% c("TRUE", "T", "1")] <- TRUE
    out[y %in% c("FALSE", "F", "0")] <- FALSE
    out
}

required_files <- c(INPUT_SPOT_WEIGHTS, INPUT_SECTION_MANIFEST, INPUT_SPOT_SECTION_MAP)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop("Required input file(s) missing:\n", paste(missing_files, collapse = "\n"), "\nRCTD input check failed; inspect this error before continuing.")
}

spot_weights_v04 <- read_csv_gz_base(INPUT_SPOT_WEIGHTS)

required_v04_cols <- c("barcode", "GSM", "sample_label", "timepoint", "spatial_region_v03", "array_x", 
    "array_y", TARGET_STATE)

missing_v04_cols <- setdiff(required_v04_cols, colnames(spot_weights_v04))

if (length(missing_v04_cols) > 0L) {
    stop("Step65 v04 normalized spot-weight file is missing required column(s): ", paste(missing_v04_cols, 
        collapse = " | "))
}

if (nrow(spot_weights_v04) != EXPECTED_STEP65_SPOTS) {
    stop("Expected ", EXPECTED_STEP65_SPOTS, " Step65 v04 spots, but found ", nrow(spot_weights_v04), 
        ".")
}

if (!setequal(unique(as.character(spot_weights_v04$sample_label)), SAMPLE_ORDER)) {
    stop("Step65 v04 sample labels do not match the expected 8 capture areas.")
}

if (any(!is.finite(spot_weights_v04$array_x)) || any(!is.finite(spot_weights_v04$array_y))) {
    stop("Step65 v04 contains non-finite array coordinates.")
}

known_metadata_cols <- c("barcode", "GSM", "sample_label", "timepoint", "spatial_region_v03", "array_x", 
    "array_y")

weight_cols <- setdiff(colnames(spot_weights_v04), known_metadata_cols)

if (length(weight_cols) != EXPECTED_REFERENCE_TYPES) {
    stop("Expected ", EXPECTED_REFERENCE_TYPES, " normalized RCTD reference-type columns, but found ", 
        length(weight_cols), ".\nDetected columns:\n", paste(weight_cols, collapse = "\n"))
}

if (!TARGET_STATE %in% weight_cols) {
    stop("Target state is not among the normalized RCTD weight columns:\n", TARGET_STATE)
}

weight_mat <- as.matrix(spot_weights_v04[, weight_cols, drop = FALSE])

storage.mode(weight_mat) <- "double"

if (any(!is.finite(weight_mat))) {
    stop("At least one normalized RCTD weight is non-finite.")
}

if (any(weight_mat < -1e-10 | weight_mat > 1 + 1e-10)) {
    stop("At least one normalized RCTD weight is outside [0,1].")
}

weight_sums <- rowSums(weight_mat)

max_abs_weight_sum_error <- max(abs(weight_sums - 1))

if (!is.finite(max_abs_weight_sum_error) || max_abs_weight_sum_error > SUM_TOL) {
    stop("Step65 v04 normalized RCTD weights no longer sum to 1.\n", "Max absolute error = ", signif(max_abs_weight_sum_error, 
        6))
}

rm(weight_mat)

gc(verbose = FALSE)

section_manifest_v01 <- read.csv(INPUT_SECTION_MANIFEST, stringsAsFactors = FALSE, check.names = FALSE)

spot_section_map_v01 <- read_csv_gz_base(INPUT_SPOT_SECTION_MAP)

required_manifest_cols <- c("sample_label", "GSM", "timepoint", "section_position_internal", "section_uid_internal", 
    "n_spots", "centroid_x", "centroid_y")

missing_manifest_cols <- setdiff(required_manifest_cols, colnames(section_manifest_v01))

if (length(missing_manifest_cols) > 0L) {
    stop("Step68 section manifest is missing required column(s): ", paste(missing_manifest_cols, collapse = " | "))
}

if (nrow(section_manifest_v01) != EXPECTED_MAJOR_SECTIONS) {
    stop("Step68 section manifest must contain exactly ", EXPECTED_MAJOR_SECTIONS, " major coordinate-defined sections.")
}

required_section_map_cols <- c("GSM", "sample_label", "timepoint", "array_x", "array_y", "is_major_component", 
    "section_position_internal", "section_uid_internal", "section_status")

missing_section_map_cols <- setdiff(required_section_map_cols, colnames(spot_section_map_v01))

if (length(missing_section_map_cols) > 0L) {
    stop("Step68 spot-to-section map is missing required column(s): ", paste(missing_section_map_cols, 
        collapse = " | "))
}

spot_section_map_v01$is_major_component <- as_logical_safe(spot_section_map_v01$is_major_component)

for (cc in c("section_position_internal", "section_uid_internal")) {
    spot_section_map_v01[[cc]] <- trimws(as.character(spot_section_map_v01[[cc]]))
    spot_section_map_v01[[cc]][!nzchar(spot_section_map_v01[[cc]])] <- NA_character_
}

if (any(is.na(spot_section_map_v01$is_major_component))) {
    stop("Could not safely parse is_major_component in Step68 map.")
}

section_key_map <- spot_section_map_v01 %>% transmute(GSM = as.character(GSM), array_x = as.numeric(array_x), 
    array_y = as.numeric(array_y), map_sample_label = as.character(sample_label), map_timepoint = as.character(timepoint), 
    is_major_component = is_major_component, section_position_internal = as.character(section_position_internal), 
    section_uid_internal = as.character(section_uid_internal), section_status = as.character(section_status))

if (anyDuplicated(section_key_map[, c("GSM", "array_x", "array_y")]) > 0L) {
    stop("Step68 map has duplicated GSM + array_x + array_y keys.")
}

spot_section_aware <- spot_weights_v04 %>% mutate(GSM = as.character(GSM), sample_label = as.character(sample_label), 
    timepoint = as.character(timepoint), array_x = as.numeric(array_x), array_y = as.numeric(array_y)) %>% 
    left_join(section_key_map, by = c("GSM", "array_x", "array_y"))

join_missing <- is.na(spot_section_aware$map_sample_label)

if (any(join_missing)) {
    stop("Section join failed for ", sum(join_missing), " Step65 v04 spots.")
}

if (any(spot_section_aware$sample_label != spot_section_aware$map_sample_label)) {
    stop("Capture-area label mismatch after Step68 join.")
}

if (any(spot_section_aware$timepoint != spot_section_aware$map_timepoint)) {
    stop("Timepoint mismatch after Step68 join.")
}

spot_section_aware$section_eligible <- spot_section_aware$is_major_component & !is.na(spot_section_aware$section_uid_internal) & 
    nzchar(spot_section_aware$section_uid_internal)

if (any(spot_section_aware$is_major_component & !spot_section_aware$section_eligible)) {
    stop("At least one Step65 spot in a major component lacks a finalized section UID.")
}

spot_section_aware$map_sample_label <- NULL
spot_section_aware$map_timepoint <- NULL

section_spot_summary <- spot_section_aware %>% filter(section_eligible) %>% group_by(sample_label, GSM, 
    timepoint, section_position_internal, section_uid_internal) %>% summarise(n_step65_spots = n(), mean_array_x = mean(array_x), 
    mean_array_y = mean(array_y), .groups = "drop")

if (nrow(section_spot_summary) != EXPECTED_MAJOR_SECTIONS) {
    stop("Expected ", EXPECTED_MAJOR_SECTIONS, " section-level groups among Step65-retained spots, but found ", 
        nrow(section_spot_summary), ".")
}

section_manifest_v05 <- section_manifest_v01 %>% select(sample_label, GSM, timepoint, section_position_internal, 
    section_uid_internal, n_spots_step68 = n_spots, centroid_x_step68 = centroid_x, centroid_y_step68 = centroid_y) %>% 
    left_join(section_spot_summary, by = c("sample_label", "GSM", "timepoint", "section_position_internal", 
        "section_uid_internal")) %>% mutate(n_removed_before_step65 = n_spots_step68 - n_step65_spots, 
    fraction_step68_spots_retained = n_step65_spots/n_spots_step68)

if (any(!is.finite(section_manifest_v05$n_step65_spots))) {
    stop("At least one finalized section has no Step65-retained spots.")
}

n_joined <- nrow(spot_section_aware)

n_section_eligible <- sum(spot_section_aware$section_eligible)

n_unassigned_tiny <- sum(!spot_section_aware$section_eligible)

n_capture_areas <- length(unique(spot_section_aware$sample_label))

n_major_sections <- length(unique(spot_section_aware$section_uid_internal[spot_section_aware$section_eligible]))

join_qc <- data.frame(metric = c("Step65_v04_input_spots", "Step68_joined_spots", "join_fraction", "section_eligible_major_spots", 
    "unassigned_tiny_component_spots", "fraction_section_eligible", "capture_areas", "major_coordinate_defined_sections", 
    "normalized_reference_types", "max_abs_RCTD_weight_sum_error"), value = c(nrow(spot_weights_v04), 
    n_joined, n_joined/nrow(spot_weights_v04), n_section_eligible, n_unassigned_tiny, n_section_eligible/n_joined, 
    n_capture_areas, n_major_sections, length(weight_cols), max_abs_weight_sum_error), stringsAsFactors = FALSE)

SPOT_OUTPUT_FILE <- file.path(STEP65_DIR, paste0("03_65_RCTD_normalized_spot_weights_section_aware_", 
    VERSION, ".csv.gz"))

write_csv_gz_base(spot_section_aware, SPOT_OUTPUT_FILE)

section_target_summary <- spot_section_aware %>% filter(section_eligible) %>% group_by(timepoint, sample_label, 
    GSM, section_position_internal, section_uid_internal) %>% summarise(n_spots = n(), mean_target_weight = mean(.data[[TARGET_STATE]], 
    na.rm = TRUE), median_target_weight = median(.data[[TARGET_STATE]], na.rm = TRUE), q25_target_weight = safe_q(.data[[TARGET_STATE]], 
    0.25), q75_target_weight = safe_q(.data[[TARGET_STATE]], 0.75), fraction_weight_gt_0_10 = mean(.data[[TARGET_STATE]] > 
    0.1, na.rm = TRUE), .groups = "drop")

capture_area_section_summary <- section_target_summary %>% group_by(timepoint, sample_label, GSM) %>% 
    summarise(n_sections = n(), mean_of_section_means = mean(mean_target_weight), median_of_section_means = median(mean_target_weight), 
        min_section_mean = min(mean_target_weight), max_section_mean = max(mean_target_weight), .groups = "drop")

capture_area_pooled_summary <- spot_section_aware %>% group_by(timepoint, sample_label, GSM) %>% summarise(pooled_all_step65_spots_mean = mean(.data[[TARGET_STATE]], 
    na.rm = TRUE), pooled_major_section_spots_mean = mean(.data[[TARGET_STATE]][section_eligible], na.rm = TRUE), 
    n_all_step65_spots = n(), n_major_section_spots = sum(section_eligible), .groups = "drop")

capture_area_sensitivity <- capture_area_section_summary %>% left_join(capture_area_pooled_summary, by = c("timepoint", 
    "sample_label", "GSM")) %>% mutate(equal_section_minus_pooled_all = mean_of_section_means - pooled_all_step65_spots_mean, 
    equal_section_minus_pooled_major = mean_of_section_means - pooled_major_section_spots_mean)

timepoint_section_summary <- section_target_summary %>% group_by(timepoint) %>% summarise(n_sections = n(), 
    n_capture_areas = n_distinct(sample_label), mean_of_section_means = mean(mean_target_weight), median_of_section_means = median(mean_target_weight), 
    sd_of_section_means = stats::sd(mean_target_weight), min_section_mean = min(mean_target_weight), 
    max_section_mean = max(mean_target_weight), .groups = "drop")

section_target_summary$timepoint <- factor(section_target_summary$timepoint, levels = TIME_LEVELS)

section_target_summary <- section_target_summary %>% arrange(timepoint, sample_label, section_position_internal)

section_target_summary$timepoint <- as.character(section_target_summary$timepoint)

timepoint_section_summary$timepoint <- factor(timepoint_section_summary$timepoint, levels = TIME_LEVELS)

timepoint_section_summary <- timepoint_section_summary %>% arrange(timepoint)

timepoint_section_summary$timepoint <- as.character(timepoint_section_summary$timepoint)

section_target_output <- bind_rows(section_target_summary %>% transmute(summary_level = "section", timepoint, 
    sample_label, GSM, section_position_internal, section_uid_internal, n_spots, n_sections = NA_integer_, 
    n_capture_areas = NA_integer_, mean_target_weight, median_target_weight, q25_target_weight, q75_target_weight, 
    fraction_weight_gt_0_10, sd_of_section_means = NA_real_, min_section_mean = NA_real_, max_section_mean = NA_real_), 
    timepoint_section_summary %>% transmute(summary_level = "timepoint_across_coordinate_defined_sections", 
        timepoint, sample_label = NA_character_, GSM = NA_character_, section_position_internal = NA_character_, 
        section_uid_internal = NA_character_, n_spots = NA_integer_, n_sections, n_capture_areas, mean_target_weight = mean_of_section_means, 
        median_target_weight = median_of_section_means, q25_target_weight = NA_real_, q75_target_weight = NA_real_, 
        fraction_weight_gt_0_10 = NA_real_, sd_of_section_means, min_section_mean, max_section_mean))

capture_area_sensitivity$timepoint <- factor(capture_area_sensitivity$timepoint, levels = TIME_LEVELS)

capture_area_sensitivity <- capture_area_sensitivity %>% arrange(timepoint, sample_label)

capture_area_sensitivity$timepoint <- as.character(capture_area_sensitivity$timepoint)

section_region_summary <- spot_section_aware %>% filter(section_eligible, !is.na(spatial_region_v03), 
    nzchar(as.character(spatial_region_v03))) %>% group_by(timepoint, sample_label, GSM, section_position_internal, 
    section_uid_internal, spatial_region_v03) %>% summarise(n_spots = n(), mean_target_weight = mean(.data[[TARGET_STATE]], 
    na.rm = TRUE), median_target_weight = median(.data[[TARGET_STATE]], na.rm = TRUE), fraction_weight_gt_0_10 = mean(.data[[TARGET_STATE]] > 
    0.1, na.rm = TRUE), .groups = "drop")

timepoint_region_summary <- section_region_summary %>% group_by(timepoint, spatial_region_v03) %>% summarise(n_sections_with_region = n(), 
    n_capture_areas_with_region = n_distinct(sample_label), total_region_spots = sum(n_spots), mean_of_section_region_means = mean(mean_target_weight), 
    median_of_section_region_means = median(mean_target_weight), sd_of_section_region_means = if (n() > 
        1L) {
        stats::sd(mean_target_weight)
    } else {
        NA_real_
    }, .groups = "drop")

section_region_output <- bind_rows(section_region_summary %>% transmute(summary_level = "section_region", 
    timepoint, sample_label, GSM, section_position_internal, section_uid_internal, spatial_region_v03, 
    n_spots, n_sections_with_region = NA_integer_, n_capture_areas_with_region = NA_integer_, total_region_spots = NA_integer_, 
    mean_target_weight, median_target_weight, fraction_weight_gt_0_10, sd_of_section_region_means = NA_real_), 
    timepoint_region_summary %>% transmute(summary_level = "timepoint_region_across_coordinate_defined_sections", 
        timepoint, sample_label = NA_character_, GSM = NA_character_, section_position_internal = NA_character_, 
        section_uid_internal = NA_character_, spatial_region_v03, n_spots = NA_integer_, n_sections_with_region, 
        n_capture_areas_with_region, total_region_spots, mean_target_weight = mean_of_section_region_means, 
        median_target_weight = median_of_section_region_means, fraction_weight_gt_0_10 = NA_real_, sd_of_section_region_means))

JOIN_QC_FILE <- file.path(STEP65_DIR, paste0("01_65_RCTD_section_join_QC_", VERSION, ".csv"))

SECTION_MANIFEST_FILE <- file.path(STEP65_DIR, paste0("02_65_RCTD_section_manifest_", VERSION, ".csv"))

SECTION_SUMMARY_FILE <- file.path(STEP65_DIR, paste0("04_65_Clec7a_section_summary_", VERSION, ".csv"))

SECTION_REGION_FILE <- file.path(STEP65_DIR, paste0("05_65_Clec7a_section_region_summary_", VERSION, 
    ".csv"))

CAPTURE_SENSITIVITY_FILE <- file.path(STEP65_DIR, paste0("07_65_Clec7a_capture_area_sensitivity_", VERSION, 
    ".csv"))

write_csv_utf8(join_qc, JOIN_QC_FILE)

write_csv_utf8(section_manifest_v05, SECTION_MANIFEST_FILE)

write_csv_utf8(section_target_output, SECTION_SUMMARY_FILE)

write_csv_utf8(section_region_output, SECTION_REGION_FILE)

write_csv_utf8(capture_area_sensitivity, CAPTURE_SENSITIVITY_FILE)

plot_df <- section_target_summary

plot_df$timepoint <- factor(plot_df$timepoint, levels = TIME_LEVELS)

plot_df$sample_label <- factor(plot_df$sample_label, levels = SAMPLE_ORDER)

p_section <- ggplot(plot_df, aes(x = timepoint, y = mean_target_weight)) + geom_boxplot(width = 0.52, 
    outlier.shape = NA, linewidth = 0.45) + geom_point(aes(shape = sample_label), position = position_jitter(width = 0.11, 
    height = 0), size = 2, alpha = 0.9) + scale_shape_manual(values = c(`3dpi_rep1` = 16, `3dpi_rep2` = 17, 
    `7dpi_rep1` = 15, `7dpi_rep2` = 3, `14dpi_rep1` = 0, `14dpi_rep2` = 8, `28dpi_rep1` = 1, `28dpi_rep2` = 4), 
    drop = FALSE) + labs(title = "Clec7a+ lipid-handling Microglia", subtitle = paste0("Section-aware normalized RCTD abundance; descriptive only; ", 
    "no animal-level independence is assumed"), x = NULL, y = "Mean normalized RCTD weight", shape = "Visium capture area") + 
    theme_classic(base_family = "Arial", base_size = 9) + theme(plot.title = element_text(size = 11), 
    plot.subtitle = element_text(size = 8), legend.position = "right", legend.text = element_text(size = 7), 
    legend.title = element_text(size = 7.5))

SECTION_PLOT_PDF <- file.path(STEP65_DIR, paste0("06_65_Clec7a_section_abundance_", VERSION, ".pdf"))

SECTION_PLOT_PNG <- file.path(STEP65_DIR, paste0("06_65_Clec7a_section_abundance_", VERSION, ".png"))

ggsave(filename = SECTION_PLOT_PDF, plot = p_section, width = 6.3, height = 4.4, units = "in", device = grDevices::cairo_pdf, 
    bg = "white")

ggsave(filename = SECTION_PLOT_PNG, plot = p_section, width = 6.3, height = 4.4, units = "in", dpi = 300, 
    bg = "white")

all_32_sections_present <- n_major_sections == EXPECTED_MAJOR_SECTIONS

all_8_capture_areas_present <- n_capture_areas == EXPECTED_CAPTURE_AREAS

four_sections_per_capture <- section_target_summary %>% count(sample_label, name = "n_sections") %>% 
    pull(n_sections)

four_sections_per_capture_ok <- all(four_sections_per_capture == 4L)

final_checks <- data.frame(check = c("Step65 v04 normalized spot-weight file has expected 15,204 spots", 
    "All 19 normalized RCTD reference types are present", "Normalized RCTD weights still sum to 1 per spot", 
    "All Step65 v04 spots match Step68 by GSM + array_x + array_y", "All 8 Visium capture areas are retained", 
    "All 32 coordinate-defined major tissue sections are represented", "Each capture area contributes exactly four coordinate-defined sections", 
    "Tiny disconnected components remain excluded from section summaries", "No RCTD computation was rerun", 
    "No RCTD weight renormalization was performed", "No animal-level independence is assumed"), pass = c(nrow(spot_weights_v04) == 
    EXPECTED_STEP65_SPOTS, length(weight_cols) == EXPECTED_REFERENCE_TYPES, max_abs_weight_sum_error <= 
    SUM_TOL, n_joined == nrow(spot_weights_v04), all_8_capture_areas_present, all_32_sections_present, 
    four_sections_per_capture_ok, all(is.na(spot_section_aware$section_uid_internal[!spot_section_aware$section_eligible]) | 
        !nzchar(spot_section_aware$section_uid_internal[!spot_section_aware$section_eligible])), TRUE, 
    TRUE, TRUE), stringsAsFactors = FALSE)

FINAL_QC_FILE <- file.path(STEP65_DIR, paste0("08_65_RCTD_section_finalize_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_checks, FINAL_QC_FILE)

REPORT_FILE <- file.path(STEP65_DIR, paste0("09_65_GSE195783_RCTD_section_finalize_REPORT_", VERSION, 
    ".txt"))

timepoint_lines <- apply(timepoint_section_summary[, c("timepoint", "n_sections", "n_capture_areas", 
    "mean_of_section_means", "median_of_section_means", "sd_of_section_means", "min_section_mean", "max_section_mean")], 
    1, function(x) {
        paste0("  ", x["timepoint"], ": n_sections=", x["n_sections"], ", n_capture_areas=", x["n_capture_areas"], 
            ", mean=", signif(as.numeric(x["mean_of_section_means"]), 5), ", median=", signif(as.numeric(x["median_of_section_means"]), 
                5), ", SD=", signif(as.numeric(x["sd_of_section_means"]), 5), ", range=", signif(as.numeric(x["min_section_mean"]), 
                5), " to ", signif(as.numeric(x["max_section_mean"]), 5))
    })

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(timepoint_section_summary)

print(capture_area_sensitivity)

print(final_checks)

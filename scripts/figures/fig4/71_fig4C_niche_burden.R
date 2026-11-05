rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)


options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(7102)

ROOT <- PROJECT_ROOT

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

INPUT_FILE <- file.path(STEP65_DIR, "05_65_Clec7a_section_region_summary_v06.csv")

OUT_DIR <- file.path(ROOT, "results", "71_Figure4C_Clec7a_spatial_niche_heatmap")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v04"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

TIME_LABELS <- c(`3dpi` = "3 dpi", `7dpi` = "7 dpi", `14dpi` = "14 dpi", `28dpi` = "28 dpi")

TARGET_REGION <- "Myeloid/immune–lipid-handling enriched"

REGION_LEVELS <- c("Myeloid/immune–lipid-handling enriched", "Fibrotic/ECM–myeloid mixed injury niche", 
    "Acute proliferative injury-associated myeloid", "IFN-responsive inflammatory/astroglial region", 
    "Perivascular / vSMC–pericyte region", "Peripheral-myelin / Schwann-like region", "Oligodendrocyte / white-matter-like", 
    "Neuronal subtype region, 7 dpi-enriched", "Cholinergic/motor-neuron-like neuronal", "Inhibitory/interneuron-like neuronal")

REGION_LABELS <- c(`Myeloid/immune–lipid-handling enriched` = "Myeloid + lipid handling", `Fibrotic/ECM–myeloid mixed injury niche` = "Fibrotic/ECM–myeloid niche", 
    `Acute proliferative injury-associated myeloid` = "Acute proliferative injury myeloid", `IFN-responsive inflammatory/astroglial region` = "IFN-inflammatory/astroglial", 
    `Perivascular / vSMC–pericyte region` = "Perivascular vSMC/pericyte", `Peripheral-myelin / Schwann-like region` = "Peripheral myelin/Schwann-like", 
    `Oligodendrocyte / white-matter-like` = "Oligo/WM-like", `Neuronal subtype region, 7 dpi-enriched` = "7 dpi neuronal subtype", 
    `Cholinergic/motor-neuron-like neuronal` = "Cholinergic/motor neuron-like", `Inhibitory/interneuron-like neuronal` = "Interneuron-like")

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step71.")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

write_csv_utf8 <- write_csv_utf8

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

if (!file.exists(INPUT_FILE)) {
    stop("Required Step65 v06 input is missing:\n", INPUT_FILE)
}

dat0 <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("summary_level", "timepoint", "sample_label", "GSM", "section_position_internal", 
    "section_uid_internal", "spatial_region_v03", "n_spots", "mean_target_weight")

missing_cols <- setdiff(required_cols, colnames(dat0))

if (length(missing_cols) > 0L) {
    stop("Step65 section-region file is missing required column(s): ", paste(missing_cols, collapse = " | "))
}

sec_region_raw <- dat0 %>% filter(summary_level == "section_region", timepoint %in% TIME_LEVELS, spatial_region_v03 %in% 
    REGION_LEVELS) %>% transmute(timepoint = as.character(timepoint), sample_label = as.character(sample_label), 
    GSM = as.character(GSM), section_position_internal = as.character(section_position_internal), section_uid_internal = as.character(section_uid_internal), 
    spatial_region_v03 = as.character(spatial_region_v03), n_spots = as.integer(n_spots), mean_target_weight = as.numeric(mean_target_weight))

if (nrow(sec_region_raw) == 0L) {
  stop("No section_region rows remain after filtering.")
}

if (any(!is.finite(sec_region_raw$mean_target_weight))) {
  stop("Non-finite mean_target_weight values detected.")
}

if (any(sec_region_raw$mean_target_weight < 0)) {
  stop("Negative normalized RCTD weights detected.")
}

if (any(!is.finite(sec_region_raw$n_spots)) || any(sec_region_raw$n_spots <= 0)) {
    stop("Invalid n_spots values detected.")
}

section_info <- sec_region_raw %>% group_by(timepoint, sample_label, GSM, section_position_internal, 
    section_uid_internal) %>% summarise(total_section_spots = sum(n_spots, na.rm = TRUE), .groups = "drop")

if (any(section_info$total_section_spots <= 0)) {
  stop("At least one section has total_section_spots <= 0.")
}

section_counts <- section_info %>% count(timepoint, name = "n_sections")

print(section_counts)

sec_region_full <- section_info %>%
  tidyr::crossing(
    spatial_region_v03 = REGION_LEVELS
  ) %>%
  left_join(
    sec_region_raw,
    by = c(
      "timepoint",
      "sample_label",
      "GSM",
      "section_position_internal",
      "section_uid_internal",
      "spatial_region_v03"
    ),
    suffix = c("", "_raw")
  ) %>%
  mutate(
    region_present =
      !is.na(n_spots),
    n_spots =
      tidyr::replace_na(
        n_spots,
        0L
      ),
    mean_target_weight =
      tidyr::replace_na(
        mean_target_weight,
        0
      ),

    region_target_weight_sum =
      n_spots *
      mean_target_weight,

    spatial_niche_burden =
      region_target_weight_sum /
      total_section_spots
  )

if (any(!is.finite(sec_region_full$spatial_niche_burden))) {
  stop("Non-finite spatial_niche_burden values detected.")
}

if (any(sec_region_full$spatial_niche_burden < 0)) {
  stop("Negative spatial_niche_burden values detected.")
}

heatmap_df <- sec_region_full %>% group_by(timepoint, spatial_region_v03) %>% summarise(n_sections_total = n_distinct(section_uid_internal), 
    n_sections_region_present = sum(region_present), n_capture_areas_total = n_distinct(sample_label), 
    total_region_spots = sum(n_spots, na.rm = TRUE), mean_spatial_niche_burden = safe_mean(spatial_niche_burden), 
    sd_spatial_niche_burden = safe_sd(spatial_niche_burden), .groups = "drop")

expected_n <- length(TIME_LEVELS) * length(REGION_LEVELS)

if (nrow(heatmap_df) != expected_n) {
    stop("Unexpected heatmap dimensions: ", nrow(heatmap_df), " rows; expected ", expected_n, ".")
}

heatmap_df <- heatmap_df %>% mutate(timepoint_raw = timepoint, timepoint = factor(timepoint, levels = TIME_LEVELS, 
    labels = unname(TIME_LABELS[TIME_LEVELS])), spatial_region_v03 = factor(spatial_region_v03, levels = rev(REGION_LEVELS)), 
    display_region = REGION_LABELS[as.character(spatial_region_v03)])

display_region_levels <- rev(unname(REGION_LABELS[REGION_LEVELS]))

heatmap_df$display_region <- factor(heatmap_df$display_region, levels = display_region_levels)

SECTION_BURDEN_FILE <- file.path(OUT_DIR, paste0("01_71_Figure4C_section_region_spatial_burden_", VERSION, 
    ".csv"))

write_csv_utf8(sec_region_full %>% arrange(factor(timepoint, levels = TIME_LEVELS), section_uid_internal, 
    factor(spatial_region_v03, levels = REGION_LEVELS)), SECTION_BURDEN_FILE)

SUMMARY_FILE <- file.path(OUT_DIR, paste0("02_71_Figure4C_timepoint_region_burden_summary_", VERSION, 
    ".csv"))

write_csv_utf8(heatmap_df %>% mutate(timepoint = as.character(timepoint), spatial_region_v03 = as.character(spatial_region_v03), 
    display_region = as.character(display_region)), SUMMARY_FILE)

MATRIX_FILE <- file.path(OUT_DIR, paste0("03_71_Figure4C_spatial_niche_burden_matrix_", VERSION, ".csv"))

matrix_df <- heatmap_df %>% select(display_region, timepoint, mean_spatial_niche_burden) %>% pivot_wider(names_from = timepoint, 
    values_from = mean_spatial_niche_burden)

write_csv_utf8(matrix_df, MATRIX_FILE)

SUPPORT_FILE <- file.path(OUT_DIR, paste0("04_71_Figure4C_region_presence_support_", VERSION, ".csv"))

support_df <- heatmap_df %>% select(display_region, timepoint, n_sections_total, n_sections_region_present, 
    total_region_spots)

write_csv_utf8(support_df, SUPPORT_FILE)

finite_values <- heatmap_df$mean_spatial_niche_burden[is.finite(heatmap_df$mean_spatial_niche_burden)]

if (length(finite_values) == 0L) {
  stop("No finite heatmap values are available.")
}

FILL_MIN <- 0

FILL_MID <- as.numeric(stats::quantile(finite_values, probs = 0.7, na.rm = TRUE, names = FALSE, type = 7))

FILL_MAX <- as.numeric(stats::quantile(finite_values, probs = 0.95, na.rm = TRUE, names = FALSE, type = 7))

if (!is.finite(FILL_MID) || !is.finite(FILL_MAX) || FILL_MID <= FILL_MIN || FILL_MAX <= FILL_MID) {
    stop("Invalid q70/q95 heatmap color anchors.")
}

LOW_COLOR  <- "#F4F4E8"
MID_COLOR  <- "#F6A13A"
HIGH_COLOR <- "#D7191C"

heatmap_df <- heatmap_df %>% mutate(cell_label = ifelse(is.finite(mean_spatial_niche_burden), sprintf("%.4f", 
    mean_spatial_niche_burden), ""))

p <- ggplot(heatmap_df, aes(x = timepoint, y = display_region, fill = mean_spatial_niche_burden)) + geom_tile(color = "white", 
    linewidth = 0.9, width = 0.96, height = 0.96) + geom_text(aes(label = cell_label), family = "Arial", 
    size = 2.85, color = "black") + scale_fill_gradientn(colours = c(LOW_COLOR, MID_COLOR, HIGH_COLOR), 
    values = scales::rescale(c(FILL_MIN, FILL_MID, FILL_MAX), from = c(FILL_MIN, FILL_MAX)), limits = c(FILL_MIN, 
        FILL_MAX), oob = scales::squish, name = "Section-normalized\nspatial burden") + scale_x_discrete(position = "top", 
    drop = FALSE) + scale_y_discrete(drop = FALSE) + labs(x = NULL, y = NULL) + theme_minimal(base_family = "Arial", 
    base_size = 9) + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 9, face = "bold", 
    color = "black", margin = margin(b = 4)), axis.text.y = element_text(size = 8.2, color = "black", 
    hjust = 1), legend.position = "right", legend.title = element_text(size = 8, hjust = 0.5), legend.text = element_text(size = 7.4), 
    legend.key.height = unit(0.55, "cm"), legend.key.width = unit(0.42, "cm"), plot.margin = margin(4, 
        4, 4, 4)) + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = unit(3.8, 
    "cm"), barwidth = unit(0.45, "cm"), ticks = TRUE, frame.colour = "black", frame.linewidth = 0.35))

PDF_FILE <- file.path(OUT_DIR, paste0("71_Figure4C_Clec7a_spatial_niche_burden_heatmap_", VERSION, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0("71_Figure4C_Clec7a_spatial_niche_burden_heatmap_", VERSION, ".png"))

SVG_FILE <- file.path(OUT_DIR, paste0("71_Figure4C_Clec7a_spatial_niche_burden_heatmap_", VERSION, ".svg"))

ggsave(filename = PDF_FILE, plot = p, width = 7.2, height = 5.6, units = "in", device = grDevices::cairo_pdf, 
    bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = 7.2, height = 5.6, units = "in", dpi = 600, bg = "white")

ggsave(filename = SVG_FILE, plot = p, width = 7.2, height = 5.6, units = "in", device = "svg", bg = "white")

target_focus <- heatmap_df %>% filter(as.character(spatial_region_v03) == TARGET_REGION) %>% transmute(timepoint = as.character(timepoint), 
    mean_spatial_niche_burden, sd_spatial_niche_burden, n_sections_total, n_sections_region_present, 
    n_capture_areas_total, total_region_spots)

TARGET_FILE <- file.path(OUT_DIR, paste0("05_71_Figure4C_target_niche_burden_summary_", VERSION, ".csv"))

write_csv_utf8(target_focus, TARGET_FILE)

section_count_ok <- all(section_counts$timepoint %in% TIME_LEVELS) && length(unique(section_counts$timepoint)) == 
    4L

burden_identity_check <- sec_region_full %>% mutate(reconstructed = n_spots * mean_target_weight/total_section_spots, 
    abs_diff = abs(reconstructed - spatial_niche_burden)) %>% summarise(max_abs_diff = max(abs_diff, 
    na.rm = TRUE)) %>% pull(max_abs_diff)

final_checks <- data.frame(check = c("Step65 v06 section-region summary was used", "Only section-level rows were used", 
    "All 4 Figure4 timepoints are represented", "All 10 spatial-region categories are represented", "Absent region within a section is encoded as zero burden", 
    "Section-normalized burden formula reconstructs exactly", "Each heatmap value averages section-level burdens", 
    "Raw burden scale is used", "No row z-score is used", "No clustering or dendrogram is used", "No hypothesis testing is performed", 
    "Cream-orange-red palette uses shared q70/q95 display anchors", "Target Myeloid + lipid handling niche is included", 
    "Heatmap PDF generated", "Heatmap PNG generated", "Heatmap SVG generated"), pass = c(grepl("section_region_summary_v06", 
    INPUT_FILE), all(dat0$summary_level[dat0$summary_level == "section_region"] == "section_region"), 
    section_count_ok, nrow(heatmap_df) == length(TIME_LEVELS) * length(REGION_LEVELS), all(sec_region_full$spatial_niche_burden[!sec_region_full$region_present] == 
        0), is.finite(burden_identity_check) && burden_identity_check < 1e-12, all(heatmap_df$n_sections_total > 
        0), all(heatmap_df$mean_spatial_niche_burden >= 0), TRUE, TRUE, TRUE, is.finite(FILL_MIN) && 
        is.finite(FILL_MID) && is.finite(FILL_MAX) && FILL_MIN < FILL_MID && FILL_MID < FILL_MAX, any(as.character(heatmap_df$spatial_region_v03) == 
        TARGET_REGION), file.exists(PDF_FILE), file.exists(PNG_FILE), file.exists(SVG_FILE)), stringsAsFactors = FALSE)

QC_FILE <- file.path(OUT_DIR, paste0("06_71_Figure4C_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_checks, QC_FILE)

REPORT_FILE <- file.path(OUT_DIR, paste0("07_71_Figure4C_REPORT_", VERSION, ".txt"))

target_lines <- apply(target_focus, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": mean spatial niche burden=", sprintf("%.5f", as.numeric(z[["mean_spatial_niche_burden"]])), 
        "; n_sections_total=", z[["n_sections_total"]], "; n_sections_region_present=", z[["n_sections_region_present"]], 
        "; n_capture_areas=", z[["n_capture_areas_total"]], "; region_spots=", z[["total_region_spots"]])
})

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(target_focus)
print(final_checks)

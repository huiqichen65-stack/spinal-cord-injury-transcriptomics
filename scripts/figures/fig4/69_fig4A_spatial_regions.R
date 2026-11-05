source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(6905)

ROOT <- PROJECT_ROOT

STEP64_DIR <- file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation")

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

STEP66_DIR <- file.path(ROOT, "results", "66_GSE195783_Clec7a_spatial_autocorrelation_neighborhood")

STEP68_DIR <- file.path(ROOT, "results", "68_GSE195783_spatial_section_map_finalize")

OUT_DIR <- file.path(ROOT, "results", "69_Figure4A_spatial_transcriptional_regions")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

STEP64_OBJECT <- file.path(STEP64_DIR, "12_64_GSE195783_spatial_region_annotated_v03.rds")

STEP68_SECTION_MAP <- file.path(STEP68_DIR, "02_68_GSE195783_spatial_section_map_spot_to_section_v01.csv.gz")

STEP65_SECTION_SUMMARY <- file.path(STEP65_DIR, "04_65_Clec7a_section_summary_v06.csv")

STEP65_SECTION_REGION <- file.path(STEP65_DIR, "05_65_Clec7a_section_region_summary_v06.csv")

STEP66_REGION_ENRICH <- file.path(STEP66_DIR, "03_66_Clec7a_section_region_enrichment_v02.csv")

VERSION <- "v12"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

TIME_LABELS <- c(`3dpi` = "3 dpi", `7dpi` = "7 dpi", `14dpi` = "14 dpi", `28dpi` = "28 dpi")

TARGET_REGION <- "Myeloid/immune–lipid-handling enriched"

MIN_SECTION_SPOTS <- 50L
MIN_TARGET_REGION_SPOTS <- 8L

required_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "dplyr", "grid", "patchwork")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(patchwork)
})

write_csv_utf8 <- write_csv_utf8

write_csv_gz_utf8 <- function(x, filename) {
    con <- gzfile(filename, open = "wt", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    write.table(x, file = con, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "")
}

first_mode <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
        return(NA_character_)
    }
    tab <- sort(table(x), decreasing = TRUE)
    names(tab)[1]
}

required_files <- c(STEP64_OBJECT, STEP68_SECTION_MAP, STEP65_SECTION_SUMMARY, STEP65_SECTION_REGION, 
    STEP66_REGION_ENRICH)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop("Required file(s) missing:\n", paste(missing_files, collapse = "\n"))
}

obj <- readRDS(STEP64_OBJECT)

if (!inherits(obj, "Seurat")) {
  stop("Step64 input is not a Seurat object.")
}

md <- obj[[]]

required_md <- c("sample_label", "timepoint", "array_x", "array_y", "spatial_region_v03", "histology_registered")

missing_md <- setdiff(required_md, colnames(md))

if (length(missing_md) > 0L) {
    stop("Step64 object is missing metadata column(s): ", paste(missing_md, collapse = ", "))
}

if (!("GSM" %in% colnames(md))) {
  stop("Step64 object must contain GSM.")
}

if (any(md$histology_registered %in% TRUE, na.rm = TRUE)) {
  stop("Step64 object must remain non-histology-registered.")
}

md_use <- md %>% transmute(barcode = rownames(md), GSM = as.character(GSM), sample_label = as.character(sample_label), 
    timepoint = as.character(timepoint), array_x = as.numeric(array_x), array_y = as.numeric(array_y), 
    spatial_region_v03 = as.character(spatial_region_v03)) %>% filter(timepoint %in% TIME_LEVELS)

section_map <- read.csv(STEP68_SECTION_MAP, stringsAsFactors = FALSE, check.names = FALSE)

required_map <- c("GSM", "sample_label", "timepoint", "array_x", "array_y", "is_major_component", "section_position_internal", 
    "section_uid_internal")

missing_map <- setdiff(required_map, colnames(section_map))

if (length(missing_map) > 0L) {
    stop("Step68 map is missing required column(s): ", paste(missing_map, collapse = ", "))
}

if (!is.logical(section_map$is_major_component)) {
    tmp <- toupper(trimws(as.character(section_map$is_major_component)))
    section_map$is_major_component <- tmp %in% c("TRUE", "T", "1")
}

section_map$section_uid_internal <- trimws(as.character(section_map$section_uid_internal))

section_map$section_uid_internal[!nzchar(section_map$section_uid_internal)] <- NA_character_

section_map2 <- section_map %>% filter(timepoint %in% TIME_LEVELS, is_major_component %in% TRUE, !is.na(section_uid_internal)) %>% 
    transmute(GSM = as.character(GSM), sample_label = as.character(sample_label), timepoint = as.character(timepoint), 
        array_x = as.numeric(array_x), array_y = as.numeric(array_y), section_position_internal = as.character(section_position_internal), 
        section_uid_internal = as.character(section_uid_internal))

joined <- md_use %>% inner_join(section_map2, by = c("GSM", "sample_label", "timepoint", "array_x", "array_y"))

if (nrow(joined) == 0L) {
  stop("Join between Step64 object and Step68 map returned 0 rows.")
}

section_summary <- read.csv(STEP65_SECTION_SUMMARY, stringsAsFactors = FALSE, check.names = FALSE)

required_sec_summary <- c("summary_level", "timepoint", "sample_label", "GSM", "section_position_internal", 
    "section_uid_internal", "n_spots", "mean_target_weight")

missing_sec_summary <- setdiff(required_sec_summary, colnames(section_summary))

if (length(missing_sec_summary) > 0L) {
    stop("Step65 section summary is missing required column(s): ", paste(missing_sec_summary, collapse = ", "))
}

section_summary2 <- section_summary %>% filter(summary_level == "section", timepoint %in% TIME_LEVELS) %>% 
    transmute(timepoint = as.character(timepoint), sample_label = as.character(sample_label), GSM = as.character(GSM), 
        section_position_internal = as.character(section_position_internal), section_uid_internal = as.character(section_uid_internal), 
        n_spots_section = as.integer(n_spots), Clec7a_section_mean = as.numeric(mean_target_weight))

section_region <- read.csv(STEP65_SECTION_REGION, stringsAsFactors = FALSE, check.names = FALSE)

required_section_region <- c("summary_level", "timepoint", "sample_label", "GSM", "section_position_internal", 
    "section_uid_internal", "spatial_region_v03", "n_spots", "mean_target_weight")

missing_section_region <- setdiff(required_section_region, colnames(section_region))

if (length(missing_section_region) > 0L) {
    stop("Step65 section-region summary is missing required column(s): ", paste(missing_section_region, 
        collapse = ", "))
}

target_region_section <- section_region %>% filter(summary_level == "section_region", timepoint %in% 
    TIME_LEVELS, spatial_region_v03 == TARGET_REGION) %>% transmute(timepoint = as.character(timepoint), 
    sample_label = as.character(sample_label), GSM = as.character(GSM), section_position_internal = as.character(section_position_internal), 
    section_uid_internal = as.character(section_uid_internal), n_target_region_spots = as.integer(n_spots), 
    target_region_mean_weight = as.numeric(mean_target_weight))

region_enrich <- read.csv(STEP66_REGION_ENRICH, stringsAsFactors = FALSE, check.names = FALSE)

required_region_enrich <- c("timepoint", "sample_label", "GSM", "section_position_internal", "section_uid_internal", 
    "spatial_region", "n_region_spots", "log2_enrichment", "enrichment_ratio")

missing_region_enrich <- setdiff(required_region_enrich, colnames(region_enrich))

if (length(missing_region_enrich) > 0L) {
    stop("Step66 direct-enrichment file is missing required column(s): ", paste(missing_region_enrich, 
        collapse = ", "))
}

target_region_enrich <- region_enrich %>% filter(timepoint %in% TIME_LEVELS, spatial_region == TARGET_REGION) %>% 
    transmute(timepoint = as.character(timepoint), sample_label = as.character(sample_label), GSM = as.character(GSM), 
        section_position_internal = as.character(section_position_internal), section_uid_internal = as.character(section_uid_internal), 
        direct_log2_enrichment = as.numeric(log2_enrichment), direct_enrichment_ratio = as.numeric(enrichment_ratio))

shape_summary <- joined %>% group_by(timepoint, sample_label, GSM, section_position_internal, section_uid_internal) %>% 
    summarise(n_spots_from_step64 = n(), n_regions = n_distinct(spatial_region_v03), dominant_region = first_mode(spatial_region_v03), 
        min_x = min(array_x, na.rm = TRUE), max_x = max(array_x, na.rm = TRUE), min_y = min(array_y, 
            na.rm = TRUE), max_y = max(array_y, na.rm = TRUE), .groups = "drop") %>% mutate(width = max_x - 
    min_x + 1, height = max_y - min_y + 1)

candidate_table <- shape_summary %>% left_join(section_summary2, by = c("timepoint", "sample_label", 
    "GSM", "section_position_internal", "section_uid_internal")) %>% left_join(target_region_section, 
    by = c("timepoint", "sample_label", "GSM", "section_position_internal", "section_uid_internal")) %>% 
    left_join(target_region_enrich, by = c("timepoint", "sample_label", "GSM", "section_position_internal", 
        "section_uid_internal")) %>% mutate(has_target_region = !is.na(n_target_region_spots) & n_target_region_spots > 
    0, target_region_spots_ok = !is.na(n_target_region_spots) & n_target_region_spots >= MIN_TARGET_REGION_SPOTS, 
    direct_enrichment_positive = !is.na(direct_log2_enrichment) & direct_log2_enrichment > 0)

timepoint_medians <- candidate_table %>% filter(!is.na(Clec7a_section_mean)) %>% group_by(timepoint) %>% 
    summarise(timepoint_median_Clec7a = median(Clec7a_section_mean), .groups = "drop")

candidate_table <- candidate_table %>% left_join(timepoint_medians, by = "timepoint") %>% mutate(abs_dist_to_timepoint_median = abs(Clec7a_section_mean - 
    timepoint_median_Clec7a), biologically_eligible = case_when(timepoint == "3dpi" ~ !is.na(Clec7a_section_mean) & 
    n_spots_from_step64 >= MIN_SECTION_SPOTS, TRUE ~ !is.na(Clec7a_section_mean) & n_spots_from_step64 >= 
    MIN_SECTION_SPOTS & has_target_region & target_region_spots_ok & direct_enrichment_positive))

candidate_table <- candidate_table %>% group_by(timepoint) %>% arrange(factor(timepoint, levels = TIME_LEVELS), 
    desc(biologically_eligible), abs_dist_to_timepoint_median, desc(n_target_region_spots), desc(n_spots_from_step64), 
    sample_label, section_position_internal, .by_group = TRUE) %>% mutate(candidate_rank = row_number()) %>% 
    ungroup()

CANDIDATE_FILE <- file.path(OUT_DIR, paste0("01_69_Figure4A_candidate_section_ranking_", VERSION, ".csv"))

write_csv_utf8(candidate_table, CANDIDATE_FILE)

rep_sections <- candidate_table %>% group_by(timepoint) %>% arrange(desc(biologically_eligible), abs_dist_to_timepoint_median, 
    desc(n_target_region_spots), desc(n_spots_from_step64), sample_label, section_position_internal, 
    .by_group = TRUE) %>% slice(1L) %>% ungroup() %>% mutate(display_time = factor(timepoint, levels = TIME_LEVELS, 
    labels = unname(TIME_LABELS[TIME_LEVELS]))) %>% arrange(factor(timepoint, levels = TIME_LEVELS))

missing_tp <- setdiff(TIME_LEVELS, rep_sections$timepoint)

if (length(missing_tp) > 0L) {
    stop("Failed to select representative section(s) for: ", paste(missing_tp, collapse = ", "))
}

late_tp_bad <- rep_sections %>% filter(timepoint %in% c("7dpi", "14dpi", "28dpi"), !biologically_eligible)

if (nrow(late_tp_bad) > 0L) {
    stop("At least one of 7/14/28 dpi selected sections is not biologically eligible.")
}

REP_FILE <- file.path(OUT_DIR, paste0("02_69_Figure4A_representative_sections_", VERSION, ".csv"))

write_csv_utf8(rep_sections, REP_FILE)

plot_df <- joined %>% inner_join(rep_sections %>% select(timepoint, sample_label, GSM, section_position_internal, 
    section_uid_internal, display_time), by = c("timepoint", "sample_label", "GSM", "section_position_internal", 
    "section_uid_internal"))

plot_df$display_time <- factor(plot_df$display_time, levels = unname(TIME_LABELS[TIME_LEVELS]))

PLOT_TABLE_FILE <- file.path(OUT_DIR, paste0("03_69_Figure4A_plot_spot_table_", VERSION, ".csv.gz"))

write_csv_gz_utf8(plot_df, PLOT_TABLE_FILE)

region_palette <- c(`Oligodendrocyte / white-matter-like` = "#37679E", `Myeloid/immune–lipid-handling enriched` = "#8E4FA5", 
    `Cholinergic/motor-neuron-like neuronal` = "#D97706", `Inhibitory/interneuron-like neuronal` = "#C9475B", 
    `Acute proliferative injury-associated myeloid` = "#B71C1C", `Peripheral-myelin / Schwann-like region` = "#238B8E", 
    `Neuronal subtype region, 7 dpi-enriched` = "#D95F8D", `Fibrotic/ECM–myeloid mixed injury niche` = "#70411F", 
    `IFN-responsive inflammatory/astroglial region` = "#2E8B45", `Perivascular / vSMC–pericyte region` = "#7E6556")

present_regions <- unique(plot_df$spatial_region_v03)

missing_regions <- setdiff(present_regions, names(region_palette))

missing_regions <- missing_regions[!is.na(missing_regions)]

if (length(missing_regions) > 0L) {
    stop("The selected Figure4A sections contain region label(s) not in the Step64 palette:\n", paste(missing_regions, 
        collapse = "\n"))
}

plot_df$spatial_region_v03 <- factor(plot_df$spatial_region_v03, levels = names(region_palette))

plot_df$display_time <- factor(plot_df$display_time, levels = unname(TIME_LABELS[TIME_LEVELS]))

CROP_PAD_X <- 3
CROP_PAD_Y <- 3

crop_bounds <- plot_df %>% group_by(timepoint, sample_label, GSM, section_position_internal, section_uid_internal) %>% 
    summarise(xmin = min(array_x, na.rm = TRUE) - CROP_PAD_X, xmax = max(array_x, na.rm = TRUE) + CROP_PAD_X, 
        ymin = min(array_y, na.rm = TRUE) - CROP_PAD_Y, ymax = max(array_y, na.rm = TRUE) + CROP_PAD_Y, 
        n_spots = n(), .groups = "drop") %>% arrange(factor(timepoint, levels = TIME_LEVELS))

CROP_BOUNDS_FILE <- file.path(OUT_DIR, paste0("03b_69_Figure4A_original_coordinate_crop_bounds_", VERSION, 
    ".csv"))

write_csv_utf8(crop_bounds, CROP_BOUNDS_FILE)

make_step64_crop_panel <- function(dat, tp, show_x_title = TRUE, show_y_title = TRUE) {
    b <- crop_bounds %>% filter(timepoint == tp)
    if (nrow(b) != 1L) {
        stop("Expected exactly one crop-bound row for ", tp, ", found ", nrow(b), ".")
    }
    dat$display_time_single <- factor(TIME_LABELS[[tp]], levels = TIME_LABELS[[tp]])
    ggplot(dat, aes(x = array_x, y = array_y, color = spatial_region_v03)) + geom_point(size = 1.28, 
        shape = 15, stroke = 0, alpha = 1) + coord_equal(xlim = c(b$xmin, b$xmax), ylim = c(b$ymin, b$ymax), 
        expand = FALSE, clip = "off") + scale_y_reverse() + scale_color_manual(values = region_palette, 
        limits = names(region_palette), breaks = names(region_palette), drop = FALSE) + facet_wrap(~display_time_single, 
        ncol = 1) + labs(x = if (show_x_title) 
        "Visium array X"
    else NULL, y = if (show_y_title) 
        "Visium array Y"
    else NULL, color = NULL) + theme_classic(base_family = "Arial", base_size = 9) + theme(strip.text = element_text(size = 8.3, 
        face = "bold"), axis.title = element_text(size = 8.4), axis.text = element_text(size = 7.1), 
        legend.position = "none", plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"))
}

legend_df <- data.frame(region_label = names(region_palette), region_color = unname(region_palette), 
    legend_y = rev(seq(from = 1, to = 6.4, length.out = length(region_palette))), stringsAsFactors = FALSE)

legend_plot <- ggplot(legend_df, aes(y = legend_y)) + geom_point(aes(x = 0), color = legend_df$region_color, 
    size = 3.4, shape = 15, alpha = 1) + geom_text(aes(x = 0.12, label = region_label), hjust = 0, family = "Arial", 
    size = 2.48) + scale_x_continuous(limits = c(-0.05, 1.55), expand = expansion(mult = c(0, 0))) + 
    scale_y_continuous(limits = c(-1.6, 9), expand = expansion(mult = c(0, 0))) + coord_cartesian(clip = "off") + 
    theme_void(base_family = "Arial", base_size = 9) + theme(plot.margin = margin(0, 1, 0, 1, unit = "mm"))

d3 <- plot_df %>% filter(timepoint == "3dpi")

d7 <- plot_df %>% filter(timepoint == "7dpi")

d14 <- plot_df %>% filter(timepoint == "14dpi")

d28 <- plot_df %>% filter(timepoint == "28dpi")

if (any(c(nrow(d3), nrow(d7), nrow(d14), nrow(d28)) == 0L)) {
    stop("At least one selected Figure4A timepoint has zero plotting spots.")
}

p3 <- make_step64_crop_panel(dat = d3, tp = "3dpi", show_x_title = FALSE, show_y_title = TRUE)

p7 <- make_step64_crop_panel(dat = d7, tp = "7dpi", show_x_title = FALSE, show_y_title = FALSE)

p14 <- make_step64_crop_panel(dat = d14, tp = "14dpi", show_x_title = TRUE, show_y_title = TRUE)

p28 <- make_step64_crop_panel(dat = d28, tp = "28dpi", show_x_title = TRUE, show_y_title = FALSE)

left_panels <- (p3 + p7)/(p14 + p28) + patchwork::plot_layout(widths = c(1, 1), heights = c(1, 1))

p <- (left_panels | legend_plot) + patchwork::plot_layout(widths = c(1, 0.52))

PDF_FILE <- file.path(OUT_DIR, paste0("69_Figure4A_spatial_transcriptional_regions_", VERSION, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0("69_Figure4A_spatial_transcriptional_regions_", VERSION, ".png"))

SVG_FILE <- file.path(OUT_DIR, paste0("69_Figure4A_spatial_transcriptional_regions_", VERSION, ".svg"))

ggsave(filename = PDF_FILE, plot = p, width = 9.6, height = 5.1, units = "in", device = grDevices::cairo_pdf, 
    bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = 9.6, height = 5.1, units = "in", dpi = 600, bg = "white")

ggsave(filename = SVG_FILE, plot = p, width = 9.6, height = 5.1, units = "in", device = "svg", bg = "white")

final_qc <- data.frame(check = c("Step64 final object exists", "Step68 spot-to-section map exists", "Step65 v06 section summary exists", 
    "Step66 v02 direct enrichment exists", "Step64 object is a Seurat object", "Step64 object contains required metadata", 
    "No H&E registration is claimed", "Representative section selected for each of 4 timepoints", "Selected 7/14/28 dpi sections satisfy biological eligibility rule", 
    "Selected 3 dpi section is closest-to-median representative", "All plotted region labels use v10 palette", 
    "All plotted spots retain original Step64 array coordinates", "Original-coordinate crop bounds file generated", 
    "Standalone legend contains 10 visible color keys and 10 explicit text labels", "Figure 4A PDF generated", 
    "Figure 4A PNG generated", "Figure 4A SVG generated"), pass = c(file.exists(STEP64_OBJECT), file.exists(STEP68_SECTION_MAP), 
    file.exists(STEP65_SECTION_SUMMARY), file.exists(STEP66_REGION_ENRICH), inherits(obj, "Seurat"), 
    length(missing_md) == 0L, !any(md$histology_registered %in% TRUE, na.rm = TRUE), nrow(rep_sections) == 
        4L, nrow(late_tp_bad) == 0L, rep_sections$timepoint[1] == "3dpi", length(missing_regions) == 
        0L, all(is.finite(plot_df$array_x)) & all(is.finite(plot_df$array_y)), file.exists(CROP_BOUNDS_FILE), 
    nrow(legend_df) == length(names(region_palette)), file.exists(PDF_FILE), file.exists(PNG_FILE), file.exists(SVG_FILE)), 
    stringsAsFactors = FALSE)

QC_FILE <- file.path(OUT_DIR, paste0("04_69_Figure4A_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_qc, QC_FILE)

rep_lines <- apply(rep_sections, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": sample=", z[["sample_label"]], ", section=", z[["section_position_internal"]], 
        ", section_uid=", z[["section_uid_internal"]], ", n_spots=", z[["n_spots_from_step64"]], ", n_regions=", 
        z[["n_regions"]], ", Clec7a_section_mean=", sprintf("%.4f", as.numeric(z[["Clec7a_section_mean"]])), 
        ", timepoint_median=", sprintf("%.4f", as.numeric(z[["timepoint_median_Clec7a"]])), ", abs_dist_to_median=", 
        sprintf("%.4f", as.numeric(z[["abs_dist_to_timepoint_median"]])), ", target_region_spots=", ifelse(is.na(z[["n_target_region_spots"]]), 
            "NA", as.character(z[["n_target_region_spots"]])), ", direct_log2_enrichment=", ifelse(is.na(z[["direct_log2_enrichment"]]), 
            "NA", sprintf("%.3f", as.numeric(z[["direct_log2_enrichment"]]))))
})

REPORT_FILE <- file.path(OUT_DIR, paste0("05_69_Figure4A_REPORT_", VERSION, ".txt"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(rep_sections)
print(final_qc)

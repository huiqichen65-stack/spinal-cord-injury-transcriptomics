rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(7001)

ROOT <- PROJECT_ROOT

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

STEP69_DIR <- file.path(ROOT, "results", "69_Figure4A_spatial_transcriptional_regions")

OUT_DIR <- file.path(ROOT, "results", "70_Figure4B_Clec7a_spatial_abundance")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REP_SECTION_FILE <- file.path(STEP69_DIR, "02_69_Figure4A_representative_sections_v12.csv")

SPOT_FILE_CANDIDATES <- c(file.path(STEP65_DIR, "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv.gz"),
  file.path(STEP65_DIR, "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv"), file.path(STEP65_DIR,
    "03_65_RCTD_normalized_spot_weights_v06.csv.gz"), file.path(STEP65_DIR, "03_65_RCTD_normalized_spot_weights_v06.csv"))

VERSION <- "v08"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

TIME_LABELS <- c(`3dpi` = "3 dpi", `7dpi` = "7 dpi", `14dpi` = "14 dpi", `28dpi` = "28 dpi")

TARGET_STATE <- "Clec7a+ lipid-handling Microglia"

COLOR_MODE <- "raw"

CROP_PAD_X <- 1
CROP_PAD_Y <- 1

required_pkgs <- c("ggplot2", "dplyr", "patchwork", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun.")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

write_csv_utf8 <- write_csv_utf8

write_csv_gz_utf8 <- function(x, filename) {
  con <- gzfile(filename, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  write.table(x, file = con, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "")
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) {
    return(NA_character_)
  }
  hit[1]
}

safe_read_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

pick_first_existing_col <- function(df, candidates, field_name) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0L) {
    stop("Could not find a usable column for '", field_name, "'. Tried:\n", paste(candidates, collapse = "\n"))
  }
  hit[1]
}

SPOT_FILE <- first_existing_file(SPOT_FILE_CANDIDATES)

if (is.na(SPOT_FILE)) {
  stop("Could not find the Step65 spot-level normalized-weight file.\n", "Tried:\n", paste(SPOT_FILE_CANDIDATES,
    collapse = "\n"))
}

required_files <- c(REP_SECTION_FILE, SPOT_FILE)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("Required file(s) missing:\n", paste(missing_files, collapse = "\n"))
}

rep_sections <- safe_read_csv(REP_SECTION_FILE)

required_rep <- c("timepoint", "sample_label", "GSM", "section_position_internal", "section_uid_internal")

missing_rep <- setdiff(required_rep, colnames(rep_sections))

if (length(missing_rep) > 0L) {
  stop("Representative-section file is missing required column(s): ", paste(missing_rep, collapse = ", "))
}

rep_sections <- rep_sections %>% filter(timepoint %in% TIME_LEVELS) %>% distinct(timepoint, sample_label, 
    GSM, section_position_internal, section_uid_internal, .keep_all = TRUE) %>% arrange(factor(timepoint, 
    levels = TIME_LEVELS)) %>% mutate(display_time = factor(timepoint, levels = TIME_LEVELS, labels = unname(TIME_LABELS[TIME_LEVELS])))

if (nrow(rep_sections) != 4L) {
  stop("Expected exactly 4 representative sections, found ", nrow(rep_sections), ".")
}

spot_df0 <- safe_read_csv(SPOT_FILE)

col_GSM <- pick_first_existing_col(spot_df0, c("GSM"), "GSM")

col_sample <- pick_first_existing_col(spot_df0, c("sample_label", "sample", "sample_id"), "sample_label")

col_time <- pick_first_existing_col(spot_df0, c("timepoint", "time"), "timepoint")

col_x <- pick_first_existing_col(spot_df0, c("array_x", "x", "spot_x"), "array_x")

col_y <- pick_first_existing_col(spot_df0, c("array_y", "y", "spot_y"), "array_y")

col_secpos <- pick_first_existing_col(spot_df0, c("section_position_internal", "section_position"), "section_position_internal")

col_secuid <- pick_first_existing_col(spot_df0, c("section_uid_internal", "section_uid"), "section_uid_internal")

if (!(TARGET_STATE %in% colnames(spot_df0))) {
  candidate_hits <- grep(pattern = "Clec7a.*lipid.*Microglia|Clec7a.*lipid", x = colnames(spot_df0),
    ignore.case = TRUE, value = TRUE)
  if (length(candidate_hits) == 0L) {
    stop("Could not find the target state column.\n", "Expected exact column: ", TARGET_STATE, "\n",
      "Available columns include:\n", paste(head(colnames(spot_df0), 60), collapse = ", "))
  }
  target_col <- candidate_hits[1]
} else {
  target_col <- TARGET_STATE
}

spot_df <- spot_df0 %>% transmute(GSM = as.character(.data[[col_GSM]]), sample_label = as.character(.data[[col_sample]]), 
    timepoint = as.character(.data[[col_time]]), array_x = as.numeric(.data[[col_x]]), array_y = as.numeric(.data[[col_y]]), 
    section_position_internal = as.character(.data[[col_secpos]]), section_uid_internal = as.character(.data[[col_secuid]]), 
    target_weight = as.numeric(.data[[target_col]])) %>% filter(timepoint %in% TIME_LEVELS)

plot_df <- spot_df %>% inner_join(rep_sections %>% dplyr::select(timepoint, sample_label, GSM, section_position_internal, 
    section_uid_internal, display_time), by = c("timepoint", "sample_label", "GSM", "section_position_internal", 
    "section_uid_internal"))

if (nrow(plot_df) == 0L) {
  stop("Join between Step65 spot-level weights and the Figure4A representative sections returned 0 rows.")
}

crop_bounds <- plot_df %>% group_by(timepoint, sample_label, GSM, section_position_internal, section_uid_internal) %>% 
    summarise(xmin = min(array_x, na.rm = TRUE) - CROP_PAD_X, xmax = max(array_x, na.rm = TRUE) + CROP_PAD_X, 
        ymin = min(array_y, na.rm = TRUE) - CROP_PAD_Y, ymax = max(array_y, na.rm = TRUE) + CROP_PAD_Y, 
        n_spots = n(), .groups = "drop") %>% arrange(factor(timepoint, levels = TIME_LEVELS))

CROP_FILE <- file.path(OUT_DIR, paste0("03_70_Figure4B_crop_bounds_", VERSION, ".csv"))

write_csv_utf8(crop_bounds, CROP_FILE)

plot_df <- plot_df %>% mutate(target_weight = ifelse(is.na(target_weight), 0, target_weight), target_weight = pmax(target_weight, 
    0))

if (COLOR_MODE != "raw") {
  stop("Figure4B v02 is locked to raw normalized RCTD weight. ", "Please keep COLOR_MODE <- \"raw\".")
}

plot_df <- plot_df %>% mutate(signal_raw = target_weight, signal_plot = target_weight)

legend_title <- "Normalized\nRCTD weight"

upper_lim <- as.numeric(quantile(plot_df$signal_plot, probs = 0.99, na.rm = TRUE, names = FALSE))

if (!is.finite(upper_lim) || upper_lim <= 0) {
  upper_lim <- max(plot_df$signal_plot, na.rm = TRUE)
}

if (!is.finite(upper_lim) || upper_lim <= 0) {
  stop("Could not determine a positive shared upper limit for the RCTD color scale.")
}

color_limits <- c(0, upper_lim)

SCALE_SUMMARY <- data.frame(metric = c("COLOR_MODE", "n_total_plotted_spots", "signal_min", "signal_max",
  "signal_q01", "signal_q50", "signal_q99", "color_limit_low", "color_limit_high"), value = c(COLOR_MODE,
  nrow(plot_df), min(plot_df$signal_plot, na.rm = TRUE), max(plot_df$signal_plot, na.rm = TRUE), as.numeric(quantile(plot_df$signal_plot,
    0.01, na.rm = TRUE)), as.numeric(quantile(plot_df$signal_plot, 0.5, na.rm = TRUE)), as.numeric(quantile(plot_df$signal_plot,
    0.99, na.rm = TRUE)), color_limits[1], color_limits[2]), stringsAsFactors = FALSE)

SCALE_FILE <- file.path(OUT_DIR, paste0("04_70_Figure4B_color_scale_summary_", VERSION, ".csv"))

write_csv_utf8(SCALE_SUMMARY, SCALE_FILE)

PLOT_TABLE_FILE <- file.path(OUT_DIR, paste0("05_70_Figure4B_plot_spot_table_", VERSION, ".csv.gz"))

write_csv_gz_utf8(plot_df, PLOT_TABLE_FILE)

make_panel <- function(dat, tp, show_x_title = TRUE, show_y_title = TRUE) {
    b <- crop_bounds %>% filter(timepoint == tp)
    if (nrow(b) != 1L) {
        stop("Expected exactly one crop-bound row for ", tp, ", found ", nrow(b), ".")
    }
    dat$display_time_single <- factor(TIME_LABELS[[tp]], levels = TIME_LABELS[[tp]])
    ggplot(dat, aes(x = array_x, y = array_y, color = signal_plot)) + geom_point(size = 1.6, shape = 16, 
        stroke = 0, alpha = 1) + coord_equal(xlim = c(b$xmin, b$xmax), ylim = c(b$ymin, b$ymax), expand = FALSE, 
        clip = "off") + scale_y_reverse() + scale_color_gradientn(colors = c("#000000", "#12002F", "#3B0A66", 
        "#6F1D70", "#A92F62", "#D94B3D", "#F47C17", "#F6B51F"), values = c(0, 0.1, 0.24, 0.4, 0.58, 0.74, 
        0.88, 1), limits = color_limits, oob = scales::squish, name = legend_title, guide = guide_colorbar(title.position = "top", 
        title.hjust = 0.5, barheight = unit(3.2, "cm"), barwidth = unit(0.42, "cm"), ticks = TRUE, frame.colour = "black", 
        frame.linewidth = 0.35)) + facet_wrap(~display_time_single, ncol = 1) + labs(x = if (show_x_title) 
        "Visium array X"
    else NULL, y = if (show_y_title) 
        "Visium array Y"
    else NULL) + theme_classic(base_family = "Arial", base_size = 9) + theme(strip.background = element_rect(fill = "white", 
        color = "black", linewidth = 0.8), strip.text = element_text(size = 8.5, face = "bold"), axis.title = element_text(size = 8.6), 
        axis.text = element_text(size = 7), legend.title = element_text(size = 8, hjust = 0.5), legend.text = element_text(size = 7.4), 
        plot.margin = margin(0.6, 0.6, 0.6, 0.6, unit = "mm"))
}

d3 <- plot_df %>% filter(timepoint == "3dpi")

d7 <- plot_df %>% filter(timepoint == "7dpi")

d14 <- plot_df %>% filter(timepoint == "14dpi")

d28 <- plot_df %>% filter(timepoint == "28dpi")

if (any(c(nrow(d3), nrow(d7), nrow(d14), nrow(d28)) == 0L)) {
  stop("At least one selected timepoint has zero plotting spots.")
}

p3 <- make_panel(dat = d3, tp = "3dpi", show_x_title = FALSE, show_y_title = TRUE)

p7 <- make_panel(dat = d7, tp = "7dpi", show_x_title = FALSE, show_y_title = FALSE)

p14 <- make_panel(dat = d14, tp = "14dpi", show_x_title = TRUE, show_y_title = TRUE)

p28 <- make_panel(dat = d28, tp = "28dpi", show_x_title = TRUE, show_y_title = FALSE)

p <- (p3 + p7)/(p14 + p28) + patchwork::plot_layout(guides = "collect", widths = c(1, 1), heights = c(1,
  1)) & theme(legend.position = "right", legend.justification = "center")

PDF_FILE <- file.path(OUT_DIR, paste0("70_Figure4B_Clec7a_spatial_abundance_", VERSION, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0("70_Figure4B_Clec7a_spatial_abundance_", VERSION, ".png"))

SVG_FILE <- file.path(OUT_DIR, paste0("70_Figure4B_Clec7a_spatial_abundance_", VERSION, ".svg"))

ggsave(filename = PDF_FILE, plot = p, width = 8.8, height = 5, units = "in", device = grDevices::cairo_pdf,
  bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = 8.8, height = 5, units = "in", dpi = 600, bg = "white")

ggsave(filename = SVG_FILE, plot = p, width = 8.8, height = 5, units = "in", device = "svg", bg = "white")

final_qc <- data.frame(check = c("Figure4A representative-section file exists", "Step65 spot-level normalized-weight file exists",
  "Exactly 4 representative sections loaded", "Representative-section join to spot table is non-empty",
  "Target-state column detected", "Raw normalized RCTD weight used as plotting signal", "One shared 0-to-q99 color scale used across all four panels",
  "Darker custom inferno-like palette is used without changing signal values", "No dummy scatter panel is used for the legend",
  "Crop-bounds file generated", "Color-scale summary generated", "Figure4B PDF generated", "Figure4B PNG generated",
  "Figure4B SVG generated"), pass = c(file.exists(REP_SECTION_FILE), file.exists(SPOT_FILE), nrow(rep_sections) ==
  4L, nrow(plot_df) > 0L, target_col %in% colnames(spot_df0), identical(COLOR_MODE, "raw") & all(abs(plot_df$signal_plot -
  plot_df$target_weight) < 1e-12), isTRUE(all.equal(color_limits[1], 0)) & is.finite(color_limits[2]) &
  color_limits[2] > 0, TRUE, TRUE, file.exists(CROP_FILE), file.exists(SCALE_FILE), file.exists(PDF_FILE),
  file.exists(PNG_FILE), file.exists(SVG_FILE)), stringsAsFactors = FALSE)

QC_FILE <- file.path(OUT_DIR, paste0("06_70_Figure4B_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_qc, QC_FILE)

rep_lines <- apply(rep_sections, 1, function(z) {
  paste0("  ", z[["timepoint"]], ": sample=", z[["sample_label"]], ", section=", z[["section_position_internal"]],
    ", section_uid=", z[["section_uid_internal"]])
})

REPORT_FILE <- file.path(OUT_DIR, paste0("07_70_Figure4B_REPORT_", VERSION, ".txt"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(rep_sections[, c("timepoint", "sample_label", "section_position_internal", "section_uid_internal")])

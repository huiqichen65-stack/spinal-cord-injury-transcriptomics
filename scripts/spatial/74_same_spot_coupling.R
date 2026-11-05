rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)


options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(7401)

ROOT <- PROJECT_ROOT

STEP73_DIR <- file.path(ROOT, "results", "73_GSE195783_Clec7a_lipid_spatial_coupling")

STEP69_DIR <- file.path(ROOT, "results", "69_Figure4A_spatial_transcriptional_regions")

INPUT_SPOT_FILE <- file.path(STEP73_DIR, "02_73_Clec7a_lipid_spot_scores_section_aware_v03.csv.gz")

REP_SECTION_FILE <- file.path(STEP69_DIR, "02_69_Figure4A_representative_sections_v12.csv")

OUT_DIR <- file.path(ROOT, "results", "74_GSE195783_Clec7a_lipid_same_spot_coupling")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v01"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

TIME_LABELS <- c(`3dpi` = "3 dpi", `7dpi` = "7 dpi", `14dpi` = "14 dpi", `28dpi` = "28 dpi")

SAMPLE_ORDER <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1", 
    "28dpi_rep2")

EXPECTED_CAPTURE_AREAS <- 8L
EXPECTED_MAJOR_SECTIONS <- 32L
MIN_SECTION_SPOTS <- 100L

CROP_PAD_X <- 0.3
CROP_PAD_Y <- 0.3

BACKGROUND_COLOR <- "#E1E1E1"
HH_LOW <- "#FEE8C8"
HH_MID <- "#F46D43"
HH_HIGH <- "#A50026"

STATE_COLORS <- c(`High–High` = "#B2182B", `High–Low` = "#EF8A62", `Low–High` = "#67A9CF", `Low–Low` = "#B8C9D9")

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "patchwork", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step74.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

write_csv_utf8 <- write_csv_utf8

write_csv_gz_base <- write_gzip_csv_utf8

read_csv_auto <- read_csv_auto

as_logical_safe <- function(x) {
    if (is.logical(x)) {
        return(x)
    }
    z <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(z))
    out[z %in% c("true", "t", "1", "yes", "y")] <- TRUE
    out[z %in% c("false", "f", "0", "no", "n")] <- FALSE
    out
}

zscore_safe <- function(x) {
    x <- as.numeric(x)
    mu <- mean(x, na.rm = TRUE)
    sig <- sd(x, na.rm = TRUE)
    if (!is.finite(sig) || sig <= 1e-12) {
        return(rep(NA_real_, length(x)))
    }
    (x - mu)/sig
}

required_files <- c(INPUT_SPOT_FILE, REP_SECTION_FILE)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop("Required input file(s) missing:\n", paste(missing_files, collapse = "\n"))
}

spot0 <- read_csv_auto(INPUT_SPOT_FILE)

required_spot_cols <- c("barcode", "GSM", "sample_label", "timepoint", "array_x", "array_y", "section_uid_internal", 
    "Clec7a_RCTDnorm", "score_Core_lipid_handling")

missing_spot_cols <- setdiff(required_spot_cols, colnames(spot0))

if (length(missing_spot_cols) > 0L) {
    stop("Step73 v03 spot file is missing required column(s): ", paste(missing_spot_cols, collapse = " | "))
}

if ("section_eligible" %in% colnames(spot0)) {
    spot0$section_eligible <- as_logical_safe(spot0$section_eligible)
    if (any(is.na(spot0$section_eligible))) {
        stop("Could not safely parse section_eligible.")
    }
} else {
    spot0$section_eligible <- !is.na(spot0$section_uid_internal) & nzchar(trimws(as.character(spot0$section_uid_internal)))
}

spot <- spot0 %>% transmute(barcode = as.character(barcode), GSM = as.character(GSM), sample_label = as.character(sample_label), 
    timepoint = as.character(timepoint), array_x = as.numeric(array_x), array_y = as.numeric(array_y), 
    section_uid_internal = as.character(section_uid_internal), section_eligible = section_eligible, Clec7a_RCTDnorm = as.numeric(Clec7a_RCTDnorm), 
    lipid_score = as.numeric(score_Core_lipid_handling)) %>% filter(section_eligible, timepoint %in% 
    TIME_LEVELS, is.finite(array_x), is.finite(array_y), is.finite(Clec7a_RCTDnorm), is.finite(lipid_score))

if (!setequal(unique(spot$sample_label), SAMPLE_ORDER)) {
    stop("Observed sample labels do not match the expected 8 capture areas.\n", "Observed: ", paste(sort(unique(spot$sample_label)), 
        collapse = " | "))
}

n_capture_areas <- n_distinct(spot$sample_label)

n_sections <- n_distinct(spot$section_uid_internal)

if (n_capture_areas != EXPECTED_CAPTURE_AREAS) {
    stop("Expected ", EXPECTED_CAPTURE_AREAS, " capture areas; found ", n_capture_areas, ".")
}

if (n_sections != EXPECTED_MAJOR_SECTIONS) {
    stop("Expected ", EXPECTED_MAJOR_SECTIONS, " major sections; found ", n_sections, ".")
}

spot$z_Clec7a_same_spot <- NA_real_
spot$z_lipid_same_spot <- NA_real_
spot$same_spot_state <- NA_character_
spot$same_spot_HH_strength <- 0

section_order <- spot %>% distinct(section_uid_internal, sample_label, GSM, timepoint) %>% arrange(factor(timepoint, 
    levels = TIME_LEVELS), factor(sample_label, levels = SAMPLE_ORDER), section_uid_internal) %>% pull(section_uid_internal)

for (sec in section_order) {
    idx <- which(spot$section_uid_internal == sec)
    d <- spot[idx, , drop = FALSE]
    if (nrow(d) < MIN_SECTION_SPOTS) {
        stop("Section ", sec, " contains only ", nrow(d), " spots (<", MIN_SECTION_SPOTS, ").")
    }
    z_clec <- zscore_safe(d$Clec7a_RCTDnorm)
    z_lipid <- zscore_safe(d$lipid_score)
    if (any(!is.finite(z_clec)) || any(!is.finite(z_lipid))) {
        stop("Non-finite within-section z-score detected in section ", sec, ".")
    }
    state <- ifelse(z_clec > 0 & z_lipid > 0, "High–High", ifelse(z_clec > 0 & z_lipid <= 0, "High–Low", 
        ifelse(z_clec <= 0 & z_lipid > 0, "Low–High", "Low–Low")))
    hh_strength <- ifelse(state == "High–High", z_clec * z_lipid, 0)
    spot$z_Clec7a_same_spot[idx] <- z_clec
    spot$z_lipid_same_spot[idx] <- z_lipid
    spot$same_spot_state[idx] <- state
    spot$same_spot_HH_strength[idx] <- hh_strength
}

section_summary <- spot %>% group_by(section_uid_internal, sample_label, GSM, timepoint) %>% summarise(n_spots = n(), 
    fraction_HH = mean(same_spot_state == "High–High"), fraction_HL = mean(same_spot_state == "High–Low"), 
    fraction_LH = mean(same_spot_state == "Low–High"), fraction_LL = mean(same_spot_state == "Low–Low"), 
    median_HH_strength = if (any(same_spot_HH_strength > 0)) {
        median(same_spot_HH_strength[same_spot_HH_strength > 0])
    } else {
        0
    }, mean_HH_strength = mean(same_spot_HH_strength), same_spot_spearman_rho = suppressWarnings(cor(Clec7a_RCTDnorm, 
        lipid_score, method = "spearman")), .groups = "drop") %>% arrange(factor(timepoint, levels = TIME_LEVELS), 
    factor(sample_label, levels = SAMPLE_ORDER), section_uid_internal)

time_summary <- section_summary %>% group_by(timepoint) %>% summarise(n_sections = n(), n_capture_areas = n_distinct(sample_label), 
    median_section_rho = median(same_spot_spearman_rho, na.rm = TRUE), q1_section_rho = quantile(same_spot_spearman_rho, 
        0.25, na.rm = TRUE), q3_section_rho = quantile(same_spot_spearman_rho, 0.75, na.rm = TRUE), median_fraction_HH = median(fraction_HH, 
        na.rm = TRUE), q1_fraction_HH = quantile(fraction_HH, 0.25, na.rm = TRUE), q3_fraction_HH = quantile(fraction_HH, 
        0.75, na.rm = TRUE), median_HH_strength = median(median_HH_strength, na.rm = TRUE), .groups = "drop") %>% 
    arrange(match(timepoint, TIME_LEVELS))

rep_sections <- read_csv_auto(REP_SECTION_FILE)

required_rep_cols <- c("timepoint", "sample_label", "GSM", "section_uid_internal")

missing_rep_cols <- setdiff(required_rep_cols, colnames(rep_sections))

if (length(missing_rep_cols) > 0L) {
    stop("Representative-section file is missing required column(s): ", paste(missing_rep_cols, collapse = " | "))
}

rep_sections <- rep_sections %>% transmute(timepoint = as.character(timepoint), sample_label = as.character(sample_label), 
    GSM = as.character(GSM), section_uid_internal = as.character(section_uid_internal)) %>% filter(timepoint %in% 
    TIME_LEVELS) %>% distinct() %>% arrange(factor(timepoint, levels = TIME_LEVELS)) %>% mutate(display_time = factor(timepoint, 
    levels = TIME_LEVELS, labels = unname(TIME_LABELS[TIME_LEVELS])))

if (nrow(rep_sections) != 4L) {
    stop("Expected exactly four Figure4A representative sections; found ", nrow(rep_sections), ".")
}

plot_df <- spot %>% inner_join(rep_sections, by = c("timepoint", "sample_label", "GSM", "section_uid_internal"))

if (nrow(plot_df) == 0L) {
    stop("Representative-section join returned zero spots.")
}

crop_bounds <- plot_df %>% group_by(timepoint, sample_label, GSM, section_uid_internal) %>% summarise(xmin = min(array_x, 
    na.rm = TRUE) - CROP_PAD_X, xmax = max(array_x, na.rm = TRUE) + CROP_PAD_X, ymin = min(array_y, na.rm = TRUE) - 
    CROP_PAD_Y, ymax = max(array_y, na.rm = TRUE) + CROP_PAD_Y, n_spots = n(), .groups = "drop") %>% 
    arrange(factor(timepoint, levels = TIME_LEVELS))

SPOT_OUTPUT_FILE <- file.path(OUT_DIR, paste0("01_74_same_spot_coupling_spot_scores_", VERSION, ".csv.gz"))

write_csv_gz_base(spot, SPOT_OUTPUT_FILE)

SECTION_FILE <- file.path(OUT_DIR, paste0("02_74_same_spot_section_summary_", VERSION, ".csv"))

write_csv_utf8(section_summary, SECTION_FILE)

TIME_FILE <- file.path(OUT_DIR, paste0("03_74_same_spot_timepoint_descriptive_", VERSION, ".csv"))

write_csv_utf8(time_summary, TIME_FILE)

REP_FILE <- file.path(OUT_DIR, paste0("04_74_same_spot_representative_sections_", VERSION, ".csv"))

write_csv_utf8(rep_sections, REP_FILE)

hh_positive <- plot_df$same_spot_HH_strength[is.finite(plot_df$same_spot_HH_strength) & plot_df$same_spot_HH_strength > 
    0]

if (length(hh_positive) == 0L) {
    stop("No High–High same-spot coupling spots were found in the four representative sections.")
}

HH_DISPLAY_MAX <- as.numeric(quantile(hh_positive, probs = 0.99, na.rm = TRUE, names = FALSE))

if (!is.finite(HH_DISPLAY_MAX) || HH_DISPLAY_MAX <= 0) {
    HH_DISPLAY_MAX <- max(hh_positive, na.rm = TRUE)
}

make_hh_panel <- function(
    dat,
    tp,
    show_x_title = TRUE,
    show_y_title = TRUE) {

  b <- crop_bounds %>%
    filter(
      timepoint ==
        tp
    )

  if (nrow(b) !=
      1L) {

    stop(
      "Expected one crop row for ",
      tp,
      "."
    )
  }

  dat$display_time_single <-
    factor(
      TIME_LABELS[[tp]],
      levels =
        TIME_LABELS[[tp]]
    )

  hh <- dat %>%
    filter(
      same_spot_state ==
        "High–High",
      same_spot_HH_strength >
        0
    )

  ggplot() +

    geom_point(
      data = dat,
      aes(
        x = array_x,
        y = array_y
      ),
      size = 1.90,
      shape = 16,
      stroke = 0,
      alpha = 1,
      color =
        BACKGROUND_COLOR
    ) +

    geom_point(
      data = hh,
      aes(
        x = array_x,
        y = array_y,
        fill =
          same_spot_HH_strength
      ),
      size = 1.90,
      shape = 21,
      stroke = 0.18,
      color = "#7F0000"
    ) +

    coord_equal(
      xlim = c(
        b$xmin,
        b$xmax
      ),
      ylim = c(
        b$ymin,
        b$ymax
      ),
      expand = FALSE,
      clip = "off"
    ) +

    scale_y_reverse() +

    scale_fill_gradientn(
      colors = c(
        HH_LOW,
        HH_MID,
        HH_HIGH
      ),
      limits = c(
        0,
        HH_DISPLAY_MAX
      ),
      oob =
        scales::squish,
      name =
        "Clec7a–lipid\ncoupling",
      guide =
        guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barheight =
            unit(
              3.0,
              "cm"
            ),
          barwidth =
            unit(
              0.40,
              "cm"
            ),
          frame.colour =
            "black",
          frame.linewidth =
            0.30,
          ticks = TRUE
        )
    ) +

    facet_wrap(
      ~ display_time_single,
      ncol = 1
    ) +

    labs(
      x =
        if (show_x_title)
          "Visium array X"
        else
          NULL,

      y =
        if (show_y_title)
          "Visium array Y"
        else
          NULL
    ) +

    theme_classic(
      base_family = "Arial",
      base_size = 9
    ) +

    theme(
      strip.background =
        element_rect(
          fill = "white",
          color = "black",
          linewidth = 0.8
        ),

      strip.text =
        element_text(
          size = 8.5,
          face = "bold"
        ),

      axis.title =
        element_text(
          size = 8.6
        ),

      axis.text =
        element_text(
          size = 7.0
        ),

      legend.title =
        element_text(
          size = 7.8,
          hjust = 0.5
        ),

      legend.text =
        element_text(
          size = 7.2
        ),

      plot.margin =
        margin(
          0.6,
          0.6,
          0.6,
          0.6,
          unit = "mm"
        )
    )
}

make_state_panel <- function(dat, tp, show_x_title = TRUE, show_y_title = TRUE) {
    b <- crop_bounds %>% filter(timepoint == tp)
    if (nrow(b) != 1L) {
        stop("Expected one crop row for ", tp, ".")
    }
    dat$display_time_single <- factor(TIME_LABELS[[tp]], levels = TIME_LABELS[[tp]])
    dat$same_spot_state_plot <- factor(dat$same_spot_state, levels = c("High–High", "High–Low", "Low–High", 
        "Low–Low"))
    ggplot(dat, aes(x = array_x, y = array_y, color = same_spot_state_plot)) + geom_point(size = 1.9, 
        shape = 16, stroke = 0, alpha = 1) + coord_equal(xlim = c(b$xmin, b$xmax), ylim = c(b$ymin, b$ymax), 
        expand = FALSE, clip = "off") + scale_y_reverse() + scale_color_manual(values = STATE_COLORS, 
        drop = FALSE, name = "Same-spot\nClec7a / lipid") + facet_wrap(~display_time_single, ncol = 1) + 
        labs(x = if (show_x_title) 
            "Visium array X"
        else NULL, y = if (show_y_title) 
            "Visium array Y"
        else NULL) + theme_classic(base_family = "Arial", base_size = 9) + theme(strip.background = element_rect(fill = "white", 
        color = "black", linewidth = 0.8), strip.text = element_text(size = 8.5, face = "bold"), axis.title = element_text(size = 8.6), 
        axis.text = element_text(size = 7), legend.title = element_text(size = 7.8), legend.text = element_text(size = 7.2), 
        plot.margin = margin(0.6, 0.6, 0.6, 0.6, unit = "mm"))
}

d3 <- plot_df %>% filter(timepoint == "3dpi")

d7 <- plot_df %>% filter(timepoint == "7dpi")

d14 <- plot_df %>% filter(timepoint == "14dpi")

d28 <- plot_df %>% filter(timepoint == "28dpi")

if (any(c(nrow(d3), nrow(d7), nrow(d14), nrow(d28)) == 0L)) {
    stop("At least one representative timepoint has zero plotting spots.")
}

hh_p3 <- make_hh_panel(d3, "3dpi", FALSE, TRUE)

hh_p7 <- make_hh_panel(d7, "7dpi", FALSE, FALSE)

hh_p14 <- make_hh_panel(d14, "14dpi", TRUE, TRUE)

hh_p28 <- make_hh_panel(d28, "28dpi", TRUE, FALSE)

p_hh <- (hh_p3 + hh_p7)/(hh_p14 + hh_p28) + patchwork::plot_layout(guides = "collect", widths = c(1, 
    1), heights = c(1, 1)) & theme(legend.position = "right", legend.justification = "center")

st_p3 <- make_state_panel(d3, "3dpi", FALSE, TRUE)

st_p7 <- make_state_panel(d7, "7dpi", FALSE, FALSE)

st_p14 <- make_state_panel(d14, "14dpi", TRUE, TRUE)

st_p28 <- make_state_panel(d28, "28dpi", TRUE, FALSE)

p_state <- (st_p3 + st_p7)/(st_p14 + st_p28) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right", 
    legend.justification = "center")

save_three_formats <- function(plot_object, stem, width = 8.8, height = 5) {
    pdf_file <- file.path(OUT_DIR, paste0(stem, ".pdf"))
    png_file <- file.path(OUT_DIR, paste0(stem, ".png"))
    svg_file <- file.path(OUT_DIR, paste0(stem, ".svg"))
    ggsave(filename = pdf_file, plot = plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    ggsave(filename = png_file, plot = plot_object, width = width, height = height, units = "in", dpi = 600, 
        bg = "white")
    ggsave(filename = svg_file, plot = plot_object, width = width, height = height, units = "in", device = "svg", 
        bg = "white")
    c(pdf = pdf_file, png = png_file, svg = svg_file)
}

HH_FILES <- save_three_formats(p_hh, paste0("05_74_Figure4E_same_spot_HH_coupling_", VERSION))

STATE_FILES <- save_three_formats(p_state, paste0("06_74_same_spot_four_state_map_", VERSION))

valid_states <- c("High–High", "High–Low", "Low–High", "Low–Low")

selected_match_ok <- setequal(rep_sections$section_uid_internal, unique(plot_df$section_uid_internal))

all_sections_big_enough <- all(section_summary$n_spots >= MIN_SECTION_SPOTS)

all_states_valid <- all(spot$same_spot_state %in% valid_states)

hh_exists_all_timepoints <- all(sapply(TIME_LEVELS, function(tp) {
    any(plot_df$timepoint == tp & plot_df$same_spot_state == "High–High")
}))

final_qc <- data.frame(check = c("Step73 v03 spot-score file was used", "Figure4A v12 representative-section file was used", 
    "All 8 capture areas were retained", "All 32 coordinate-defined major sections were analyzed", "All analyzed sections contain at least 100 spots", 
    "Section-wise Clec7a z-scores are finite", "Section-wise lipid-program z-scores are finite", "All same-spot state labels are valid", 
    "High–High means z_Clec7a > 0 AND z_lipid > 0 in the same spot", "The same four representative sections as Figure4A/B were used", 
    "Each representative timepoint contains at least one HH spot", "No spot-level P values were computed", 
    "No H&E-defined lesion compartment is inferred", "Main HH PDF generated", "Main HH PNG generated", 
    "Main HH SVG generated", "Four-state PDF generated", "Four-state PNG generated", "Four-state SVG generated"), 
    pass = c(identical(basename(INPUT_SPOT_FILE), "02_73_Clec7a_lipid_spot_scores_section_aware_v03.csv.gz"), 
        identical(basename(REP_SECTION_FILE), "02_69_Figure4A_representative_sections_v12.csv"), n_capture_areas == 
            EXPECTED_CAPTURE_AREAS, n_sections == EXPECTED_MAJOR_SECTIONS, all_sections_big_enough, all(is.finite(spot$z_Clec7a_same_spot)), 
        all(is.finite(spot$z_lipid_same_spot)), all_states_valid, TRUE, selected_match_ok, hh_exists_all_timepoints, 
        TRUE, TRUE, file.exists(HH_FILES[["pdf"]]), file.exists(HH_FILES[["png"]]), file.exists(HH_FILES[["svg"]]), 
        file.exists(STATE_FILES[["pdf"]]), file.exists(STATE_FILES[["png"]]), file.exists(STATE_FILES[["svg"]])), 
    stringsAsFactors = FALSE)

QC_FILE <- file.path(OUT_DIR, paste0("07_74_same_spot_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_qc, QC_FILE)

REPORT_FILE <- file.path(OUT_DIR, paste0("08_74_same_spot_REPORT_", VERSION, ".txt"))

time_lines <- apply(time_summary, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": n_sections=", z[["n_sections"]], "; median section rho=", sprintf("%.3f", 
        as.numeric(z[["median_section_rho"]])), "; median HH fraction=", sprintf("%.3f", as.numeric(z[["median_fraction_HH"]])), 
        "; median HH strength=", sprintf("%.3f", as.numeric(z[["median_HH_strength"]])))
})

rep_lines <- apply(rep_sections, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": ", z[["section_uid_internal"]], " [", z[["sample_label"]], "]")
})

writeLines("Numerical results are provided in the accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(time_summary)

print(rep_sections)

print(final_qc)

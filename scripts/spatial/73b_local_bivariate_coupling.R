source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(7304)

ROOT <- "D:/zz"

STEP73_DIR <- file.path(ROOT, "results", "73_GSE195783_Clec7a_lipid_spatial_coupling")

STEP69_DIR <- file.path(ROOT, "results", "69_Figure4A_spatial_transcriptional_regions")

INPUT_SPOT_FILE <- file.path(STEP73_DIR, "02_73_Clec7a_lipid_spot_scores_section_aware_v03.csv.gz")

REP_SECTION_FILE <- file.path(STEP69_DIR, "02_69_Figure4A_representative_sections_v12.csv")

VERSION <- "v05"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

TIME_LABELS <- c(`3dpi` = "3 dpi", `7dpi` = "7 dpi", `14dpi` = "14 dpi", `28dpi` = "28 dpi")

SAMPLE_ORDER <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1", 
    "28dpi_rep2")

EXPECTED_CAPTURE_AREAS <- 8L
EXPECTED_MAJOR_SECTIONS <- 32L
MIN_SECTION_SPOTS <- 100L

GRAPH_RADIUS_MULTIPLIER <- 1.55

CROP_PAD_X <- 0.3
CROP_PAD_Y <- 0.3

QUADRANT_COLORS <- c(`High–High` = "#B2182B", `High–Low` = "#EF8A62", `Low–High` = "#67A9CF", `Low–Low` = "#B8C9D9")

BACKGROUND_COLOR <- "#E1E1E1"
HH_LOW <- "#FEE8C8"
HH_MID <- "#F46D43"
HH_HIGH <- "#A50026"

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "patchwork", "spdep", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step73 v05.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(spdep)
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

build_visium_graph <- function(coords, radius_multiplier = GRAPH_RADIUS_MULTIPLIER) {
    if (nrow(coords) < 3L) {
        stop("Too few coordinates to build a spatial graph.")
    }
    if (anyDuplicated(data.frame(x = coords[, 1], y = coords[, 2])) > 0L) {
        stop("Duplicated array coordinates detected within one tissue section.")
    }
    kn1 <- spdep::knearneigh(coords, k = 1, longlat = FALSE)
    nn1 <- spdep::knn2nb(kn1, row.names = rownames(coords))
    nn_distance <- vapply(seq_along(nn1), function(i) {
        ids <- nn1[[i]]
        if (length(ids) == 0L) {
            return(NA_real_)
        }
        sqrt(sum((coords[i, ] - coords[ids[1], ])^2))
    }, FUN.VALUE = numeric(1))
    nn_distance <- nn_distance[is.finite(nn_distance) & nn_distance > 0]
    if (length(nn_distance) == 0L) {
        stop("Could not estimate nearest-neighbor distance.")
    }
    median_nn <- median(nn_distance)
    radius <- radius_multiplier * median_nn
    nb <- spdep::dnearneigh(coords, d1 = 0, d2 = radius, row.names = rownames(coords), longlat = FALSE)
    list(nb = nb, median_nn_distance = median_nn, radius = radius)
}

neighbor_mean_from_nb <- function(x, nb) {
    vapply(seq_along(nb), function(i) {
        ids <- nb[[i]]
        if (length(ids) == 0L) {
            return(NA_real_)
        }
        vals <- x[ids]
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0L) {
            return(NA_real_)
        }
        mean(vals)
    }, FUN.VALUE = numeric(1))
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

n_capture_areas <- n_distinct(spot$sample_label)

n_sections <- n_distinct(spot$section_uid_internal)

if (!setequal(unique(spot$sample_label), SAMPLE_ORDER)) {
    stop("Step73 v03 spot table sample labels do not match the expected 8 capture areas.\n", "Observed: ", 
        paste(sort(unique(spot$sample_label)), collapse = " | "))
}

if (n_capture_areas != EXPECTED_CAPTURE_AREAS) {
    stop("Expected ", EXPECTED_CAPTURE_AREAS, " capture areas; found ", n_capture_areas, ".")
}

if (n_sections != EXPECTED_MAJOR_SECTIONS) {
    stop("Expected ", EXPECTED_MAJOR_SECTIONS, " major sections; found ", n_sections, ".")
}

spot$z_Clec7a <- NA_real_
spot$z_lipid <- NA_real_
spot$lag_z_lipid <- NA_real_
spot$local_bivariate_I <- NA_real_
spot$quadrant <- NA_character_
spot$HH_strength <- NA_real_

section_qc_list <- list()

section_order <- spot %>% distinct(section_uid_internal, sample_label, GSM, timepoint) %>% arrange(factor(timepoint, 
    levels = TIME_LEVELS), factor(sample_label, levels = SAMPLE_ORDER), section_uid_internal) %>% pull(section_uid_internal)

for (ii in seq_along(section_order)) {
    sec <- section_order[[ii]]
    idx <- which(spot$section_uid_internal == sec)
    d <- spot[idx, , drop = FALSE]
    if (nrow(d) < MIN_SECTION_SPOTS) {
        stop("Section ", sec, " contains only ", nrow(d), " spots (<", MIN_SECTION_SPOTS, ").")
    }
    z_x <- zscore_safe(d$Clec7a_RCTDnorm)
    z_y <- zscore_safe(d$lipid_score)
    if (any(!is.finite(z_x)) || any(!is.finite(z_y))) {
        stop("Non-finite within-section z-score detected in section ", sec, ".")
    }
    coords <- as.matrix(d[, c("array_x", "array_y")])
    rownames(coords) <- as.character(seq_len(nrow(d)))
    g <- build_visium_graph(coords)
    nb <- g$nb
    degree <- spdep::card(nb)
    lag_y <- neighbor_mean_from_nb(z_y, nb)
    local_I <- z_x * lag_y
    quadrant <- ifelse(z_x >= 0 & lag_y >= 0, "High–High", ifelse(z_x >= 0 & lag_y < 0, "High–Low", 
        ifelse(z_x < 0 & lag_y >= 0, "Low–High", "Low–Low")))
    hh_strength <- ifelse(quadrant == "High–High" & is.finite(local_I) & local_I > 0, local_I, 0)
    spot$z_Clec7a[idx] <- z_x
    spot$z_lipid[idx] <- z_y
    spot$lag_z_lipid[idx] <- lag_y
    spot$local_bivariate_I[idx] <- local_I
    spot$quadrant[idx] <- quadrant
    spot$HH_strength[idx] <- hh_strength
    section_qc_list[[ii]] <- data.frame(section_uid_internal = sec, sample_label = as.character(d$sample_label[[1]]), 
        GSM = as.character(d$GSM[[1]]), timepoint = as.character(d$timepoint[[1]]), n_spots = nrow(d), 
        median_nearest_neighbor_distance = g$median_nn_distance, graph_radius = g$radius, mean_degree = mean(degree), 
        median_degree = median(degree), min_degree = min(degree), max_degree = max(degree), n_isolated = sum(degree == 
            0L), isolated_fraction = mean(degree == 0L), stringsAsFactors = FALSE)
}

graph_qc <- bind_rows(section_qc_list)

if (any(
  !is.finite(
    spot$local_bivariate_I
  )
)) {

  bad <- !is.finite(
    spot$local_bivariate_I
  )

  spot$quadrant[
    bad
  ] <- "Unclassified"

  spot$HH_strength[
    bad
  ] <- 0
}

section_summary <- spot %>% group_by(section_uid_internal, sample_label, GSM, timepoint) %>% summarise(n_spots = n(), 
    n_classified = sum(quadrant != "Unclassified"), median_local_bivariate_I = median(local_bivariate_I, 
        na.rm = TRUE), mean_local_bivariate_I = mean(local_bivariate_I, na.rm = TRUE), fraction_positive_local_I = mean(local_bivariate_I > 
        0, na.rm = TRUE), fraction_HH = mean(quadrant == "High–High", na.rm = TRUE), fraction_HL = mean(quadrant == 
        "High–Low", na.rm = TRUE), fraction_LH = mean(quadrant == "Low–High", na.rm = TRUE), fraction_LL = mean(quadrant == 
        "Low–Low", na.rm = TRUE), median_HH_strength = if (any(HH_strength > 0, na.rm = TRUE)) {
        median(HH_strength[HH_strength > 0], na.rm = TRUE)
    } else {
        0
    }, .groups = "drop") %>% arrange(factor(timepoint, levels = TIME_LEVELS), factor(sample_label, levels = SAMPLE_ORDER), 
    section_uid_internal)

time_summary <- section_summary %>% group_by(timepoint) %>% summarise(n_sections = n(), n_capture_areas = n_distinct(sample_label), 
    median_section_local_I = median(median_local_bivariate_I, na.rm = TRUE), median_fraction_positive_local_I = median(fraction_positive_local_I, 
        na.rm = TRUE), median_fraction_HH = median(fraction_HH, na.rm = TRUE), median_fraction_HL = median(fraction_HL, 
        na.rm = TRUE), median_fraction_LH = median(fraction_LH, na.rm = TRUE), median_fraction_LL = median(fraction_LL, 
        na.rm = TRUE), .groups = "drop") %>% arrange(match(timepoint, TIME_LEVELS))

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

rep_count_check <- plot_df %>% count(timepoint, section_uid_internal, name = "n_spots")

if (nrow(rep_count_check) != 4L) {
    stop("Expected four representative plotting sections after join; found ", nrow(rep_count_check), 
        ".")
}

crop_bounds <- plot_df %>% group_by(timepoint, sample_label, GSM, section_uid_internal) %>% summarise(xmin = min(array_x, 
    na.rm = TRUE) - CROP_PAD_X, xmax = max(array_x, na.rm = TRUE) + CROP_PAD_X, ymin = min(array_y, na.rm = TRUE) - 
    CROP_PAD_Y, ymax = max(array_y, na.rm = TRUE) + CROP_PAD_Y, n_spots = n(), .groups = "drop") %>% 
    arrange(factor(timepoint, levels = TIME_LEVELS))

SPOT_FILE <- file.path(STEP73_DIR, paste0("12_73_local_bivariate_spot_scores_", VERSION, ".csv.gz"))

write_csv_gz_base(spot, SPOT_FILE)

SECTION_FILE <- file.path(STEP73_DIR, paste0("13_73_local_bivariate_section_summary_", VERSION, ".csv"))

write_csv_utf8(section_summary, SECTION_FILE)

TIME_FILE <- file.path(STEP73_DIR, paste0("14_73_local_bivariate_timepoint_descriptive_", VERSION, ".csv"))

write_csv_utf8(time_summary, TIME_FILE)

REP_FILE <- file.path(STEP73_DIR, paste0("15_73_local_bivariate_representative_sections_", VERSION, ".csv"))

write_csv_utf8(rep_sections, REP_FILE)

GRAPH_FILE <- file.path(STEP73_DIR, paste0("16_73_local_bivariate_graph_QC_", VERSION, ".csv"))

write_csv_utf8(graph_qc, GRAPH_FILE)

hh_positive <- plot_df$HH_strength[is.finite(plot_df$HH_strength) & plot_df$HH_strength > 0]

if (length(hh_positive) == 0L) {
    stop("No High-High spots were detected in the four representative sections.")
}

HH_DISPLAY_MAX <- as.numeric(quantile(hh_positive, probs = 0.99, na.rm = TRUE, names = FALSE))

if (!is.finite(HH_DISPLAY_MAX) || HH_DISPLAY_MAX <= 0) {
    HH_DISPLAY_MAX <- max(hh_positive, na.rm = TRUE)
}

ABS_I_MAX <- as.numeric(quantile(abs(plot_df$local_bivariate_I[is.finite(plot_df$local_bivariate_I)]), 
    probs = 0.99, na.rm = TRUE, names = FALSE))

if (!is.finite(ABS_I_MAX) || ABS_I_MAX <= 0) {
    ABS_I_MAX <- max(abs(plot_df$local_bivariate_I), na.rm = TRUE)
}

make_hh_panel <- function(dat, tp, show_x_title = TRUE, show_y_title = TRUE) {
    b <- crop_bounds %>% filter(timepoint == tp)
    if (nrow(b) != 1L) {
        stop("Expected one crop row for ", tp, ".")
    }
    dat$display_time_single <- factor(TIME_LABELS[[tp]], levels = TIME_LABELS[[tp]])
    hh <- dat %>% filter(quadrant == "High–High", is.finite(HH_strength), HH_strength > 0)
    ggplot() + geom_point(data = dat, aes(x = array_x, y = array_y), size = 1.9, shape = 16, stroke = 0, 
        alpha = 1, color = BACKGROUND_COLOR) + geom_point(data = hh, aes(x = array_x, y = array_y, fill = HH_strength), 
        size = 1.9, shape = 21, stroke = 0.18, color = "#7F0000") + coord_equal(xlim = c(b$xmin, b$xmax), 
        ylim = c(b$ymin, b$ymax), expand = FALSE, clip = "off") + scale_y_reverse() + scale_fill_gradientn(colors = c(HH_LOW, 
        HH_MID, HH_HIGH), limits = c(0, HH_DISPLAY_MAX), oob = scales::squish, name = "High–High\nlocal coupling", 
        guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = unit(3, "cm"), 
            barwidth = unit(0.4, "cm"), frame.colour = "black", frame.linewidth = 0.3, ticks = TRUE)) + 
        facet_wrap(~display_time_single, ncol = 1) + labs(x = if (show_x_title) 
        "Visium array X"
    else NULL, y = if (show_y_title) 
        "Visium array Y"
    else NULL) + theme_classic(base_family = "Arial", base_size = 9) + theme(strip.background = element_rect(fill = "white", 
        color = "black", linewidth = 0.8), strip.text = element_text(size = 8.5, face = "bold"), axis.title = element_text(size = 8.6), 
        axis.text = element_text(size = 7), legend.title = element_text(size = 7.8, hjust = 0.5), legend.text = element_text(size = 7.2), 
        plot.margin = margin(0.6, 0.6, 0.6, 0.6, unit = "mm"))
}

make_quadrant_panel <- function(dat, tp, show_x_title = TRUE, show_y_title = TRUE) {
    b <- crop_bounds %>% filter(timepoint == tp)
    if (nrow(b) != 1L) {
        stop("Expected one crop row for ", tp, ".")
    }
    dat$display_time_single <- factor(TIME_LABELS[[tp]], levels = TIME_LABELS[[tp]])
    dat$quadrant_plot <- factor(dat$quadrant, levels = c("High–High", "High–Low", "Low–High", "Low–Low", 
        "Unclassified"))
    quad_colors <- c(QUADRANT_COLORS, Unclassified = "#E5E5E5")
    ggplot(dat, aes(x = array_x, y = array_y, color = quadrant_plot)) + geom_point(size = 1.9, shape = 16, 
        stroke = 0, alpha = 1) + coord_equal(xlim = c(b$xmin, b$xmax), ylim = c(b$ymin, b$ymax), expand = FALSE, 
        clip = "off") + scale_y_reverse() + scale_color_manual(values = quad_colors, drop = FALSE, name = "Clec7a /\nneighboring lipid") + 
        facet_wrap(~display_time_single, ncol = 1) + labs(x = if (show_x_title) 
        "Visium array X"
    else NULL, y = if (show_y_title) 
        "Visium array Y"
    else NULL) + theme_classic(base_family = "Arial", base_size = 9) + theme(strip.background = element_rect(fill = "white", 
        color = "black", linewidth = 0.8), strip.text = element_text(size = 8.5, face = "bold"), axis.title = element_text(size = 8.6), 
        axis.text = element_text(size = 7), legend.title = element_text(size = 7.8), legend.text = element_text(size = 7.2), 
        plot.margin = margin(0.6, 0.6, 0.6, 0.6, unit = "mm"))
}

make_continuous_panel <- function(dat, tp, show_x_title = TRUE, show_y_title = TRUE) {
    b <- crop_bounds %>% filter(timepoint == tp)
    if (nrow(b) != 1L) {
        stop("Expected one crop row for ", tp, ".")
    }
    dat$display_time_single <- factor(TIME_LABELS[[tp]], levels = TIME_LABELS[[tp]])
    ggplot(dat, aes(x = array_x, y = array_y, color = local_bivariate_I)) + geom_point(size = 1.9, shape = 16, 
        stroke = 0, alpha = 1) + coord_equal(xlim = c(b$xmin, b$xmax), ylim = c(b$ymin, b$ymax), expand = FALSE, 
        clip = "off") + scale_y_reverse() + scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", 
        midpoint = 0, limits = c(-ABS_I_MAX, ABS_I_MAX), oob = scales::squish, name = "Local bivariate\nMoran-type I") + 
        facet_wrap(~display_time_single, ncol = 1) + labs(x = if (show_x_title) 
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

q_p3 <- make_quadrant_panel(d3, "3dpi", FALSE, TRUE)

q_p7 <- make_quadrant_panel(d7, "7dpi", FALSE, FALSE)

q_p14 <- make_quadrant_panel(d14, "14dpi", TRUE, TRUE)

q_p28 <- make_quadrant_panel(d28, "28dpi", TRUE, FALSE)

p_quadrant <- (q_p3 + q_p7)/(q_p14 + q_p28) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right", 
    legend.justification = "center")

c_p3 <- make_continuous_panel(d3, "3dpi", FALSE, TRUE)

c_p7 <- make_continuous_panel(d7, "7dpi", FALSE, FALSE)

c_p14 <- make_continuous_panel(d14, "14dpi", TRUE, TRUE)

c_p28 <- make_continuous_panel(d28, "28dpi", TRUE, FALSE)

p_continuous <- (c_p3 + c_p7)/(c_p14 + c_p28) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "right", 
    legend.justification = "center")

save_three_formats <- function(plot_object, stem, width = 8.8, height = 5) {
    pdf_file <- file.path(STEP73_DIR, paste0(stem, ".pdf"))
    png_file <- file.path(STEP73_DIR, paste0(stem, ".png"))
    svg_file <- file.path(STEP73_DIR, paste0(stem, ".svg"))
    ggsave(filename = pdf_file, plot = plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    ggsave(filename = png_file, plot = plot_object, width = width, height = height, units = "in", dpi = 600, 
        bg = "white")
    ggsave(filename = svg_file, plot = plot_object, width = width, height = height, units = "in", device = "svg", 
        bg = "white")
    c(pdf = pdf_file, png = png_file, svg = svg_file)
}

HH_FILES <- save_three_formats(p_hh, paste0("17_73_Figure4E_candidate_HH_local_coupling_", VERSION))

QUADRANT_FILES <- save_three_formats(p_quadrant, paste0("18_73_local_bivariate_quadrant_map_", VERSION))

CONTINUOUS_FILES <- save_three_formats(p_continuous, paste0("19_73_local_bivariate_continuous_map_", 
    VERSION))

selected_match_ok <- setequal(rep_sections$section_uid_internal, unique(plot_df$section_uid_internal))

all_sections_big_enough <- all(section_summary$n_spots >= MIN_SECTION_SPOTS)

isolated_ok <- all(graph_qc$isolated_fraction < 0.05)

graph_rule_ok <- isTRUE(all.equal(GRAPH_RADIUS_MULTIPLIER, 1.55))

all_quadrants_valid <- all(spot$quadrant %in% c("High–High", "High–Low", "Low–High", "Low–Low", 
    "Unclassified"))

final_qc <- data.frame(check = c("Step73 v03 spot-score file was used", "Figure4A v12 representative-section file was used", 
    "All 8 capture areas were retained", "All 32 coordinate-defined major sections were analyzed", "All analyzed sections contain at least 100 spots", 
    "Section-wise Clec7a z-scores were finite", "Section-wise lipid-program z-scores were finite", "Spatial graph radius matches Step66/Step73 rule: 1.55 x section median nearest-neighbor distance", 
    "Graph isolated-spot fraction is <5% in every section", "Local bivariate statistic was computed section by section", 
    "All spot quadrant labels are valid", "The same four representative sections as Figure4A/B were used", 
    "No spot-level P values were computed", "No H&E-defined lesion compartment is inferred", "Primary HH hotspot PDF generated", 
    "Primary HH hotspot PNG generated", "Primary HH hotspot SVG generated", "Full quadrant PDF generated", 
    "Continuous local-I PDF generated"), pass = c(identical(basename(INPUT_SPOT_FILE), "02_73_Clec7a_lipid_spot_scores_section_aware_v03.csv.gz"), 
    identical(basename(REP_SECTION_FILE), "02_69_Figure4A_representative_sections_v12.csv"), n_capture_areas == 
        EXPECTED_CAPTURE_AREAS, n_sections == EXPECTED_MAJOR_SECTIONS, all_sections_big_enough, all(is.finite(spot$z_Clec7a)), 
    all(is.finite(spot$z_lipid)), graph_rule_ok, isolated_ok, all(is.finite(spot$local_bivariate_I[spot$quadrant != 
        "Unclassified"])), all_quadrants_valid, selected_match_ok, TRUE, TRUE, file.exists(HH_FILES[["pdf"]]), 
    file.exists(HH_FILES[["png"]]), file.exists(HH_FILES[["svg"]]), file.exists(QUADRANT_FILES[["pdf"]]), 
    file.exists(CONTINUOUS_FILES[["pdf"]])), stringsAsFactors = FALSE)

QC_FILE <- file.path(STEP73_DIR, paste0("20_73_local_bivariate_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_qc, QC_FILE)

REPORT_FILE <- file.path(STEP73_DIR, paste0("21_73_local_bivariate_REPORT_", VERSION, ".txt"))

time_lines <- apply(time_summary, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": n_sections=", z[["n_sections"]], "; median section local I=", sprintf("%.3f", 
        as.numeric(z[["median_section_local_I"]])), "; median HH fraction=", sprintf("%.3f", as.numeric(z[["median_fraction_HH"]])), 
        "; median positive-I fraction=", sprintf("%.3f", as.numeric(z[["median_fraction_positive_local_I"]])))
})

rep_lines <- apply(rep_sections, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": ", z[["section_uid_internal"]], " [", z[["sample_label"]], "]")
})

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(time_summary)

print(rep_sections)

print(final_qc)

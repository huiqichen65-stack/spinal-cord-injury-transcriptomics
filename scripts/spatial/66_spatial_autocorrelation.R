rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)


options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(6602)

ROOT <- PROJECT_ROOT

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

INPUT_FILE <- file.path(STEP65_DIR, "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv.gz")

OUT_DIR <- file.path(ROOT, "results", "66_GSE195783_Clec7a_spatial_autocorrelation_neighborhood")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v02"

TARGET_STATE <- "Clec7a+ lipid-handling Microglia"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

SAMPLE_ORDER <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1", 
    "28dpi_rep2")

SECTION_POSITION_LEVELS <- c("upper_left", "upper_right", "lower_left", "lower_right")

FOCUS_REGIONS <- c("Myeloid/immune–lipid-handling enriched", "Fibrotic/ECM–myeloid mixed injury niche", 
    "IFN-responsive inflammatory/astroglial region")

EXPECTED_STEP65_SPOTS <- 15204L
EXPECTED_CAPTURE_AREAS <- 8L
EXPECTED_MAJOR_SECTIONS <- 32L

MIN_SECTION_SPOTS <- 100L
MIN_REGION_SPOTS <- 10L
MIN_OTHER_SPOTS <- 10L

MORAN_NSIM <- 499L
GRAPH_RADIUS_MULTIPLIER <- 1.55

EPS <- 1e-8

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "spdep", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step66 v02.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(spdep)
  library(scales)
})

write_csv_utf8 <- write_csv_utf8

read_csv_gz_base <- read_gzip_csv

write_csv_gz_base <- write_gzip_csv

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

safe_wilcox_greater_zero <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 4L) {
        return(NA_real_)
    }
    if (all(abs(x) < .Machine$double.eps^0.5)) {
        return(1)
    }
    suppressWarnings(stats::wilcox.test(x, mu = 0, alternative = "greater", exact = FALSE, conf.int = FALSE)$p.value)
}

safe_fraction_positive <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    mean(x > 0)
}

safe_median <- finite_median

build_visium_graph <- function(coords, radius_multiplier = GRAPH_RADIUS_MULTIPLIER) {
    if (nrow(coords) < 3L) {
        stop("Too few coordinates to build a spatial graph.")
    }
    if (anyDuplicated(data.frame(x = coords[, 1], y = coords[, 2])) > 0L) {
        stop("Duplicated array coordinates detected within a tissue section.")
    }
    kn1 <- spdep::knearneigh(coords, k = 1, longlat = FALSE)
    nn_index <- as.integer(kn1$nn[, 1])
    nn_distance <- sqrt(rowSums((coords - coords[nn_index, , drop = FALSE])^2))
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
        mean(x[ids], na.rm = TRUE)
    }, numeric(1))
}

summarize_section_focus <- function(effect_df, effect_col, ratio_col, analysis_label) {
    rows <- list()
    for (reg in FOCUS_REGIONS) {
        z <- effect_df[effect_df$spatial_region == reg, , drop = FALSE]
        x <- z[[effect_col]]
        rows[[reg]] <- data.frame(analysis_level = "section_level_exploratory", analysis = analysis_label, 
            spatial_region = reg, n_sections_tested = sum(is.finite(x)), n_capture_areas_represented = length(unique(z$sample_label[is.finite(x)])), 
            median_log2_effect = safe_median(x), median_enrichment_ratio = ifelse(is.finite(safe_median(x)), 
                2^safe_median(x), NA_real_), fraction_sections_effect_gt_1 = safe_fraction_positive(x), 
            exploratory_wilcox_p_greater_zero = safe_wilcox_greater_zero(x), exploratory_FDR_BH = NA_real_, 
            stringsAsFactors = FALSE)
    }
    out <- bind_rows(rows)
    out$exploratory_FDR_BH <- p.adjust(out$exploratory_wilcox_p_greater_zero, method = "BH")
    out
}

build_blocked_effects <- function(effect_df, effect_col, analysis_label) {
    z <- effect_df %>% filter(spatial_region %in% FOCUS_REGIONS) %>% group_by(sample_label, timepoint, 
        spatial_region) %>% summarise(analysis = analysis_label, n_sections_contributing = sum(is.finite(.data[[effect_col]])), 
        blocked_log2_effect = if (sum(is.finite(.data[[effect_col]])) > 0L) {
            median(.data[[effect_col]][is.finite(.data[[effect_col]])])
        }
        else {
            NA_real_
        }, .groups = "drop")
    z
}

summarize_blocked_focus <- function(blocked_df, analysis_label) {
    rows <- list()
    for (reg in FOCUS_REGIONS) {
        z <- blocked_df[blocked_df$spatial_region == reg, , drop = FALSE]
        x <- z$blocked_log2_effect
        rows[[reg]] <- data.frame(analysis_level = "capture_area_blocked_sensitivity", analysis = analysis_label, 
            spatial_region = reg, n_capture_areas_tested = sum(is.finite(x)), n_sections_total_contributing = sum(z$n_sections_contributing[is.finite(x)], 
                na.rm = TRUE), median_blocked_log2_effect = safe_median(x), median_blocked_enrichment_ratio = ifelse(is.finite(safe_median(x)), 
                2^safe_median(x), NA_real_), fraction_capture_areas_effect_gt_1 = safe_fraction_positive(x), 
            wilcox_p_greater_zero = safe_wilcox_greater_zero(x), FDR_BH = NA_real_, stringsAsFactors = FALSE)
    }
    out <- bind_rows(rows)
    out$FDR_BH <- p.adjust(out$wilcox_p_greater_zero, method = "BH")
    out
}

if (!file.exists(INPUT_FILE)) {
    stop("Step65 v06 section-aware spot-weight file not found:\n", INPUT_FILE)
}

md_all <- read_csv_gz_base(INPUT_FILE)

required_meta <- c("barcode", "GSM", "sample_label", "timepoint", "spatial_region_v03", "array_x", "array_y", 
    TARGET_STATE, "section_eligible", "section_position_internal", "section_uid_internal")

missing_meta <- setdiff(required_meta, colnames(md_all))

if (length(missing_meta) > 0L) {
    stop("Step65 v06 input is missing required column(s): ", paste(missing_meta, collapse = " | "))
}

if (nrow(md_all) != EXPECTED_STEP65_SPOTS) {
    stop("Expected ", EXPECTED_STEP65_SPOTS, " Step65 spots, but found ", nrow(md_all), ".")
}

md_all$section_eligible <- as_logical_safe(md_all$section_eligible)

if (any(is.na(md_all$section_eligible))) {
    stop("Could not safely parse section_eligible.")
}

md <- md_all %>% filter(section_eligible)

if (any(is.na(md$section_uid_internal)) | any(!nzchar(trimws(as.character(md$section_uid_internal))))) {
    stop("At least one section-eligible spot lacks section_uid_internal.")
}

if (any(!is.finite(md$array_x)) | any(!is.finite(md$array_y))) {
    stop("Section-eligible spots contain non-finite array coordinates.")
}

target_weight <- as.numeric(md[[TARGET_STATE]])

if (any(!is.finite(target_weight))) {
    stop("Target normalized RCTD weight contains non-finite values.")
}

if (any(target_weight < -1e-10 | target_weight > 1 + 1e-10)) {
    stop("Target normalized RCTD weight is outside [0,1].")
}

md$Clec7a_RCTDnorm <- target_weight

n_capture_areas <- length(unique(md$sample_label))

n_sections <- length(unique(md$section_uid_internal))

if (n_capture_areas != EXPECTED_CAPTURE_AREAS) {
    stop("Expected ", EXPECTED_CAPTURE_AREAS, " capture areas; found ", n_capture_areas, ".")
}

if (n_sections != EXPECTED_MAJOR_SECTIONS) {
    stop("Expected ", EXPECTED_MAJOR_SECTIONS, " coordinate-defined tissue sections; found ", n_sections, 
        ".")
}

sections_per_capture <- md %>% distinct(sample_label, section_uid_internal) %>% count(sample_label, name = "n_sections")

if (any(sections_per_capture$n_sections != 4L)) {
    stop("Each capture area must contribute exactly four coordinate-defined sections.")
}

section_info <- md %>% distinct(section_uid_internal, sample_label, GSM, timepoint, section_position_internal) %>% 
    mutate(timepoint = factor(timepoint, levels = TIME_LEVELS), sample_label = factor(sample_label, levels = SAMPLE_ORDER), 
        section_position_internal = factor(section_position_internal, levels = SECTION_POSITION_LEVELS)) %>% 
    arrange(timepoint, sample_label, section_position_internal)

SECTION_ORDER <- as.character(section_info$section_uid_internal)

neighbor_score_all <- rep(NA_real_, nrow(md))

names(neighbor_score_all) <- as.character(seq_len(nrow(md)))

graph_qc_list <- list()
moran_list <- list()

for (i in seq_along(SECTION_ORDER)) {
    sec <- SECTION_ORDER[i]
    idx <- which(md$section_uid_internal == sec)
    d <- md[idx, , drop = FALSE]
    if (nrow(d) < MIN_SECTION_SPOTS) {
        stop(sec, ": fewer than ", MIN_SECTION_SPOTS, " spots.")
    }
    coords <- as.matrix(d[, c("array_x", "array_y"), drop = FALSE])
    rownames(coords) <- as.character(idx)
    g <- build_visium_graph(coords)
    nb <- g$nb
    degree <- spdep::card(nb)
    listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    x <- d$Clec7a_RCTDnorm
    if (!is.finite(stats::sd(x)) || stats::sd(x) <= 0) {
        stop(sec, ": target weight has zero/non-finite variance.")
    }
    set.seed(660200L + i)
    moran_mc <- spdep::moran.mc(x, listw = listw, nsim = MORAN_NSIM, alternative = "greater", zero.policy = TRUE, 
        na.action = na.fail)
    neigh_mean <- neighbor_mean_from_nb(x, nb)
    neighbor_score_all[idx] <- neigh_mean
    graph_qc_list[[sec]] <- data.frame(section_uid_internal = sec, sample_label = as.character(unique(d$sample_label)[1]), 
        GSM = as.character(unique(d$GSM)[1]), timepoint = as.character(unique(d$timepoint)[1]), section_position_internal = as.character(unique(d$section_position_internal)[1]), 
        n_spots = nrow(d), median_nearest_neighbor_distance = g$median_nn_distance, graph_radius = g$radius, 
        mean_degree = mean(degree), median_degree = median(degree), min_degree = min(degree), max_degree = max(degree), 
        n_isolated = sum(degree == 0L), fraction_isolated = mean(degree == 0L), stringsAsFactors = FALSE)
    moran_list[[sec]] <- data.frame(section_uid_internal = sec, sample_label = as.character(unique(d$sample_label)[1]), 
        GSM = as.character(unique(d$GSM)[1]), timepoint = as.character(unique(d$timepoint)[1]), section_position_internal = as.character(unique(d$section_position_internal)[1]), 
        n_spots = nrow(d), Moran_I = as.numeric(moran_mc$statistic), permutation_p = as.numeric(moran_mc$p.value), 
        nsim = MORAN_NSIM, alternative = "greater", stringsAsFactors = FALSE)
}

graph_qc <- bind_rows(graph_qc_list)

moran_results <- bind_rows(moran_list)

moran_results$FDR_BH <- p.adjust(moran_results$permutation_p, method = "BH")

graph_qc$timepoint <- factor(graph_qc$timepoint, levels = TIME_LEVELS)

graph_qc$sample_label <- factor(graph_qc$sample_label, levels = SAMPLE_ORDER)

graph_qc$section_position_internal <- factor(graph_qc$section_position_internal, levels = SECTION_POSITION_LEVELS)

graph_qc <- graph_qc %>% arrange(timepoint, sample_label, section_position_internal)

graph_qc$timepoint <- as.character(graph_qc$timepoint)

graph_qc$sample_label <- as.character(graph_qc$sample_label)

graph_qc$section_position_internal <- as.character(graph_qc$section_position_internal)

moran_results$timepoint <- factor(moran_results$timepoint, levels = TIME_LEVELS)

moran_results$sample_label <- factor(moran_results$sample_label, levels = SAMPLE_ORDER)

moran_results$section_position_internal <- factor(moran_results$section_position_internal, levels = SECTION_POSITION_LEVELS)

moran_results <- moran_results %>% arrange(timepoint, sample_label, section_position_internal)

moran_results$timepoint <- as.character(moran_results$timepoint)

moran_results$sample_label <- as.character(moran_results$sample_label)

moran_results$section_position_internal <- as.character(moran_results$section_position_internal)

md$Clec7a_neighbor_mean_v02 <- neighbor_score_all

moran_time_summary <- moran_results %>% group_by(timepoint) %>% summarise(n_sections = n(), n_capture_areas = n_distinct(sample_label), 
    median_Moran_I = median(Moran_I), min_Moran_I = min(Moran_I), max_Moran_I = max(Moran_I), n_FDR_lt_0_05 = sum(FDR_BH < 
        0.05, na.rm = TRUE), .groups = "drop")

moran_time_summary$timepoint <- factor(moran_time_summary$timepoint, levels = TIME_LEVELS)

moran_time_summary <- moran_time_summary %>% arrange(timepoint)

moran_time_summary$timepoint <- as.character(moran_time_summary$timepoint)

region_section_rows <- list()
counter <- 1L

for (sec in SECTION_ORDER) {
    d <- md[md$section_uid_internal == sec, , drop = FALSE]
    regions <- sort(unique(as.character(d$spatial_region_v03)))
    regions <- regions[!is.na(regions) & nzchar(regions)]
    for (reg in regions) {
        in_idx <- as.character(d$spatial_region_v03) == reg
        out_idx <- !in_idx
        n_in <- sum(in_idx)
        n_out <- sum(out_idx)
        if (n_in < MIN_REGION_SPOTS || n_out < MIN_OTHER_SPOTS) {
            next
        }
        mean_in <- mean(d$Clec7a_RCTDnorm[in_idx])
        mean_out <- mean(d$Clec7a_RCTDnorm[out_idx])
        region_section_rows[[counter]] <- data.frame(section_uid_internal = sec, sample_label = as.character(unique(d$sample_label)[1]), 
            GSM = as.character(unique(d$GSM)[1]), timepoint = as.character(unique(d$timepoint)[1]), section_position_internal = as.character(unique(d$section_position_internal)[1]), 
            spatial_region = reg, n_region_spots = n_in, n_other_spots = n_out, mean_target_region = mean_in, 
            mean_target_other = mean_out, difference = mean_in - mean_out, enrichment_ratio = (mean_in + 
                EPS)/(mean_out + EPS), log2_enrichment = log2((mean_in + EPS)/(mean_out + EPS)), stringsAsFactors = FALSE)
        counter <- counter + 1L
    }
}

region_section <- bind_rows(region_section_rows)

direct_section_focus <- summarize_section_focus(effect_df = region_section, effect_col = "log2_enrichment", 
    ratio_col = "enrichment_ratio", analysis_label = "direct_region_enrichment")

direct_blocked_effects <- build_blocked_effects(effect_df = region_section, effect_col = "log2_enrichment", 
    analysis_label = "direct_region_enrichment")

direct_blocked_summary <- summarize_blocked_focus(blocked_df = direct_blocked_effects, analysis_label = "direct_region_enrichment")

neigh_section_rows <- list()
counter <- 1L

for (sec in SECTION_ORDER) {
    d <- md[md$section_uid_internal == sec, , drop = FALSE]
    regions <- sort(unique(as.character(d$spatial_region_v03)))
    regions <- regions[!is.na(regions) & nzchar(regions)]
    valid_neighbor <- is.finite(d$Clec7a_neighbor_mean_v02)
    for (reg in regions) {
        in_idx <- as.character(d$spatial_region_v03) == reg
        out_idx <- !in_idx
        in_valid <- in_idx & valid_neighbor
        out_valid <- out_idx & valid_neighbor
        n_in <- sum(in_valid)
        n_out <- sum(out_valid)
        if (n_in < MIN_REGION_SPOTS || n_out < MIN_OTHER_SPOTS) {
            next
        }
        mean_in <- mean(d$Clec7a_neighbor_mean_v02[in_valid])
        mean_out <- mean(d$Clec7a_neighbor_mean_v02[out_valid])
        neigh_section_rows[[counter]] <- data.frame(section_uid_internal = sec, sample_label = as.character(unique(d$sample_label)[1]), 
            GSM = as.character(unique(d$GSM)[1]), timepoint = as.character(unique(d$timepoint)[1]), section_position_internal = as.character(unique(d$section_position_internal)[1]), 
            spatial_region = reg, n_region_spots_with_neighbors = n_in, n_other_spots_with_neighbors = n_out, 
            mean_neighbor_target_region = mean_in, mean_neighbor_target_other = mean_out, neighbor_difference = mean_in - 
                mean_out, neighbor_enrichment_ratio = (mean_in + EPS)/(mean_out + EPS), log2_neighbor_enrichment = log2((mean_in + 
                EPS)/(mean_out + EPS)), stringsAsFactors = FALSE)
        counter <- counter + 1L
    }
}

neigh_section <- bind_rows(neigh_section_rows)

neigh_section_focus <- summarize_section_focus(effect_df = neigh_section, effect_col = "log2_neighbor_enrichment", 
    ratio_col = "neighbor_enrichment_ratio", analysis_label = "neighborhood_enrichment")

neigh_blocked_effects <- build_blocked_effects(effect_df = neigh_section, effect_col = "log2_neighbor_enrichment", 
    analysis_label = "neighborhood_enrichment")

neigh_blocked_summary <- summarize_blocked_focus(blocked_df = neigh_blocked_effects, analysis_label = "neighborhood_enrichment")

GRAPH_QC_FILE <- file.path(OUT_DIR, paste0("01_66_section_spatial_graph_QC_", VERSION, ".csv"))

MORAN_FILE <- file.path(OUT_DIR, paste0("02_66_Clec7a_section_global_MoransI_", VERSION, ".csv"))

MORAN_TIME_FILE <- file.path(OUT_DIR, paste0("02b_66_Clec7a_section_MoransI_timepoint_summary_", VERSION, 
    ".csv"))

DIRECT_SECTION_FILE <- file.path(OUT_DIR, paste0("03_66_Clec7a_section_region_enrichment_", VERSION, 
    ".csv"))

DIRECT_FOCUS_FILE <- file.path(OUT_DIR, paste0("04_66_Clec7a_region_enrichment_focus_tests_", VERSION, 
    ".csv"))

NEIGH_SECTION_FILE <- file.path(OUT_DIR, paste0("05_66_Clec7a_section_neighborhood_enrichment_", VERSION, 
    ".csv"))

NEIGH_FOCUS_FILE <- file.path(OUT_DIR, paste0("06_66_Clec7a_neighborhood_enrichment_focus_tests_", VERSION, 
    ".csv"))

SPOT_NEIGHBOR_FILE <- file.path(OUT_DIR, paste0("07_66_Clec7a_section_neighbor_scores_", VERSION, ".csv.gz"))

write_csv_utf8(graph_qc, GRAPH_QC_FILE)

write_csv_utf8(moran_results, MORAN_FILE)

write_csv_utf8(moran_time_summary, MORAN_TIME_FILE)

write_csv_utf8(region_section, DIRECT_SECTION_FILE)

write_csv_utf8(bind_rows(direct_section_focus, direct_blocked_summary), DIRECT_FOCUS_FILE)

write_csv_utf8(neigh_section, NEIGH_SECTION_FILE)

write_csv_utf8(bind_rows(neigh_section_focus, neigh_blocked_summary), NEIGH_FOCUS_FILE)

spot_neighbor_output <- md %>% transmute(barcode, GSM, sample_label, timepoint, section_position_internal, 
    section_uid_internal, spatial_region_v03, array_x, array_y, Clec7a_RCTDnorm, Clec7a_neighbor_mean_v02)

write_csv_gz_base(spot_neighbor_output, SPOT_NEIGHBOR_FILE)

shape_values <- c(`3dpi_rep1` = 16, `3dpi_rep2` = 17, `7dpi_rep1` = 15, `7dpi_rep2` = 3, `14dpi_rep1` = 0, 
    `14dpi_rep2` = 8, `28dpi_rep1` = 1, `28dpi_rep2` = 4)

plot_moran <- moran_results

plot_moran$timepoint <- factor(plot_moran$timepoint, levels = TIME_LEVELS)

plot_moran$sample_label <- factor(plot_moran$sample_label, levels = SAMPLE_ORDER)

p1 <- ggplot(plot_moran, aes(x = timepoint, y = Moran_I)) + geom_hline(yintercept = 0, linewidth = 0.35, 
    linetype = 2) + geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.45) + geom_point(aes(shape = sample_label), 
    position = position_jitter(width = 0.11, height = 0), size = 2, alpha = 0.9) + scale_shape_manual(values = shape_values, 
    drop = FALSE) + labs(title = "Section-level spatial autocorrelation of Clec7a+ normalized RCTD weight", 
    subtitle = paste0("Moran's I; ", MORAN_NSIM, " within-section permutations"), x = NULL, y = "Moran's I", 
    shape = "Visium capture area") + theme_classic(base_family = "Arial", base_size = 9) + theme(legend.position = "right", 
    legend.text = element_text(size = 7), legend.title = element_text(size = 7.5))

plot_region <- region_section %>% filter(spatial_region %in% FOCUS_REGIONS)

plot_region$timepoint <- factor(plot_region$timepoint, levels = TIME_LEVELS)

plot_region$sample_label <- factor(plot_region$sample_label, levels = SAMPLE_ORDER)

plot_region$spatial_region <- factor(plot_region$spatial_region, levels = FOCUS_REGIONS)

p2 <- ggplot(plot_region, aes(x = timepoint, y = log2_enrichment)) + geom_hline(yintercept = 0, linewidth = 0.35, 
    linetype = 2) + geom_point(aes(shape = sample_label), position = position_jitter(width = 0.11, height = 0), 
    size = 1.8, alpha = 0.9) + scale_shape_manual(values = shape_values, drop = FALSE) + facet_wrap(~spatial_region, 
    ncol = 1, scales = "free_y") + labs(title = "Section-level direct Clec7a+ enrichment in predefined spatial niches", 
    subtitle = "Each point is one within-section log2(region / other) effect size", x = NULL, y = "log2(region / other)", 
    shape = "Visium capture area") + theme_classic(base_family = "Arial", base_size = 8.5) + theme(strip.text = element_text(size = 7.5), 
    legend.position = "right", legend.text = element_text(size = 7), legend.title = element_text(size = 7.5))

plot_neigh <- neigh_section %>% filter(spatial_region %in% FOCUS_REGIONS)

plot_neigh$timepoint <- factor(plot_neigh$timepoint, levels = TIME_LEVELS)

plot_neigh$sample_label <- factor(plot_neigh$sample_label, levels = SAMPLE_ORDER)

plot_neigh$spatial_region <- factor(plot_neigh$spatial_region, levels = FOCUS_REGIONS)

p3 <- ggplot(plot_neigh, aes(x = timepoint, y = log2_neighbor_enrichment)) + geom_hline(yintercept = 0, 
    linewidth = 0.35, linetype = 2) + geom_point(aes(shape = sample_label), position = position_jitter(width = 0.11, 
    height = 0), size = 1.8, alpha = 0.9) + scale_shape_manual(values = shape_values, drop = FALSE) + 
    facet_wrap(~spatial_region, ncol = 1, scales = "free_y") + labs(title = "Section-level Clec7a+ neighborhood enrichment", 
    subtitle = "Each point is one within-section log2(neighbor signal in region / other) effect size", 
    x = NULL, y = "log2(neighborhood region / other)", shape = "Visium capture area") + theme_classic(base_family = "Arial", 
    base_size = 8.5) + theme(strip.text = element_text(size = 7.5), legend.position = "right", legend.text = element_text(size = 7), 
    legend.title = element_text(size = 7.5))

blocked_plot <- bind_rows(direct_blocked_effects %>% transmute(sample_label, timepoint, spatial_region, 
    analysis = "Direct", blocked_log2_effect), neigh_blocked_effects %>% transmute(sample_label, timepoint, 
    spatial_region, analysis = "Neighborhood", blocked_log2_effect))

blocked_plot$spatial_region <- factor(blocked_plot$spatial_region, levels = FOCUS_REGIONS)

blocked_plot$sample_label <- factor(blocked_plot$sample_label, levels = SAMPLE_ORDER)

p4 <- ggplot(blocked_plot, aes(x = analysis, y = blocked_log2_effect)) + geom_hline(yintercept = 0, linewidth = 0.35, 
    linetype = 2) + geom_point(aes(shape = sample_label), position = position_jitter(width = 0.1, height = 0), 
    size = 2, alpha = 0.9) + scale_shape_manual(values = shape_values, drop = FALSE) + facet_wrap(~spatial_region, 
    ncol = 1, scales = "free_y") + labs(title = "Capture-area-blocked sensitivity analysis", subtitle = "Each point is the median of valid section-level effects within one Visium capture area", 
    x = NULL, y = "Blocked log2 enrichment", shape = "Visium capture area") + theme_classic(base_family = "Arial", 
    base_size = 8.5) + theme(strip.text = element_text(size = 7.5), legend.position = "right", legend.text = element_text(size = 7), 
    legend.title = element_text(size = 7.5))

plot_spot <- spot_neighbor_output

plot_spot$sample_label <- factor(plot_spot$sample_label, levels = SAMPLE_ORDER)

p5 <- ggplot(plot_spot, aes(x = array_x, y = array_y, color = Clec7a_neighbor_mean_v02)) + geom_point(size = 0.48, 
    alpha = 0.95) + coord_equal() + scale_y_reverse() + scale_color_viridis_c(option = "magma", limits = c(0, 
    1), oob = scales::squish, na.value = "grey85") + facet_wrap(~sample_label, ncol = 4) + labs(title = "Immediate-neighborhood Clec7a+ normalized RCTD signal", 
    subtitle = "Neighbor means are calculated independently within each coordinate-defined tissue section", 
    x = "Visium array X", y = "Visium array Y", color = "Neighbor\nmean") + theme_classic(base_family = "Arial", 
    base_size = 8.5) + theme(strip.text = element_text(size = 7.5))

DIAG_PDF <- file.path(OUT_DIR, paste0("08_66_Clec7a_section_spatial_statistics_diagnostics_", VERSION, 
    ".pdf"))

grDevices::cairo_pdf(DIAG_PDF, width = 10.5, height = 7.8, family = "Arial")

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

dev.off()

BLOCKED_PNG <- file.path(OUT_DIR, paste0("08b_66_Clec7a_capture_area_blocked_focus_", VERSION, ".png"))

ggsave(filename = BLOCKED_PNG, plot = p4, width = 7.4, height = 7, units = "in", dpi = 300, bg = "white")

graph_reasonable <- all(graph_qc$mean_degree >= 2 & graph_qc$mean_degree <= 8)

isolates_reasonable <- all(graph_qc$fraction_isolated <= 0.05)

all_moran_finite <- all(is.finite(moran_results$Moran_I))

focus_direct_present <- all(FOCUS_REGIONS %in% direct_section_focus$spatial_region)

focus_neigh_present <- all(FOCUS_REGIONS %in% neigh_section_focus$spatial_region)

blocked_direct_present <- all(FOCUS_REGIONS %in% direct_blocked_summary$spatial_region)

blocked_neigh_present <- all(FOCUS_REGIONS %in% neigh_blocked_summary$spatial_region)

section_counts_ok <- all(sections_per_capture$n_sections == 4L)

final_checks <- data.frame(check = c("Step65 v06 section-aware spot-weight file was used", "All 8 Visium capture areas are retained", 
    "All 32 coordinate-defined major tissue sections are analyzed", "Each capture area contributes exactly four coordinate-defined sections", 
    "Spatial graph is rebuilt independently within each tissue section", "Spatial graph mean degree is between 2 and 8 in every section", 
    "Spatial graph isolated-spot fraction is <=5% in every section", "Moran's I is finite for all 32 sections", 
    "All 3 prespecified regions are represented in direct-enrichment summaries", "All 3 prespecified regions are represented in neighborhood summaries", 
    "Capture-area-blocked summaries are available for all 3 prespecified regions", "Insufficient blocked replication is reported as NA rather than forced significance", 
    "No timepoint-specific hypothesis testing is performed", "No per-spot significance testing is performed", 
    "No animal-level independence is assumed", "No H&E-defined lesion-core/border claim is made"), pass = c(grepl("section_aware_v06", 
    INPUT_FILE), n_capture_areas == EXPECTED_CAPTURE_AREAS, n_sections == EXPECTED_MAJOR_SECTIONS, section_counts_ok, 
    TRUE, graph_reasonable, isolates_reasonable, all_moran_finite, focus_direct_present, focus_neigh_present, 
    blocked_direct_present & blocked_neigh_present, TRUE, TRUE, TRUE, TRUE, TRUE), stringsAsFactors = FALSE)

FINAL_QC_FILE <- file.path(OUT_DIR, paste0("09_66_GSE195783_Clec7a_spatial_statistics_FINAL_QC_", VERSION, 
    ".csv"))

write_csv_utf8(final_checks, FINAL_QC_FILE)

fmt_p <- function(x) {
    if (!is.finite(x)) {
        return("NA")
    }
    format(x, scientific = TRUE, digits = 3)
}

moran_time_lines <- apply(moran_time_summary, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": n_sections=", z[["n_sections"]], "; median Moran I=", sprintf("%.4f", 
        as.numeric(z[["median_Moran_I"]])), "; range=", sprintf("%.4f", as.numeric(z[["min_Moran_I"]])), 
        " to ", sprintf("%.4f", as.numeric(z[["max_Moran_I"]])), "; FDR<0.05 sections=", z[["n_FDR_lt_0_05"]])
})

direct_section_lines <- apply(direct_section_focus, 1, function(z) {
    paste0("  ", z[["spatial_region"]], ": n_sections=", z[["n_sections_tested"]], "; median ratio=", 
        sprintf("%.3f", as.numeric(z[["median_enrichment_ratio"]])), "; fraction sections >1=", sprintf("%.3f", 
            as.numeric(z[["fraction_sections_effect_gt_1"]])), "; exploratory section-level P=", fmt_p(as.numeric(z[["exploratory_wilcox_p_greater_zero"]])), 
        "; exploratory FDR=", fmt_p(as.numeric(z[["exploratory_FDR_BH"]])))
})

direct_blocked_lines <- apply(direct_blocked_summary, 1, function(z) {
    paste0("  ", z[["spatial_region"]], ": n_capture_areas=", z[["n_capture_areas_tested"]], "; median blocked ratio=", 
        sprintf("%.3f", as.numeric(z[["median_blocked_enrichment_ratio"]])), "; fraction capture areas >1=", 
        sprintf("%.3f", as.numeric(z[["fraction_capture_areas_effect_gt_1"]])), "; blocked Wilcoxon P=", 
        fmt_p(as.numeric(z[["wilcox_p_greater_zero"]])), "; blocked FDR=", fmt_p(as.numeric(z[["FDR_BH"]])))
})

neigh_section_lines <- apply(neigh_section_focus, 1, function(z) {
    paste0("  ", z[["spatial_region"]], ": n_sections=", z[["n_sections_tested"]], "; median ratio=", 
        sprintf("%.3f", as.numeric(z[["median_enrichment_ratio"]])), "; fraction sections >1=", sprintf("%.3f", 
            as.numeric(z[["fraction_sections_effect_gt_1"]])), "; exploratory section-level P=", fmt_p(as.numeric(z[["exploratory_wilcox_p_greater_zero"]])), 
        "; exploratory FDR=", fmt_p(as.numeric(z[["exploratory_FDR_BH"]])))
})

neigh_blocked_lines <- apply(neigh_blocked_summary, 1, function(z) {
    paste0("  ", z[["spatial_region"]], ": n_capture_areas=", z[["n_capture_areas_tested"]], "; median blocked ratio=", 
        sprintf("%.3f", as.numeric(z[["median_blocked_enrichment_ratio"]])), "; fraction capture areas >1=", 
        sprintf("%.3f", as.numeric(z[["fraction_capture_areas_effect_gt_1"]])), "; blocked Wilcoxon P=", 
        fmt_p(as.numeric(z[["wilcox_p_greater_zero"]])), "; blocked FDR=", fmt_p(as.numeric(z[["FDR_BH"]])))
})

REPORT_FILE <- file.path(OUT_DIR, paste0("10_66_GSE195783_Clec7a_spatial_statistics_REPORT_", VERSION, 
    ".txt"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(moran_time_summary)

print(direct_section_focus)

print(direct_blocked_summary)

print(neigh_section_focus)

print(neigh_blocked_summary)

print(final_checks)

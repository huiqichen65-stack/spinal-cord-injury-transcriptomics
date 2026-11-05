source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(1234)

ROOT <- "D:/zz"

INPUT_FILE <- file.path(ROOT, "results", "60_GSE195783_Visium_import_coordinate_reconstruction", "03_60_GSE195783_final_spot_metadata_v05.csv.gz")

OUT_DIR <- file.path(ROOT, "results", "67_GSE195783_spatial_section_identity_audit")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v01"
STEP <- 67L

EXPECTED_TOTAL_SPOTS <- 15449L

EXPECTED_SAMPLES <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", 
    "28dpi_rep1", "28dpi_rep2")

EXPECTED_MAJOR_COMPONENTS <- 4L

PRIMARY_ADJ_FACTOR <- 1.60

SENSITIVITY_FACTORS <- c(1.4, 1.6, 1.8, 2)

MAJOR_MIN_SPOTS_ABSOLUTE <- 50L
MAJOR_MIN_FRACTION_SAMPLE <- 0.05

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

get_connected_components <- function(dmat, threshold) {
    n <- nrow(dmat)
    if (n == 0L) {
        return(integer(0))
    }
    if (n == 1L) {
        return(1L)
    }
    adj <- (dmat <= threshold & dmat > 0)
    component <- integer(n)
    current_component <- 0L
    for (start in seq_len(n)) {
        if (component[start] != 0L) {
            next
        }
        current_component <- current_component + 1L
        queue <- start
        component[start] <- current_component
        q_index <- 1L
        while (q_index <= length(queue)) {
            v <- queue[q_index]
            neigh <- which(adj[v, ] & component == 0L)
            if (length(neigh) > 0L) {
                component[neigh] <- current_component
                queue <- c(queue, neigh)
            }
            q_index <- q_index + 1L
        }
    }
    component
}

audit_one_sample <- function(
    df,
    adj_factor,
    major_min_abs,
    major_min_fraction) {

  coords <- as.matrix(
    df[
      ,
      c(
        "array_x",
        "array_y"
      )
    ]
  )

  storage.mode(coords) <- "double"

  if (nrow(coords) < 2L) {
    stop(
      "Too few spots in sample: ",
      unique(df$sample_label)
    )
  }

  dmat <- as.matrix(
    stats::dist(
      coords,
      method = "euclidean"
    )
  )

  diag(dmat) <- Inf

  nearest_neighbor_distance <- apply(
    dmat,
    1,
    min,
    na.rm = TRUE
  )

  base_nn <- stats::median(
    nearest_neighbor_distance,
    na.rm = TRUE
  )

  if (
    !is.finite(base_nn) ||
      base_nn <= 0
  ) {
    stop(
      "Could not determine nearest-neighbor spacing for sample: ",
      unique(df$sample_label)
    )
  }

  threshold <- base_nn * adj_factor

  component_id <- get_connected_components(
    dmat = dmat,
    threshold = threshold
  )

  comp_tab <- as.data.frame(
    table(component_id),
    stringsAsFactors = FALSE
  )

  colnames(comp_tab) <- c(
    "component_id",
    "n_spots"
  )

  comp_tab$component_id <- as.integer(
    as.character(comp_tab$component_id)
  )

  n_sample <- nrow(df)

  major_cutoff <- max(
    major_min_abs,
    ceiling(
      n_sample *
        major_min_fraction
    )
  )

  comp_tab$is_major <- (
    comp_tab$n_spots >=
      major_cutoff
  )

  major_ids <- comp_tab %>%
    filter(is_major) %>%
    arrange(desc(n_spots)) %>%
    pull(component_id)

  major_rank_lookup <- setNames(
    seq_along(major_ids),
    major_ids
  )

  major_rank <- rep(
    NA_integer_,
    length(component_id)
  )

  idx_major <- component_id %in% major_ids

  major_rank[idx_major] <- unname(
    major_rank_lookup[
      as.character(
        component_id[idx_major]
      )
    ]
  )

  list(
    component_id = component_id,
    major_rank = major_rank,
    base_nn = base_nn,
    threshold = threshold,
    major_cutoff = major_cutoff,
    component_table = comp_tab,
    n_raw_components = nrow(comp_tab),
    n_major_components = sum(comp_tab$is_major),
    n_tiny_components = sum(!comp_tab$is_major),
    n_spots_major = sum(
      comp_tab$n_spots[
        comp_tab$is_major
      ]
    ),
    n_spots_tiny = sum(
      comp_tab$n_spots[
        !comp_tab$is_major
      ]
    )
  )
}

if (!file.exists(INPUT_FILE)) {
    stop("Input metadata file not found:\n", INPUT_FILE)
}

con <- gzfile(INPUT_FILE, open = "rt")

spot_md <- read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)

close(con)

required_cols <- c("barcode", "GSM", "timepoint", "replicate", "sample_label", "array_x", "array_y", 
    "histology_registered")

missing_cols <- setdiff(required_cols, colnames(spot_md))

if (length(missing_cols) > 0L) {
    stop("Input metadata is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

if (nrow(spot_md) != EXPECTED_TOTAL_SPOTS) {
    warning("Expected ", EXPECTED_TOTAL_SPOTS, " Step60 spots, but found ", nrow(spot_md), ". Audit will continue, but inspect this discrepancy.")
}

if (any(spot_md$histology_registered %in% TRUE, na.rm = TRUE)) {
    stop("histology_registered must remain FALSE.")
}

observed_samples <- unique(as.character(spot_md$sample_label))

missing_expected_samples <- setdiff(EXPECTED_SAMPLES, observed_samples)

unexpected_samples <- setdiff(observed_samples, EXPECTED_SAMPLES)

if (length(missing_expected_samples) > 0L) {
    stop("Missing expected sample(s): ", paste(missing_expected_samples, collapse = ", "))
}

if (length(unexpected_samples) > 0L) {
    stop("Unexpected sample(s): ", paste(unexpected_samples, collapse = ", "))
}

spot_md$sample_label <- factor(spot_md$sample_label, levels = EXPECTED_SAMPLES)

sensitivity_rows <- list()

sensitivity_counter <- 1L

for (samp in EXPECTED_SAMPLES) {
    df_s <- spot_md %>% filter(as.character(sample_label) == samp)
    for (fac in SENSITIVITY_FACTORS) {
        aud <- audit_one_sample(df = df_s, adj_factor = fac, major_min_abs = MAJOR_MIN_SPOTS_ABSOLUTE, 
            major_min_fraction = MAJOR_MIN_FRACTION_SAMPLE)
        sensitivity_rows[[sensitivity_counter]] <- data.frame(sample_label = samp, timepoint = unique(as.character(df_s$timepoint))[1], 
            replicate = unique(as.character(df_s$replicate))[1], GSM = unique(as.character(df_s$GSM))[1], 
            adjacency_factor = fac, median_nearest_neighbor_distance = aud$base_nn, adjacency_threshold = aud$threshold, 
            major_component_min_spots = aud$major_cutoff, n_spots_total = nrow(df_s), n_raw_components = aud$n_raw_components, 
            n_major_components = aud$n_major_components, n_tiny_components = aud$n_tiny_components, n_spots_in_major_components = aud$n_spots_major, 
            n_spots_in_tiny_components = aud$n_spots_tiny, fraction_spots_in_major_components = aud$n_spots_major/nrow(df_s), 
            expected_four_major_components = aud$n_major_components == EXPECTED_MAJOR_COMPONENTS, stringsAsFactors = FALSE)
        sensitivity_counter <- sensitivity_counter + 1L
    }
}

sensitivity_df <- bind_rows(sensitivity_rows)

primary_rows <- list()

component_rows <- list()

spot_rows <- list()

pcounter <- 1L
ccounter <- 1L
scounter <- 1L

for (
  samp in EXPECTED_SAMPLES
) {

  df_s <- spot_md %>%
    filter(
      as.character(sample_label) ==
        samp
    )

  aud <- audit_one_sample(
    df = df_s,
    adj_factor =
      PRIMARY_ADJ_FACTOR,
    major_min_abs =
      MAJOR_MIN_SPOTS_ABSOLUTE,
    major_min_fraction =
      MAJOR_MIN_FRACTION_SAMPLE
  )

  primary_rows[[pcounter]] <- data.frame(

    sample_label = samp,

    timepoint = unique(
      as.character(
        df_s$timepoint
      )
    )[1],

    replicate = unique(
      as.character(
        df_s$replicate
      )
    )[1],

    GSM = unique(
      as.character(
        df_s$GSM
      )
    )[1],

    n_spots_total =
      nrow(df_s),

    median_nearest_neighbor_distance =
      aud$base_nn,

    adjacency_factor =
      PRIMARY_ADJ_FACTOR,

    adjacency_threshold =
      aud$threshold,

    major_component_min_spots =
      aud$major_cutoff,

    n_raw_components =
      aud$n_raw_components,

    n_major_components =
      aud$n_major_components,

    n_tiny_components =
      aud$n_tiny_components,

    n_spots_in_major_components =
      aud$n_spots_major,

    n_spots_in_tiny_components =
      aud$n_spots_tiny,

    fraction_spots_in_major_components =
      aud$n_spots_major /
      nrow(df_s),

    expected_four_major_components =
      aud$n_major_components ==
      EXPECTED_MAJOR_COMPONENTS,

    stringsAsFactors = FALSE
  )

  comp_table <- aud$component_table

  coords_s <- data.frame(
    component_id =
      aud$component_id,
    array_x =
      as.numeric(
        df_s$array_x
      ),
    array_y =
      as.numeric(
        df_s$array_y
      ),
    stringsAsFactors = FALSE
  )

  centroid_table <- coords_s %>%
    group_by(
      component_id
    ) %>%
    summarise(
      centroid_x =
        mean(array_x),
      centroid_y =
        mean(array_y),
      min_x =
        min(array_x),
      max_x =
        max(array_x),
      min_y =
        min(array_y),
      max_y =
        max(array_y),
      .groups = "drop"
    )

  comp_table <- comp_table %>%
    left_join(
      centroid_table,
      by = "component_id"
    ) %>%
    mutate(
      sample_label = samp,
      timepoint = unique(
        as.character(
          df_s$timepoint
        )
      )[1],
      replicate = unique(
        as.character(
          df_s$replicate
        )
      )[1],
      GSM = unique(
        as.character(
          df_s$GSM
        )
      )[1],
      fraction_of_sample =
        n_spots /
        nrow(df_s)
    ) %>%
    select(
      sample_label,
      timepoint,
      replicate,
      GSM,
      component_id,
      is_major,
      n_spots,
      fraction_of_sample,
      centroid_x,
      centroid_y,
      min_x,
      max_x,
      min_y,
      max_y
    ) %>%
    arrange(
      desc(is_major),
      desc(n_spots)
    )

  component_rows[[ccounter]] <-
    comp_table

  spot_tmp <- df_s %>%
    transmute(
      barcode =
        as.character(barcode),
      GSM =
        as.character(GSM),
      timepoint =
        as.character(timepoint),
      replicate =
        as.character(replicate),
      sample_label =
        as.character(sample_label),
      array_x =
        as.numeric(array_x),
      array_y =
        as.numeric(array_y)
    )

  spot_tmp$component_id_internal <-
    aud$component_id

  major_component_ids <-
    aud$component_table$component_id[
      aud$component_table$is_major
    ]

  spot_tmp$is_major_component <-
    spot_tmp$component_id_internal %in%
      major_component_ids

  spot_tmp$major_component_rank_internal <-
    aud$major_rank

  spot_rows[[scounter]] <-
    spot_tmp

  pcounter <- pcounter + 1L
  ccounter <- ccounter + 1L
  scounter <- scounter + 1L
}

primary_df <- bind_rows(primary_rows)

component_df <- bind_rows(component_rows)

spot_component_df <- bind_rows(spot_rows)

stability_df <- sensitivity_df %>% group_by(sample_label, timepoint, replicate, GSM) %>% summarise(n_thresholds_tested = n(), 
    n_thresholds_with_four_major_components = sum(expected_four_major_components), all_thresholds_support_four = all(expected_four_major_components), 
    major_component_count_range = paste0(min(n_major_components), "-", max(n_major_components)), max_fraction_spots_tiny = max(1 - 
        fraction_spots_in_major_components), .groups = "drop")

SAMPLE_SUMMARY_FILE <- file.path(OUT_DIR, paste0("01_67_GSE195783_spatial_section_identity_audit_sample_summary_", 
    VERSION, ".csv"))

COMPONENT_SUMMARY_FILE <- file.path(OUT_DIR, paste0("02_67_GSE195783_spatial_section_identity_audit_component_summary_", 
    VERSION, ".csv"))

SENSITIVITY_FILE <- file.path(OUT_DIR, paste0("03_67_GSE195783_spatial_section_identity_audit_threshold_sensitivity_", 
    VERSION, ".csv"))

STABILITY_FILE <- file.path(OUT_DIR, paste0("04_67_GSE195783_spatial_section_identity_audit_stability_summary_", 
    VERSION, ".csv"))

SPOT_COMPONENT_FILE <- file.path(OUT_DIR, paste0("05_67_GSE195783_spatial_section_identity_audit_spot_component_map_", 
    VERSION, ".csv.gz"))

write_csv_utf8(primary_df, SAMPLE_SUMMARY_FILE)

write_csv_utf8(component_df, COMPONENT_SUMMARY_FILE)

write_csv_utf8(sensitivity_df, SENSITIVITY_FILE)

write_csv_utf8(stability_df, STABILITY_FILE)

gz_con <- gzfile(SPOT_COMPONENT_FILE, open = "wt")

write.csv(spot_component_df, file = gz_con, row.names = FALSE, fileEncoding = "UTF-8", na = "")

close(gz_con)

plot_df <- spot_component_df

plot_df$display_group <- ifelse(plot_df$is_major_component, paste0("major_", plot_df$major_component_rank_internal), 
    "tiny")

audit_palette <- c(major_1 = "#4C78A8", major_2 = "#F58518", major_3 = "#54A24B", major_4 = "#B279A2", 
    tiny = "#BDBDBD")

plot_df$sample_label <- factor(plot_df$sample_label, levels = EXPECTED_SAMPLES)

p_audit <- ggplot(plot_df, aes(x = array_x, y = array_y, color = display_group)) + geom_point(size = 0.55, 
    alpha = 0.95) + coord_equal() + scale_y_reverse() + scale_color_manual(values = audit_palette, drop = FALSE) + 
    facet_wrap(~sample_label, ncol = 4) + labs(title = "GSE195783 spatial tissue-component audit", subtitle = paste0("Major disconnected tissue pieces detected from Visium array coordinates; ", 
    "no H&E/pixel registration"), x = "Visium array X", y = "Visium array Y") + theme_classic(base_family = "Arial", 
    base_size = 9) + theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5), 
    strip.text = element_text(size = 8), legend.position = "none")

PDF_FILE <- file.path(OUT_DIR, paste0("06_67_GSE195783_spatial_section_identity_audit_component_map_", 
    VERSION, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0("06_67_GSE195783_spatial_section_identity_audit_component_map_", 
    VERSION, ".png"))

ggsave(filename = PDF_FILE, plot = p_audit, width = 11.5, height = 6.8, units = "in", device = grDevices::cairo_pdf, 
    bg = "white")

ggsave(filename = PNG_FILE, plot = p_audit, width = 11.5, height = 6.8, units = "in", dpi = 300, bg = "white")

all_primary_four <- all(primary_df$n_major_components == EXPECTED_MAJOR_COMPONENTS)

all_stable_four <- all(stability_df$all_thresholds_support_four)

n_samples_primary_four <- sum(primary_df$n_major_components == EXPECTED_MAJOR_COMPONENTS)

n_samples_stable_four <- sum(stability_df$all_thresholds_support_four)

final_qc <- data.frame(check = c("Input contains 8 expected Visium capture areas", "No H&E registration is claimed", 
    "Primary threshold identifies four major tissue pieces in every capture area", "Four-major-piece result is stable across all tested thresholds", 
    "No object was modified"), pass = c(length(unique(as.character(spot_md$sample_label))) == 8L, all(spot_md$histology_registered == 
    FALSE), all_primary_four, all_stable_four, TRUE), stringsAsFactors = FALSE)

QC_FILE <- file.path(OUT_DIR, paste0("07_67_GSE195783_spatial_section_identity_audit_FINAL_QC_", VERSION, 
    ".csv"))

write_csv_utf8(final_qc, QC_FILE)

REPORT_FILE <- file.path(OUT_DIR, paste0("08_67_GSE195783_spatial_section_identity_audit_REPORT_", VERSION, 
    ".txt"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(primary_df)

print(stability_df)

print(final_qc)

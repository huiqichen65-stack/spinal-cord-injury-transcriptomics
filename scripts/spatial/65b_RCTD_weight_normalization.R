source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(1234)

ROOT <- "D:/zz"

STEP64_DIR <- file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation")

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

CHECKPOINT_DIR <- file.path(STEP65_DIR, "checkpoints_v02")

SPATIAL_FILE <- file.path(STEP64_DIR, "12_64_GSE195783_spatial_region_annotated_v03.rds")

REF_MANIFEST_FILE <- file.path(STEP65_DIR, "01_65_reference_manifest_v03.csv")

LABEL_MAP_FILE <- file.path(STEP65_DIR, "01a_65_RCTD_internal_label_map_v02.csv")

VERSION <- "v04"

TARGET_STATE <- "Clec7a+ lipid-handling Microglia"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

SAMPLE_ORDER <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1",
  "28dpi_rep2")

MYELOID_STATE_DISPLAY_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-handling Microglia", "Activated inflammatory Microglia",
  "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage",
  "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

TARGET_SANITY_GENES <- c("Lpl", "Apoe", "Trem2", "Abca1", "Cst7", "Apoc1")

SUM_TOL <- 1e-08

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr", "SummarizedExperiment",
  "SpatialExperiment", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step65 v04.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(scales)
})

write_csv_utf8 <- write_csv_utf8

save_pdf_png <- function(p, stem, width, height) {
    ggsave(filename = file.path(STEP65_DIR, paste0(stem, "_", VERSION, ".pdf")), plot = p, width = width, 
        height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
    ggsave(filename = file.path(STEP65_DIR, paste0(stem, "_", VERSION, ".png")), plot = p, width = width, 
        height = height, units = "in", dpi = 300, bg = "white")
}

get_data_layer <- function(obj, assay = "Spatial") {
    obj2 <- obj
    if (inherits(obj2[[assay]], "Assay5")) {
        obj2 <- JoinLayers(object = obj2, assay = assay)
    }
    ans <- tryCatch({
        SeuratObject::LayerData(object = obj2, assay = assay, layer = "data")
    }, error = function(e) {
        NULL
    })
    if (is.null(ans)) {
        ans <- Seurat::GetAssayData(object = obj2, assay = assay, slot = "data")
    }
    ans
}

safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10L) {
        return(NA_real_)
    }
    suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

checkpoint_files <- file.path(CHECKPOINT_DIR, paste0(SAMPLE_ORDER, "_RCTD_full_v02.rds"))

required_files <- c(SPATIAL_FILE, REF_MANIFEST_FILE, LABEL_MAP_FILE, checkpoint_files)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("Required saved file(s) missing:\n", paste(missing_files, collapse = "\n"), "\nRCTD input check failed; inspect this error before continuing.")
}

spatial_obj <- readRDS(SPATIAL_FILE)

if (!inherits(spatial_obj, "Seurat")) {
  stop("Spatial input is not a Seurat object.")
}

spatial_md <- spatial_obj[[]]

required_spatial_meta <- c("GSM", "sample_label", "timepoint", "array_x", "array_y", "spatial_region_v03",
  "histology_registered")

missing_spatial_meta <- setdiff(required_spatial_meta, colnames(spatial_md))

if (length(missing_spatial_meta) > 0L) {
  stop("Spatial object missing metadata column(s): ", paste(missing_spatial_meta, collapse = ", "))
}

if (any(spatial_obj$histology_registered)) {
  stop("H&E registration must remain FALSE.")
}

reference_manifest <- read.csv(REF_MANIFEST_FILE, stringsAsFactors = FALSE, check.names = FALSE)

label_map <- read.csv(LABEL_MAP_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_map_cols <- c("display_label", "RCTD_internal_label")

missing_map_cols <- setdiff(required_map_cols, colnames(label_map))

if (length(missing_map_cols) > 0L) {
  stop("Label map is missing column(s): ", paste(missing_map_cols, collapse = ", "))
}

internal_to_display <- setNames(label_map$display_label, label_map$RCTD_internal_label)

write_csv_utf8(reference_manifest, file.path(STEP65_DIR, paste0("01_65_reference_manifest_", VERSION,
  ".csv")))

weights_list <- list()
sample_qc_list <- list()

for (i in seq_along(SAMPLE_ORDER)) {

  samp <- SAMPLE_ORDER[i]

  result_spe <- readRDS(checkpoint_files[i])

  if (!"weights" %in% SummarizedExperiment::assayNames(result_spe)) {
    stop(samp, ": checkpoint lacks assay 'weights'.")
  }

  weights_raw <- SummarizedExperiment::assay(result_spe, "weights")

  original_weight_names <- rownames(weights_raw)

  if (all(original_weight_names %in% names(internal_to_display))) {

    rownames(weights_raw) <- unname(internal_to_display[original_weight_names])

  } else if (all(original_weight_names %in% label_map$display_label)) {

    rownames(weights_raw) <- original_weight_names

  } else {

    bad_names <- setdiff(original_weight_names, c(names(internal_to_display), label_map$display_label))

    stop(samp, ": checkpoint contains unmapped RCTD type(s): ", paste(bad_names, collapse = " | "))
  }

  if (anyDuplicated(rownames(weights_raw)) > 0L) {
    stop(samp, ": mapped RCTD display labels are not unique.")
  }

  if (!TARGET_STATE %in% rownames(weights_raw)) {
    stop(samp, ": target state missing after label mapping:\n", TARGET_STATE)
  }

  raw_sums <- Matrix::colSums(weights_raw)

  if (any(!is.finite(raw_sums) | raw_sums <= 0)) {
    stop(samp, ": at least one spot has a non-positive/non-finite raw weight sum.")
  }

  weights_norm <- sweep(weights_raw, MARGIN = 2, STATS = raw_sums, FUN = "/")

  norm_sums <- Matrix::colSums(weights_norm)

  max_abs_sum_error <- max(abs(norm_sums - 1))

  if (!is.finite(max_abs_sum_error) || max_abs_sum_error > SUM_TOL) {
    stop(samp, ": normalized weight sums are not 1 within tolerance. ", "Max absolute error = ",
      signif(max_abs_sum_error, 6))
  }

  cells_input <- rownames(spatial_md)[as.character(spatial_md$sample_label) == samp]

  spot_common <- intersect(colnames(weights_norm), cells_input)

  if (length(spot_common) == 0L) {
    stop(samp, ": no checkpoint spot names match spatial object.")
  }

  weights_raw <- weights_raw[, spot_common, drop = FALSE]

  weights_norm <- weights_norm[, spot_common, drop = FALSE]

  raw_sums <- raw_sums[spot_common]

  norm_sums <- norm_sums[spot_common]

  md_out <- spatial_md[spot_common, , drop = FALSE]

  sample_qc_list[[samp]] <- data.frame(sample_label = samp, timepoint = as.character(unique(md_out$timepoint)[1]),
    n_spots_input = length(cells_input), n_spots_RCTD = ncol(weights_norm), fraction_spots_returned = ncol(weights_norm)/length(cells_input),
    n_reference_types = nrow(weights_norm), raw_mean_weight_sum = mean(raw_sums), raw_median_weight_sum = median(raw_sums),
    raw_min_weight_sum = min(raw_sums), raw_max_weight_sum = max(raw_sums), normalized_mean_weight_sum = mean(norm_sums),
    normalized_median_weight_sum = median(norm_sums), normalized_min_weight_sum = min(norm_sums),
    normalized_max_weight_sum = max(norm_sums), max_abs_normalized_sum_error = max(abs(norm_sums -
      1)), checkpoint_file = checkpoint_files[i], stringsAsFactors = FALSE)

  weights_df <- as.data.frame(t(as.matrix(weights_norm)), check.names = FALSE, stringsAsFactors = FALSE)

  weights_df$barcode <- rownames(weights_df)

  weights_df$GSM <- as.character(md_out$GSM)

  weights_df$sample_label <- as.character(md_out$sample_label)

  weights_df$timepoint <- as.character(md_out$timepoint)

  weights_df$spatial_region_v03 <- as.character(md_out$spatial_region_v03)

  weights_df$array_x <- md_out$array_x

  weights_df$array_y <- md_out$array_y

  weights_list[[samp]] <- weights_df

  rm(result_spe, weights_raw, weights_norm, weights_df)

  gc(verbose = FALSE)
}

sample_qc <- bind_rows(sample_qc_list)

sample_qc$timepoint <- factor(sample_qc$timepoint, levels = TIME_LEVELS)

sample_qc <- sample_qc %>% arrange(timepoint, sample_label)

write_csv_utf8(sample_qc, file.path(STEP65_DIR, paste0("02_65_RCTD_sample_QC_", VERSION, ".csv")))

spot_weights <- bind_rows(weights_list)

spot_weights$timepoint <- factor(spot_weights$timepoint, levels = TIME_LEVELS)

spot_weights$sample_label <- factor(spot_weights$sample_label, levels = SAMPLE_ORDER)

SPOT_WEIGHT_FILE <- file.path(STEP65_DIR, paste0("03_65_RCTD_normalized_spot_weights_", VERSION, ".csv.gz"))

con_weights <- gzfile(SPOT_WEIGHT_FILE, open = "wt")

write.csv(spot_weights, file = con_weights, row.names = FALSE, na = "")

close(con_weights)

target_sample_rows <- list()

for (samp in SAMPLE_ORDER) {
    d <- spot_weights[as.character(spot_weights$sample_label) == samp, , drop = FALSE]
    x <- d[[TARGET_STATE]]
    target_sample_rows[[samp]] <- data.frame(summary_level = "sample", timepoint = as.character(unique(d$timepoint)[1]), 
        sample_label = samp, n_spots = nrow(d), n_replicates = NA_integer_, mean_weight = mean(x, na.rm = TRUE), 
        median_weight = median(x, na.rm = TRUE), q25_target_weight = as.numeric(quantile(x, 0.25, na.rm = TRUE)), 
        q75_target_weight = as.numeric(quantile(x, 0.75, na.rm = TRUE)), fraction_weight_gt_0_10 = mean(x > 
            0.1, na.rm = TRUE), sd_of_sample_means = NA_real_, min_sample_mean = NA_real_, max_sample_mean = NA_real_, 
        stringsAsFactors = FALSE)
}

target_sample_summary <- bind_rows(target_sample_rows)

target_time_rows <- list()

for (tp in TIME_LEVELS) {
    vals <- target_sample_summary$mean_weight[target_sample_summary$timepoint == tp]
    target_time_rows[[tp]] <- data.frame(summary_level = "timepoint_across_replicates", timepoint = tp, 
        sample_label = NA_character_, n_spots = NA_integer_, n_replicates = length(vals), mean_weight = mean(vals, 
            na.rm = TRUE), median_weight = NA_real_, q25_target_weight = NA_real_, q75_target_weight = NA_real_, 
        fraction_weight_gt_0_10 = NA_real_, sd_of_sample_means = if (length(vals) > 1L) {
            stats::sd(vals, na.rm = TRUE)
        }
        else {
            NA_real_
        }, min_sample_mean = min(vals, na.rm = TRUE), max_sample_mean = max(vals, na.rm = TRUE), stringsAsFactors = FALSE)
}

target_time_summary <- bind_rows(target_time_rows)

target_time_output <- bind_rows(target_sample_summary, target_time_summary)

target_time_output$timepoint <- factor(target_time_output$timepoint, levels = TIME_LEVELS)

target_time_output <- target_time_output %>% arrange(timepoint, summary_level, sample_label)

write_csv_utf8(target_time_output, file.path(STEP65_DIR, paste0("04_65_Clec7a_timepoint_summary_", VERSION,
  ".csv")))

target_region_sample_rows <- list()
counter <- 1L

for (samp in SAMPLE_ORDER) {
    d_samp <- spot_weights[as.character(spot_weights$sample_label) == samp, , drop = FALSE]
    regions <- unique(d_samp$spatial_region_v03)
    regions <- regions[!is.na(regions) & nzchar(regions)]
    for (reg in regions) {
        d <- d_samp[d_samp$spatial_region_v03 == reg, , drop = FALSE]
        x <- d[[TARGET_STATE]]
        target_region_sample_rows[[counter]] <- data.frame(summary_level = "sample_region", timepoint = as.character(unique(d$timepoint)[1]), 
            sample_label = samp, spatial_region_v03 = reg, n_spots = nrow(d), n_replicates = NA_integer_, 
            mean_target_weight = mean(x, na.rm = TRUE), median_target_weight = median(x, na.rm = TRUE), 
            fraction_weight_gt_0_10 = mean(x > 0.1, na.rm = TRUE), sd_of_sample_means = NA_real_, min_sample_mean = NA_real_, 
            max_sample_mean = NA_real_, stringsAsFactors = FALSE)
        counter <- counter + 1L
    }
}

target_region_sample <- bind_rows(target_region_sample_rows)

region_keys <- unique(target_region_sample[, c("timepoint", "spatial_region_v03")])

target_region_time_rows <- list()

for (i in seq_len(nrow(region_keys))) {
    tp <- region_keys$timepoint[i]
    reg <- region_keys$spatial_region_v03[i]
    d <- target_region_sample[target_region_sample$timepoint == tp & target_region_sample$spatial_region_v03 == 
        reg, , drop = FALSE]
    vals <- d$mean_target_weight
    target_region_time_rows[[i]] <- data.frame(summary_level = "timepoint_region_across_replicates", 
        timepoint = tp, sample_label = NA_character_, spatial_region_v03 = reg, n_spots = NA_integer_, 
        n_replicates = length(vals), mean_target_weight = mean(vals, na.rm = TRUE), median_target_weight = NA_real_, 
        fraction_weight_gt_0_10 = NA_real_, sd_of_sample_means = if (length(vals) > 1L) {
            stats::sd(vals, na.rm = TRUE)
        }
        else {
            NA_real_
        }, min_sample_mean = min(vals, na.rm = TRUE), max_sample_mean = max(vals, na.rm = TRUE), stringsAsFactors = FALSE)
}

target_region_time <- bind_rows(target_region_time_rows)

target_region_output <- bind_rows(target_region_sample, target_region_time)

target_region_output$timepoint <- factor(target_region_output$timepoint, levels = TIME_LEVELS)

target_region_output <- target_region_output %>%
  arrange(timepoint, spatial_region_v03, summary_level, sample_label)

write_csv_utf8(target_region_output, file.path(STEP65_DIR, paste0("05_65_Clec7a_region_summary_", VERSION,
  ".csv")))

plot_df <- spot_weights

plot_df$sample_label <- factor(plot_df$sample_label, levels = SAMPLE_ORDER)

p_target_map <- ggplot(plot_df, aes(x = array_x, y = array_y, color = .data[[TARGET_STATE]])) + geom_point(size = 0.55,
  alpha = 0.95) + coord_equal() + scale_y_reverse() + scale_color_viridis_c(option = "magma", limits = c(0,
  1), oob = scales::squish) + facet_wrap(~sample_label, ncol = 4) + labs(title = "Normalized RCTD: Clec7a+ lipid-handling Microglia",
  subtitle = "Per-spot weights normalized to sum to 1; Visium V1 array coordinates", x = "Visium array X",
  y = "Visium array Y", color = "Normalized\nRCTD weight") + theme_classic(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5), strip.text = element_text(size = 8),
    legend.position = "right")

save_pdf_png(p = p_target_map, stem = "06_65_Clec7a_RCTD_spatial_maps", width = 11.7, height = 7.5)

region_plot_df <- target_region_sample

region_plot_df$timepoint <- factor(region_plot_df$timepoint, levels = TIME_LEVELS)

p_target_region <- ggplot(region_plot_df, aes(x = timepoint, y = mean_target_weight, group = sample_label)) +
  geom_point(size = 2.2) + geom_line(linewidth = 0.45, alpha = 0.65) + facet_wrap(~spatial_region_v03,
  scales = "free_y", ncol = 3) + labs(title = "Normalized Clec7a+ RCTD weight by spatial region", subtitle = "Each point is one biological spatial sample; descriptive only",
  x = NULL, y = "Mean normalized RCTD weight") + theme_classic(base_family = "Arial", base_size = 8.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(size = 7.5), plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 8))

save_pdf_png(p = p_target_region, stem = "07_65_Clec7a_region_time_diagnostic", width = 11, height = 9)

myeloid_cols_present <- intersect(MYELOID_STATE_DISPLAY_LEVELS, colnames(spot_weights))

if (length(myeloid_cols_present) != length(MYELOID_STATE_DISPLAY_LEVELS)) {
  stop("Not all 9 finalized myeloid states are present in normalized weights.")
}

myeloid_long <- spot_weights %>% select(sample_label, timepoint, all_of(myeloid_cols_present)) %>% pivot_longer(cols = all_of(myeloid_cols_present), 
    names_to = "myeloid_state", values_to = "weight")

myeloid_sample_mean <- myeloid_long %>% group_by(sample_label, timepoint, myeloid_state) %>% summarise(mean_weight = mean(weight, 
    na.rm = TRUE), .groups = "drop")

myeloid_time_mean <- myeloid_sample_mean %>% group_by(timepoint, myeloid_state) %>% summarise(mean_weight = mean(mean_weight, 
    na.rm = TRUE), .groups = "drop")

myeloid_time_mean$timepoint <- factor(myeloid_time_mean$timepoint, levels = TIME_LEVELS)

myeloid_time_mean$myeloid_state <- factor(myeloid_time_mean$myeloid_state, levels = rev(MYELOID_STATE_DISPLAY_LEVELS))

p_heat <- ggplot(myeloid_time_mean, aes(x = timepoint, y = myeloid_state, fill = mean_weight)) + geom_tile(linewidth = 0.35) +
  scale_fill_viridis_c(option = "magma") + labs(title = "Normalized RCTD weights of the 9 myeloid states",
  subtitle = "Timepoint values are means of the two spatial-sample means", x = NULL, y = NULL, fill = "Mean\nweight") +
  theme_classic(base_family = "Arial", base_size = 9) + theme(axis.text.y = element_text(size = 7.5),
  plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8))

save_pdf_png(p = p_heat, stem = "08_65_all_myeloid_state_timepoint_heatmap", width = 7.8, height = 5.5)

spatial_data <- get_data_layer(spatial_obj, assay = "Spatial")

sanity_genes <- intersect(TARGET_SANITY_GENES, rownames(spatial_data))

if (length(sanity_genes) >= 3L) {
    sanity_score <- Matrix::colMeans(spatial_data[sanity_genes, , drop = FALSE])
    target_weight_named <- setNames(spot_weights[[TARGET_STATE]], spot_weights$barcode)
    sanity_rows <- list()
    for (samp in SAMPLE_ORDER) {
        cells <- rownames(spatial_md)[as.character(spatial_md$sample_label) == samp]
        cells <- intersect(cells, names(target_weight_named))
        sanity_rows[[samp]] <- data.frame(sample_label = samp, timepoint = as.character(unique(spatial_md[cells, 
            "timepoint"])[1]), n_spots = length(cells), n_signature_genes = length(sanity_genes), signature_genes = paste(sanity_genes, 
            collapse = "; "), spearman_RCTD_vs_marker_score = safe_cor(target_weight_named[cells], sanity_score[cells]), 
            stringsAsFactors = FALSE)
    }
    sanity_cor <- bind_rows(sanity_rows)
    write_csv_utf8(sanity_cor, file.path(STEP65_DIR, paste0("08b_65_target_marker_sanity_correlation_", 
        VERSION, ".csv")))
}

weight_type_cols <- setdiff(colnames(spot_weights), c("barcode", "GSM", "sample_label", "timepoint",
  "spatial_region_v03", "array_x", "array_y"))

weight_column_map <- data.frame(RCTD_reference_type = weight_type_cols, Seurat_metadata_column = paste0("RCTDnorm_",
  make.names(weight_type_cols, unique = TRUE)), stringsAsFactors = FALSE)

write_csv_utf8(weight_column_map, file.path(STEP65_DIR, paste0("03b_65_RCTD_normalized_weight_column_map_",
  VERSION, ".csv")))

for (j in seq_len(nrow(weight_column_map))) {
    type_name <- weight_column_map$RCTD_reference_type[j]
    meta_name <- weight_column_map$Seurat_metadata_column[j]
    vals <- setNames(spot_weights[[type_name]], spot_weights$barcode)
    spatial_obj[[meta_name]] <- vals[Cells(spatial_obj)]
}

target_meta_col <- weight_column_map$Seurat_metadata_column[weight_column_map$RCTD_reference_type ==
  TARGET_STATE]

if (length(target_meta_col) != 1L) {
  stop("Could not uniquely identify normalized target metadata column.")
}

spatial_obj@misc$GSE195783_step65 <- list(step = 65L, final_version = VERSION, RCTD_computation_version = "v02",
  method = "Bioconductor spacexr RCTD", rctd_mode = "full", weight_normalization = "per-spot divide by sum across all 19 reference types",
  reference_design = "10 broad non-myeloid classes + 9 finalized myeloid states", target_state = TARGET_STATE,
  target_metadata_column = target_meta_col, RCTD_rerun_in_v04 = FALSE, per_spot_significance_testing = FALSE,
  histology_registered = FALSE)

ANNOTATED_OBJECT_FILE <- file.path(STEP65_DIR, paste0("11_65_GSE195783_RCTD_annotated_spatial_", VERSION,
  ".rds"))

saveRDS(spatial_obj, ANNOTATED_OBJECT_FILE, compress = TRUE)

normalized_sum_ok <- all(abs(sample_qc$normalized_mean_weight_sum - 1) < SUM_TOL) && all(abs(sample_qc$normalized_median_weight_sum -
  1) < SUM_TOL) && all(sample_qc$max_abs_normalized_sum_error < SUM_TOL)

target_nonzero <- any(spot_weights[[TARGET_STATE]] > 0)

target_values_in_object <- spatial_obj[[target_meta_col]][, 1]

final_checks <- data.frame(check = c("All 8 v02 RCTD checkpoint files exist", "All 8 checkpoint objects contain RCTD weights",
  "Reference manifest contains 19 reference types", "All 9 finalized myeloid states are present", "Every sample retains at least 95% of input spots",
  "Every normalized spot weight vector sums to 1 within tolerance", "Clec7a+ lipid-handling Microglia normalized weight is non-zero",
  "Normalized RCTD weights were added to the clean Step64 object", "No RCTD computation was rerun in v04",
  "No per-spot significance testing was performed", "No H&E registration is claimed"), pass = c(all(file.exists(checkpoint_files)),
  length(weights_list) == 8L, nrow(reference_manifest) == 19L, all(MYELOID_STATE_DISPLAY_LEVELS %in%
    colnames(spot_weights)), all(sample_qc$fraction_spots_returned >= 0.95), normalized_sum_ok, target_nonzero,
  sum(!is.na(target_values_in_object)) == nrow(spot_weights), TRUE, TRUE, all(spatial_obj$histology_registered ==
    FALSE)), stringsAsFactors = FALSE)

write_csv_utf8(final_checks, file.path(STEP65_DIR, paste0("09_65_GSE195783_RCTD_FINAL_QC_", VERSION,
  ".csv")))

sample_lines <- apply(sample_qc, 1, function(z) {
    paste0("  ", z[["sample_label"]], ": ", z[["n_spots_RCTD"]], "/", z[["n_spots_input"]], " spots; raw mean sum=", 
        sprintf("%.4f", as.numeric(z[["raw_mean_weight_sum"]])), "; normalized mean sum=", sprintf("%.8f", 
            as.numeric(z[["normalized_mean_weight_sum"]])), "; max abs normalized sum error=", format(as.numeric(z[["max_abs_normalized_sum_error"]]), 
            scientific = TRUE, digits = 3))
})

target_time_report <- target_time_output[target_time_output$summary_level == "timepoint_across_replicates",
  , drop = FALSE]

target_time_lines <- apply(target_time_report, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": normalized mean of sample means=", sprintf("%.4f", as.numeric(z[["mean_weight"]])), 
        " (range ", sprintf("%.4f", as.numeric(z[["min_sample_mean"]])), "-", sprintf("%.4f", as.numeric(z[["max_sample_mean"]])), 
        ")")
})

REPORT_FILE <- file.path(STEP65_DIR, paste0("10_65_GSE195783_RCTD_REPORT_", VERSION, ".txt"))

report_con <- file(REPORT_FILE, open = "wt", encoding = "UTF-8")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = report_con)

close(report_con)

print(target_time_report)
print(final_checks)

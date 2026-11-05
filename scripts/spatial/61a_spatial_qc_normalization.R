rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 1200)
options(future.globals.maxSize = 8 * 1024^3)

set.seed(1234)

ROOT <- PROJECT_ROOT

STEP60_DIR <- file.path(ROOT, "results", "60_GSE195783_Visium_import_coordinate_reconstruction")

INPUT_RDS <- file.path(STEP60_DIR, "05_60_GSE195783_spatial_Seurat_list_v05.rds")

OUT_DIR <- file.path(ROOT, "results", "61_GSE195783_spatial_QC_normalization")

FILTERED_OBJECT_DIR <- file.path(OUT_DIR, "objects_filtered_raw_v01")

NORMALIZED_OBJECT_DIR <- file.path(OUT_DIR, "objects_normalized_v01")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(FILTERED_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(NORMALIZED_OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v01"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

ABS_MIN_NCOUNT <- 1000
ABS_MIN_NFEATURE <- 500
MT_FLOOR <- 25
MT_CEILING <- 35
MAD_MULTIPLIER <- 3

MIN_RETENTION_FOR_AUTO_PROCEED <- 0.6
LOGNORMALIZE_SCALE_FACTOR <- 10000
N_VARIABLE_FEATURES <- 3000

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "patchwork", "sctransform")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step 61.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(sctransform)
})

write_csv_utf8 <- write_csv_utf8

safe_median <- finite_median

safe_mad <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) {
        return(0)
    }
    ans <- stats::mad(x, center = stats::median(x), constant = 1.4826, na.rm = TRUE)
    if (!is.finite(ans) || ans <= 0) {
        iqr_val <- stats::IQR(x, na.rm = TRUE)
        ans <- iqr_val/1.349
    }
    if (!is.finite(ans) || ans < 0) {
        ans <- 0
    }
    ans
}

get_counts_matrix <- function(seu) {
    ans <- tryCatch({
        SeuratObject::LayerData(object = seu, assay = "Spatial", layer = "counts")
    }, error = function(e) {
        NULL
    })
    if (is.null(ans)) {
        ans <- Seurat::GetAssayData(object = seu, assay = "Spatial", slot = "counts")
    }
    ans
}

calc_percent_mt <- function(counts) {
    mt_idx <- grepl("^mt-", rownames(counts), ignore.case = TRUE)
    ncount <- Matrix::colSums(counts)
    if (!any(mt_idx)) {
        return(rep(NA_real_, ncol(counts)))
    }
    100 * Matrix::colSums(counts[mt_idx, , drop = FALSE])/pmax(ncount, 1)
}

calc_thresholds <- function(ncount, nfeature, percent_mt) {
    med_count <- safe_median(ncount)
    mad_count <- safe_mad(ncount)
    med_feature <- safe_median(nfeature)
    mad_feature <- safe_mad(nfeature)
    min_count <- max(ABS_MIN_NCOUNT, med_count - MAD_MULTIPLIER * mad_count)
    min_feature <- max(ABS_MIN_NFEATURE, med_feature - MAD_MULTIPLIER * mad_feature)
    if (all(!is.finite(percent_mt))) {
        med_mt <- NA_real_
        mad_mt <- NA_real_
        max_mt <- Inf
    }
    else {
        med_mt <- safe_median(percent_mt)
        mad_mt <- safe_mad(percent_mt)
        mt_dynamic <- med_mt + MAD_MULTIPLIER * mad_mt
        max_mt <- min(MT_CEILING, max(MT_FLOOR, mt_dynamic))
    }
    data.frame(median_nCount = med_count, MAD_nCount = mad_count, min_nCount = min_count, median_nFeature = med_feature, 
        MAD_nFeature = mad_feature, min_nFeature = min_feature, median_percent_mt = med_mt, MAD_percent_mt = mad_mt, 
        max_percent_mt = max_mt, stringsAsFactors = FALSE)
}

if (!file.exists(INPUT_RDS)) {
  stop("Step 60 object not found:\n", INPUT_RDS)
}

seurat_list_raw <- readRDS(INPUT_RDS)

if (!is.list(seurat_list_raw)) {
  stop("Step 60 input is not a list of Seurat objects.")
}

if (length(seurat_list_raw) != 8L) {
  stop("Expected 8 Seurat objects, found ", length(seurat_list_raw), ".")
}

if (is.null(names(seurat_list_raw)) || any(!nzchar(names(seurat_list_raw)))) {
  stop("Step 60 Seurat list must be named by GSM.")
}

input_validation <- list()

for (gsm in names(seurat_list_raw)) {
    seu <- seurat_list_raw[[gsm]]
    if (!inherits(seu, "Seurat")) {
        stop(gsm, " is not a Seurat object.")
    }
    if (!"Spatial" %in% names(seu@assays)) {
        stop(gsm, ": Spatial assay is missing.")
    }
    md <- seu[[]]
    required_md <- c("GSM", "timepoint", "replicate", "sample_label", "coordinate_system", "histology_registered", 
        "array_x", "array_y")
    missing_md <- setdiff(required_md, colnames(md))
    if (length(missing_md) > 0L) {
        stop(gsm, ": missing metadata column(s): ", paste(missing_md, collapse = ", "))
    }
    if (!all(md$coordinate_system == "Visium_V1_array")) {
        stop(gsm, ": coordinate_system is not uniformly Visium_V1_array.")
    }
    if (any(md$histology_registered)) {
        stop(gsm, ": Step 61 must not assume H&E registration.")
    }
    input_validation[[gsm]] <- data.frame(GSM = gsm, n_spots = ncol(seu), n_features = nrow(seu), timepoint = unique(as.character(md$timepoint))[1], 
        sample_label = unique(as.character(md$sample_label))[1], coordinate_system = unique(as.character(md$coordinate_system))[1], 
        stringsAsFactors = FALSE)
}

input_validation <- bind_rows(input_validation)

input_validation$timepoint <- factor(input_validation$timepoint, levels = TIME_LEVELS)

input_validation <- input_validation %>% arrange(timepoint, sample_label)

write_csv_utf8(input_validation, file.path(OUT_DIR, paste0("00_61_GSE195783_input_validation_", VERSION,
  ".csv")))

spot_qc_list <- list()
threshold_list <- list()

for (gsm in names(seurat_list_raw)) {
    seu <- seurat_list_raw[[gsm]]
    md <- seu[[]]
    counts <- get_counts_matrix(seu)
    if (ncol(counts) != ncol(seu)) {
        stop(gsm, ": count matrix and Seurat object have different spot counts.")
    }
    ncount <- Matrix::colSums(counts)
    nfeature <- Matrix::colSums(counts > 0)
    percent_mt <- calc_percent_mt(counts)
    if (!identical(colnames(counts), colnames(seu))) {
        idx <- match(colnames(seu), colnames(counts))
        if (any(is.na(idx))) {
            stop(gsm, ": count barcode alignment failed.")
        }
        ncount <- ncount[idx]
        nfeature <- nfeature[idx]
        percent_mt <- percent_mt[idx]
    }
    thresholds <- calc_thresholds(ncount = as.numeric(ncount), nfeature = as.numeric(nfeature), percent_mt = as.numeric(percent_mt))
    thresholds$GSM <- gsm
    thresholds$timepoint <- unique(as.character(md$timepoint))[1]
    thresholds$sample_label <- unique(as.character(md$sample_label))[1]
    thresholds <- thresholds[, c("GSM", "timepoint", "sample_label", "median_nCount", "MAD_nCount", "min_nCount", 
        "median_nFeature", "MAD_nFeature", "min_nFeature", "median_percent_mt", "MAD_percent_mt", "max_percent_mt")]
    threshold_list[[gsm]] <- thresholds
    fail_nCount <- ncount < thresholds$min_nCount
    fail_nFeature <- nfeature < thresholds$min_nFeature
    if (is.finite(thresholds$max_percent_mt)) {
        fail_percent_mt <- percent_mt > thresholds$max_percent_mt
        fail_percent_mt[!is.finite(percent_mt)] <- FALSE
    }
    else {
        fail_percent_mt <- rep(FALSE, length(percent_mt))
    }
    pass_qc <- !(fail_nCount | fail_nFeature | fail_percent_mt)
    spot_qc_list[[gsm]] <- data.frame(barcode = colnames(seu), GSM = gsm, timepoint = unique(as.character(md$timepoint))[1], 
        replicate = unique(md$replicate)[1], sample_label = unique(as.character(md$sample_label))[1], 
        array_x = md$array_x, array_y = md$array_y, nCount_Spatial = as.numeric(ncount), nFeature_Spatial = as.numeric(nfeature), 
        percent_mt = as.numeric(percent_mt), min_nCount_threshold = thresholds$min_nCount, min_nFeature_threshold = thresholds$min_nFeature, 
        max_percent_mt_threshold = thresholds$max_percent_mt, fail_nCount = fail_nCount, fail_nFeature = fail_nFeature, 
        fail_percent_mt = fail_percent_mt, pass_qc = pass_qc, stringsAsFactors = FALSE)
}

spot_qc <- bind_rows(spot_qc_list)

threshold_table <- bind_rows(threshold_list)

threshold_table$timepoint <- factor(threshold_table$timepoint, levels = TIME_LEVELS)

threshold_table <- threshold_table %>% arrange(timepoint, sample_label)

spot_qc$timepoint <- factor(spot_qc$timepoint, levels = TIME_LEVELS)

sample_order <- threshold_table$sample_label

spot_qc$sample_label <- factor(spot_qc$sample_label, levels = sample_order)

write_csv_utf8(threshold_table, file.path(OUT_DIR, paste0("01_61_GSE195783_QC_thresholds_", VERSION,
  ".csv")))

spot_qc_file <- file.path(OUT_DIR, paste0("02_61_GSE195783_spot_QC_decisions_", VERSION, ".csv.gz"))

spot_qc_con <- gzfile(spot_qc_file, open = "wt")

write.csv(spot_qc, file = spot_qc_con, row.names = FALSE, na = "")

close(spot_qc_con)

retention_summary <- spot_qc %>% group_by(GSM, timepoint, replicate, sample_label) %>% summarise(n_spots_before = n(), 
    n_pass = sum(pass_qc), n_fail = sum(!pass_qc), retention_rate = n_pass/n_spots_before, fail_nCount = sum(fail_nCount), 
    fail_nFeature = sum(fail_nFeature), fail_percent_mt = sum(fail_percent_mt), .groups = "drop") %>% 
    arrange(timepoint, replicate)

write_csv_utf8(retention_summary, file.path(OUT_DIR, paste0("03_61_GSE195783_QC_retention_summary_",
  VERSION, ".csv")))

plot_long <- bind_rows(data.frame(sample_label = spot_qc$sample_label, metric = "nCount_Spatial", value = spot_qc$nCount_Spatial),
  data.frame(sample_label = spot_qc$sample_label, metric = "nFeature_Spatial", value = spot_qc$nFeature_Spatial),
  data.frame(sample_label = spot_qc$sample_label, metric = "percent_mt", value = spot_qc$percent_mt))

plot_long$sample_label <- factor(plot_long$sample_label, levels = sample_order)

plot_long$metric <- factor(plot_long$metric, levels = c("nCount_Spatial", "nFeature_Spatial", "percent_mt"),
  labels = c("UMI counts", "Detected genes", "Mitochondrial (%)"))

p_pre_violin <- ggplot(plot_long, aes(x = sample_label, y = value)) + geom_violin(scale = "width", trim = TRUE,
  linewidth = 0.25, fill = "grey85") + geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25,
  fill = "white") + facet_wrap(~metric, scales = "free_y", ncol = 1) + labs(title = "GSE195783 spatial spot QC before filtering",
  x = NULL, y = NULL) + theme_classic(base_family = "Arial", base_size = 10) + theme(plot.title = element_text(size = 13,
  family = "Arial"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), strip.text = element_text(size = 9,
  family = "Arial"))

ggsave(filename = file.path(OUT_DIR, paste0("04_61_preQC_distributions_", VERSION, ".pdf")), plot = p_pre_violin,
  width = 8.8, height = 9, units = "in", device = grDevices::cairo_pdf, bg = "white")

ggsave(filename = file.path(OUT_DIR, paste0("04_61_preQC_distributions_", VERSION, ".png")), plot = p_pre_violin,
  width = 8.8, height = 9, units = "in", dpi = 300, bg = "white")

spot_qc$qc_status <- ifelse(spot_qc$pass_qc, "PASS", "FAIL")

spot_qc$qc_status <- factor(spot_qc$qc_status, levels = c("PASS", "FAIL"))

p_spatial_qc <- ggplot(spot_qc, aes(x = array_x, y = array_y, color = qc_status)) + geom_point(size = 0.52,
  alpha = 0.9) + coord_equal() + scale_y_reverse() + scale_color_manual(values = c(PASS = "grey25",
  FAIL = "#D73027")) + facet_wrap(~sample_label, ncol = 4) + labs(title = "GSE195783 spatial distribution of QC decisions",
  subtitle = "Array coordinates only; no H&E/pixel registration", x = "Visium array X", y = "Visium array Y",
  color = NULL) + theme_classic(base_family = "Arial", base_size = 10) + theme(plot.title = element_text(size = 13,
  family = "Arial"), plot.subtitle = element_text(size = 9, family = "Arial"), strip.text = element_text(size = 8.5,
  family = "Arial"), legend.position = "top")

ggsave(filename = file.path(OUT_DIR, paste0("05_61_spatial_QC_pass_fail_", VERSION, ".pdf")), plot = p_spatial_qc,
  width = 11, height = 7.2, units = "in", device = grDevices::cairo_pdf, bg = "white")

ggsave(filename = file.path(OUT_DIR, paste0("05_61_spatial_QC_pass_fail_", VERSION, ".png")), plot = p_spatial_qc,
  width = 11, height = 7.2, units = "in", dpi = 300, bg = "white")

low_retention <- retention_summary[retention_summary$retention_rate < MIN_RETENTION_FOR_AUTO_PROCEED,
  , drop = FALSE]

if (nrow(low_retention) > 0L) {
    stop("QC safety gate triggered: at least one sample retains <", 100 * MIN_RETENTION_FOR_AUTO_PROCEED, 
        "% of spots.\n", "Review files 01-05 before normalization.\nAffected sample(s):\n", paste(paste0(low_retention$sample_label, 
            " = ", sprintf("%.1f%%", 100 * low_retention$retention_rate)), collapse = "\n"))
}

filtered_list <- vector("list", length(seurat_list_raw))

names(filtered_list) <- names(seurat_list_raw)

for (gsm in names(seurat_list_raw)) {

  seu <- seurat_list_raw[[gsm]]

  q <- spot_qc[as.character(spot_qc$GSM) == gsm, , drop = FALSE]

  keep_barcodes <- q$barcode[q$pass_qc]

  if (length(keep_barcodes) < 200L) {
    stop(gsm, ": fewer than 200 spots remain after filtering.")
  }

  seu_filtered <- subset(x = seu, cells = keep_barcodes)

  q_keep <- q[match(colnames(seu_filtered), q$barcode), , drop = FALSE]

  if (any(is.na(q_keep$barcode))) {
    stop(gsm, ": failed to align retained QC metadata.")
  }

  seu_filtered$nCount_Spatial_recomputed <- q_keep$nCount_Spatial

  seu_filtered$nFeature_Spatial_recomputed <- q_keep$nFeature_Spatial

  seu_filtered$percent_mt <- q_keep$percent_mt

  seu_filtered$QC_pass <- TRUE

  seu_filtered$QC_min_nCount <- q_keep$min_nCount_threshold

  seu_filtered$QC_min_nFeature <- q_keep$min_nFeature_threshold

  seu_filtered$QC_max_percent_mt <- q_keep$max_percent_mt_threshold

  seu_filtered@misc$GSE195783_step61 <- list(step = 61L, version = VERSION, qc_method = paste0("Per-sample robust QC: lower nCount/nFeature median-",
    MAD_MULTIPLIER, "*MAD with absolute floors; percent.mt upper median+", MAD_MULTIPLIER, "*MAD with 25-35% bounds"),
    normalization = "not yet applied in filtered raw object", histology_registered = FALSE)

  filtered_list[[gsm]] <- seu_filtered

  saveRDS(seu_filtered, file = file.path(FILTERED_OBJECT_DIR, paste0(gsm, "_filtered_raw_", VERSION,
    ".rds")), compress = TRUE)
}

filtered_list_file <- file.path(OUT_DIR, paste0("06_61_GSE195783_filtered_raw_Seurat_list_", VERSION,
  ".rds"))

saveRDS(filtered_list, file = filtered_list_file, compress = TRUE)

normalized_list <- vector("list", length(filtered_list))

names(normalized_list) <- names(filtered_list)

normalization_audit_list <- list()

for (gsm in names(filtered_list)) {

  seu <- filtered_list[[gsm]]

  seu <- Seurat::NormalizeData(object = seu, assay = "Spatial", normalization.method = "LogNormalize",
    scale.factor = LOGNORMALIZE_SCALE_FACTOR, verbose = FALSE)

  seu <- Seurat::FindVariableFeatures(object = seu, assay = "Spatial", selection.method = "vst", nfeatures = N_VARIABLE_FEATURES,
    verbose = FALSE)

  seu <- Seurat::SCTransform(object = seu, assay = "Spatial", new.assay.name = "SCT", vst.flavor = "v2",
    variable.features.n = N_VARIABLE_FEATURES, return.only.var.genes = TRUE, conserve.memory = TRUE,
    verbose = FALSE)

  Seurat::DefaultAssay(seu) <- "SCT"

  seu@misc$GSE195783_step61 <- list(step = 61L, version = VERSION, QC_filtered = TRUE, Spatial_LogNormalize = TRUE,
    Spatial_scale_factor = LOGNORMALIZE_SCALE_FACTOR, SCT_normalization = TRUE, SCT_vst_flavor = "v2",
    SCT_variable_features = N_VARIABLE_FEATURES, integrated = FALSE, clustered = FALSE, histology_registered = FALSE)

  normalized_list[[gsm]] <- seu

  saveRDS(seu, file = file.path(NORMALIZED_OBJECT_DIR, paste0(gsm, "_QCfiltered_normalized_", VERSION,
    ".rds")), compress = TRUE)

  normalization_audit_list[[gsm]] <- data.frame(GSM = gsm, timepoint = unique(as.character(seu$timepoint))[1],
    sample_label = unique(as.character(seu$sample_label))[1], n_spots = ncol(seu), Spatial_assay_present = "Spatial" %in%
      names(seu@assays), SCT_assay_present = "SCT" %in% names(seu@assays), n_Spatial_features = nrow(seu[["Spatial"]]),
    n_SCT_features = nrow(seu[["SCT"]]), n_Spatial_variable_features = length(Seurat::VariableFeatures(seu,
      assay = "Spatial")), n_SCT_variable_features = length(Seurat::VariableFeatures(seu, assay = "SCT")),
    stringsAsFactors = FALSE)
}

normalized_list_file <- file.path(OUT_DIR, paste0("07_61_GSE195783_QCfiltered_normalized_Seurat_list_",
  VERSION, ".rds"))

saveRDS(normalized_list, file = normalized_list_file, compress = TRUE)

normalization_audit <- bind_rows(normalization_audit_list)

normalization_audit$timepoint <- factor(normalization_audit$timepoint, levels = TIME_LEVELS)

normalization_audit <- normalization_audit %>% arrange(timepoint, sample_label)

write_csv_utf8(normalization_audit, file.path(OUT_DIR, paste0("08_61_GSE195783_normalization_audit_",
  VERSION, ".csv")))

post_qc <- spot_qc[spot_qc$pass_qc, , drop = FALSE]

post_long <- bind_rows(data.frame(sample_label = post_qc$sample_label, metric = "nCount_Spatial", value = post_qc$nCount_Spatial),
  data.frame(sample_label = post_qc$sample_label, metric = "nFeature_Spatial", value = post_qc$nFeature_Spatial),
  data.frame(sample_label = post_qc$sample_label, metric = "percent_mt", value = post_qc$percent_mt))

post_long$sample_label <- factor(post_long$sample_label, levels = sample_order)

post_long$metric <- factor(post_long$metric, levels = c("nCount_Spatial", "nFeature_Spatial", "percent_mt"),
  labels = c("UMI counts", "Detected genes", "Mitochondrial (%)"))

p_post_violin <- ggplot(post_long, aes(x = sample_label, y = value)) + geom_violin(scale = "width", trim = TRUE,
  linewidth = 0.25, fill = "grey85") + geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25,
  fill = "white") + facet_wrap(~metric, scales = "free_y", ncol = 1) + labs(title = "GSE195783 spatial spot QC after filtering",
  x = NULL, y = NULL) + theme_classic(base_family = "Arial", base_size = 10) + theme(plot.title = element_text(size = 13,
  family = "Arial"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), strip.text = element_text(size = 9,
  family = "Arial"))

ggsave(filename = file.path(OUT_DIR, paste0("09_61_postQC_distributions_", VERSION, ".pdf")), plot = p_post_violin,
  width = 8.8, height = 9, units = "in", device = grDevices::cairo_pdf, bg = "white")

ggsave(filename = file.path(OUT_DIR, paste0("09_61_postQC_distributions_", VERSION, ".png")), plot = p_post_violin,
  width = 8.8, height = 9, units = "in", dpi = 300, bg = "white")

final_checks <- data.frame(check = c("Step 60 input contains 8 Seurat objects", "All input objects use Visium V1 array coordinates",
  "No input object claims H&E registration", "All 8 samples retain at least 60% of spots", "All 8 filtered objects contain at least 200 spots",
  "Filtered raw Seurat list contains 8 objects", "Normalized Seurat list contains 8 objects", "Spatial assay retained in all normalized objects",
  "SCT assay created in all normalized objects", "No integration performed in Step 61", "No clustering performed in Step 61"),
  pass = c(length(seurat_list_raw) == 8L, all(input_validation$coordinate_system == "Visium_V1_array"),
    all(vapply(seurat_list_raw, function(x) {
      all(x$histology_registered == FALSE)
    }, logical(1))), all(retention_summary$retention_rate >= MIN_RETENTION_FOR_AUTO_PROCEED), all(vapply(filtered_list,
      ncol, integer(1)) >= 200L), length(filtered_list) == 8L, length(normalized_list) == 8L, all(normalization_audit$Spatial_assay_present),
    all(normalization_audit$SCT_assay_present), TRUE, TRUE), stringsAsFactors = FALSE)

write_csv_utf8(final_checks, file.path(OUT_DIR, paste0("10_61_GSE195783_FINAL_QC_", VERSION, ".csv")))

report_file <- file.path(OUT_DIR, paste0("11_61_GSE195783_QC_normalization_REPORT_", VERSION, ".txt"))

report_con <- file(report_file, open = "wt", encoding = "UTF-8")

threshold_lines <- apply(threshold_table, 1, function(z) {
  paste0("  ", z[["sample_label"]], ": min_nCount=", round(as.numeric(z[["min_nCount"]]), 1), ", min_nFeature=",
    round(as.numeric(z[["min_nFeature"]]), 1), ", max_percent_mt=", round(as.numeric(z[["max_percent_mt"]]),
      2))
})

retention_lines <- apply(retention_summary, 1, function(z) {
  paste0("  ", z[["sample_label"]], ": ", z[["n_pass"]], "/", z[["n_spots_before"]], " retained (",
    sprintf("%.2f%%", 100 * as.numeric(z[["retention_rate"]])), ")")
})

writeLines("Numerical results are provided in the accompanying output tables.", 
    con = report_con)

close(report_con)

print(threshold_table)

print(retention_summary)

print(normalization_audit)

print(final_checks)

message("  05_61_spatial_QC_pass_fail_", VERSION, ".pdf")

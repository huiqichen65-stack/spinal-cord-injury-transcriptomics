source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 1200)

set.seed(1234)

ROOT <- "D:/zz"

STEP61_DIR <- file.path(ROOT, "results", "61_GSE195783_spatial_QC_normalization")

INPUT_RDS <- file.path(STEP61_DIR, "07_61_GSE195783_QCfiltered_normalized_Seurat_list_v01.rds")

OUT_DIR <- file.path(ROOT, "results", "62_GSE195783_joint_embedding_batch_assessment_clustering")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v01"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

N_VARIABLE_FEATURES <- 3000
N_PCS <- 50
DIMS_USE <- 1:30

UMAP_N_NEIGHBORS <- 30
UMAP_MIN_DIST <- 0.3

CLUSTER_RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8)

WORKING_RESOLUTION <- 0.4

detected_cores <- parallel::detectCores(logical = FALSE)

if (is.na(detected_cores) || detected_cores < 1) {
  detected_cores <- 2
}

N_THREADS <- max(1L, min(4L, as.integer(detected_cores - 1L)))

if (N_THREADS < 1L) {
  N_THREADS <- 1L
}

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "patchwork", "irlba", "uwot",
  "RcppAnnoy")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step 62.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(irlba)
  library(uwot)
  library(RcppAnnoy)
})

write_csv_utf8 <- write_csv_utf8

save_pdf_png <- function(p, stem, width, height) {
    pdf_file <- file.path(OUT_DIR, paste0(stem, "_", VERSION, ".pdf"))
    png_file <- file.path(OUT_DIR, paste0(stem, "_", VERSION, ".png"))
    ggsave(filename = pdf_file, plot = p, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    ggsave(filename = png_file, plot = p, width = width, height = height, units = "in", dpi = 300, bg = "white")
}

r2_from_lm <- function(y, x) {
    dat <- data.frame(y = as.numeric(y), x = factor(x))
    dat <- dat[is.finite(dat$y) & !is.na(dat$x), , drop = FALSE]
    if (nrow(dat) < 5L || nlevels(dat$x) < 2L) {
        return(NA_real_)
    }
    fit <- stats::lm(y ~ x, data = dat)
    summary(fit)$r.squared
}

if (!file.exists(INPUT_RDS)) {
  stop("Step 61 normalized object not found:\n", INPUT_RDS)
}

obj_list <- readRDS(INPUT_RDS)

if (!is.list(obj_list)) {
  stop("Input is not a list.")
}

if (length(obj_list) != 8L) {
  stop("Expected 8 Seurat objects, found ", length(obj_list), ".")
}

if (is.null(names(obj_list)) || any(!nzchar(names(obj_list)))) {
  stop("Input list must be named by GSM.")
}

validation_list <- vector("list", length(obj_list))

names(validation_list) <- names(obj_list)

for (gsm in names(obj_list)) {
    seu <- obj_list[[gsm]]
    if (!inherits(seu, "Seurat")) {
        stop(gsm, " is not a Seurat object.")
    }
    if (!"Spatial" %in% names(seu@assays)) {
        stop(gsm, ": Spatial assay is missing.")
    }
    if (!"SCT" %in% names(seu@assays)) {
        stop(gsm, ": Step 61 SCT assay is missing.")
    }
    md <- seu[[]]
    required_md <- c("GSM", "timepoint", "sample_label", "array_x", "array_y", "percent_mt", "coordinate_system", 
        "histology_registered")
    missing_md <- setdiff(required_md, colnames(md))
    if (length(missing_md) > 0L) {
        stop(gsm, ": missing metadata column(s): ", paste(missing_md, collapse = ", "))
    }
    if (!all(md$coordinate_system == "Visium_V1_array")) {
        stop(gsm, ": coordinate system is not Visium_V1_array.")
    }
    if (any(md$histology_registered)) {
        stop(gsm, ": H&E registration must remain FALSE.")
    }
    validation_list[[gsm]] <- data.frame(GSM = gsm, timepoint = unique(as.character(md$timepoint))[1], 
        sample_label = unique(as.character(md$sample_label))[1], n_spots = ncol(seu), n_features = nrow(seu[["Spatial"]]), 
        Spatial_present = TRUE, SCT_present = TRUE, stringsAsFactors = FALSE)
}

validation <- bind_rows(validation_list)

validation$timepoint <- factor(validation$timepoint, levels = TIME_LEVELS)

validation <- validation %>% arrange(timepoint, sample_label)

write_csv_utf8(validation, file.path(OUT_DIR, paste0("00_62_input_validation_", VERSION, ".csv")))

merge_list <- vector("list", length(obj_list))

names(merge_list) <- names(obj_list)

for (gsm in names(obj_list)) {

  seu <- obj_list[[gsm]]

  DefaultAssay(seu) <- "Spatial"

  seu[["SCT"]] <- NULL

  merge_list[[gsm]] <- seu
}

combined <- merge(x = merge_list[[1]], y = merge_list[-1], add.cell.ids = names(merge_list), merge.data = TRUE,
  merge.dr = FALSE)

DefaultAssay(combined) <- "Spatial"

if (inherits(combined[["Spatial"]], "Assay5")) {
    combined <- JoinLayers(object = combined, assay = "Spatial")
}

data_check <- tryCatch({
  LayerData(object = combined, assay = "Spatial", layer = "data")
}, error = function(e) {
  NULL
})

if (is.null(data_check) || nrow(data_check) == 0L || ncol(data_check) == 0L) {
  stop("Merged Spatial normalized data layer is unavailable. ", "Normalization check failed; inspect this error before continuing.")
}

rm(data_check)

gc()

combined$timepoint <- factor(as.character(combined$timepoint), levels = TIME_LEVELS)

sample_levels <- validation$sample_label

combined$sample_label <- factor(as.character(combined$sample_label), levels = sample_levels)

if (any(is.na(combined$timepoint))) {
  stop("Merged object contains unresolved timepoints.")
}

if (any(is.na(combined$sample_label))) {
  stop("Merged object contains unresolved sample_label values.")
}

combined <- FindVariableFeatures(object = combined, assay = "Spatial", selection.method = "vst", nfeatures = N_VARIABLE_FEATURES,
  verbose = FALSE)

if (length(VariableFeatures(combined, assay = "Spatial")) < 1000L) {
  stop("Too few joint variable features were identified.")
}

combined <- ScaleData(object = combined, assay = "Spatial", features = VariableFeatures(combined, assay = "Spatial"),
  verbose = FALSE)

combined <- RunPCA(object = combined, assay = "Spatial", features = VariableFeatures(combined, assay = "Spatial"),
  npcs = N_PCS, approx = TRUE, seed.use = 1234, verbose = FALSE)

pca_stdev <- Stdev(object = combined, reduction = "pca")

pca_var <- pca_stdev^2

pca_pct <- 100 * pca_var/sum(pca_var)

pca_summary <- data.frame(PC = seq_along(pca_pct), variance_percent = pca_pct, cumulative_percent = cumsum(pca_pct),
  stringsAsFactors = FALSE)

write_csv_utf8(pca_summary, file.path(OUT_DIR, paste0("01_62_PCA_variance_", VERSION, ".csv")))

emb <- Embeddings(combined, reduction = "pca")

md <- combined[[]]

pc_batch_list <- vector("list", min(N_PCS, ncol(emb)))

for (i in seq_along(pc_batch_list)) {
    pc_name <- paste0("PC_", i)
    y <- emb[, i]
    r2_time <- r2_from_lm(y = y, x = md$timepoint)
    r2_sample <- r2_from_lm(y = y, x = md$sample_label)
    extra_sample <- r2_sample - r2_time
    if (is.finite(extra_sample) && extra_sample < 0) {
        extra_sample <- 0
    }
    pc_batch_list[[i]] <- data.frame(PC = pc_name, PC_number = i, variance_percent = pca_pct[i], R2_timepoint = r2_time, 
        R2_sample_label = r2_sample, extra_R2_sample_beyond_timepoint = extra_sample, stringsAsFactors = FALSE)
}

pc_batch <- bind_rows(pc_batch_list)

pc_batch$weighted_extra_sample_structure <- pc_batch$variance_percent * pc_batch$extra_R2_sample_beyond_timepoint/100

write_csv_utf8(pc_batch, file.path(OUT_DIR, paste0("02_62_PC_timepoint_vs_sample_structure_", VERSION,
  ".csv")))

combined <- RunUMAP(object = combined, reduction = "pca", dims = DIMS_USE, n.neighbors = UMAP_N_NEIGHBORS,
  min.dist = UMAP_MIN_DIST, metric = "cosine", umap.method = "uwot", seed.use = 1234, n.threads = N_THREADS,
  verbose = FALSE)

combined <- FindNeighbors(object = combined, reduction = "pca", dims = DIMS_USE, nn.method = "annoy",
  annoy.metric = "euclidean", n.trees = 50, verbose = FALSE)

combined <- FindClusters(object = combined, resolution = CLUSTER_RESOLUTIONS, algorithm = 1, random.seed = 1234,
  verbose = FALSE)

working_col <- paste0("Spatial_snn_res.", WORKING_RESOLUTION)

if (!working_col %in% colnames(combined[[]])) {
    possible_cols <- grep(paste0("res\\.", WORKING_RESOLUTION, "$"), colnames(combined[[]]), value = TRUE)
    if (length(possible_cols) != 1L) {
        stop("Could not uniquely identify working clustering column for resolution ", WORKING_RESOLUTION, 
            ".")
    }
    working_col <- possible_cols[1]
}

combined$prelim_cluster_r0.4 <- combined[[working_col]][, 1]

Idents(combined) <- "prelim_cluster_r0.4"

cluster_sample_counts <- as.data.frame(table(cluster = combined$prelim_cluster_r0.4, sample_label = combined$sample_label),
  stringsAsFactors = FALSE)

cluster_sample_counts <- cluster_sample_counts[cluster_sample_counts$Freq > 0, , drop = FALSE]

cluster_sample_props <- cluster_sample_counts %>% group_by(cluster) %>% mutate(proportion_within_cluster = Freq/sum(Freq)) %>% 
    ungroup()

write_csv_utf8(cluster_sample_props, file.path(OUT_DIR, paste0("03_62_cluster_by_sample_composition_r0.4_",
  VERSION, ".csv")))

cluster_time_counts <- as.data.frame(table(cluster = combined$prelim_cluster_r0.4, timepoint = combined$timepoint),
  stringsAsFactors = FALSE)

cluster_time_counts <- cluster_time_counts[cluster_time_counts$Freq > 0, , drop = FALSE]

cluster_time_props <- cluster_time_counts %>% group_by(cluster) %>% mutate(proportion_within_cluster = Freq/sum(Freq)) %>% 
    ungroup()

write_csv_utf8(cluster_time_props, file.path(OUT_DIR, paste0("04_62_cluster_by_timepoint_composition_r0.4_",
  VERSION, ".csv")))

sample_cluster_counts <- as.data.frame(table(sample_label = combined$sample_label, cluster = combined$prelim_cluster_r0.4),
  stringsAsFactors = FALSE)

sample_cluster_counts <- sample_cluster_counts[sample_cluster_counts$Freq > 0, , drop = FALSE]

sample_cluster_props <- sample_cluster_counts %>% group_by(sample_label) %>% mutate(proportion_within_sample = Freq/sum(Freq)) %>% 
    ungroup()

write_csv_utf8(sample_cluster_props, file.path(OUT_DIR, paste0("05_62_sample_cluster_proportions_r0.4_",
  VERSION, ".csv")))

p_elbow <- ElbowPlot(object = combined, ndims = min(N_PCS, 50)) + ggtitle("GSE195783 joint PCA") + theme_classic(base_family = "Arial",
  base_size = 10)

save_pdf_png(p = p_elbow, stem = "06_62_PCA_elbow", width = 5.8, height = 4.6)

p_pca_sample <- DimPlot(object = combined, reduction = "pca", group.by = "sample_label", pt.size = 0.25,
  raster = TRUE) + ggtitle("PCA by sample") + theme_classic(base_family = "Arial", base_size = 10)

p_pca_time <- DimPlot(object = combined, reduction = "pca", group.by = "timepoint", pt.size = 0.25, raster = TRUE) +
  ggtitle("PCA by timepoint") + theme_classic(base_family = "Arial", base_size = 10)

save_pdf_png(p = p_pca_sample + p_pca_time, stem = "07_62_PCA_sample_timepoint", width = 11, height = 5.2)

p_umap_sample <- DimPlot(object = combined, reduction = "umap", group.by = "sample_label", pt.size = 0.28,
  raster = TRUE) + ggtitle("UMAP by sample") + theme_classic(base_family = "Arial", base_size = 10)

p_umap_time <- DimPlot(object = combined, reduction = "umap", group.by = "timepoint", pt.size = 0.28,
  raster = TRUE) + ggtitle("UMAP by timepoint") + theme_classic(base_family = "Arial", base_size = 10)

save_pdf_png(p = p_umap_sample + p_umap_time, stem = "08_62_UMAP_sample_timepoint", width = 11, height = 5.2)

p_umap_cluster <- DimPlot(object = combined, reduction = "umap", group.by = "prelim_cluster_r0.4", label = TRUE,
  repel = TRUE, pt.size = 0.3, raster = TRUE) + ggtitle("Preliminary spatial clusters (resolution 0.4)") +
  theme_classic(base_family = "Arial", base_size = 10)

save_pdf_png(p = p_umap_cluster, stem = "09_62_UMAP_preliminary_clusters_r0.4", width = 6.8, height = 5.5)

p_umap_qc <- FeaturePlot(object = combined, reduction = "umap", features = c("nCount_Spatial", "nFeature_Spatial",
  "percent_mt"), pt.size = 0.25, ncol = 3, raster = TRUE) & theme_classic(base_family = "Arial", base_size = 9)

save_pdf_png(p = p_umap_qc, stem = "10_62_UMAP_QC_covariates", width = 12, height = 4.3)

cluster_sample_props$cluster <- factor(cluster_sample_props$cluster, levels = sort(unique(cluster_sample_props$cluster)))

cluster_sample_props$sample_label <- factor(cluster_sample_props$sample_label, levels = sample_levels)

p_cluster_sample <- ggplot(cluster_sample_props, aes(x = cluster, y = proportion_within_cluster, fill = sample_label)) +
  geom_col(width = 0.85) + scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,
  0)) + labs(title = "Sample composition of preliminary clusters", x = "Preliminary cluster", y = "Composition within cluster",
  fill = "Sample") + theme_classic(base_family = "Arial", base_size = 10) + theme(legend.position = "right")

save_pdf_png(p = p_cluster_sample, stem = "11_62_cluster_sample_composition_r0.4", width = 8.4, height = 5.3)

pc_batch_plot_df <- pc_batch[pc_batch$PC_number <= 30, , drop = FALSE]

pc_batch_plot_df$PC <- factor(pc_batch_plot_df$PC, levels = paste0("PC_", 1:30))

pc_long <- bind_rows(data.frame(PC = pc_batch_plot_df$PC, structure = "Timepoint R2", value = pc_batch_plot_df$R2_timepoint),
  data.frame(PC = pc_batch_plot_df$PC, structure = "Extra sample R2 beyond timepoint", value = pc_batch_plot_df$extra_R2_sample_beyond_timepoint))

p_batch_structure <- ggplot(pc_long, aes(x = PC, y = value, fill = structure)) + geom_col(position = "dodge",
  width = 0.8) + labs(title = "PC structure: biological timepoint vs additional sample effect", x = NULL,
  y = "R²", fill = NULL) + theme_classic(base_family = "Arial", base_size = 10) + theme(axis.text.x = element_text(angle = 60,
  hjust = 1, size = 7), legend.position = "top")

save_pdf_png(p = p_batch_structure, stem = "12_62_PC_batch_structure", width = 9.5, height = 5)

combined@misc$GSE195783_step62 <- list(step = 62L, version = VERSION, input = INPUT_RDS, joint_assay = "Spatial",
  joint_normalization = "Reused Step61 LogNormalize data; no renormalization", integrated = FALSE,
  batch_corrected = FALSE, variable_features = N_VARIABLE_FEATURES, PCA_npcs = N_PCS, UMAP_dims = DIMS_USE,
  clustering_resolutions = CLUSTER_RESOLUTIONS, working_resolution = WORKING_RESOLUTION, working_cluster_column = "prelim_cluster_r0.4",
  histology_registered = FALSE)

OBJECT_FILE <- file.path(OUT_DIR, paste0("13_62_GSE195783_joint_unintegrated_clustered_", VERSION, ".rds"))

saveRDS(combined, file = OBJECT_FILE, compress = TRUE)

res_cols <- grep("res\\.", colnames(combined[[]]), value = TRUE)

final_checks <- data.frame(check = c("Input contains 8 Step61 objects", "Merged object contains all retained Step61 spots",
  "Merged Spatial normalized data layer is available", "Joint variable features >=1000", "PCA contains at least 30 PCs",
  "UMAP embedding exists", "Annoy-based neighbor graph exists", "Four requested clustering resolutions were generated",
  "Working resolution 0.4 cluster column exists", "No batch integration/correction was performed",
  "No H&E registration is claimed"), pass = c(length(obj_list) == 8L, ncol(combined) == sum(validation$n_spots),
  TRUE, length(VariableFeatures(combined, assay = "Spatial")) >= 1000L, ncol(Embeddings(combined, "pca")) >=
    30L, "umap" %in% names(combined@reductions), any(grepl("nn", names(combined@neighbors), ignore.case = TRUE)),
  length(res_cols) >= length(CLUSTER_RESOLUTIONS), "prelim_cluster_r0.4" %in% colnames(combined[[]]),
  isFALSE(combined@misc$GSE195783_step62$integrated), all(combined$histology_registered == FALSE)),
  stringsAsFactors = FALSE)

write_csv_utf8(final_checks, file.path(OUT_DIR, paste0("14_62_GSE195783_FINAL_QC_", VERSION, ".csv")))

n_clusters <- length(unique(combined$prelim_cluster_r0.4))

top_extra <- pc_batch %>% arrange(desc(weighted_extra_sample_structure)) %>% slice_head(n = 10)

top_extra_lines <- apply(top_extra, 1, function(z) {
    paste0("  ", z[["PC"]], ": variance=", sprintf("%.2f%%", as.numeric(z[["variance_percent"]])), ", R2_timepoint=", 
        sprintf("%.3f", as.numeric(z[["R2_timepoint"]])), ", extra_sample_R2=", sprintf("%.3f", as.numeric(z[["extra_R2_sample_beyond_timepoint"]])))
})

report_file <- file.path(OUT_DIR, paste0("15_62_GSE195783_joint_embedding_batch_assessment_REPORT_",
  VERSION, ".txt"))

report_con <- file(report_file, open = "wt", encoding = "UTF-8")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = report_con)

close(report_con)

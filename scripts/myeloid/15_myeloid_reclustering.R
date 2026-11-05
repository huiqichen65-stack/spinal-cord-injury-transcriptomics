source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234

N_HVG <- 3000
N_PCS_COMPUTE <- 50
N_PCS_HARMONY <- 30
MAX_PLOT_CELLS <- 50000

MYELOID_TYPES <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid")

DATASET_LEVELS <- c("17702045", "GSE162610", "GSE196928")

TIME_LEVELS <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "harmony", "ggplot2", "patchwork", "dplyr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\nPlease install them before rerunning this script.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

find_input <- function(expected_path, pattern) {
    if (file.exists(expected_path)) {
        return(normalizePath(expected_path, winslash = "/", mustWork = TRUE))
    }
    hits <- list.files(project_dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Input file not found:\n", basename(expected_path))
    }
    normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

in_file <- find_input(expected_input, "^09_Figure1A_major_celltype_UMAP_object_v04\\.rds$")

out_dir <- file.path(project_dir, "results", "15_myeloid_extract_reintegrate")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "15_myeloid_extract_reintegrate_v01.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

fmt <- format_number

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

md <- obj[[]]

required_meta <- c("major_celltype", "source_dataset")

missing_meta <- setdiff(required_meta, colnames(md))

if (length(missing_meta) > 0) {
  stop("Input atlas is missing required metadata field(s): ", paste(missing_meta, collapse = ", "))
}

capture.output(print(sort(table(obj$major_celltype), decreasing = TRUE)), file = log_file, append = TRUE)

keep_cells <- Cells(obj)[as.character(obj$major_celltype) %in% MYELOID_TYPES]

if (length(keep_cells) == 0) {
  stop("No cells matched the predefined myeloid labels.")
}

myeloid <- subset(obj, cells = keep_cells)

rm(obj, keep_cells, md)
gc(verbose = FALSE)

DefaultAssay(myeloid) <- "RNA"

myeloid$parent_major_celltype <- factor(as.character(myeloid$major_celltype), levels = MYELOID_TYPES)

myeloid$source_dataset <- factor(as.character(myeloid$source_dataset), levels = DATASET_LEVELS)

if (anyNA(myeloid$source_dataset)) {
    bad_ds <- unique(as.character(myeloid$source_dataset))
    stop("Unrecognized source_dataset value(s) detected. ", "Please inspect metadata before continuing.")
}

if ("timepoint_main" %in% colnames(myeloid[[]])) {
    myeloid$timepoint_main <- factor(as.character(myeloid$timepoint_main), levels = TIME_LEVELS)
} else if ("time" %in% colnames(myeloid[[]])) {
    time_raw <- tolower(trimws(as.character(myeloid$time)))
    tp <- rep(NA_character_, length(time_raw))
    tp[grepl("uninj|sham|control|ctrl|^0+d$|^00\\s*d$", time_raw)] <- "Control"
    tp[grepl("^4\\s*h|4hr|4 hour", time_raw)] <- "4 h"
    tp[grepl("^0?1\\s*d|1dpi|1 day", time_raw)] <- "1 dpi"
    tp[grepl("^0?3\\s*d|3dpi|3 day", time_raw)] <- "3 dpi"
    tp[grepl("^0?7\\s*d|7dpi|7 day", time_raw)] <- "7 dpi"
    tp[grepl("^14\\s*d|14dpi|14 day", time_raw)] <- "14 dpi"
    tp[grepl("^28\\s*d|28dpi|28 day", time_raw)] <- "28 dpi"
    tp[grepl("^38\\s*d|38dpi|38 day", time_raw)] <- "38 dpi"
    myeloid$timepoint_main <- factor(tp, levels = TIME_LEVELS)
    rm(time_raw, tp)
} else {
    warning("No timepoint_main/time field found. ", "Step 15 can continue, but time-point diagnostics will be skipped.")
    myeloid$timepoint_main <- factor(rep(NA_character_, ncol(myeloid)), levels = TIME_LEVELS)
}

capture.output(print(table(myeloid$parent_major_celltype, useNA = "ifany")), file = log_file, append = TRUE)

capture.output(print(table(myeloid$source_dataset, useNA = "ifany")), file = log_file, append = TRUE)

capture.output(print(table(myeloid$timepoint_main, useNA = "ifany")), file = log_file, append = TRUE)

remaining_types <- unique(as.character(myeloid$major_celltype))

unexpected_types <- setdiff(remaining_types, MYELOID_TYPES)

unexpected_types <- unexpected_types[!is.na(unexpected_types)]

if (length(unexpected_types) > 0) {
  stop("Unexpected non-target cell type(s) remain: ", paste(unexpected_types, collapse = ", "))
}

myeloid@reductions <- list()
myeloid@graphs <- list()
myeloid@neighbors <- list()

myeloid <- tryCatch(JoinLayers(myeloid, assay = "RNA"), error = function(e) {
  log_msg("JoinLayers skipped/failed safely: ", conditionMessage(e))
  myeloid
})

DefaultAssay(myeloid) <- "RNA"

counts_test <- tryCatch(GetAssayData(myeloid, assay = "RNA", layer = "counts"), error = function(e) NULL)

if (is.null(counts_test)) {
  stop("RNA counts layer could not be accessed after subsetting.")
}

rm(counts_test)

myeloid <- NormalizeData(myeloid, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
  verbose = TRUE)

myeloid <- FindVariableFeatures(myeloid, assay = "RNA", selection.method = "vst", nfeatures = N_HVG,
  verbose = TRUE)

hvg <- VariableFeatures(myeloid)

if (length(hvg) < 1000) {
  stop("Unexpectedly few HVGs were identified: ", length(hvg))
}

myeloid <- ScaleData(myeloid, assay = "RNA", features = hvg, verbose = TRUE)

set.seed(SEED)

myeloid <- RunPCA(myeloid, assay = "RNA", features = hvg, npcs = N_PCS_COMPUTE, reduction.name = "pca_myeloid",
  reduction.key = "mPC_", seed.use = SEED, verbose = TRUE)

if (!"pca_myeloid" %in% names(myeloid@reductions)) {
  stop("RunPCA finished but pca_myeloid was not found.")
}

pca_sd <- Stdev(myeloid, reduction = "pca_myeloid")

pca_var <- pca_sd^2
pca_pct <- pca_var/sum(pca_var) * 100

pca_df <- data.frame(PC = seq_along(pca_sd), stdev = pca_sd, variance = pca_var, variance_percent = pca_pct,
  cumulative_variance_percent = cumsum(pca_pct), used_for_Harmony = seq_along(pca_sd) <= N_PCS_HARMONY,
  stringsAsFactors = FALSE)

write.csv(pca_df, file.path(out_dir, "15_myeloid_PCA_variance_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

set.seed(SEED)

myeloid <- RunHarmony(object = myeloid, group.by.vars = "source_dataset", reduction.use = "pca_myeloid",
  dims.use = 1:N_PCS_HARMONY, reduction.save = "harmony_myeloid30", max.iter.harmony = 20, plot_convergence = FALSE,
  verbose = TRUE)

if (!"harmony_myeloid30" %in% names(myeloid@reductions)) {
  stop("RunHarmony finished but harmony_myeloid30 was not found.")
}

set.seed(SEED)

myeloid <- RunUMAP(object = myeloid, reduction = "harmony_myeloid30", dims = 1:N_PCS_HARMONY, reduction.name = "umap_myeloid_harmony30",
  reduction.key = "mUMAP_", seed.use = SEED, verbose = TRUE)

if (!"umap_myeloid_harmony30" %in% names(myeloid@reductions)) {
  stop("RunUMAP finished but umap_myeloid_harmony30 was not found.")
}

md_m <- myeloid[[]]

summary_total <- data.frame(summary_type = "total", group1 = "myeloid_total", group2 = NA_character_,
  n_cells = ncol(myeloid), stringsAsFactors = FALSE)

summary_dataset <- md_m %>% count(source_dataset, name = "n_cells") %>% transmute(summary_type = "source_dataset", 
    group1 = as.character(source_dataset), group2 = NA_character_, n_cells = n_cells)

summary_parent <- md_m %>% count(parent_major_celltype, name = "n_cells") %>% transmute(summary_type = "parent_major_celltype", 
    group1 = as.character(parent_major_celltype), group2 = NA_character_, n_cells = n_cells)

summary_time <- md_m %>% filter(!is.na(timepoint_main)) %>% count(timepoint_main, name = "n_cells", .drop = FALSE) %>% 
    transmute(summary_type = "timepoint", group1 = as.character(timepoint_main), group2 = NA_character_, 
        n_cells = n_cells)

summary_time_dataset <- md_m %>% filter(!is.na(timepoint_main)) %>% count(timepoint_main, source_dataset, 
    name = "n_cells") %>% transmute(summary_type = "timepoint_by_dataset", group1 = as.character(timepoint_main), 
    group2 = as.character(source_dataset), n_cells = n_cells)

summary_df <- bind_rows(summary_total, summary_dataset, summary_parent, summary_time, summary_time_dataset)

write.csv(summary_df, file.path(out_dir, "15_myeloid_extract_reintegrate_summary_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

set.seed(SEED)

plot_cells <- Cells(myeloid)

if (length(plot_cells) > MAX_PLOT_CELLS) {
  plot_cells <- sample(plot_cells, MAX_PLOT_CELLS)
}

obj_plot <- subset(myeloid, cells = plot_cells)

diagnostic_pdf <- file.path(out_dir, "15_myeloid_extract_reintegrate_diagnostics_v01.pdf")

pdf(diagnostic_pdf, width = 11, height = 8.5, onefile = TRUE)

count_ds <- md_m %>% count(source_dataset, parent_major_celltype, name = "n_cells")

p_count <- ggplot(count_ds, aes(x = source_dataset, y = n_cells, fill = parent_major_celltype)) + geom_col(width = 0.72) +
  labs(title = "Myeloid cells retained from the global atlas", x = NULL, y = "Cell count", fill = "Parent annotation") +
  theme_classic(base_size = 11)

print(p_count)

p_elbow <- ElbowPlot(myeloid, reduction = "pca_myeloid", ndims = N_PCS_COMPUTE) + geom_vline(xintercept = N_PCS_HARMONY,
  linetype = 2) + labs(title = paste0("Myeloid PCA | dashed line = PC ", N_PCS_HARMONY)) + theme_classic(base_size = 11)

p_pca <- DimPlot(obj_plot, reduction = "pca_myeloid", group.by = "source_dataset", pt.size = 0.15, raster = TRUE) +
  ggtitle("Before Harmony: myeloid PCA") + theme_classic(base_size = 11)

print(p_elbow + p_pca + plot_layout(ncol = 2))

p_dataset <- DimPlot(obj_plot, reduction = "umap_myeloid_harmony30", group.by = "source_dataset", pt.size = 0.15,
  raster = TRUE) + ggtitle("Myeloid Harmony UMAP by source dataset") + theme_classic(base_size = 11)

print(p_dataset)

p_parent <- DimPlot(obj_plot, reduction = "umap_myeloid_harmony30", group.by = "parent_major_celltype",
  pt.size = 0.15, raster = TRUE) + ggtitle("Myeloid Harmony UMAP by parent annotation") + theme_classic(base_size = 11)

print(p_parent)

time_cells <- Cells(obj_plot)[!is.na(obj_plot$timepoint_main)]

if (length(time_cells) > 0) {
    obj_time <- subset(obj_plot, cells = time_cells)
    p_time <- DimPlot(obj_time, reduction = "umap_myeloid_harmony30", group.by = "timepoint_main", pt.size = 0.15, 
        raster = TRUE) + ggtitle("Myeloid Harmony UMAP by time point") + theme_classic(base_size = 11)
    print(p_time)
    rm(obj_time)
}

dev.off()

rm(obj_plot, count_ds, md_m)

gc(verbose = FALSE)

out_rds <- file.path(out_dir, "15_myeloid_extract_reintegrate_3scRNA_v01.rds")

saveRDS(myeloid, out_rds, compress = FALSE)

log_msg("HVGs = ", fmt(length(VariableFeatures(myeloid))))

capture.output(sessionInfo(), file = log_file, append = TRUE)

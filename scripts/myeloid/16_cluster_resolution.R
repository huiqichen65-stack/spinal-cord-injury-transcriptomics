source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
RESOLUTIONS <- c(0.1, 0.15, 0.2)
MAX_PLOT_CELLS <- 50000

SNN_GRAPH <- "myeloid30_snn"
UMAP_REDUCTION <- "umap_myeloid_harmony30"
HARMONY_REDUCTION <- "harmony_myeloid30"

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "dplyr", "tidyr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\nPlease install them before rerunning.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected_input <- file.path(project_dir, "results", "16_myeloid_cluster_resolution_assessment", "16_myeloid_cluster_resolution_assessment_3scRNA_v02.rds")

if (file.exists(expected_input)) {
    in_file <- normalizePath(expected_input, winslash = "/", mustWork = TRUE)
} else {
    hits <- list.files(project_dir, pattern = "^16_myeloid_cluster_resolution_assessment_3scRNA_v02\\.rds$", 
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) {
        stop("Input file not found:\n", "16_myeloid_cluster_resolution_assessment_3scRNA_v02.rds")
    }
    in_file <- normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(project_dir, "results", "16_myeloid_cluster_resolution_assessment")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "16_myeloid_cluster_resolution_assessment_v03.log")

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- write_log

res_label <- function(res) {
  if (abs(res - 0.1) < 1e-08)
    return("0.10")
  if (abs(res - 0.15) < 1e-08)
    return("0.15")
  if (abs(res - 0.2) < 1e-08)
    return("0.20")
  format(res, trim = TRUE, scientific = FALSE)
}

cluster_col_name <- function(res) {

  paste0("myeloid_res.", format(res, trim = TRUE, scientific = FALSE))
}

numeric_cluster_levels <- function(x) {
  as.character(sort(unique(as.integer(as.character(x)))))
}

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

if (!SNN_GRAPH %in% names(obj@graphs)) {
  stop("Required SNN graph is missing: ", SNN_GRAPH, "\nPlease use the Step 16 v02 object.")
}

required_reductions <- c(HARMONY_REDUCTION, UMAP_REDUCTION)

missing_reductions <- setdiff(required_reductions, names(obj@reductions))

if (length(missing_reductions) > 0) {
  stop("Missing reduction(s): ", paste(missing_reductions, collapse = ", "))
}

md0 <- obj[[]]

required_meta <- c("source_dataset", "parent_major_celltype")

missing_meta <- setdiff(required_meta, colnames(md0))

if (length(missing_meta) > 0) {
  stop("Missing metadata field(s): ", paste(missing_meta, collapse = ", "))
}

cluster_cols <- character(length(RESOLUTIONS))

for (i in seq_along(RESOLUTIONS)) {
    res <- RESOLUTIONS[i]
    cluster_col <- cluster_col_name(res)
    cluster_cols[i] <- cluster_col
    if (cluster_col %in% colnames(obj[[]])) {
        log_msg("resolution = ", res_label(res), " | existing column reused: ", cluster_col)
    }
    else {
        set.seed(SEED)
        obj <- FindClusters(object = obj, graph.name = SNN_GRAPH, resolution = res, algorithm = 1, random.seed = SEED, 
            cluster.name = cluster_col, verbose = TRUE)
        log_msg("resolution = ", res_label(res), " | new clustering created: ", cluster_col)
    }
    tb <- table(obj[[cluster_col]][, 1])
}

summary_list <- vector("list", length(RESOLUTIONS))

for (i in seq_along(RESOLUTIONS)) {
    res <- RESOLUTIONS[i]
    cluster_col <- cluster_cols[i]
    cl <- as.character(obj[[cluster_col]][, 1])
    tb <- table(cl)
    sizes <- as.numeric(tb)
    summary_list[[i]] <- data.frame(resolution = res, resolution_label = res_label(res), cluster_column = cluster_col, 
        n_clusters = length(tb), min_cluster_cells = min(sizes), q1_cluster_cells = as.numeric(quantile(sizes, 
            0.25, names = FALSE)), median_cluster_cells = median(sizes), mean_cluster_cells = mean(sizes), 
        q3_cluster_cells = as.numeric(quantile(sizes, 0.75, names = FALSE)), max_cluster_cells = max(sizes), 
        n_clusters_lt100 = sum(sizes < 100), n_clusters_lt200 = sum(sizes < 200), n_clusters_lt500 = sum(sizes < 
            500), n_clusters_lt1000 = sum(sizes < 1000), stringsAsFactors = FALSE)
}

summary_df <- bind_rows(summary_list)

write.csv(summary_df, file.path(out_dir, "16_myeloid_cluster_resolution_summary_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

capture.output(print(summary_df), file = log_file, append = TRUE)

composition_list <- list()
cluster_qc_list <- list()

for (i in seq_along(RESOLUTIONS)) {

  res <- RESOLUTIONS[i]
  cluster_col <- cluster_cols[i]

  md_i <- obj[[]]

  md_i$cluster_tmp <- as.character(md_i[[cluster_col]])

  cluster_size <- md_i %>%
    count(cluster_tmp, name = "cluster_n")

  ds <- md_i %>%
    count(cluster_tmp, source_dataset, name = "n_cells") %>%
    left_join(cluster_size, by = "cluster_tmp") %>%
    mutate(fraction = n_cells/cluster_n, resolution = res, resolution_label = res_label(res), variable = "source_dataset",
      category = as.character(source_dataset)) %>%
    select(resolution, resolution_label, cluster = cluster_tmp, cluster_n, variable, category, n_cells,
      fraction)

  composition_list[[paste0("dataset_", res_label(res))]] <- ds

  parent <- md_i %>%
    count(cluster_tmp, parent_major_celltype, name = "n_cells") %>%
    left_join(cluster_size, by = "cluster_tmp") %>%
    mutate(fraction = n_cells/cluster_n, resolution = res, resolution_label = res_label(res), variable = "parent_major_celltype",
      category = as.character(parent_major_celltype)) %>%
    select(resolution, resolution_label, cluster = cluster_tmp, cluster_n, variable, category, n_cells,
      fraction)

  composition_list[[paste0("parent_", res_label(res))]] <- parent

  has_timepoint <- ("timepoint_main" %in% colnames(md_i) && any(!is.na(md_i$timepoint_main)))

  if (has_timepoint) {

    time <- md_i %>%
      filter(!is.na(timepoint_main)) %>%
      count(cluster_tmp, timepoint_main, name = "n_cells") %>%
      left_join(cluster_size, by = "cluster_tmp") %>%
      mutate(fraction = n_cells/cluster_n, resolution = res, resolution_label = res_label(res),
        variable = "timepoint_main", category = as.character(timepoint_main)) %>%
      select(resolution, resolution_label, cluster = cluster_tmp, cluster_n, variable, category,
        n_cells, fraction)

    composition_list[[paste0("time_", res_label(res))]] <- time
  }

  ds_qc <- ds %>%
    group_by(cluster) %>%
    slice_max(order_by = fraction, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(cluster, dominant_dataset = category, dominant_dataset_fraction = fraction)

  parent_qc <- parent %>%
    group_by(cluster) %>%
    slice_max(order_by = fraction, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(cluster, dominant_parent = category, dominant_parent_fraction = fraction)

  if (has_timepoint) {

    time_qc <- time %>%
      group_by(cluster) %>%
      slice_max(order_by = fraction, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(cluster, dominant_timepoint = category, dominant_timepoint_fraction = fraction)

  } else {

    time_qc <- data.frame(cluster = cluster_size$cluster_tmp, dominant_timepoint = NA_character_,
      dominant_timepoint_fraction = NA_real_, stringsAsFactors = FALSE)
  }

  qc_i <- cluster_size %>%
    rename(cluster = cluster_tmp) %>%
    left_join(ds_qc, by = "cluster") %>%
    left_join(parent_qc, by = "cluster") %>%
    left_join(time_qc, by = "cluster") %>%
    mutate(resolution = res, resolution_label = res_label(res)) %>%
    select(resolution, resolution_label, cluster, cluster_n, dominant_parent, dominant_parent_fraction,
      dominant_dataset, dominant_dataset_fraction, dominant_timepoint, dominant_timepoint_fraction)

  cluster_qc_list[[i]] <- qc_i
}

composition_df <- bind_rows(composition_list)

cluster_qc_df <- bind_rows(cluster_qc_list)

write.csv(composition_df, file.path(out_dir, "16_myeloid_cluster_composition_v03.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

write.csv(cluster_qc_df, file.path(out_dir, "16_myeloid_cluster_QC_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

set.seed(SEED)

plot_cells <- Cells(obj)

if (length(plot_cells) > MAX_PLOT_CELLS) {
  plot_cells <- sample(plot_cells, MAX_PLOT_CELLS)
}

obj_plot <- subset(obj, cells = plot_cells)

diagnostic_pdf <- file.path(out_dir, "16_myeloid_cluster_resolution_diagnostics_v03.pdf")

pdf(diagnostic_pdf, width = 12, height = 8.5, onefile = TRUE)

umap_plots <- vector("list", length(RESOLUTIONS))

for (i in seq_along(RESOLUTIONS)) {
    res <- RESOLUTIONS[i]
    cluster_col <- cluster_cols[i]
    umap_plots[[i]] <- DimPlot(obj_plot, reduction = UMAP_REDUCTION, group.by = cluster_col, label = TRUE, 
        repel = TRUE, pt.size = 0.1, raster = TRUE) + ggtitle(paste0("Resolution ", res_label(res))) + 
        theme_classic(base_size = 10) + theme(legend.position = "none")
}

print(wrap_plots(umap_plots, ncol = 3) + plot_annotation(title = "Low-resolution myeloid clustering comparison in the same UMAP space"))

size_plot_df <- cluster_qc_df %>% group_by(resolution_label) %>% mutate(cluster = factor(cluster, levels = numeric_cluster_levels(cluster))) %>% 
    ungroup()

p_size <- ggplot(size_plot_df, aes(x = cluster, y = cluster_n)) + geom_col(width = 0.75) + facet_wrap(~resolution_label,
  scales = "free_x", ncol = 1) + labs(title = "Cluster sizes across low-resolution settings", x = "Cluster",
  y = "Cell count") + theme_classic(base_size = 11)

print(p_size)

summary_long <- summary_df %>%
  select(resolution_label, n_clusters, n_clusters_lt100, n_clusters_lt200, n_clusters_lt500, n_clusters_lt1000) %>%
  pivot_longer(cols = -resolution_label, names_to = "metric", values_to = "value")

p_summary <- ggplot(summary_long, aes(x = factor(resolution_label, levels = c("0.10", "0.15", "0.20")),
  y = value, group = metric)) + geom_line() + geom_point(size = 2.6) + facet_wrap(~metric, scales = "free_y",
  ncol = 2) + labs(title = "Resolution-level cluster-number and small-cluster diagnostics", x = "Resolution",
  y = NULL) + theme_classic(base_size = 11)

print(p_summary)

ds_plot_df <- composition_df %>% filter(variable == "source_dataset") %>% group_by(resolution_label) %>% 
    mutate(cluster = factor(cluster, levels = numeric_cluster_levels(cluster))) %>% ungroup()

p_dataset <- ggplot(ds_plot_df, aes(x = cluster, y = fraction, fill = category)) + geom_col(width = 0.82) +
  facet_wrap(~resolution_label, scales = "free_x", ncol = 1) + scale_y_continuous(limits = c(0, 1),
  labels = function(x) {
    paste0(round(x * 100), "%")
  }, expand = c(0, 0)) + labs(title = "Source-dataset composition of each myeloid cluster", x = "Cluster",
  y = "Within-cluster fraction", fill = "Dataset") + theme_classic(base_size = 10)

print(p_dataset)

parent_plot_df <- composition_df %>% filter(variable == "parent_major_celltype") %>% group_by(resolution_label) %>% 
    mutate(cluster = factor(cluster, levels = numeric_cluster_levels(cluster))) %>% ungroup()

p_parent <- ggplot(parent_plot_df, aes(x = cluster, y = fraction, fill = category)) + geom_col(width = 0.82) +
  facet_wrap(~resolution_label, scales = "free_x", ncol = 1) + scale_y_continuous(limits = c(0, 1),
  labels = function(x) {
    paste0(round(x * 100), "%")
  }, expand = c(0, 0)) + labs(title = "Parent major-cell-type composition of each cluster", x = "Cluster",
  y = "Within-cluster fraction", fill = "Parent annotation") + theme_classic(base_size = 10)

print(p_parent)

time_plot_df <- composition_df %>% filter(variable == "timepoint_main")

if (nrow(time_plot_df) > 0) {
    time_plot_df <- time_plot_df %>% group_by(resolution_label) %>% mutate(cluster = factor(cluster, 
        levels = numeric_cluster_levels(cluster))) %>% ungroup()
    p_time <- ggplot(time_plot_df, aes(x = cluster, y = fraction, fill = category)) + geom_col(width = 0.82) + 
        facet_wrap(~resolution_label, scales = "free_x", ncol = 1) + scale_y_continuous(limits = c(0, 
        1), labels = function(x) {
        paste0(round(x * 100), "%")
    }, expand = c(0, 0)) + labs(title = "Time-point composition of each myeloid cluster", x = "Cluster", 
        y = "Within-cluster fraction", fill = "Time point") + theme_classic(base_size = 10)
    print(p_time)
}

dev.off()

rm(obj_plot, plot_cells)

gc(verbose = FALSE)

out_rds <- file.path(out_dir, "16_myeloid_cluster_resolution_assessment_3scRNA_v03.rds")

saveRDS(obj, out_rds, compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)

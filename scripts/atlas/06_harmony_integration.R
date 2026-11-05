source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
N_PCS <- 30
RESOLUTIONS <- c(0.4, 0.6, 0.8)
MAX_PLOT_CELLS <- 50000

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "harmony", "ggplot2", "patchwork", "data.table")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("缺少R包：", paste(missing_pkgs, collapse = ", "), "\n请先安装后再运行。")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

in_file <- file.path(project_dir, "results", "05_parameter_assessment", "05_parameter_assessment_3scRNA_merged_v01.rds")

out_dir <- file.path(project_dir, "results", "06_Harmony_cluster_resolution")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "06_Harmony_cluster_resolution_3scRNA_v01.log")

if (file.exists(log_file)) file.remove(log_file)

if (!file.exists(in_file)) {
  stop("未找到输入文件：\n", in_file, "\n请确认05步骤已经完成。")
}

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

if (!"pca3000" %in% names(obj@reductions)) {
  stop("对象中没有pca3000 reduction。", "\n请使用05_parameter_assessment_3scRNA_merged_v01.rds。")
}

pca_emb <- Embeddings(obj, reduction = "pca3000")

if (ncol(pca_emb) < N_PCS) {
  stop("pca3000只有", ncol(pca_emb), "个PC，少于计划使用的", N_PCS, "个。")
}

md <- obj[[]]

if (!"source_dataset" %in% colnames(md)) {
    cn <- Cells(obj)
    inferred <- ifelse(grepl("^F177_", cn), "17702045", ifelse(grepl("^GSE162610_", cn), "GSE162610", 
        ifelse(grepl("^GSE196928_", cn), "GSE196928", NA_character_)))
    if (anyNA(inferred)) {
        stop("缺少source_dataset字段，且无法从cell ID完整恢复数据集来源。")
    }
    obj$source_dataset <- inferred
    md <- obj[[]]
}

capture.output(print(table(obj$source_dataset)), file = log_file, append = TRUE)

set.seed(SEED)

obj <- RunHarmony(object = obj, group.by.vars = "source_dataset", reduction.use = "pca3000", dims.use = 1:N_PCS,
  reduction.save = "harmony30", max.iter.harmony = 20, plot_convergence = FALSE, verbose = TRUE)

if (!"harmony30" %in% names(obj@reductions)) {
  stop("RunHarmony完成后未找到harmony30 reduction。")
}

obj <- FindNeighbors(object = obj, reduction = "harmony30", dims = 1:N_PCS, graph.name = c("harmony30_nn",
  "harmony30_snn"), verbose = TRUE)

cluster_cols <- character()

for (res in RESOLUTIONS) {
    res_chr <- format(res, trim = TRUE, scientific = FALSE)
    cluster_name <- paste0("harmony30_res.", res_chr)
    set.seed(SEED)
    obj <- FindClusters(object = obj, graph.name = "harmony30_snn", resolution = res, algorithm = 1, 
        random.seed = SEED, cluster.name = cluster_name, verbose = TRUE)
    cluster_cols <- c(cluster_cols, cluster_name)
    tb <- table(obj[[cluster_name]][, 1])
}

set.seed(SEED)

obj <- RunUMAP(object = obj, reduction = "harmony30", dims = 1:N_PCS, reduction.name = "umap_harmony30",
  reduction.key = "hUMAP_", seed.use = SEED, verbose = TRUE)

summary_list <- list()

for (i in seq_along(RESOLUTIONS)) {
    res <- RESOLUTIONS[i]
    cluster_name <- cluster_cols[i]
    cl <- as.character(obj[[cluster_name]][, 1])
    tb <- table(cl)
    summary_list[[i]] <- data.frame(resolution = res, cluster_column = cluster_name, n_clusters = length(tb), 
        min_cluster_cells = min(tb), median_cluster_cells = as.numeric(median(tb)), max_cluster_cells = max(tb), 
        n_clusters_lt100 = sum(tb < 100), n_clusters_lt200 = sum(tb < 200), n_clusters_lt500 = sum(tb < 
            500), stringsAsFactors = FALSE)
}

summary_df <- do.call(rbind, summary_list)

write.csv(summary_df, file.path(out_dir, "06_Harmony_cluster_resolution_3scRNA_summary_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

set.seed(SEED)

plot_cells <- Cells(obj)

if (length(plot_cells) > MAX_PLOT_CELLS) {
  plot_cells <- sample(plot_cells, MAX_PLOT_CELLS)
}

obj_plot <- subset(obj, cells = plot_cells)

pdf_file <- file.path(out_dir, "06_Harmony_cluster_resolution_3scRNA_diagnostics_v01.pdf")

pdf(pdf_file, width = 11, height = 8.5, onefile = TRUE)

p_before <- DimPlot(obj_plot, reduction = "pca3000", group.by = "source_dataset", pt.size = 0.15, raster = TRUE) +
  ggtitle("Before Harmony: PCA3000") + theme_classic(base_size = 11)

harmony_df <- as.data.frame(Embeddings(obj_plot, reduction = "harmony30")[, 1:2, drop = FALSE])

colnames(harmony_df) <- c("Harmony1", "Harmony2")

harmony_df$source_dataset <- obj_plot$source_dataset

p_after <- ggplot(harmony_df, aes(x = Harmony1, y = Harmony2, color = source_dataset)) + geom_point(size = 0.15,
  alpha = 0.6) + labs(title = "After Harmony: first 2 Harmony dimensions", color = "Dataset") + theme_classic(base_size = 11)

print(p_before + p_after + plot_layout(ncol = 2))

p_dataset <- DimPlot(obj_plot, reduction = "umap_harmony30", group.by = "source_dataset", pt.size = 0.15,
  raster = TRUE) + ggtitle("Harmony UMAP by source dataset") + theme_classic(base_size = 11)

print(p_dataset)

for (i in seq_along(RESOLUTIONS)) {
    res <- RESOLUTIONS[i]
    cluster_name <- cluster_cols[i]
    p <- DimPlot(obj_plot, reduction = "umap_harmony30", group.by = cluster_name, label = TRUE, repel = TRUE, 
        pt.size = 0.12, raster = TRUE) + ggtitle(paste0("Harmony UMAP | resolution = ", res)) + theme_classic(base_size = 11) + 
        theme(legend.position = "none")
    print(p)
}

p_ncluster <- ggplot(summary_df, aes(x = factor(resolution), y = n_clusters, group = 1)) + geom_line(linewidth = 0.7) +
  geom_point(size = 3) + geom_text(aes(label = n_clusters), vjust = -0.7, size = 4) + labs(x = "Resolution",
  y = "Number of clusters", title = "Cluster number across tested resolutions") + theme_classic(base_size = 12)

print(p_ncluster)

dev.off()

rm(obj_plot)
gc(verbose = FALSE)

saveRDS(obj, file.path(out_dir, "06_Harmony_cluster_resolution_3scRNA_v01.rds"), compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)

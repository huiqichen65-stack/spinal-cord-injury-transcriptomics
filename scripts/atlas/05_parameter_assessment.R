source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
N_HVG_MAX <- 3000
N_HVG_COMPARE <- 2000
N_PCS_DIAGNOSTIC <- 50

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "patchwork")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("缺少R包：", paste(missing_pkgs, collapse = ", "), "\n请先安装后再运行。")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
in_dir <- file.path(project_dir, "results", "04_doublet")
out_dir <- file.path(project_dir, "results", "05_parameter_assessment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "05_parameter_assessment_3scRNA_v01.log")
if (file.exists(log_file)) file.remove(log_file)

paths <- c(F177 = file.path(in_dir, "04_doublet_17702045_singlet_v03.rds"), GSE162610 = file.path(in_dir,
  "04_doublet_GSE162610_singlet_v03.rds"), GSE196928 = file.path(in_dir, "04_doublet_GSE196928_Control_singlet_v03.rds"))

if (any(!file.exists(paths))) {
  stop("缺少04双细胞去除后的RDS：\n", paste(paths[!file.exists(paths)], collapse = "\n"))
}

obj177 <- readRDS(paths["F177"])
obj162 <- readRDS(paths["GSE162610"])
obj196 <- readRDS(paths["GSE196928"])

obj177$source_dataset <- "17702045"
obj162$source_dataset <- "GSE162610"
obj196$source_dataset <- "GSE196928"

n177 <- ncol(obj177)
n162 <- ncol(obj162)
n196 <- ncol(obj196)

DefaultAssay(obj177) <- "RNA"
DefaultAssay(obj162) <- "RNA"
DefaultAssay(obj196) <- "RNA"

merged <- merge(x = obj177, y = list(obj162, obj196), add.cell.ids = c("F177", "GSE162610", "GSE196928"),
  merge.data = FALSE)

rm(obj177, obj162, obj196)
gc(verbose = FALSE)

merged <- JoinLayers(merged, assay = "RNA")

DefaultAssay(merged) <- "RNA"

merged@reductions <- list()
merged@graphs <- list()
merged@neighbors <- list()

merged <- NormalizeData(merged, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
  verbose = TRUE)

merged <- FindVariableFeatures(merged, assay = "RNA", selection.method = "vst", nfeatures = N_HVG_MAX,
  verbose = TRUE)

hvg3000 <- VariableFeatures(merged)
hvg2000 <- hvg3000[seq_len(min(N_HVG_COMPARE, length(hvg3000)))]

if (length(hvg3000) < N_HVG_MAX) {
  warning("实际得到的HVG少于", N_HVG_MAX, "：", length(hvg3000))
}

hvg_df <- data.frame(gene = hvg3000, rank = seq_along(hvg3000), in_top2000 = seq_along(hvg3000) <= length(hvg2000),
  stringsAsFactors = FALSE)

write.csv(hvg_df, file.path(out_dir, "05_parameter_assessment_3scRNA_HVG_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

merged <- ScaleData(merged, assay = "RNA", features = hvg2000, verbose = TRUE)

merged <- RunPCA(merged, assay = "RNA", features = hvg2000, npcs = N_PCS_DIAGNOSTIC, reduction.name = "pca2000",
  reduction.key = "PC2000_", seed.use = SEED, verbose = TRUE)

merged <- ScaleData(merged, assay = "RNA", features = hvg3000, verbose = TRUE)

merged <- RunPCA(merged, assay = "RNA", features = hvg3000, npcs = N_PCS_DIAGNOSTIC, reduction.name = "pca3000",
  reduction.key = "PC3000_", seed.use = SEED, verbose = TRUE)

get_pca_stats <- function(obj, reduction_name, hvg_n) {
    sdv <- obj[[reduction_name]]@stdev
    variance <- sdv^2
    pct <- variance/sum(variance) * 100
    cum <- cumsum(pct)
    data.frame(reduction = reduction_name, n_HVG = hvg_n, PC = seq_along(sdv), stdev = sdv, variance_pct_within_first50 = pct, 
        cumulative_pct_within_first50 = cum, stringsAsFactors = FALSE)
}

pca_stats <- rbind(get_pca_stats(merged, "pca2000", length(hvg2000)), get_pca_stats(merged, "pca3000",
  length(hvg3000)))

write.csv(pca_stats, file.path(out_dir, "05_parameter_assessment_3scRNA_PCA_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

md <- merged[[]]

dataset_summary <- as.data.frame(table(md$source_dataset), stringsAsFactors = FALSE)
colnames(dataset_summary) <- c("source_dataset", "n_cells")
dataset_summary$percent <- (dataset_summary$n_cells/sum(dataset_summary$n_cells) * 100)

time_candidates <- c("timepoint", "time", "Time", "group_time", "injury_time", "day", "dpi")
time_field <- time_candidates[time_candidates %in% colnames(md)][1]

if (length(time_field) == 0 || is.na(time_field)) {
  time_field <- NA_character_
} else {
}

summary_rows <- dataset_summary
summary_rows$category <- "dataset"
summary_rows$group <- summary_rows$source_dataset
summary_rows$source_dataset <- NULL

summary_rows <- summary_rows[, c("category", "group", "n_cells", "percent")]

if (!is.na(time_field)) {
    tt <- as.data.frame(table(dataset = md$source_dataset, time = md[[time_field]]), stringsAsFactors = FALSE)
    tt <- tt[tt$Freq > 0, , drop = FALSE]
    tt2 <- data.frame(category = "dataset_timepoint", group = paste(tt$dataset, tt$time, sep = " | "), 
        n_cells = tt$Freq, percent = tt$Freq/nrow(md) * 100, stringsAsFactors = FALSE)
    summary_rows <- rbind(summary_rows, tt2)
}

write.csv(summary_rows, file.path(out_dir, "05_parameter_assessment_3scRNA_summary_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

pdf_file <- file.path(out_dir, "05_parameter_assessment_3scRNA_diagnostics_v01.pdf")

pdf(pdf_file, width = 11, height = 8.5, onefile = TRUE)

p_count <- ggplot(dataset_summary, aes(x = source_dataset, y = n_cells)) + geom_col() + geom_text(aes(label = scales::comma(n_cells)),
  vjust = -0.4, size = 4) + labs(x = NULL, y = "Singlet cell count", title = "Main-atlas cells after DecontX and doublet removal") +
  theme_classic(base_size = 12)

print(p_count)

p_elbow <- ggplot(pca_stats, aes(x = PC, y = stdev, linetype = factor(n_HVG))) + geom_line(linewidth = 0.7) +
  geom_point(size = 1.1) + scale_x_continuous(breaks = seq(0, N_PCS_DIAGNOSTIC, 5)) + labs(x = "PC",
  y = "Standard deviation", linetype = "HVG number", title = "PCA elbow: 2000 vs 3000 HVGs") + theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(p_elbow)

p_var <- ggplot(pca_stats, aes(x = PC, y = variance_pct_within_first50, linetype = factor(n_HVG))) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) + scale_x_continuous(breaks = seq(0, N_PCS_DIAGNOSTIC,
  5)) + labs(x = "PC", y = "Variance contribution within first 50 PCs (%)", linetype = "HVG number",
  title = "Variance contribution by PC") + theme_classic(base_size = 12) + theme(legend.position = "top")

print(p_var)

p_cum <- ggplot(pca_stats, aes(x = PC, y = cumulative_pct_within_first50, linetype = factor(n_HVG))) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) + scale_x_continuous(breaks = seq(0, N_PCS_DIAGNOSTIC,
  5)) + labs(x = "PC", y = "Cumulative contribution within first 50 PCs (%)", linetype = "HVG number",
  title = "Cumulative PCA contribution") + theme_classic(base_size = 12) + theme(legend.position = "top")

print(p_cum)

emb2000 <- Embeddings(merged, reduction = "pca2000")
df2000 <- data.frame(PC1 = emb2000[, 1], PC2 = emb2000[, 2], source_dataset = md[colnames(merged), "source_dataset"],
  stringsAsFactors = FALSE)

set.seed(SEED)
if (nrow(df2000) > 40000) {
  idx <- sample(seq_len(nrow(df2000)), 40000)
  df2000_plot <- df2000[idx, , drop = FALSE]
} else {
  df2000_plot <- df2000
}

p_pca2000 <- ggplot(df2000_plot, aes(x = PC1, y = PC2, color = source_dataset)) + geom_point(size = 0.25,
  alpha = 0.5) + labs(title = "PCA before Harmony: Top 2000 HVGs", color = "Dataset") + theme_classic(base_size = 11) +
  theme(legend.position = "top")

print(p_pca2000)

emb3000 <- Embeddings(merged, reduction = "pca3000")
df3000 <- data.frame(PC1 = emb3000[, 1], PC2 = emb3000[, 2], source_dataset = md[colnames(merged), "source_dataset"],
  stringsAsFactors = FALSE)

set.seed(SEED)
if (nrow(df3000) > 40000) {
  idx <- sample(seq_len(nrow(df3000)), 40000)
  df3000_plot <- df3000[idx, , drop = FALSE]
} else {
  df3000_plot <- df3000
}

p_pca3000 <- ggplot(df3000_plot, aes(x = PC1, y = PC2, color = source_dataset)) + geom_point(size = 0.25,
  alpha = 0.5) + labs(title = "PCA before Harmony: Top 3000 HVGs", color = "Dataset") + theme_classic(base_size = 11) +
  theme(legend.position = "top")

print(p_pca3000)

dev.off()

saveRDS(merged, file.path(out_dir, "05_parameter_assessment_3scRNA_merged_v01.rds"), compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)

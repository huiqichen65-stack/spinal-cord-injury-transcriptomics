source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234
MIN_FEATURES_G196 <- 300
ROUND_DECONTX_COUNTS <- TRUE

F177_BATCH_MODE <- "whole_dataset"

set.seed(SEED)

cran_pkgs <- c("Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "patchwork")

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_cran) > 0) {
  stop("缺少R包：", paste(missing_cran, collapse = ", "), "\n请先安装：install.packages(c(",
    paste(sprintf("\"%s\"", missing_cran), collapse = ", "), "))")
}

if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  stop("缺少SingleCellExperiment。\n", "请先运行：\n", "if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")\n",
    "BiocManager::install(\"SingleCellExperiment\")")
}

has_decontX <- requireNamespace("decontX", quietly = TRUE)
has_celda <- requireNamespace("celda", quietly = TRUE)

if (!has_decontX && !has_celda) {
  stop("未安装DecontX实现。\n", "推荐安装当前Bioconductor版本：\n", "if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")\n",
    "BiocManager::install(\"decontX\")\n", "若你的环境使用旧版celda，也可安装：BiocManager::install(\"celda\")")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(SingleCellExperiment)
})

decontx_backend <- if (has_decontX) "decontX" else "celda"

run_decontx_core <- function(sce, seed = SEED) {
  if (decontx_backend == "decontX") {
    decontX::decontX(sce, assayName = "counts", seed = seed, verbose = TRUE)
  } else {
    celda::decontX(sce, assayName = "counts", seed = seed, verbose = TRUE)
  }
}

get_decontx_counts <- function(sce) {
  if (decontx_backend == "decontX") {
    decontX::decontXcounts(sce)
  } else {
    celda::decontXcounts(sce)
  }
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_dir, "data")
qc_dir <- file.path(project_dir, "results", "02_QC")
out_dir <- file.path(project_dir, "results", "03_DecontX")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "03_DecontX_3scRNA_v01.log")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- write_log

if (decontx_backend == "decontX") {
  log_msg("decontX version = ", as.character(utils::packageVersion("decontX")))
} else {
  log_msg("celda version = ", as.character(utils::packageVersion("celda")))
}
log_msg("Seurat version = ", as.character(utils::packageVersion("Seurat")))

get_rna_counts <- function(obj) {
    DefaultAssay(obj) <- "RNA"
    x <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
    if (is.null(x)) {
        x <- tryCatch(GetAssayData(obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
    }
    if (is.null(x)) {
        obj2 <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) NULL)
        if (!is.null(obj2)) {
            x <- tryCatch(GetAssayData(obj2, assay = "RNA", layer = "counts"), error = function(e) NULL)
        }
    }
    if (is.null(x)) 
        stop("无法从Seurat对象提取RNA counts。")
    as(x, "dgCMatrix")
}

extract_16bp <- function(x) {
  x <- as.character(x)
  m <- regexpr("[ACGT]{16}", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- m > 0
  out[ok] <- regmatches(x, m)[ok]
  out
}

read_mtx_triplet <- function(matrix_file, barcode_file, feature_file) {
    mat <- Matrix::readMM(gzfile(matrix_file))
    mat <- as(mat, "dgCMatrix")
    bc <- data.table::fread(barcode_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
    ft <- data.table::fread(feature_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
    barcodes <- as.character(bc[[1]])
    if (ncol(ft) >= 2) {
        genes <- as.character(ft[[2]])
        bad <- is.na(genes) | genes == ""
        genes[bad] <- as.character(ft[[1]])[bad]
    }
    else {
        genes <- as.character(ft[[1]])
    }
    genes <- make.unique(genes)
    stopifnot(nrow(mat) == length(genes), ncol(mat) == length(barcodes))
    rownames(mat) <- genes
    colnames(mat) <- barcodes
    mat
}

contam_to_pct <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0)
    return(x)

  if (max(x, na.rm = TRUE) <= 1.000001) {
    x * 100
  } else {
    x
  }
}

run_decontx_one <- function(obj, dataset_name, sample_name) {
  counts <- get_rna_counts(obj)

  keep_gene <- Matrix::rowSums(counts) > 0
  counts_run <- counts[keep_gene, , drop = FALSE]

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_run))

  sce <- run_decontx_core(sce, seed = SEED)

  corrected <- get_decontx_counts(sce)
  corrected <- as(corrected, "dgCMatrix")

  if (ROUND_DECONTX_COUNTS) {
    corrected@x <- round(corrected@x)
    corrected <- Matrix::drop0(corrected)
  }

  contamination <- as.numeric(SummarizedExperiment::colData(sce)$decontX_contamination)

  clusters <- NULL
  if ("decontX_clusters" %in% colnames(SummarizedExperiment::colData(sce))) {
    clusters <- as.character(SummarizedExperiment::colData(sce)$decontX_clusters)
  }

  md <- obj[[]]
  md <- md[colnames(counts_run), , drop = FALSE]

  if ("nCount_RNA" %in% colnames(md)) {
    md$nCount_RNA_before_DecontX <- md$nCount_RNA
    md$nCount_RNA <- NULL
  } else {
    md$nCount_RNA_before_DecontX <- Matrix::colSums(counts_run)
  }

  if ("nFeature_RNA" %in% colnames(md)) {
    md$nFeature_RNA_before_DecontX <- md$nFeature_RNA
    md$nFeature_RNA <- NULL
  } else {
    md$nFeature_RNA_before_DecontX <- Matrix::colSums(counts_run > 0)
  }

  md$decontX_contamination <- contamination
  md$decontX_contamination_pct <- contam_to_pct(contamination)
  md$decontX_run_unit <- sample_name

  if (!is.null(clusters)) {
    md$decontX_cluster <- clusters
  }

  raw_umi <- Matrix::colSums(counts_run)
  corrected_umi <- Matrix::colSums(corrected)

  md$decontX_removed_UMI <- pmax(raw_umi - corrected_umi, 0)
  md$decontX_removed_UMI_pct <- ifelse(raw_umi > 0, md$decontX_removed_UMI/raw_umi * 100, 0)

  new_obj <- CreateSeuratObject(counts = corrected, meta.data = md, project = paste0(dataset_name,
    "_DecontX"), min.cells = 0, min.features = 0)

  rm(sce, counts, counts_run, corrected)
  gc(verbose = FALSE)

  new_obj
}

split_run_merge <- function(obj, dataset_name, sample_field) {
    md <- obj[[]]
    if (!sample_field %in% colnames(md)) {
        stop(dataset_name, " 缺少sample字段：", sample_field)
    }
    samples <- unique(as.character(md[[sample_field]]))
    samples <- samples[!is.na(samples) & nzchar(samples)]
    if (length(samples) == 0) {
        stop(dataset_name, " 没有可用sample标签。")
    }
    res <- vector("list", length(samples))
    names(res) <- samples
    for (i in seq_along(samples)) {
        s <- samples[i]
        cells <- rownames(md)[as.character(md[[sample_field]]) == s]
        one <- subset(obj, cells = cells)
        res[[s]] <- run_decontx_one(one, dataset_name = dataset_name, sample_name = s)
        rm(one)
        gc(verbose = FALSE)
    }
    if (length(res) == 1) 
        return(res[[1]])
    merged <- merge(x = res[[1]], y = res[-1], merge.data = FALSE)
    rm(res)
    gc(verbose = FALSE)
    merged
}

summarize_decontx <- function(obj) {
    md <- obj[[]]
    needed <- c("dataset", "sample_id", "decontX_contamination_pct", "decontX_removed_UMI_pct")
    for (nm in needed) {
        if (!nm %in% colnames(md)) {
            if (nm == "dataset") 
                md$dataset <- "unknown"
            if (nm == "sample_id") 
                md$sample_id <- "unknown"
        }
    }
    key <- interaction(md$dataset, md$sample_id, drop = TRUE, lex.order = TRUE)
    ids <- split(seq_len(nrow(md)), key)
    do.call(rbind, lapply(ids, function(i) {
        raw_umi <- md$nCount_RNA_before_DecontX[i]
        after_umi <- md$nCount_RNA[i]
        data.frame(dataset = as.character(md$dataset[i[1]]), sample_id = as.character(md$sample_id[i[1]]), 
            time = if ("time" %in% colnames(md)) {
                as.character(md$time[i[1]])
            }
            else {
                NA_character_
            }, treatment = if ("treatment" %in% colnames(md)) {
                as.character(md$treatment[i[1]])
            }
            else {
                NA_character_
            }, n_cells = length(i), contamination_median_pct = median(md$decontX_contamination_pct[i], 
                na.rm = TRUE), contamination_p95_pct = unname(quantile(md$decontX_contamination_pct[i], 
                0.95, na.rm = TRUE)), contamination_max_pct = max(md$decontX_contamination_pct[i], na.rm = TRUE), 
            raw_total_UMI = sum(raw_umi, na.rm = TRUE), corrected_total_UMI = sum(after_umi, na.rm = TRUE), 
            removed_UMI_pct = (sum(raw_umi, na.rm = TRUE) - sum(after_umi, na.rm = TRUE))/max(sum(raw_umi, 
                na.rm = TRUE), 1) * 100, stringsAsFactors = FALSE)
    }))
}

f177 <- file.path(qc_dir, "02_QC_17702045_v01.rds")
if (!file.exists(f177)) {
  stop("未找到QC对象：", f177)
}

obj177_qc <- readRDS(f177)

if (F177_BATCH_MODE == "timepoint") {
  warning("17702045当前没有可靠真实library字段；你手动选择了timepoint模式。", "请仅在确认timepoint对应独立library时使用。")

  if (!"time" %in% colnames(obj177_qc[[]])) {
    stop("17702045对象中没有time字段。")
  }

  obj177_qc$decontX_batch <- paste0("F177_", obj177_qc$time)

  obj177_dc <- split_run_merge(obj177_qc, dataset_name = "17702045", sample_field = "decontX_batch")
} else {

  obj177_dc <- run_decontx_one(obj177_qc, dataset_name = "17702045", sample_name = "whole_dataset")
}

rm(obj177_qc)
gc(verbose = FALSE)

f162 <- file.path(qc_dir, "02_QC_GSE162610_v01.rds")
if (!file.exists(f162)) {
  stop("未找到QC对象：", f162)
}

obj162_qc <- readRDS(f162)

if (!"sample_id" %in% colnames(obj162_qc[[]])) {
  stop("GSE162610 QC对象缺少sample_id。")
}

obj162_dc <- split_run_merge(obj162_qc, dataset_name = "GSE162610", sample_field = "sample_id")

rm(obj162_qc)
gc(verbose = FALSE)

d196 <- file.path(data_dir, "GSE196928")
meta_file196 <- file.path(d196, "GSE196928_metadata.csv.gz")
raw_tar196 <- file.path(d196, "GSE196928_RAW.tar")
raw_dir196 <- file.path(d196, "RAW_extracted")

if (!file.exists(meta_file196)) {
  stop("未找到：", meta_file196)
}

if (!dir.exists(raw_dir196)) {
  if (!file.exists(raw_tar196)) {
    stop("未找到：", raw_tar196)
  }
  dir.create(raw_dir196, recursive = TRUE, showWarnings = FALSE)
  untar(raw_tar196, exdir = raw_dir196)
}

meta196 <- data.table::fread(meta_file196, data.table = FALSE, check.names = FALSE, showProgress = FALSE)

for (nm in c("V1", "orig.ident")) {
  if (!nm %in% colnames(meta196)) {
    stop("GSE196928 metadata缺少字段：", nm)
  }
}

meta196$barcode16 <- extract_16bp(meta196$V1)

if (anyNA(meta196$barcode16)) {
  stop("GSE196928 metadata有 ", sum(is.na(meta196$barcode16)), " 个细胞无法提取16bp barcode。")
}

if ("nCount_RNA" %in% colnames(meta196)) {
  meta196$source_nCount_RNA <- meta196$nCount_RNA
  meta196$nCount_RNA <- NULL
}
if ("nFeature_RNA" %in% colnames(meta196)) {
  meta196$source_nFeature_RNA <- meta196$nFeature_RNA
  meta196$nFeature_RNA <- NULL
} else {
  stop("GSE196928 processed metadata中没有nFeature_RNA，", "无法按>=300完成本项目的可验证QC。")
}

sample196 <- data.frame(gsm = c("GSM5904825", "GSM5904826", "GSM5904827", "GSM5904828", "GSM5904829",
  "GSM5904830"), raw_id = c("B1", "B2", "B3", "B4", "B7", "B8"), orig_ident = c("00 d Control", "00 d PLX5622",
  "07 d Control", "07 d PLX5622", "28 d Control", "28 d PLX5622"), time = c("Uninjured", "Uninjured",
  "7d", "7d", "28d", "28d"), treatment = c("Control", "PLX5622", "Control", "PLX5622", "Control", "PLX5622"),
  stringsAsFactors = FALSE)

obj196_list <- vector("list", nrow(sample196))
names(obj196_list) <- sample196$raw_id

for (i in seq_len(nrow(sample196))) {
  gsm <- sample196$gsm[i]
  bid <- sample196$raw_id[i]
  group_name <- sample196$orig_ident[i]
  prefix <- paste0(gsm, "_", bid)

  matrix_file <- file.path(raw_dir196, paste0(prefix, "_matrix.mtx.gz"))
  barcode_file <- file.path(raw_dir196, paste0(prefix, "_barcodes.tsv.gz"))
  feature_file <- file.path(raw_dir196, paste0(prefix, "_features.tsv.gz"))

  for (f in c(matrix_file, barcode_file, feature_file)) {
    if (!file.exists(f))
      stop("未找到：", f)
  }

  mat <- read_mtx_triplet(matrix_file, barcode_file, feature_file)

  raw_bc16 <- extract_16bp(colnames(mat))

  target <- meta196[as.character(meta196$orig.ident) == group_name, , drop = FALSE]

  if (nrow(target) == 0) {
    stop("metadata没有找到组：", group_name)
  }

  if (anyDuplicated(target$barcode16)) {
    stop("GSE196928组内出现重复16bp barcode：", group_name, "。不能安全映射。")
  }

  if (anyDuplicated(raw_bc16)) {
    stop("RAW样本内出现重复16bp barcode：", bid, "。不能安全映射。")
  }

  idx <- match(target$barcode16, raw_bc16)

  if (anyNA(idx)) {
    stop(bid, " 有 ", sum(is.na(idx)), " 个processed metadata细胞无法在RAW MTX中找到。")
  }

  mat <- mat[, idx, drop = FALSE]

  new_cell_id <- paste0("G196_", bid, "_", target$barcode16)
  colnames(mat) <- new_cell_id

  md <- target
  rownames(md) <- new_cell_id

  md$dataset <- "GSE196928"
  md$sample_id <- bid
  md$sample_id_source <- "RAW_filename+processed_metadata"
  md$time <- sample196$time[i]
  md$condition <- ifelse(md$time == "Uninjured", "Uninjured", "SCI")
  md$treatment <- sample196$treatment[i]
  md$cell_id_original <- as.character(target$V1)

  keep <- as.numeric(md$source_nFeature_RNA) >= MIN_FEATURES_G196

  mat <- mat[, keep, drop = FALSE]
  md <- md[keep, , drop = FALSE]

  obj196_list[[bid]] <- CreateSeuratObject(counts = mat, meta.data = md, project = paste0("G196_",
    bid), min.cells = 0, min.features = 0)

  rm(mat, md, target, idx)
  gc(verbose = FALSE)
}

obj196_dc_list <- vector("list", length(obj196_list))
names(obj196_dc_list) <- names(obj196_list)

for (bid in names(obj196_list)) {
  obj196_dc_list[[bid]] <- run_decontx_one(obj196_list[[bid]], dataset_name = "GSE196928", sample_name = bid)
}

control_ids <- c("B1", "B3", "B7")
plx_ids <- c("B2", "B4", "B8")

obj196_control_dc <- merge(obj196_dc_list[[control_ids[1]]], y = obj196_dc_list[control_ids[-1]], merge.data = FALSE)

obj196_plx_dc <- merge(obj196_dc_list[[plx_ids[1]]], y = obj196_dc_list[plx_ids[-1]], merge.data = FALSE)

rm(obj196_list, obj196_dc_list, meta196)
gc(verbose = FALSE)

summary_df <- rbind(summarize_decontx(obj177_dc), summarize_decontx(obj162_dc), summarize_decontx(obj196_control_dc),
  summarize_decontx(obj196_plx_dc))

summary_df <- summary_df[order(summary_df$dataset, summary_df$treatment, summary_df$time, summary_df$sample_id),
  , drop = FALSE]

write.csv(summary_df, file.path(out_dir, "03_DecontX_3scRNA_summary_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

print(summary_df)

collect_plot_md <- function(obj) {
  md <- obj[[]]
  keep <- intersect(c("dataset", "sample_id", "time", "treatment", "decontX_contamination_pct", "decontX_removed_UMI_pct"),
    colnames(md))
  md[, keep, drop = FALSE]
}

plot_md <- rbind(collect_plot_md(obj177_dc), collect_plot_md(obj162_dc), collect_plot_md(obj196_control_dc),
  collect_plot_md(obj196_plx_dc))

plot_md$dataset <- factor(plot_md$dataset, levels = c("17702045", "GSE162610", "GSE196928"))

p1 <- ggplot(plot_md, aes(x = decontX_contamination_pct, fill = dataset)) + geom_histogram(bins = 60,
  alpha = 0.55, position = "identity") + facet_wrap(~dataset, scales = "free_y", ncol = 1) + labs(x = "DecontX estimated contamination (%)",
  y = "Cell count", title = "Ambient RNA contamination across datasets") + theme_classic(base_size = 11) +
  theme(legend.position = "none")

plot_md$sample_label <- paste(plot_md$dataset, plot_md$sample_id, sep = " | ")

p2 <- ggplot(plot_md, aes(x = sample_label, y = decontX_contamination_pct)) + geom_boxplot(outlier.size = 0.2,
  linewidth = 0.35) + labs(x = NULL, y = "DecontX estimated contamination (%)", title = "DecontX contamination by sample/library") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

p3 <- ggplot(plot_md, aes(x = sample_label, y = decontX_removed_UMI_pct)) + geom_boxplot(outlier.size = 0.2,
  linewidth = 0.35) + labs(x = NULL, y = "UMI removed by DecontX (%)", title = "Fraction of counts removed by DecontX") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

pdf(file.path(out_dir, "03_DecontX_3scRNA_contamination_v01.pdf"), width = 11, height = 8.5)
print(p1)
print(p2)
print(p3)
dev.off()

saveRDS(obj177_dc, file.path(out_dir, "03_DecontX_17702045_v01.rds"), compress = FALSE)

saveRDS(obj162_dc, file.path(out_dir, "03_DecontX_GSE162610_v01.rds"), compress = FALSE)

saveRDS(obj196_control_dc, file.path(out_dir, "03_DecontX_GSE196928_Control_v01.rds"), compress = FALSE)

saveRDS(obj196_plx_dc, file.path(out_dir, "03_DecontX_GSE196928_PLX5622_v01.rds"), compress = FALSE)

capture.output(sessionInfo(), file = log_file, append = TRUE)

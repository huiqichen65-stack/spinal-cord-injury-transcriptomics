source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

MIN_FEATURES <- 300
MAX_MT <- 10
SEED <- 1234
set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "patchwork")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("缺少R包：", paste(missing_pkgs, collapse = ", "), "\n请先运行：\ninstall.packages(c(",
    paste(sprintf("\"%s\"", missing_pkgs), collapse = ", "), "))")
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
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "02_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(data_dir)) {
  stop("未找到data文件夹。当前目录：", project_dir, "\n请打开zz.Rproj后再运行本脚本。")
}

log_file <- file.path(out_dir, "02_QC_3scRNA_basic_v01.log")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- write_log

standardize_time <- function(x) {
    z <- tolower(trimws(as.character(x)))
    z <- gsub("\\s+", "", z)
    z <- gsub("dpi", "d", z)
    z <- gsub("days?", "d", z)
    z <- gsub("hours?", "h", z)
    z <- gsub("hrs?", "h", z)
    out <- rep(NA_character_, length(z))
    out[z %in% c("uninjured", "uninj", "naive", "sham", "control", "0d", "00d", "0", "day0", "d0")] <- "Uninjured"
    out[z %in% c("4h", "4hr")] <- "4h"
    out[z %in% c("1d", "01d", "day1", "d1")] <- "1d"
    out[z %in% c("3d", "03d", "day3", "d3")] <- "3d"
    out[z %in% c("7d", "07d", "day7", "d7")] <- "7d"
    out[z %in% c("14d", "day14", "d14")] <- "14d"
    out[z %in% c("28d", "day28", "d28")] <- "28d"
    out[z %in% c("38d", "day38", "d38")] <- "38d"
    out[is.na(out)] <- as.character(x)[is.na(out)]
    out
}

time_to_numeric <- function(x) {
  x <- standardize_time(x)
  map <- c(Uninjured = 0, `4h` = 4/24, `1d` = 1, `3d` = 3, `7d` = 7, `14d` = 14, `28d` = 28, `38d` = 38)
  unname(map[x])
}

get_rna_counts <- function(obj) {
    DefaultAssay(obj) <- "RNA"
    counts <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
    if (is.null(counts)) {
        counts <- tryCatch(GetAssayData(obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
    }
    if (is.null(counts)) {
        obj2 <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) NULL)
        if (!is.null(obj2)) {
            counts <- tryCatch(GetAssayData(obj2, assay = "RNA", layer = "counts"), error = function(e) NULL)
        }
    }
    if (is.null(counts)) {
        stop("无法从Seurat对象提取RNA counts。")
    }
    counts
}

read_mtx_triplet <- function(matrix_file, barcode_file, feature_file) {
    mat <- Matrix::readMM(gzfile(matrix_file))
    mat <- as(mat, "dgCMatrix")
    bc <- fread(barcode_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
    ft <- fread(feature_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
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
    if (nrow(mat) != length(genes)) {
        stop("matrix行数与gene/features数量不一致。")
    }
    if (ncol(mat) != length(barcodes)) {
        stop("matrix列数与barcode数量不一致。")
    }
    rownames(mat) <- genes
    colnames(mat) <- barcodes
    mat
}

prefix_cells <- function(obj, prefix) {
  RenameCells(obj, new.names = paste0(prefix, "_", Cells(obj)))
}

add_qc_metrics <- function(obj) {
  DefaultAssay(obj) <- "RNA"

  mt_features <- grep("^mt-|^MT-", rownames(obj), value = TRUE)

  if (length(mt_features) == 0) {
    warning("未识别到线粒体基因，请核对gene symbol。")
    obj$percent.mt <- 0
  } else {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_features)
  }

  hb_features <- grep("^Hba|^Hbb|^HBA|^HBB", rownames(obj), value = TRUE)

  if (length(hb_features) > 0) {
    obj[["percent.HB"]] <- PercentageFeatureSet(obj, features = hb_features)
  } else {
    obj$percent.HB <- 0
  }

  obj
}

apply_basic_qc <- function(obj) {
  subset(obj, subset = nFeature_RNA >= MIN_FEATURES & percent.mt < MAX_MT)
}

qc_count_summary <- function(obj, dataset_name, stage) {
    md <- obj[[]]
    if (!"sample_id" %in% colnames(md)) {
        md$sample_id <- "NA"
    }
    if (!"time" %in% colnames(md)) {
        md$time <- "NA"
    }
    if (!"treatment" %in% colnames(md)) {
        md$treatment <- "Control"
    }
    key <- interaction(md$sample_id, md$time, md$treatment, drop = TRUE, lex.order = TRUE)
    idx <- split(seq_len(nrow(md)), key)
    do.call(rbind, lapply(idx, function(i) {
        data.frame(dataset = dataset_name, stage = stage, sample_id = as.character(md$sample_id[i[1]]), 
            time = as.character(md$time[i[1]]), treatment = as.character(md$treatment[i[1]]), n_cells = length(i), 
            median_nFeature_RNA = median(md$nFeature_RNA[i], na.rm = TRUE), median_nCount_RNA = median(md$nCount_RNA[i], 
                na.rm = TRUE), median_percent_mt = median(md$percent.mt[i], na.rm = TRUE), median_percent_HB = median(md$percent.HB[i], 
                na.rm = TRUE), stringsAsFactors = FALSE)
    }))
}

qc_quantiles <- function(obj, dataset_name, stage) {
    probs <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.5, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 1)
    md <- obj[[]]
    vars <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
    do.call(rbind, lapply(vars, function(v) {
        q <- quantile(md[[v]], probs = probs, na.rm = TRUE, names = FALSE)
        data.frame(dataset = dataset_name, stage = stage, metric = v, quantile = probs, value = as.numeric(q), 
            stringsAsFactors = FALSE)
    }))
}

save_qc_violin <- function(obj_before, obj_after, group_by, dataset_name) {
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
    make_panel <- function(obj, title_prefix) {
        plots <- lapply(features, function(f) {
            VlnPlot(obj, features = f, group.by = group_by, pt.size = 0) + NoLegend() + theme(axis.title.x = element_blank(), 
                axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) + ggtitle(paste0(title_prefix, 
                ": ", f))
        })
        wrap_plots(plots, nrow = 2)
    }
    p_before <- make_panel(obj_before, "Before QC")
    p_after <- make_panel(obj_after, "After QC")
    pdf(file.path(out_dir, paste0("02_QC_", dataset_name, "_violin_before_after_v01.pdf")), width = 15, 
        height = 12)
    print(p_before)
    print(p_after)
    dev.off()
}

save_qc_scatter <- function(obj, group_by, dataset_name) {
    p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = group_by)
    p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = group_by)
    p <- p1 + p2
    ggsave(filename = file.path(out_dir, paste0("02_QC_", dataset_name, "_scatter_before_v01.pdf")), 
        plot = p, width = 12, height = 5.5)
}

f177 <- file.path(data_dir, "17702045", "sci_scrna.rds")

if (!file.exists(f177)) {
  stop("未找到：", f177)
}

old177 <- readRDS(f177)

if (!inherits(old177, "Seurat")) {
  stop("sci_scrna.rds不是Seurat对象。")
}

counts177 <- get_rna_counts(old177)
meta177 <- old177[[]]

if (!"time" %in% colnames(meta177)) {
  stop("17702045 metadata中没有time字段。")
}

meta177$dataset <- "17702045"
meta177$time_raw <- as.character(meta177$time)
meta177$time <- standardize_time(meta177$time)
meta177$time_numeric <- time_to_numeric(meta177$time)
meta177$condition <- ifelse(meta177$time == "Uninjured", "Uninjured", "SCI")
meta177$treatment <- "Control"

sample_candidates <- c("sample_id", "sample", "orig.ident", "library", "replicate", "animal", "mouse")
sample_candidates <- sample_candidates[sample_candidates %in% colnames(meta177)]

if (length(sample_candidates) > 0) {
  sample_field177 <- sample_candidates[1]
  meta177$sample_id <- as.character(meta177[[sample_field177]])
  meta177$sample_id_source <- paste0("metadata:", sample_field177)
} else {
  meta177$sample_id <- paste0("F177_", meta177$time)
  meta177$sample_id_source <- "derived_from_timepoint"
}

meta177$cell_id_original <- rownames(meta177)

obj177 <- CreateSeuratObject(counts = counts177, meta.data = meta177, project = "F177", min.cells = 0,
  min.features = 0)

obj177 <- prefix_cells(obj177, "F177")
obj177 <- add_qc_metrics(obj177)

before177 <- ncol(obj177)

save_qc_violin(obj177, apply_basic_qc(obj177), group_by = "time", dataset_name = "17702045")

save_qc_scatter(obj177, group_by = "time", dataset_name = "17702045")

obj177_qc <- apply_basic_qc(obj177)
after177 <- ncol(obj177_qc)

d162 <- file.path(data_dir, "GSE162610")

f162_mtx <- file.path(d162, "GSE162610_sci_mat.mtx.gz")
f162_bc <- file.path(d162, "GSE162610_barcodes.tsv.gz")
f162_gene <- file.path(d162, "GSE162610_genes.tsv.gz")
f162_meta <- file.path(d162, "GSE162610_barcode_metadata.tsv.gz")

for (f in c(f162_mtx, f162_bc, f162_gene, f162_meta)) {
  if (!file.exists(f)) {
    stop("未找到：", f)
  }
}

mat162 <- read_mtx_triplet(f162_mtx, f162_bc, f162_gene)

meta162 <- fread(f162_meta, data.table = FALSE, check.names = FALSE, showProgress = FALSE)

if (!"V1" %in% colnames(meta162)) {
  stop("GSE162610 metadata中未找到barcode字段V1。")
}

idx162 <- match(colnames(mat162), as.character(meta162$V1))

if (anyNA(idx162)) {
  stop("GSE162610有 ", sum(is.na(idx162)), " 个barcode无法匹配metadata。")
}

meta162 <- meta162[idx162, , drop = FALSE]
rownames(meta162) <- colnames(mat162)

if (!"sample_id" %in% colnames(meta162)) {
  if ("orig.ident" %in% colnames(meta162)) {
    meta162$sample_id <- as.character(meta162$orig.ident)
  } else {
    stop("GSE162610缺少sample_id/orig.ident。")
  }
}

if (!"time" %in% colnames(meta162)) {
  stop("GSE162610 metadata中未找到time字段。")
}

meta162$dataset <- "GSE162610"
meta162$time_raw <- as.character(meta162$time)
meta162$time <- standardize_time(meta162$time)
meta162$time_numeric <- time_to_numeric(meta162$time)
meta162$condition <- ifelse(meta162$time == "Uninjured", "Uninjured", "SCI")
meta162$treatment <- "Control"
meta162$sample_id_source <- "metadata:sample_id"
meta162$cell_id_original <- colnames(mat162)

obj162 <- CreateSeuratObject(counts = mat162, meta.data = meta162, project = "GSE162610", min.cells = 0,
  min.features = 0)

obj162 <- prefix_cells(obj162, "G162")
obj162 <- add_qc_metrics(obj162)

before162 <- ncol(obj162)

save_qc_violin(obj162, apply_basic_qc(obj162), group_by = "sample_id", dataset_name = "GSE162610")

save_qc_scatter(obj162, group_by = "sample_id", dataset_name = "GSE162610")

obj162_qc <- apply_basic_qc(obj162)
after162 <- ncol(obj162_qc)

d196 <- file.path(data_dir, "GSE196928")
raw_tar <- file.path(d196, "GSE196928_RAW.tar")
raw_dir <- file.path(d196, "RAW_extracted")

if (!file.exists(raw_tar)) {
  stop("未找到：", raw_tar)
}

if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
    untar(raw_tar, exdir = raw_dir)
} else {
    log_msg("RAW_extracted已存在，跳过重复解压。")
}

sample196 <- data.frame(gsm = c("GSM5904825", "GSM5904826", "GSM5904827", "GSM5904828", "GSM5904829",
  "GSM5904830"), raw_id = c("B1", "B2", "B3", "B4", "B7", "B8"), sample_id = c("B1", "B2", "B3", "B4",
  "B7", "B8"), time = c("Uninjured", "Uninjured", "7d", "7d", "28d", "28d"), treatment = c("Control",
  "PLX5622", "Control", "PLX5622", "Control", "PLX5622"), stringsAsFactors = FALSE)

sample196$time_numeric <- time_to_numeric(sample196$time)

sample196$condition <- ifelse(sample196$time == "Uninjured", "Uninjured", "SCI")

obj196_list <- list()

for (i in seq_len(nrow(sample196))) {
    gsm <- sample196$gsm[i]
    bid <- sample196$raw_id[i]
    prefix <- paste0(gsm, "_", bid)
    matrix_file <- file.path(raw_dir, paste0(prefix, "_matrix.mtx.gz"))
    barcode_file <- file.path(raw_dir, paste0(prefix, "_barcodes.tsv.gz"))
    feature_file <- file.path(raw_dir, paste0(prefix, "_features.tsv.gz"))
    for (f in c(matrix_file, barcode_file, feature_file)) {
        if (!file.exists(f)) {
            stop("解压后未找到：", f)
        }
    }
    mat <- read_mtx_triplet(matrix_file, barcode_file, feature_file)
    md <- data.frame(dataset = "GSE196928", sample_id = bid, sample_id_source = "RAW_filename", time_raw = sample196$time[i], 
        time = sample196$time[i], time_numeric = sample196$time_numeric[i], condition = sample196$condition[i], 
        treatment = sample196$treatment[i], cell_id_original = colnames(mat), row.names = colnames(mat), 
        stringsAsFactors = FALSE)
    obj <- CreateSeuratObject(counts = mat, meta.data = md, project = paste0("G196_", bid), min.cells = 0, 
        min.features = 0)
    obj <- prefix_cells(obj, paste0("G196_", bid))
    obj <- add_qc_metrics(obj)
    obj196_list[[bid]] <- obj
    rm(mat, md, obj)
    gc(verbose = FALSE)
}

obj196_all <- merge(obj196_list[["B1"]], y = list(obj196_list[["B2"]], obj196_list[["B3"]], obj196_list[["B4"]],
  obj196_list[["B7"]], obj196_list[["B8"]]), merge.data = FALSE)

before196 <- ncol(obj196_all)

save_qc_violin(obj196_all, apply_basic_qc(obj196_all), group_by = "sample_id", dataset_name = "GSE196928")

save_qc_scatter(obj196_all, group_by = "sample_id", dataset_name = "GSE196928")

obj196_all_qc <- apply_basic_qc(obj196_all)

after196 <- ncol(obj196_all_qc)

obj196_control_qc <- subset(obj196_all_qc, subset = treatment == "Control")

obj196_plx_qc <- subset(obj196_all_qc, subset = treatment == "PLX5622")

summary_counts <- rbind(qc_count_summary(obj177, "17702045", "before_QC"), qc_count_summary(obj177_qc,
  "17702045", "after_QC"), qc_count_summary(obj162, "GSE162610", "before_QC"), qc_count_summary(obj162_qc,
  "GSE162610", "after_QC"), qc_count_summary(obj196_all, "GSE196928", "before_QC"), qc_count_summary(obj196_all_qc,
  "GSE196928", "after_QC"))

write.csv(summary_counts, file.path(out_dir, "02_QC_3scRNA_cell_counts_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

quantile_table <- rbind(qc_quantiles(obj177, "17702045", "before_QC"), qc_quantiles(obj177_qc, "17702045",
  "after_QC"), qc_quantiles(obj162, "GSE162610", "before_QC"), qc_quantiles(obj162_qc, "GSE162610",
  "after_QC"), qc_quantiles(obj196_all, "GSE196928", "before_QC"), qc_quantiles(obj196_all_qc, "GSE196928",
  "after_QC"))

write.csv(quantile_table, file.path(out_dir, "02_QC_3scRNA_quantiles_v01.csv"), row.names = FALSE, fileEncoding = "UTF-8")

dataset_summary <- data.frame(dataset = c("17702045", "GSE162610", "GSE196928"), before_QC = c(before177,
  before162, before196), after_QC = c(after177, after162, after196), removed = c(before177 - after177,
  before162 - after162, before196 - after196), retained_pct = round(c(after177/before177, after162/before162,
  after196/before196) * 100, 2))

write.csv(dataset_summary, file.path(out_dir, "02_QC_3scRNA_dataset_summary_v01.csv"), row.names = FALSE,
  fileEncoding = "UTF-8")

print(dataset_summary)

saveRDS(obj177_qc, file.path(out_dir, "02_QC_17702045_v01.rds"), compress = FALSE)

saveRDS(obj162_qc, file.path(out_dir, "02_QC_GSE162610_v01.rds"), compress = FALSE)

saveRDS(obj196_control_qc, file.path(out_dir, "02_QC_GSE196928_Control_v01.rds"), compress = FALSE)

saveRDS(obj196_plx_qc, file.path(out_dir, "02_QC_GSE196928_PLX5622_v01.rds"), compress = FALSE)

capture.output(sessionInfo(), file = file.path(out_dir, "02_QC_sessionInfo_v01.txt"))

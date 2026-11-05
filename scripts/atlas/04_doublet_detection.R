source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

SEED <- 1234

D177_FIXED_DBR <- 0.07

set.seed(SEED)

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
  "scDblFinder", "ggplot2", "data.table")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop("缺少R包：", paste(missing_pkgs, collapse = ", "), "\n请先安装后再运行。")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scDblFinder)
  library(ggplot2)
  library(data.table)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
in_dir <- file.path(project_dir, "results", "03_DecontX")
out_dir <- file.path(project_dir, "results", "04_doublet")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "04_doublet_scDblFinder_3scRNA_v03.log")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- write_log

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
        stop("无法从Seurat对象读取RNA counts。")
    as(x, "dgCMatrix")
}

run_scDblFinder_one <- function(obj, dataset_name, unit_name, dbr_value = NULL) {
  counts <- get_rna_counts(obj)

  dbr_text <- if (is.null(dbr_value)) {
    "AUTO"
  } else {
    paste0(round(dbr_value * 100, 2), "% fixed")
  }

  sce <- SingleCellExperiment(assays = list(counts = counts))

  set.seed(SEED)

  if (is.null(dbr_value)) {
    sce <- scDblFinder(sce, dbr = NULL)
  } else {
    sce <- scDblFinder(sce, dbr = dbr_value)
  }

  cd <- as.data.frame(colData(sce))

  if (!all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(cd))) {
    stop("scDblFinder未返回预期的score/class字段。")
  }

  score <- as.numeric(cd$scDblFinder.score)
  cls <- as.character(cd$scDblFinder.class)

  if (length(score) != ncol(obj) || length(cls) != ncol(obj)) {
    stop("scDblFinder结果长度与细胞数不一致。")
  }

  obj$scDblFinder_score <- score
  obj$scDblFinder_class <- cls
  obj$doublet_unit <- unit_name
  obj$doublet_dbr_mode <- if (is.null(dbr_value))
    "automatic" else "fixed"
  obj$doublet_dbr_input <- if (is.null(dbr_value))
    NA_real_ else dbr_value

  is_doublet <- cls == "doublet"

  n_before <- ncol(obj)
  n_doublet <- sum(is_doublet, na.rm = TRUE)
  n_after <- n_before - n_doublet

  observed_rate <- n_doublet/n_before

  summary_row <- data.frame(dataset = dataset_name, unit = unit_name, n_cells_before = n_before, dbr_mode = if (is.null(dbr_value))
    "automatic" else "fixed", dbr_input = if (is.null(dbr_value))
    NA_real_ else dbr_value, n_doublets = n_doublet, doublet_pct = observed_rate * 100, n_singlets = n_after,
    score_median = median(score, na.rm = TRUE), score_p95 = unname(quantile(score, 0.95, na.rm = TRUE)),
    score_max = max(score, na.rm = TRUE), stringsAsFactors = FALSE)

  score_df <- data.frame(dataset = dataset_name, unit = unit_name, score = score, class = cls, dbr_mode = summary_row$dbr_mode,
    stringsAsFactors = FALSE)

  singlet_cells <- Cells(obj)[!is_doublet]
  obj_singlet <- subset(obj, cells = singlet_cells)

  rm(counts, sce, cd)
  gc(verbose = FALSE)

  list(singlet = obj_singlet, summary = summary_row, score_df = score_df)
}

run_by_unit <- function(obj, dataset_name, split_field = NULL, dbr_value = NULL, whole_unit_name = "whole_dataset") {
    if (is.null(split_field)) {
        res <- run_scDblFinder_one(obj, dataset_name = dataset_name, unit_name = whole_unit_name, dbr_value = dbr_value)
        return(list(singlet = res$singlet, summary = res$summary, score_df = res$score_df))
    }
    md <- obj[[]]
    if (!split_field %in% colnames(md)) {
        stop(dataset_name, " 缺少字段：", split_field)
    }
    units <- unique(as.character(md[[split_field]]))
    units <- units[!is.na(units) & nzchar(units)]
    singlet_list <- list()
    summary_list <- list()
    score_list <- list()
    for (u in units) {
        cells <- rownames(md)[as.character(md[[split_field]]) == u]
        one <- subset(obj, cells = cells)
        res <- run_scDblFinder_one(one, dataset_name = dataset_name, unit_name = u, dbr_value = dbr_value)
        singlet_list[[u]] <- res$singlet
        summary_list[[u]] <- res$summary
        score_list[[u]] <- res$score_df
        rm(one, res)
        gc(verbose = FALSE)
    }
    if (length(singlet_list) == 1) {
        merged <- singlet_list[[1]]
    }
    else {
        merged <- merge(singlet_list[[1]], y = singlet_list[-1], merge.data = FALSE)
    }
    list(singlet = merged, summary = do.call(rbind, summary_list), score_df = do.call(rbind, score_list))
}

paths <- c(f177 = file.path(in_dir, "03_DecontX_17702045_v01.rds"), f162 = file.path(in_dir, "03_DecontX_GSE162610_v01.rds"),
  f196c = file.path(in_dir, "03_DecontX_GSE196928_Control_v01.rds"), f196p = file.path(in_dir, "03_DecontX_GSE196928_PLX5622_v01.rds"))

if (any(!file.exists(paths))) {
  stop("缺少DecontX RDS：\n", paste(paths[!file.exists(paths)], collapse = "\n"))
}

obj177 <- readRDS(paths["f177"])
obj162 <- readRDS(paths["f162"])
obj196c <- readRDS(paths["f196c"])
obj196p <- readRDS(paths["f196p"])

res177 <- run_by_unit(obj177, dataset_name = "17702045", split_field = NULL, dbr_value = D177_FIXED_DBR,
  whole_unit_name = "whole_dataset")

rm(obj177)
gc(verbose = FALSE)

res162 <- run_by_unit(obj162, dataset_name = "GSE162610", split_field = "sample_id", dbr_value = NULL)

rm(obj162)
gc(verbose = FALSE)

res196c <- run_by_unit(obj196c, dataset_name = "GSE196928_Control", split_field = "sample_id", dbr_value = NULL)

rm(obj196c)
gc(verbose = FALSE)

res196p <- run_by_unit(obj196p, dataset_name = "GSE196928_PLX5622", split_field = "sample_id", dbr_value = NULL)

rm(obj196p)
gc(verbose = FALSE)

summary_df <- rbind(res177$summary, res162$summary, res196c$summary, res196p$summary)

write.csv(summary_df, file.path(out_dir, "04_doublet_3scRNA_summary_v03.csv"), row.names = FALSE, fileEncoding = "UTF-8")

score_df <- rbind(res177$score_df, res162$score_df, res196c$score_df, res196p$score_df)

score_df$label <- paste(score_df$dataset, score_df$unit, sep = " | ")

pdf(file.path(out_dir, "04_doublet_3scRNA_diagnostics_v03.pdf"), width = 12, height = 8.5)

p1 <- ggplot(summary_df, aes(x = interaction(dataset, unit, sep = " | "), y = doublet_pct, fill = dbr_mode)) +
  geom_col() + labs(x = NULL, y = "scDblFinder doublets (%)", title = "Doublet proportion by sample/library",
  fill = "dbr mode") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 60,
  hjust = 1, vjust = 1))

print(p1)

p2 <- ggplot(score_df, aes(x = label, y = score)) + geom_boxplot(outlier.size = 0.2, linewidth = 0.35) +
  labs(x = NULL, y = "scDblFinder score", title = "scDblFinder score distribution") + theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

print(p2)

p3 <- ggplot(score_df, aes(x = class, y = score)) + geom_boxplot(outlier.size = 0.2) + facet_wrap(~dataset,
  scales = "free_y") + labs(x = NULL, y = "scDblFinder score", title = "scDblFinder classification") +
  theme_classic(base_size = 11)

print(p3)

dev.off()

saveRDS(res177$singlet, file.path(out_dir, "04_doublet_17702045_singlet_v03.rds"), compress = FALSE)

saveRDS(res162$singlet, file.path(out_dir, "04_doublet_GSE162610_singlet_v03.rds"), compress = FALSE)

saveRDS(res196c$singlet, file.path(out_dir, "04_doublet_GSE196928_Control_singlet_v03.rds"), compress = FALSE)

saveRDS(res196p$singlet, file.path(out_dir, "04_doublet_GSE196928_PLX5622_singlet_v03.rds"), compress = FALSE)

main_idx <- summary_df$dataset %in% c("17702045", "GSE162610", "GSE196928_Control")

main_before <- sum(summary_df$n_cells_before[main_idx])
main_doublets <- sum(summary_df$n_doublets[main_idx])
main_after <- sum(summary_df$n_singlets[main_idx])

log_msg("2. 17702045因缺少可靠library ID，暂固定dbr=7%。")
capture.output(sessionInfo(), file = log_file, append = TRUE)

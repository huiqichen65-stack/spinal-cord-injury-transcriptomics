options(stringsAsFactors = FALSE)
options(warn = 1)

required_pkgs <- c("Matrix", "data.table", "SeuratObject")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop("缺少R包：", paste(missing_pkgs, collapse = ", "), "\n请先安装： install.packages(c(",
    paste(sprintf("\"%s\"", missing_pkgs), collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SeuratObject)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "01_QC")

if (!dir.exists(data_dir)) {
  stop("未找到 data 文件夹。当前工作目录：", project_dir, "\n请先打开 zz.Rproj，使 getwd() 指向项目根目录。")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "01_QC_3scRNA_check_v01.log")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

safe_nunique <- function(x) length(unique(x[!is.na(x)]))
collapse_unique <- function(x, max_n = 20) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0)
    return("")
  if (length(x) > max_n)
    return(paste0(paste(head(x, max_n), collapse = " | "), " | ..."))
  paste(x, collapse = " | ")
}

find_candidate_fields <- function(meta) {
  if (is.null(meta) || ncol(meta) == 0)
    return(list(sample = character(), time = character(), group = character()))
  nms <- colnames(meta)
  low <- tolower(nms)
  list(sample = nms[grepl("sample|orig\\.ident|library|replicate|animal|mouse", low)], time = nms[grepl("time|dpi|day|stage",
    low)], group = nms[grepl("group|condition|treat|vehicle|plx|injury|status", low)])
}

print_candidate_fields <- function(meta, candidates) {
  for (nm in names(candidates)) {
    fields <- candidates[[nm]]
    if (length(fields) == 0) {
    } else {
      for (f in fields) log_msg("    - ", f, ": ", collapse_unique(meta[[f]]))
    }
  }
}

meta_profile <- function(meta, dataset, outfile) {
  if (is.null(meta) || ncol(meta) == 0)
    return(invisible(NULL))
  ans <- data.frame(dataset = dataset, field = colnames(meta), class = vapply(meta, function(x) paste(class(x),
    collapse = "/"), character(1)), n_unique = vapply(meta, safe_nunique, integer(1)), example_values = vapply(meta,
    collapse_unique, character(1)), stringsAsFactors = FALSE)
  write.csv(ans, outfile, row.names = FALSE, fileEncoding = "UTF-8")
}

best_overlap_with_meta <- function(ids, meta) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) == 0 || is.null(meta) || nrow(meta) == 0) 
        return(list(field = NA_character_, matched = 0L, ratio = 0))
    candidates <- list()
    rn <- rownames(meta)
    if (!is.null(rn) && !all(rn == as.character(seq_len(nrow(meta))))) 
        candidates[["<rownames>"]] <- rn
    for (nm in colnames(meta)) {
        x <- meta[[nm]]
        if (is.character(x) || is.factor(x)) 
            candidates[[nm]] <- as.character(x)
    }
    if (length(candidates) == 0) 
        return(list(field = NA_character_, matched = 0L, ratio = 0))
    matched <- vapply(candidates, function(x) sum(ids %in% x), integer(1))
    best <- which.max(matched)
    list(field = names(candidates)[best], matched = unname(matched[best]), ratio = unname(matched[best])/length(ids))
}

file_ok <- function(path) file.exists(path) && file.info(path)$size > 0

summary_rows <- list()
add_summary <- function(dataset, status, n_genes = NA, n_cells = NA, n_samples = NA, time_fields = "",
  group_fields = "", notes = "") {
  summary_rows[[length(summary_rows) + 1]] <<- data.frame(dataset = dataset, status = status, n_genes = n_genes,
    n_cells = n_cells, n_samples = n_samples, time_fields = time_fields, group_fields = group_fields,
    notes = notes, stringsAsFactors = FALSE)
}

f177 <- file.path(data_dir, "17702045", "sci_scrna.rds")

if (!file_ok(f177)) {
  log_msg("[FAIL] 未找到或文件为空：", f177)
  add_summary("17702045", "FAIL", notes = "sci_scrna.rds missing/empty")
} else {
  tryCatch({
    obj177 <- readRDS(f177)
    if (inherits(obj177, "Seurat")) {
      meta177 <- obj177[[]]
      id_identical <- identical(colnames(obj177), rownames(meta177))
      candidates <- find_candidate_fields(meta177)
      print_candidate_fields(meta177, candidates)
      meta_profile(meta177, "17702045", file.path(out_dir, "01_QC_17702045_metadata_fields_v01.csv"))
      sf <- if ("orig.ident" %in% candidates$sample)
        "orig.ident" else if (length(candidates$sample) > 0)
        candidates$sample[1] else NA_character_
      ns <- if (!is.na(sf))
        safe_nunique(meta177[[sf]]) else NA_integer_
      add_summary("17702045", "PASS", nrow(obj177), ncol(obj177), ns, paste(candidates$time, collapse = ";"),
        paste(candidates$group, collapse = ";"), "Seurat object readable")
    } else {
      log_msg("[WARN] RDS可读，但不是Seurat对象")
      add_summary("17702045", "WARN", notes = paste0("RDS class: ", paste(class(obj177), collapse = "/")))
    }
    rm(obj177)
    gc(verbose = FALSE)
  }, error = function(e) {
    log_msg("[FAIL] 17702045读取失败：", conditionMessage(e))
    add_summary("17702045", "FAIL", notes = conditionMessage(e))
  })
}

d162 <- file.path(data_dir, "GSE162610")
f162_mtx <- file.path(d162, "GSE162610_sci_mat.mtx.gz")
f162_bc <- file.path(d162, "GSE162610_barcodes.tsv.gz")
f162_gene <- file.path(d162, "GSE162610_genes.tsv.gz")
f162_bcmet <- file.path(d162, "GSE162610_barcode_metadata.tsv.gz")
f162_gmet <- file.path(d162, "GSE162610_gene_metadata.tsv.gz")
req162 <- c(matrix = f162_mtx, barcodes = f162_bc, genes = f162_gene, barcode_metadata = f162_bcmet,
  gene_metadata = f162_gmet)
missing162 <- names(req162)[!vapply(req162, file_ok, logical(1))]

if (length(missing162) > 0) {
  log_msg("[FAIL] 缺少/空文件：", paste(missing162, collapse = ", "))
  add_summary("GSE162610", "FAIL", notes = paste0("missing/empty: ", paste(missing162, collapse = ";")))
} else {
  tryCatch({
    mat162 <- Matrix::readMM(gzfile(f162_mtx))
    bc162 <- data.table::fread(f162_bc, header = FALSE, data.table = FALSE)
    gn162 <- data.table::fread(f162_gene, header = FALSE, data.table = FALSE)
    bm162 <- data.table::fread(f162_bcmet, data.table = FALSE)
    gm162 <- data.table::fread(f162_gmet, data.table = FALSE)

    gene_dim_ok <- nrow(mat162) == nrow(gn162)
    bc_dim_ok <- ncol(mat162) == nrow(bc162)

    gene_ids <- as.character(gn162[[1]])
    barcodes <- as.character(bc162[[1]])
    log_msg(if (sum(duplicated(gene_ids)) == 0)
      "[PASS] gene ID无重复" else "[WARN] gene ID存在重复")
    log_msg(if (sum(duplicated(barcodes)) == 0)
      "[PASS] barcode无重复" else "[WARN] barcode存在重复")

    bc_match <- best_overlap_with_meta(barcodes, bm162)
    gene_match <- best_overlap_with_meta(gene_ids, gm162)

    candidates <- find_candidate_fields(bm162)
    print_candidate_fields(bm162, candidates)
    meta_profile(bm162, "GSE162610", file.path(out_dir, "01_QC_GSE162610_barcode_metadata_fields_v01.csv"))
    meta_profile(gm162, "GSE162610_gene_metadata", file.path(out_dir, "01_QC_GSE162610_gene_metadata_fields_v01.csv"))

    sf <- if (length(candidates$sample) > 0)
      candidates$sample[1] else NA_character_
    ns <- if (!is.na(sf))
      safe_nunique(bm162[[sf]]) else NA_integer_
    status162 <- if (gene_dim_ok && bc_dim_ok)
      "PASS" else "FAIL"
    add_summary("GSE162610", status162, nrow(mat162), ncol(mat162), ns, paste(candidates$time, collapse = ";"),
      paste(candidates$group, collapse = ";"), paste0("barcode-meta=", bc_match$field, ":", round(bc_match$ratio *
        100, 2), "%; gene-meta=", gene_match$field, ":", round(gene_match$ratio * 100, 2), "%"))
    rm(mat162, bc162, gn162, bm162, gm162)
    gc(verbose = FALSE)
  }, error = function(e) {
    log_msg("[FAIL] GSE162610读取/检查失败：", conditionMessage(e))
    add_summary("GSE162610", "FAIL", notes = conditionMessage(e))
  })
}

d196 <- file.path(data_dir, "GSE196928")
f196_count <- file.path(d196, "GSE196928_counts.csv.gz")
f196_meta <- file.path(d196, "GSE196928_metadata.csv.gz")
f196_raw <- file.path(d196, "GSE196928_RAW.tar")
req196 <- c(counts = f196_count, metadata = f196_meta, RAWtar = f196_raw)
missing196 <- names(req196)[!vapply(req196, file_ok, logical(1))]

if (length(missing196) > 0) {
  log_msg("[FAIL] 缺少/空文件：", paste(missing196, collapse = ", "))
  add_summary("GSE196928", "FAIL", notes = paste0("missing/empty: ", paste(missing196, collapse = ";")))
} else {
  tryCatch({
    meta196 <- data.table::fread(f196_meta, data.table = FALSE)

    cnt196 <- data.table::fread(f196_count, data.table = FALSE, check.names = FALSE, showProgress = TRUE)

    first_col_vals <- as.character(cnt196[[1]])
    count_col_ids <- colnames(cnt196)[-1]
    col_match <- best_overlap_with_meta(count_col_ids, meta196)
    row_match <- best_overlap_with_meta(first_col_vals, meta196)

    if (col_match$ratio >= 0.5 && col_match$ratio >= row_match$ratio) {
      orientation <- "genes x cells"
      n_genes196 <- nrow(cnt196)
      n_cells196 <- ncol(cnt196) - 1
      best_match <- col_match
    } else if (row_match$ratio >= 0.5 && row_match$ratio > col_match$ratio) {
      orientation <- "cells x genes"
      n_cells196 <- nrow(cnt196)
      n_genes196 <- ncol(cnt196) - 1
      best_match <- row_match
    } else {
      orientation <- "uncertain"
      n_genes196 <- NA_integer_
      n_cells196 <- NA_integer_
      best_match <- if (col_match$ratio >= row_match$ratio)
        col_match else row_match
      log_msg("[WARN] 暂无法可靠判断counts方向")
    }

    tar_members <- tryCatch(utils::untar(f196_raw, list = TRUE), error = function(e) e)
    if (inherits(tar_members, "error")) {
      log_msg("[WARN] RAW.tar目录读取失败：", conditionMessage(tar_members))
      tar_note <- "RAW.tar listing failed"
    } else {
      tar_note <- paste0("RAW.tar members=", length(tar_members))
    }

    candidates <- find_candidate_fields(meta196)
    print_candidate_fields(meta196, candidates)
    meta_profile(meta196, "GSE196928", file.path(out_dir, "01_QC_GSE196928_metadata_fields_v01.csv"))

    sf <- if (length(candidates$sample) > 0)
      candidates$sample[1] else NA_character_
    ns <- if (!is.na(sf))
      safe_nunique(meta196[[sf]]) else NA_integer_
    status196 <- if (orientation != "uncertain" && best_match$ratio >= 0.5)
      "PASS" else "WARN"
    add_summary("GSE196928", status196, n_genes196, n_cells196, ns, paste(candidates$time, collapse = ";"),
      paste(candidates$group, collapse = ";"), paste0("orientation=", orientation, "; best metadata match=",
        best_match$field, ":", round(best_match$ratio * 100, 2), "%; ", tar_note))
    rm(meta196, cnt196)
    gc(verbose = FALSE)
  }, error = function(e) {
    log_msg("[FAIL] GSE196928读取/检查失败：", conditionMessage(e))
    add_summary("GSE196928", "FAIL", notes = conditionMessage(e))
  })
}

summary_df <- do.call(rbind, summary_rows)
summary_file <- file.path(out_dir, "01_QC_3scRNA_summary_v01.csv")
write.csv(summary_df, summary_file, row.names = FALSE, fileEncoding = "UTF-8")
print(summary_df, row.names = FALSE)

if (all(summary_df$status == "PASS")) {
} else if (any(summary_df$status == "FAIL")) {
  log_msg("[总体结论] 至少1套数据存在FAIL，请先修复后再进行整合。")
} else {
  log_msg("[总体结论] 数据可读取，但存在WARN，请先核对metadata字段/ID匹配。")
}

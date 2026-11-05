source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

ROOT <- PROJECT_ROOT
IN_DIR <- file.path(ROOT, "data", "GSE281205")
OUT_DIR <- file.path(ROOT, "results", "116_GSE281205_cell_QC")
FIG_DIR <- file.path(OUT_DIR, "figures")
RDS_DIR <- file.path(OUT_DIR, "RDS")

MIN_FEATURES <- 300L
MAX_MT_PERCENT <- 10
SAVE_RDS <- TRUE
SAVE_BARCODE_QC <- TRUE
DRAW_FIGURES <- TRUE

SAMPLES <- data.frame(gsm = paste0("GSM", 8614662:8614667), intended_group = rep(c("Young-14dpi", "Young-28dpi"),
  each = 3L), dpi = rep(c(14L, 28L), each = 3L), replicate = rep(1:3, times = 2L), stringsAsFactors = FALSE)

assert <- function(condition, message) {
  if (!isTRUE(condition))
    stop(message, call. = FALSE)
}
assert(dir.exists(IN_DIR), paste("Missing input directory:", IN_DIR))
assert(requireNamespace("Matrix", quietly = TRUE), "Missing package 'Matrix'. Please run install.packages('Matrix') first.")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "00_116_run_log_v02.log")
cat("", file = LOG_FILE)
log_msg <- function(...) {
  message_text <- paste0(..., collapse = "")
  cat(message_text, "\n", sep = "")
  cat(message_text, "\n", file = LOG_FILE, append = TRUE, sep = "")
}
save_csv <- function(df, filename) {
  assert(is.data.frame(df), paste("Expected data.frame:", filename))
  target <- file.path(OUT_DIR, filename)
  utils::write.csv(df, target, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(target)
}
read_gz_table <- function(path, ncols) {
  assert(file.exists(path), paste("Missing file:", path))
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  ans <- utils::read.delim(con, header = FALSE, sep = "\t", quote = "", comment.char = "", colClasses = "character",
    stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE)
  assert(ncol(ans) >= ncols, paste("Unexpected TSV structure:", basename(path)))
  ans
}
resolve_triplet <- function(gsm, suffix) {
  all_paths <- list.files(IN_DIR, full.names = TRUE, recursive = FALSE)
  nm <- basename(all_paths)
  matched <- all_paths[startsWith(nm, paste0(gsm, "_")) & endsWith(nm, paste0("_", suffix))]
  assert(length(matched) == 1L, paste("Expected exactly one", suffix, "for", gsm, "but found", length(matched)))
  matched[[1L]]
}
read_geo_soft <- function(path) {
  assert(file.exists(path), paste("Missing GEO SOFT file:", path))
  con <- gzfile(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE[[:space:]]*=[[:space:]]*GSM", lines)
  assert(length(starts) >= 6L, "GEO SOFT has fewer than six SAMPLE blocks")
  ends <- c(starts[-1L] - 1L, length(lines))
  result <- lapply(seq_along(starts), function(i) {
    b <- lines[starts[[i]]:ends[[i]]]
    field <- function(prefix) {
      z <- b[startsWith(b, prefix)]
      if (!length(z))
        return(NA_character_)
      trimws(substring(z[[1L]], nchar(prefix) + 1L))
    }
    ch <- b[startsWith(b, "!Sample_characteristics_ch1 = ")]
    ch <- substring(ch, nchar("!Sample_characteristics_ch1 = ") + 1L)
    trt <- ch[grepl("^treatment[[:space:]]*:", ch, ignore.case = TRUE)]
    data.frame(gsm = sub("^\\^SAMPLE[[:space:]]*=[[:space:]]*", "", b[[1L]]), geo_title = field("!Sample_title = "),
      geo_treatment = if (length(trt)) {
        trimws(sub("^[^:]+:[[:space:]]*", "", trt[[1L]]))
      } else NA_character_, stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}
median_safe <- function(x) {
  if (!length(x) || all(is.na(x)))
    return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}
quantile_safe <- function(x, p = 0.95) {
  if (!length(x) || all(is.na(x)))
    return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

plot_sample_qc <- function(md, gsm) {
    f <- file.path(FIG_DIR, paste0("Fig_116_", gsm, "_cell_QC_v02.png"))
    grDevices::png(f, width = 1800, height = 1300, res = 160)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3.2, 1), oma = c(0, 0, 2, 0))
    graphics::hist(md$nFeature_RNA, breaks = 65, col = "#8DA9C4", border = "white", main = "Detected genes: before QC", 
        xlab = "nFeature_RNA")
    graphics::abline(v = MIN_FEATURES, col = "#C3423F", lwd = 2.5)
    graphics::hist(md$percent.mt, breaks = 65, col = "#8DA9C4", border = "white", main = "Mitochondrial percentage: before QC", 
        xlab = "percent.mt (%)")
    graphics::abline(v = MAX_MT_PERCENT, col = "#C3423F", lwd = 2.5)
    graphics::plot(md$nCount_RNA + 1, md$nFeature_RNA + 1, log = "xy", pch = 16, cex = 0.23, col = grDevices::adjustcolor(ifelse(md$QC_pass, 
        "#13866B", "#C3423F"), alpha.f = 0.4), xlab = "UMI counts + 1 (log scale)", ylab = "Detected genes + 1 (log scale)", 
        main = "Retained (green) / excluded (red)")
    graphics::abline(h = MIN_FEATURES + 1, lty = 2, col = "#333333")
    kept <- md$nFeature_RNA[md$QC_pass]
    graphics::hist(kept, breaks = 65, col = "#13866B", border = "white", main = "Detected genes: after basic QC", 
        xlab = "nFeature_RNA")
    graphics::mtext(paste0(gsm, " | retained ", sum(md$QC_pass), "/", nrow(md), " barcodes"), side = 3, 
        outer = TRUE, line = 0.5, font = 2, cex = 1.1)
    invisible(f)
}

plot_overview <- function(sumtab) {
  f <- file.path(FIG_DIR, "Fig_116_all_samples_cell_QC_overview_v02.png")
  grDevices::png(f, width = 1800, height = 930, res = 155)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(6, 4.2, 3.5, 1))
  lab <- paste0(sub("GSM", "", sumtab$gsm), "\n", sumtab$dpi, " dpi")
  x <- rbind(sumtab$raw_barcodes, sumtab$qc_retained)
  bp <- graphics::barplot(x, beside = TRUE, col = c("#9BB0C5", "#13866B"), names.arg = lab, las = 1,
    cex.names = 0.9, main = "Counts before / after cell QC", ylab = "Barcodes")
  graphics::legend("topleft", c("Input", "Retained"), fill = c("#9BB0C5", "#13866B"), bty = "n")
  graphics::barplot(100 * sumtab$qc_retained/sumtab$raw_barcodes, names.arg = lab, col = "#13866B",
    ylim = c(0, 100), main = "Basic QC retention per sample", ylab = "Retained (%)", las = 1, cex.names = 0.9)
  invisible(f)
}

main <- function() {
    geo <- read_geo_soft(file.path(IN_DIR, "GSE281205_family.soft.gz"))
    ii <- match(SAMPLES$gsm, geo$gsm)
    assert(!anyNA(ii), "Six planned GSMs not found in GEO SOFT metadata")
    sample_meta <- cbind(SAMPLES, geo[ii, c("geo_title", "geo_treatment"), drop = FALSE])
    rownames(sample_meta) <- NULL
    sample_meta$metadata_conflict <- grepl("Young", sample_meta$geo_title, ignore.case = TRUE) & grepl("Old", 
        sample_meta$geo_treatment, ignore.case = TRUE)
    sample_meta$age_assignment_status <- ifelse(sample_meta$metadata_conflict, "PROVISIONAL_YOUNG_28DPI__GEO_CONFLICT", 
        "GEO_TITLE_TREATMENT_CONSISTENT")
    save_csv(sample_meta, "01_116_sample_metadata_audit_v02.csv")
    if (any(sample_meta$metadata_conflict)) {
        log_msg("WARNING: Young/Old GEO metadata conflict in ", paste(sample_meta$gsm[sample_meta$metadata_conflict], 
            collapse = ", "), ". Preserve for manuscript provenance.")
    }
    all_sample_stats <- vector("list", nrow(SAMPLES))
    all_barcode_data <- vector("list", nrow(SAMPLES))
    feature_reference <- NULL
    feature_concordance <- vector("list", nrow(SAMPLES))
    rds_paths <- character()
    for (i in seq_len(nrow(SAMPLES))) {
        s <- sample_meta[i, , drop = FALSE]
        gsm <- s$gsm[[1L]]
        file_b <- resolve_triplet(gsm, "barcodes.tsv.gz")
        file_f <- resolve_triplet(gsm, "features.tsv.gz")
        file_m <- resolve_triplet(gsm, "matrix.mtx.gz")
        bar <- read_gz_table(file_b, 1L)
        feat <- read_gz_table(file_f, 2L)
        assert(!anyDuplicated(bar[[1L]]), paste("Duplicate barcodes:", gsm))
        assert(!anyDuplicated(feat[[1L]]), paste("Duplicate feature IDs:", gsm))
        feature_types <- if (ncol(feat) >= 3L) 
            trimws(feat[[3L]])
        else {
            rep("Gene Expression", nrow(feat))
        }
        gex <- which(feature_types == "Gene Expression")
        assert(length(gex) >= 1000L, paste("Too few Gene Expression features:", gsm))
        con <- gzfile(file_m, open = "rt")
        mm <- tryCatch(Matrix::readMM(con), finally = close(con))
        assert(nrow(mm) == nrow(feat) && ncol(mm) == nrow(bar), paste("Matrix/TSV dimension mismatch:", 
            gsm))
        x <- as(mm[gex, , drop = FALSE], "dgCMatrix")
        rm(mm)
        gene_ids <- trimws(as.character(feat[[1L]][gex]))
        gene_symbols <- trimws(as.character(feat[[2L]][gex]))
        assert(all(nzchar(gene_ids)) && !anyDuplicated(gene_ids), paste("Invalid gene IDs:", gsm))
        assert(all(is.finite(x@x)) && all(x@x >= 0) && all(abs(x@x - round(x@x)) < 1e-08), paste("Non-finite, negative, or non-integer matrix counts:", 
            gsm))
        rownames(x) <- gene_ids
        raw_barcode <- as.character(bar[[1L]])
        colnames(x) <- paste0(gsm, "_", raw_barcode)
        gene_map <- data.frame(gene_id = gene_ids, gene_symbol = gene_symbols, feature_type = "Gene Expression")
        if (is.null(feature_reference)) {
            feature_reference <- gene_map
            save_csv(gene_map, "06_116_gene_feature_mapping_reference_v02.csv")
            is_same <- TRUE
        }
        else {
            is_same <- identical(gene_map, feature_reference)
            if (!is_same) {
                save_csv(gene_map, paste0("06_116_", gsm, "_feature_mapping_exception_v02.csv"))
                log_msg("WARNING: Gene feature mapping differs from first sample: ", gsm)
            }
        }
        feature_concordance[[i]] <- data.frame(gsm = gsm, feature_map_identical_to_first_sample = is_same, 
            n_gene_expression_features = nrow(gene_map), n_duplicated_symbols = sum(duplicated(gene_symbols[nzchar(gene_symbols)])), 
            stringsAsFactors = FALSE)
        mt_rows <- grep("^mt-", gene_symbols, ignore.case = TRUE)
        assert(length(mt_rows) >= 5L, paste("Too few mouse mitochondrial genes; inspect symbols:", gsm))
        nCount <- as.numeric(Matrix::colSums(x))
        nFeature <- as.numeric(Matrix::colSums(x > 0))
        mtCount <- as.numeric(Matrix::colSums(x[mt_rows, , drop = FALSE]))
        pctMt <- ifelse(nCount > 0, 100 * mtCount/nCount, 0)
        fail_low <- nFeature < MIN_FEATURES
        fail_mt <- pctMt >= MAX_MT_PERCENT
        pass <- !fail_low & !fail_mt
        reason <- ifelse(pass, "PASS", ifelse(fail_low & fail_mt, "LOW_FEATURES_AND_HIGH_MT", ifelse(fail_low, 
            "LOW_FEATURES", "HIGH_MT")))
        md <- data.frame(cell_id = colnames(x), gsm = gsm, dpi = s$dpi[[1L]], replicate = s$replicate[[1L]], 
            intended_group = s$intended_group[[1L]], age_assignment_status = s$age_assignment_status[[1L]], 
            original_barcode = raw_barcode, nCount_RNA = nCount, nFeature_RNA = nFeature, percent.mt = pctMt, 
            mt_UMI = mtCount, fail_low_features = fail_low, fail_high_mt = fail_mt, QC_pass = pass, QC_reason = reason, 
            stringsAsFactors = FALSE)
        assert(sum(pass) > 0L, paste("No cells passed primary QC:", gsm))
        all_barcode_data[[i]] <- md
        stat <- data.frame(gsm = gsm, dpi = s$dpi[[1L]], replicate = s$replicate[[1L]], intended_group = s$intended_group[[1L]], 
            age_assignment_status = s$age_assignment_status[[1L]], raw_barcodes = nrow(md), qc_retained = sum(pass), 
            qc_excluded = sum(!pass), qc_retention_percent = 100 * mean(pass), failed_low_features_only = sum(fail_low & 
                !fail_mt), failed_high_mt_only = sum(!fail_low & fail_mt), failed_both = sum(fail_low & 
                fail_mt), median_nFeature_before = median_safe(nFeature), median_nCount_before = median_safe(nCount), 
            median_percent_mt_before = median_safe(pctMt), median_nFeature_after = median_safe(nFeature[pass]), 
            median_nCount_after = median_safe(nCount[pass]), median_percent_mt_after = median_safe(pctMt[pass]), 
            p95_nFeature_after = quantile_safe(nFeature[pass]), p95_nCount_after = quantile_safe(nCount[pass]), 
            n_mitochondrial_genes = length(mt_rows), gene_features = nrow(x), stringsAsFactors = FALSE)
        all_sample_stats[[i]] <- stat
        log_msg("Counts: input ", nrow(md), ", pass ", sum(pass), " (", round(100 * mean(pass), 2), "%), failed ", 
            sum(!pass))
        log_msg("  Median genes/UMI/mt pre: ", median_safe(nFeature), " / ", median_safe(nCount), " / ", 
            round(median_safe(pctMt), 2), "%")
        log_msg("  Median genes/UMI/mt post: ", median_safe(nFeature[pass]), " / ", median_safe(nCount[pass]), 
            " / ", round(median_safe(pctMt[pass]), 2), "%")
        if (DRAW_FIGURES) {
            tryCatch(plot_sample_qc(md, gsm), error = function(e) {
                log_msg("WARNING: plot failed for ", gsm, ": ", conditionMessage(e))
            })
        }
        if (SAVE_RDS) {
            rds_file <- file.path(RDS_DIR, paste0("116_", gsm, "_QC_filtered_raw_counts_v02.rds"))
            result <- list(counts = x[, pass, drop = FALSE], cell_metadata = md[pass, , drop = FALSE], 
                feature_metadata = gene_map, provenance = list(dataset = "GSE281205", gsm = gsm, intended_group = s$intended_group[[1L]], 
                  geo_title = s$geo_title[[1L]], geo_treatment = s$geo_treatment[[1L]], age_assignment_status = s$age_assignment_status[[1L]], 
                  note = "CELL QC ONLY: raw GEX UMI counts; neither decontaminated nor doublet-filtered", 
                  thresholds = list(min_features = MIN_FEATURES, mitochondrial_percent_strictly_less_than = MAX_MT_PERCENT), 
                  input_files = basename(c(file_b, file_f, file_m))))
            saveRDS(result, file = rds_file, compress = "gzip")
            rds_paths <- c(rds_paths, rds_file)
            rm(result)
        }
        rm(x, md, bar, feat, gene_map)
        invisible(gc(verbose = FALSE))
    }
    sumtab <- do.call(rbind, all_sample_stats)
    barcode <- do.call(rbind, all_barcode_data)
    fmap <- do.call(rbind, feature_concordance)
    save_csv(sumtab, "03_116_per_sample_cell_QC_summary_v02.csv")
    if (SAVE_BARCODE_QC) {
        save_csv(barcode, "02_116_per_barcode_cell_QC_v02.csv")
    }
    save_csv(fmap, "05_116_feature_mapping_concordance_v02.csv")
    failtab <- as.data.frame(table(barcode$gsm, barcode$QC_reason), stringsAsFactors = FALSE)
    names(failtab) <- c("gsm", "QC_reason", "n_barcodes")
    failtab <- failtab[failtab$n_barcodes > 0L, , drop = FALSE]
    save_csv(failtab, "04_116_QC_failure_reason_v02.csv")
    stagelist <- split(sumtab, sumtab$dpi)
    stages <- do.call(rbind, lapply(stagelist, function(t) {
        data.frame(dpi = t$dpi[[1L]], n_gsm_libraries = nrow(t), raw_barcodes = sum(t$raw_barcodes), 
            qc_retained = sum(t$qc_retained), qc_excluded = sum(t$qc_excluded), qc_retention_percent = 100 * 
                sum(t$qc_retained)/sum(t$raw_barcodes), min_cells_per_gsm = min(t$qc_retained), max_cells_per_gsm = max(t$qc_retained), 
            age_metadata_conflict_samples = sum(grepl("PROVISIONAL", t$age_assignment_status)), stringsAsFactors = FALSE)
    }))
    rownames(stages) <- NULL
    save_csv(stages, "07_116_stage_cell_QC_summary_v02.csv")
    if (DRAW_FIGURES) {
        tryCatch(plot_overview(sumtab), error = function(e) {
            log_msg("WARNING: overview plot failed: ", conditionMessage(e))
        })
    }
    tabtxt <- capture.output(print(sumtab[, c("gsm", "dpi", "raw_barcodes", "qc_retained", "qc_retention_percent", 
        "median_nFeature_after", "median_percent_mt_after")], row.names = FALSE))
    writeLines("Numerical results are provided in the accompanying output tables.", 
        con = file.path(OUT_DIR, "08_116_REPORT_v02.txt"), useBytes = TRUE)
    writeLines(capture.output(utils::sessionInfo()), con = file.path(OUT_DIR, "09_116_sessionInfo_v02.txt"), 
        useBytes = TRUE)
    index <- data.frame(relative_path = c("00_116_run_log_v02.log", "01_116_sample_metadata_audit_v02.csv", 
        if (SAVE_BARCODE_QC) "02_116_per_barcode_cell_QC_v02.csv", "03_116_per_sample_cell_QC_summary_v02.csv", 
        "04_116_QC_failure_reason_v02.csv", "05_116_feature_mapping_concordance_v02.csv", "06_116_gene_feature_mapping_reference_v02.csv", 
        "07_116_stage_cell_QC_summary_v02.csv", "08_116_REPORT_v02.txt", "09_116_sessionInfo_v02.txt", 
        if (DRAW_FIGURES) file.path("figures", list.files(FIG_DIR, pattern = "^Fig_116_.*_v02\\.png$")), 
        if (SAVE_RDS) file.path("RDS", basename(rds_paths))), stringsAsFactors = FALSE)
    save_csv(index, "10_116_OUTPUT_INDEX_v02.csv")
    invisible(sumtab)
}

tryCatch(main(), error = function(e) {
  log_msg("\nERROR: ", conditionMessage(e))
  stop(e)
})

source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

ROOT <- PROJECT_ROOT
IN116_DIR <- file.path(ROOT, "results", "116_GSE281205_cell_QC")
IN118_DIR <- file.path(ROOT, "results", "118_GSE281205_DecontX")
OUT_DIR <- file.path(ROOT, "results", "119_GSE281205_scDblFinder")
RDS_DIR <- file.path(OUT_DIR, "RDS")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
SAMPLE_IDS <- paste0("GSM", 8614662:8614667)
EXPECTED_DPI <- rep(c(14L, 28L), each = 3L)
RUN_VERSION <- "v02"
SEED_BASE <- 11900L
DBR_PER_1K <- 0.008
SAVE_SINGLET_RAW_RDS <- TRUE
SAVE_SINGLET_DECONTX_RDS <- TRUE
MAKE_FIGURES <- TRUE

check <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}
for (p in c("Matrix", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder", "BiocParallel")) {
    check(requireNamespace(p, quietly = TRUE), paste0("Required package missing: ", p, ". Install BiocManager then BiocManager::install('", 
        p, "')."))
}
check(dir.exists(file.path(IN116_DIR, "RDS")),
      paste("Cannot locate 116 QC RDS directory:", file.path(IN116_DIR, "RDS")))
check(dir.exists(IN118_DIR), paste("Cannot locate 118 directory:", IN118_DIR))
for (d in c(OUT_DIR, RDS_DIR, FIG_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
LOG_FILE <- file.path(OUT_DIR, "00_119_run_log_v02.log")
cat("", file = LOG_FILE)
log_msg <- function(...) {
  line <- paste0(..., collapse = "")
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG_FILE, append = TRUE)
  invisible(line)
}
save_table <- function(df, file) {
    check(is.data.frame(df), paste("Not a data.frame:", file))
    utils::write.csv(df, file.path(OUT_DIR, file), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
med <- function(x) {
  if (length(x) > 0L && !all(is.na(x))) {
    return(as.numeric(stats::median(x, na.rm = TRUE)))
  }
  NA_real_
}
read_input_list <- function(file, gsm, corrected = FALSE) {
  check(file.exists(file), paste("Input RDS not found:", file))
  obj <- readRDS(file)
  req <- if (corrected) {
    c("counts_decontX", "cell_metadata", "feature_metadata", "provenance")
  } else {
    c("counts", "cell_metadata", "feature_metadata", "provenance")
  }
  check(is.list(obj) && all(req %in% names(obj)),
        paste("Input RDS structure mismatch:", basename(file)))
  x <- if (corrected) obj$counts_decontX else obj$counts
  check(inherits(x, "sparseMatrix"),
        paste("Matrix is not sparse:", basename(file)))
  x <- as(x, "dgCMatrix")
  check(is.data.frame(obj$cell_metadata) && is.data.frame(obj$feature_metadata),
        paste("Metadata is not a data.frame:", basename(file)))
  check(ncol(x) == nrow(obj$cell_metadata) &&
          nrow(x) == nrow(obj$feature_metadata),
        paste("Matrix/metadata dimensions disagree:", basename(file)))
  check(identical(colnames(x), as.character(obj$cell_metadata$cell_id)),
        paste("Cell IDs/metadata order disagree:", basename(file)))
  check(identical(rownames(x), as.character(obj$feature_metadata$gene_id)),
        paste("Gene IDs/feature-metadata order disagree:", basename(file)))
  check(all(as.character(obj$cell_metadata$gsm) == gsm),
        paste("Mixed/incorrect GSM identifiers:", basename(file)))
  if (!corrected) {
    check(all(is.finite(x@x)) && all(x@x >= 0) &&
            all(abs(x@x - round(x@x)) < 1e-8),
          paste("116 counts must be nonnegative INTEGER UMI:", gsm))
    check("QC_pass" %in% names(obj$cell_metadata) &&
            all(obj$cell_metadata$QC_pass),
          paste("116 RDS contains cells failing step 116 QC:", gsm))
  } else {
    check(all(is.finite(x@x)) && all(x@x >= 0),
          paste("Invalid corrected DecontX values:", gsm))
  }
  obj[[if (corrected) "counts_decontX" else "counts"]] <- x
  obj
}

plot_gsm <- function(md, gsm) {
    f <- file.path(FIG_DIR, paste0("Fig_119_", gsm, "_scDblFinder_QC_v02.png"))
    grDevices::png(f, width = 1750, height = 1320, res = 155)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3.3, 1.1), oma = c(0, 0, 1.6, 0))
    scores <- md$scDblFinder_score
    graphics::hist(scores, breaks = 55, col = "#7FAAC2", border = "white", main = "scDblFinder score distribution", 
        xlab = "scDblFinder.score")
    grp <- factor(md$scDblFinder_class, levels = c("singlet", "doublet"))
    graphics::boxplot(log10(md$nCount_RNA + 1) ~ grp, col = c("#67A98E", "#D89078"), outline = FALSE, 
        main = "Original UMI depth by call", xlab = "Predicted class", ylab = "log10(original UMI + 1)")
    cols <- ifelse(md$scDblFinder_class == "doublet", "#CE7458", "#3C8C81")
    graphics::plot(md$nCount_RNA + 1, md$nFeature_RNA + 1, log = "xy", pch = 16, cex = 0.3, col = grDevices::adjustcolor(cols, 
        alpha.f = 0.3), main = "Original UMI vs detected genes", xlab = "Original UMI + 1 (log)", ylab = "Detected genes + 1 (log)")
    contamination_group <- cut(md$decontX_contamination_percent, breaks = c(-Inf, 20, 50, Inf), right = FALSE, 
        labels = c("<20%", "20-<50%", ">=50%"))
    tt <- table(contamination_group, grp)
    rate <- ifelse(rowSums(tt) > 0, 100 * tt[, "doublet"]/rowSums(tt), NA_real_)
    graphics::barplot(rate, names.arg = names(rate), col = "#658DB2", ylim = c(0, max(5, rate, na.rm = TRUE) * 
        1.2), main = "Doublet calls by DecontX group", ylab = "Predicted doublets (%)")
    graphics::mtext(paste0(gsm, " | input ", nrow(md), "; singlet ", sum(md$scDblFinder_class == "singlet"), 
        "; doublet ", sum(md$scDblFinder_class == "doublet")), side = 3, outer = TRUE, line = 0.25, cex = 1.05, 
        font = 2)
    invisible(f)
}
plot_overview <- function(tab) {
    f <- file.path(FIG_DIR, "Fig_119_all_samples_scDblFinder_overview_v02.png")
    grDevices::png(f, width = 1620, height = 820, res = 145)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(1, 2), mar = c(6.2, 4.6, 3.5, 0.6))
    labels <- paste0(sub("GSM", "", tab$gsm), "\n", tab$dpi, " dpi")
    graphics::barplot(rbind(tab$n_input_QC_barcodes, tab$n_predicted_singlets), beside = TRUE, names.arg = labels, 
        las = 1, col = c("#9CB5CE", "#3A9B82"), ylab = "Barcodes", main = "116 QC input and retained singlets")
    graphics::legend("topleft", legend = c("Input", "Singlets"), fill = c("#9CB5CE", "#3A9B82"), bty = "n")
    graphics::barplot(tab$predicted_doublet_percent, names.arg = labels, las = 1, col = "#CA785B", ylim = c(0, 
        max(5, tab$predicted_doublet_percent, na.rm = TRUE) * 1.15), ylab = "Called doublets (%)", main = "scDblFinder by GSM")
    invisible(f)
}

main <- function() {

  contamination_file <- file.path(IN118_DIR,
                                   "02_118_per_cell_contamination_v01.csv")
  check(file.exists(contamination_file),
        paste("Missing 118 contamination table:", contamination_file))
  contamination <- utils::read.csv(contamination_file,
                                   check.names = FALSE, stringsAsFactors = FALSE)
  check(all(c("cell_id", "gsm", "dpi", "decontX_contamination_fraction",
              "decontX_contamination_percent") %in% names(contamination)),
        "118 per-cell table is missing mandatory matching fields")
  check(!anyNA(contamination$cell_id) &&
          anyDuplicated(contamination$cell_id) == 0L,
        "118 per-cell cell_id contains NA or duplicates")
  check(setequal(sort(unique(contamination$gsm)), sort(SAMPLE_IDS)),
        "118 contamination data does not contain exactly the expected six GSMs")

  manifest <- vector("list", length(SAMPLE_IDS))
  per_cell <- vector("list", length(SAMPLE_IDS))
  sample_summary <- vector("list", length(SAMPLE_IDS))
  cross_tables <- vector("list", length(SAMPLE_IDS))
  out_files <- character(0)

  for (i in seq_along(SAMPLE_IDS)) {
    gsm <- SAMPLE_IDS[[i]]
    expected_dpi <- EXPECTED_DPI[[i]]
    fin116 <- file.path(IN116_DIR, "RDS",
                         paste0("116_", gsm, "_QC_filtered_raw_counts_v02.rds"))
    fin118 <- file.path(IN118_DIR, "RDS",
                         paste0("118_", gsm, "_DecontX_corrected_counts_v01.rds"))
    original <- read_input_list(fin116, gsm, corrected = FALSE)
    raw <- original$counts
    md <- original$cell_metadata
    feat <- original$feature_metadata
    check(identical(unique(as.integer(md$dpi)), expected_dpi),
          paste("Unexpected dpi in input:", gsm))
    check(all(md$nCount_RNA == as.numeric(Matrix::colSums(raw))),
          paste("Mismatch between 116 original UMI and QC metadata:", gsm))
    check(!anyDuplicated(colnames(raw)) && !anyDuplicated(rownames(raw)),
          paste("Duplicate cell/gene identifiers:", gsm))
    record <- contamination[contamination$gsm == gsm, , drop = FALSE]
    ix <- match(colnames(raw), record$cell_id)
    check(nrow(record) == ncol(raw) && !anyNA(ix) &&
            setequal(record$cell_id, colnames(raw)),
          paste("118 contamination rows do not match 116 cells 1:1:", gsm))
    record <- record[ix, , drop = FALSE]
    check(all(record$gsm == gsm) &&
            all(as.integer(record$dpi) == expected_dpi),
          paste("118 sample/stage mismatch:", gsm))
    contamination_frac <- as.numeric(record$decontX_contamination_fraction)
    check(all(is.finite(contamination_frac)) &&
            all(contamination_frac >= 0 & contamination_frac <= 1),
          paste("118 invalid contamination estimate:", gsm))
    check(all(abs(100 * contamination_frac -
                    as.numeric(record$decontX_contamination_percent)) < 1e-5),
          paste("118 fraction/percentage disagreement:", gsm))

    usable <- which(Matrix::rowSums(raw) > 0)
    check(length(usable) > 1000L,
          paste("Too few expressed features for scDblFinder:", gsm))
    classifier_counts <- raw[usable, , drop = FALSE]
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = classifier_counts))
    check(identical(colnames(sce), colnames(raw)),
          paste("SCE cell order changed:", gsm))
    seed <- as.integer(SEED_BASE + i)
    set.seed(seed)
    time_begin <- Sys.time()
    modeled <- tryCatch(
      scDblFinder::scDblFinder(
        sce, clusters = NULL, samples = NULL,
        dbr = NULL, dbr.per1k = DBR_PER_1K,
        BPPARAM = BiocParallel::SerialParam(), verbose = TRUE),
      error = function(e) {
        log_msg("ERROR scDblFinder failed for ", gsm, ": ", conditionMessage(e))
        stop(e)
      }
    )
    elapsed_min <- as.numeric(difftime(Sys.time(), time_begin, units = "mins"))
    rescol <- SummarizedExperiment::colData(modeled)
    check(all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(rescol)),
          paste("Missing scDblFinder result columns:", gsm))
    check(identical(colnames(modeled), colnames(raw)),
          paste("scDblFinder changed cell order:", gsm))
    calls <- as.character(rescol$scDblFinder.class)
    score <- as.numeric(rescol$scDblFinder.score)
    check(!anyNA(calls) && all(calls %in% c("singlet", "doublet")),
          paste("Unknown scDblFinder class:", gsm))
    check(all(is.finite(score)), paste("Invalid doublet scores:", gsm))
    keep <- calls == "singlet"
    check(sum(keep) > 0L, paste("Zero predicted singlets:", gsm))

    new_md <- md
    new_md$scDblFinder_score <- score
    new_md$scDblFinder_class <- calls
    new_md$predicted_singlet <- keep
    new_md$decontX_contamination_fraction <- contamination_frac
    new_md$decontX_contamination_percent <- 100 * contamination_frac
    new_md$decontX_contamination_ge20 <- contamination_frac >= 0.2
    new_md$decontX_contamination_ge50 <- contamination_frac >= 0.5
    new_md$decontX_contamination_ge80 <- contamination_frac >= 0.8
    new_md$scDblFinder_seed <- seed
    per_cell[[i]] <- new_md[, c(
      "cell_id", "gsm", "dpi", "replicate", "intended_group",
      "age_assignment_status", "nCount_RNA", "nFeature_RNA", "percent.mt",
      "scDblFinder_score", "scDblFinder_class", "predicted_singlet",
      "decontX_contamination_fraction", "decontX_contamination_percent",
      "decontX_contamination_ge20", "decontX_contamination_ge50",
      "decontX_contamination_ge80", "scDblFinder_seed"), drop = FALSE]

    hgrp <- cut(100 * contamination_frac, breaks = c(-Inf, 20, 50, Inf),
                right = FALSE, labels = c("<20%", "20-<50%", ">=50%"))
    tab <- as.data.frame(table(contamination_group = hgrp,
                               scDblFinder_class = factor(calls,
                                  levels = c("singlet", "doublet"))))
    tab$gsm <- gsm
    tab$dpi <- expected_dpi
    tab$n_class_group <- as.integer(tab$Freq)
    tab$Freq <- NULL
    cross_tables[[i]] <- tab[, c("gsm", "dpi", "contamination_group",
                                "scDblFinder_class", "n_class_group")]

    sample_summary[[i]] <- data.frame(
      gsm = gsm, dpi = expected_dpi,
      replicate = unique(md$replicate)[[1]],
      intended_group = unique(md$intended_group)[[1]],
      age_assignment_status = unique(md$age_assignment_status)[[1]],
      n_input_QC_barcodes = ncol(raw),
      n_predicted_singlets = sum(keep),
      n_predicted_doublets = sum(!keep),
      predicted_doublet_percent = 100 * mean(!keep),
      nominal_10x_expected_percent_approx =
        100 * min(1, DBR_PER_1K * ncol(raw) / 1000),
      score_median_all = med(score),
      score_median_singlet = med(score[keep]),
      score_median_doublet = med(score[!keep]),
      n_contamination_ge20 = sum(contamination_frac >= 0.2),
      n_contamination_ge50 = sum(contamination_frac >= 0.5),
      n_contamination_ge80 = sum(contamination_frac >= 0.8),
      n_doublets_and_contamination_ge20 = sum(!keep & contamination_frac >= 0.2),
      n_doublets_and_contamination_ge50 = sum(!keep & contamination_frac >= 0.5),
      n_singlets_and_contamination_ge50 = sum(keep & contamination_frac >= 0.5),
      median_contamination_singlets_percent = med(100 * contamination_frac[keep]),
      median_contamination_doublets_percent = med(100 * contamination_frac[!keep]),
      raw_UMI_median_singlet = med(new_md$nCount_RNA[keep]),
      raw_UMI_median_doublet = med(new_md$nCount_RNA[!keep]),
      classifier_genes_expressed = length(usable),
      run_minutes = elapsed_min,
      stringsAsFactors = FALSE)

    rm(sce, modeled, classifier_counts, rescol)
    invisible(gc(verbose = FALSE))

    if (SAVE_SINGLET_RAW_RDS) {
      raw_name <- paste0("119_", gsm, "_singlet_raw_counts_v02.rds")
      raw_outfile <- file.path(RDS_DIR, raw_name)
      saved <- list(
        counts = raw[, keep, drop = FALSE],
        cell_metadata = new_md[keep, , drop = FALSE],
        feature_metadata = feat,
        provenance = list(
          input_116 = fin116, input_118_contamination = contamination_file,
          original_116_provenance = original$provenance,
          method = "scDblFinder on per-GSM original integer UMI counts",
          scDblFinder_version = as.character(utils::packageVersion("scDblFinder")),
          seed = seed, clusters = NULL, dbr = NULL,
          dbr_per1k = DBR_PER_1K, n_input = ncol(raw),
          n_predicted_singlets = sum(keep),
          note = "Single-GSM singlet candidates only; no cell-type annotation. Original counts preserved."
        )
      )
      saveRDS(saved, raw_outfile, compress = "gzip")
      check(file.exists(raw_outfile), paste("Failed saving:", raw_outfile))
      out_files <- c(out_files, raw_outfile)
      rm(saved)
      invisible(gc(verbose = FALSE))
    }

    if (SAVE_SINGLET_DECONTX_RDS) {
      corrected <- read_input_list(fin118, gsm, corrected = TRUE)
      corx <- corrected$counts_decontX
      check(identical(rownames(raw), rownames(corx)) &&
              identical(colnames(raw), colnames(corx)),
            paste("118 corrected matrix differs from 116 in gene/cell order:", gsm))
      check(all(as.character(corrected$cell_metadata$cell_id) == colnames(raw)),
            paste("118 corrected metadata has mismatched cells:", gsm))
      if ("decontX_contamination_fraction" %in%
          names(corrected$cell_metadata)) {
        check(all(abs(as.numeric(corrected$cell_metadata$decontX_contamination_fraction) -
                        contamination_frac) < 1e-6),
              paste("118 RDS and 118 CSV disagree on contamination:", gsm))
      }
      corr_name <- paste0("119_", gsm, "_singlet_DecontX_counts_v02.rds")
      corr_outfile <- file.path(RDS_DIR, corr_name)
      saved_corr <- list(
        counts_decontX = corx[, keep, drop = FALSE],
        cell_metadata = new_md[keep, , drop = FALSE],
        feature_metadata = feat,
        provenance = list(
          input_118_corrected = fin118, input_116_raw = fin116,
          original_118_provenance = corrected$provenance,
          doublet_detection_method = "scDblFinder on original UMI (116), applied to 118 corrected counts",
          cell_alignment_verified = TRUE,
          note = "DecontX estimated (possibly fractional) counts for annotation only; original UMI saved separately."
        )
      )
      saveRDS(saved_corr, corr_outfile, compress = "gzip")
      check(file.exists(corr_outfile), paste("Failed saving:", corr_outfile))
      out_files <- c(out_files, corr_outfile)
      rm(corrected, corx, saved_corr)
      invisible(gc(verbose = FALSE))
    }

    if (MAKE_FIGURES) {
      tryCatch(plot_gsm(new_md, gsm), error = function(e) {
        log_msg("WARNING: plot failed ", gsm, ": ", conditionMessage(e))
      })
    }
    manifest[[i]] <- data.frame(
      gsm = gsm, dpi = expected_dpi,
      n_cells_input = ncol(raw), n_singlets = sum(keep),
      source_116_raw_RDS = fin116,
      source_118_corrected_RDS = if (SAVE_SINGLET_DECONTX_RDS) fin118 else NA_character_,
      source_118_per_cell_csv = contamination_file,
      output_singlet_raw_RDS = if (SAVE_SINGLET_RAW_RDS) {
        file.path(RDS_DIR, paste0("119_", gsm, "_singlet_raw_counts_v02.rds"))
      } else {
        NA_character_
      },
      output_singlet_DecontX_RDS = if (SAVE_SINGLET_DECONTX_RDS) {
        file.path(RDS_DIR, paste0("119_", gsm, "_singlet_DecontX_counts_v02.rds"))
      } else {
        NA_character_
      },
      seed = seed, source_classification = "original_integer_UMI",
      stringsAsFactors = FALSE)
    writeLines(c(
      paste("GSM:", gsm), paste("Timestamp:", Sys.time()),
      paste("scDblFinder version:", as.character(utils::packageVersion("scDblFinder"))),
      paste("Original UMI cells:", ncol(raw)),
      paste("Usable features:", length(usable)),
      paste("Seed:", seed),
      paste("Default automatic doublet rate; dbr.per1k:", DBR_PER_1K),
      paste("Called doublets:", sum(!keep)),
      paste("Retained singlets:", sum(keep)),
      paste("Doublet percent:", round(100 * mean(!keep), 3)),
      paste("High contamination (>=50%) remaining singlets:",
            sum(keep & contamination_frac >= 0.5)),
      "Contamination percentage is not a deletion criterion."
    ), file.path(LOG_DIR, paste0("119_", gsm, "_scDblFinder_model_v02.log")),
    useBytes = TRUE)
    rm(original, raw, md, feat, new_md)
    invisible(gc(verbose = FALSE))
  }

  man <- do.call(rbind, manifest)
  cells <- do.call(rbind, per_cell)
  stats <- do.call(rbind, sample_summary)
  ct <- do.call(rbind, cross_tables)
  save_table(man, "01_119_input_output_manifest_v02.csv")
  save_table(cells, "02_119_per_cell_doublet_calls_v02.csv")
  save_table(stats, "03_119_per_sample_doublet_summary_v02.csv")
  save_table(ct, "04_119_doublet_vs_contamination_v02.csv")
  stage <- do.call(rbind, lapply(split(stats, stats$dpi), function(d) {
    data.frame(
      dpi = d$dpi[[1L]], n_gsm = nrow(d),
      qc_input_barcodes = sum(d$n_input_QC_barcodes),
      predicted_singlets = sum(d$n_predicted_singlets),
      predicted_doublets = sum(d$n_predicted_doublets),
      doublet_call_percent = 100 * sum(d$n_predicted_doublets) /
        sum(d$n_input_QC_barcodes),
      remaining_singlets_contamination_ge50 =
        sum(d$n_singlets_and_contamination_ge50),
      n_provisional_age_gsm = sum(grepl("PROVISIONAL", d$age_assignment_status)),
      stringsAsFactors = FALSE)
  }))
  rownames(stage) <- NULL
  save_table(stage, "05_119_stage_doublet_summary_v02.csv")
  if (MAKE_FIGURES) {
    tryCatch(plot_overview(stats), error = function(e) {
      log_msg("WARNING: overall plot failed: ", conditionMessage(e))
    })
  }

  total_input <- sum(stats$n_input_QC_barcodes)
  total_singlets <- sum(stats$n_predicted_singlets)
  total_doublets <- sum(stats$n_predicted_doublets)
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, "06_119_REPORT_v02.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()),
             file.path(OUT_DIR, "07_119_sessionInfo_v02.txt"), useBytes = TRUE)
  files <- c(
    "00_119_run_log_v02.log", "01_119_input_output_manifest_v02.csv",
    "02_119_per_cell_doublet_calls_v02.csv",
    "03_119_per_sample_doublet_summary_v02.csv",
    "04_119_doublet_vs_contamination_v02.csv",
    "05_119_stage_doublet_summary_v02.csv", "06_119_REPORT_v02.txt",
    "07_119_sessionInfo_v02.txt",
    file.path("RDS", basename(out_files)),
    file.path("logs", list.files(LOG_DIR, pattern = "^119_GSM.*_v02\\.log$")),
    if (MAKE_FIGURES) file.path("figures",
      list.files(FIG_DIR, pattern = "^Fig_119_.*_v02\\.png$"))
  )
  index <- data.frame(relative_path = files,
                      file_exists = file.exists(file.path(OUT_DIR, files)))
  save_table(index, "08_119_OUTPUT_INDEX_v02.csv")
  invisible(stats)
}

tryCatch(main(), error = function(e) {
  log_msg("\nSTOPPED: ", conditionMessage(e))
  stop(e)
})

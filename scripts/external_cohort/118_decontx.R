source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

ROOT <- PROJECT_ROOT
SOURCE_DIR <- file.path(ROOT, "results", "116_GSE281205_cell_QC")
IN_RDS_DIR <- file.path(SOURCE_DIR, "RDS")
OUT_DIR <- file.path(ROOT, "results", "118_GSE281205_DecontX")
OUT_RDS_DIR <- file.path(OUT_DIR, "RDS")
FIG_DIR <- file.path(OUT_DIR, "figures")
MODEL_LOG_DIR <- file.path(OUT_DIR, "logs")

SAMPLE_IDS <- paste0("GSM", 8614662:8614667)
SEED_BASE <- 11800L
RUN_VERSION <- "v01"
SAVE_PER_CELL <- TRUE
SAVE_PLOTS <- TRUE

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}
assert(dir.exists(IN_RDS_DIR), paste0("Missing 116 RDS folder: ", IN_RDS_DIR))
assert(requireNamespace("Matrix", quietly = TRUE), "Package Matrix is required: install.packages('Matrix')")
if (requireNamespace("decontX", quietly = TRUE)) {
    DECONTX_FUNCTION <- decontX::decontX
    DECONTX_PACKAGE <- "decontX"
} else if (requireNamespace("celda", quietly = TRUE)) {
    DECONTX_FUNCTION <- celda::decontX
    DECONTX_PACKAGE <- "celda"
} else {
    stop(paste0("Missing DecontX implementation. Install BiocManager and then run: ", "BiocManager::install('decontX') [or 'celda' for an older release]."), 
        call. = FALSE)
}

for (d in c(OUT_DIR, OUT_RDS_DIR, FIG_DIR, MODEL_LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOG_FILE <- file.path(OUT_DIR, "00_118_run_log_v01.log")
cat("", file = LOG_FILE)
log_msg <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n", sep = "")
  cat(s, "\n", file = LOG_FILE, append = TRUE, sep = "")
}
write_table <- function(x, filename) {
    assert(is.data.frame(x), paste0("Expected data frame: ", filename))
    utils::write.csv(x, file.path(OUT_DIR, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
med <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else as.numeric(stats::median(x, na.rm = TRUE))
qtl <- function(x, p) if (!length(x) || all(is.na(x))) NA_real_ else
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

plot_gsm <- function(gsm, md, raw_umi, corrected_umi) {
  fout <- file.path(FIG_DIR, paste0("Fig_118_", gsm, "_DecontX_v01.png"))
  grDevices::png(fout, width = 1760, height = 1280, res = 155)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.6, 3.2, 1.2),
                oma = c(0, 0, 2, 0))
  pct <- md$decontX_contamination_percent
  graphics::hist(pct, breaks = 60, col = "#86B5CA", border = "white",
                 xlab = "Estimated contamination (%)",
                 main = "DecontX contamination estimate")
  graphics::abline(v = med(pct), col = "#BC4935", lwd = 2)
  graphics::plot(raw_umi + 1, pct, log = "x", pch = 16, cex = 0.30,
                 col = grDevices::adjustcolor("#336F9D", alpha.f = 0.27),
                 xlab = "Raw UMI + 1 (log scale)",
                 ylab = "Estimated contamination (%)",
                 main = "Contamination versus original depth")
  graphics::plot(raw_umi + 1, corrected_umi + 1,
                 log = "xy", pch = 16, cex = 0.30,
                 col = grDevices::adjustcolor("#2C8C70", alpha.f = 0.27),
                 xlab = "Original UMI + 1", ylab = "DecontX UMI + 1",
                 main = "Original vs corrected counts")
  ul <- max(c(raw_umi, corrected_umi), na.rm = TRUE) + 1
  graphics::abline(a = 0, b = 1, col = "#C14B41", lty = 2, lwd = 2)

  graphics::hist(100 * corrected_umi / pmax(raw_umi, 1), breaks = 60,
                 col = "#70A88E", border = "white",
                 xlab = "Corrected / original UMI (%)",
                 main = "Retained expression per barcode")
  graphics::mtext(paste0(gsm, " | cells = ", nrow(md),
                          " | DecontX (no empty-droplet background)"),
                  side = 3, outer = TRUE, line = 0.3,
                  font = 2, cex = 1.05)
  invisible(fout)
}

plot_overall <- function(sample_stats) {
    f <- file.path(FIG_DIR, "Fig_118_all_samples_DecontX_overview_v01.png")
    grDevices::png(f, width = 1580, height = 820, res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(1, 2), mar = c(6.5, 4.3, 3, 1))
    labs <- paste0(sub("GSM", "", sample_stats$gsm), "\n", sample_stats$dpi, " dpi")
    graphics::barplot(sample_stats$median_contamination_percent, names.arg = labs, las = 1, col = "#6EA0BF", 
        main = "Median estimated contamination", ylab = "Percent of transcripts (%)")
    graphics::barplot(sample_stats$corrected_over_raw_UMI_percent, names.arg = labs, las = 1, col = "#3A8C75", 
        ylim = c(0, max(100, sample_stats$corrected_over_raw_UMI_percent, na.rm = TRUE) * 1.03), main = "Total UMI retained after DecontX", 
        ylab = "Corrected / original UMI (%)")
    invisible(f)
}

MARKERS <- c("Clec7a", "P2ry12", "Tmem119", "Cx3cr1", "Trem2", "Apoe", "Lpl", "Gpnmb", "Lyz2", "Csf1r", 
    "Tyrobp", "S100a8", "S100a9", "Mbp", "Plp1", "Aqp4", "Pecam1")
marker_audit <- function(gsm, stage, raw, corr, feature_metadata) {
  out <- lapply(MARKERS, function(sym) {
    idx <- which(as.character(feature_metadata$gene_symbol) == sym)
    if (!length(idx)) {
      return(data.frame(gsm = gsm, dpi = stage, gene_symbol = sym,
                        n_matched_gene_ids = 0L, raw_detection_pct = NA_real_,
                        corrected_rounded_detection_pct = NA_real_,
                        raw_total_UMI = NA_real_,
                        corrected_total_estimated_UMI = NA_real_))
    }

    a <- as.numeric(Matrix::colSums(raw[idx, , drop = FALSE]))
    b <- as.numeric(Matrix::colSums(corr[idx, , drop = FALSE]))
    data.frame(gsm = gsm, dpi = stage, gene_symbol = sym,
               n_matched_gene_ids = length(idx),
               raw_detection_pct = 100 * mean(a >= 1),
               corrected_rounded_detection_pct = 100 * mean(round(b) >= 1),
               raw_total_UMI = sum(a),
               corrected_total_estimated_UMI = sum(b))
  })
  do.call(rbind, out)
}

main <- function() {

  manifest <- vector("list", length(SAMPLE_IDS))
  per_cell <- vector("list", length(SAMPLE_IDS))
  sample_stats <- vector("list", length(SAMPLE_IDS))
  markers <- vector("list", length(SAMPLE_IDS))
  output_rds <- character(length(SAMPLE_IDS))

  qc_summary_file <- file.path(SOURCE_DIR,
                               "03_116_per_sample_cell_QC_summary_v02.csv")
  assert(file.exists(qc_summary_file),
         paste0("Missing 116 sample QC summary: ", qc_summary_file))
  qc_summary <- utils::read.csv(qc_summary_file, check.names = FALSE)
  assert(all(SAMPLE_IDS %in% qc_summary$gsm),
         "116 QC summary does not contain all six GSM identifiers")

  for (i in seq_along(SAMPLE_IDS)) {
    gsm <- SAMPLE_IDS[[i]]
    fin <- file.path(IN_RDS_DIR,
                      paste0("116_", gsm, "_QC_filtered_raw_counts_v02.rds"))
    assert(file.exists(fin), paste0("116 input RDS missing: ", fin))
    input <- readRDS(fin)
    assert(is.list(input) &&
             all(c("counts", "cell_metadata", "feature_metadata", "provenance") %in%
                   names(input)), paste0("116 RDS structure mismatch: ", gsm))
    raw <- input$counts
    assert(inherits(raw, "sparseMatrix"),
           paste0("Original UMI matrix is not sparse: ", gsm))
    raw <- as(raw, "dgCMatrix")
    md <- input$cell_metadata
    feature_md <- input$feature_metadata
    assert(is.data.frame(md) && is.data.frame(feature_md),
           paste0("Missing feature or cell metadata: ", gsm))
    assert(nrow(raw) == nrow(feature_md),
           paste0("Feature metadata/matrix row mismatch: ", gsm))
    assert(ncol(raw) == nrow(md),
           paste0("Cell metadata/matrix column mismatch: ", gsm))
    assert(identical(rownames(raw), as.character(feature_md$gene_id)),
           paste0("Gene IDs not aligned with feature metadata: ", gsm))
    assert(!anyDuplicated(rownames(raw)) && !anyDuplicated(colnames(raw)),
           paste0("Duplicate gene IDs or cell IDs: ", gsm))
    assert("cell_id" %in% names(md),
           paste0("cell_id not in 116 cell_metadata: ", gsm))
    assert(identical(colnames(raw), as.character(md$cell_id)),
           paste0("Cell IDs not aligned with metadata: ", gsm))
    assert(all(as.character(md$gsm) == gsm),
           paste0("Mixed sample identities in RDS: ", gsm))
    assert(all(is.finite(raw@x)) && all(raw@x >= 0) &&
             all(abs(raw@x - round(raw@x)) < 1e-8),
           paste0("Original matrix contains invalid UMI counts: ", gsm))
    assert("QC_pass" %in% names(md) && all(md$QC_pass),
           paste0("116 filtered RDS includes QC-failed cells: ", gsm))

    qc_expected <- qc_summary$qc_retained[match(gsm, qc_summary$gsm)]
    assert(ncol(raw) == qc_expected,
           paste0("Cell count disagrees with 116 QC report: ", gsm))
    stage <- unique(as.integer(md$dpi))
    assert(length(stage) == 1L && stage %in% c(14L, 28L),
           paste0("Unexpected or mixed dpi: ", gsm))
    ageflag <- unique(as.character(md$age_assignment_status))
    assert(length(ageflag) == 1L,
           paste0("Conflicting within-library age audit labels: ", gsm))
    if (stage == 28L) {
      log_msg("WARNING: provisional age assignment, original GEO label conflict: ",
              ageflag)
    }

    expressed <- which(Matrix::rowSums(raw) > 0)
    assert(length(expressed) >= 1000L,
           paste0("Too few expressed genes for DecontX: ", gsm))

    model_input <- raw[expressed, , drop = FALSE]
    model_log <- file.path(MODEL_LOG_DIR,
                           paste0("118_", gsm, "_DecontX_model_v01.log"))
    seed <- as.integer(SEED_BASE + i)
    set.seed(seed)
    result <- tryCatch(
      DECONTX_FUNCTION(x = model_input, z = NULL,
                       seed = seed, logfile = model_log, verbose = TRUE),
      error = function(e) {
        log_msg("ERROR: DecontX failed for ", gsm, ": ", conditionMessage(e))
        stop(e)
      }
    )
    assert(is.list(result) && all(c("decontXcounts", "contamination") %in%
                                      names(result)),
           paste0("Unexpected DecontX return structure: ", gsm))
    corrected_fit <- as(result$decontXcounts, "dgCMatrix")
    assert(identical(dim(corrected_fit), dim(model_input)),
           paste0("DecontX changed the model matrix dimensions: ", gsm))
    assert(all(is.finite(corrected_fit@x)) && all(corrected_fit@x >= -1e-8),
           paste0("Corrected values are non-finite/negative: ", gsm))
    corrected_fit@x[corrected_fit@x < 0] <- 0
    corrected_fit <- Matrix::drop0(corrected_fit)
    if (is.null(rownames(corrected_fit))) rownames(corrected_fit) <- rownames(model_input)
    if (is.null(colnames(corrected_fit))) colnames(corrected_fit) <- colnames(model_input)
    assert(identical(dimnames(corrected_fit), dimnames(model_input)),
           paste0("Gene/barcode order changed inside DecontX: ", gsm))

    corrected <- as(Matrix::Matrix(0, nrow = nrow(raw), ncol = ncol(raw),
                                   sparse = TRUE), "dgCMatrix")
    dimnames(corrected) <- dimnames(raw)
    corrected[expressed, ] <- corrected_fit
    corrected <- Matrix::drop0(corrected)
    assert(identical(dimnames(corrected), dimnames(raw)),
           paste0("Unable to restore complete gene IDs: ", gsm))

    contamination <- as.numeric(result$contamination)
    assert(length(contamination) == ncol(raw) &&
             all(is.finite(contamination)) &&
             all(contamination >= -1e-8) &&
             all(contamination <= 1 + 1e-8),
           paste0("Invalid DecontX contamination fractions: ", gsm))
    contamination <- pmax(0, pmin(contamination, 1))
    provisional_clusters <- if (!is.null(result$z) &&
                                length(result$z) == ncol(raw)) {
      as.character(result$z)
    } else rep(NA_character_, ncol(raw))

    raw_umi <- as.numeric(Matrix::colSums(raw))
    corrected_umi <- as.numeric(Matrix::colSums(corrected))
    assert(all(is.finite(corrected_umi)) && all(corrected_umi >= 0),
           paste0("Invalid corrected total counts: ", gsm))
    updated_md <- md
    updated_md$decontX_contamination_fraction <- contamination
    updated_md$decontX_contamination_percent <- 100 * contamination
    updated_md$decontX_heuristic_cluster <- provisional_clusters
    updated_md$decontX_raw_UMI <- raw_umi
    updated_md$decontX_corrected_estimated_UMI <- corrected_umi
    updated_md$decontX_retained_UMI_percent <-
      100 * corrected_umi / pmax(raw_umi, 1)

    outname <- paste0("118_", gsm, "_DecontX_corrected_counts_v01.rds")
    outfile <- file.path(OUT_RDS_DIR, outname)
    saved <- list(
      counts_decontX = corrected,
      cell_metadata = updated_md,
      feature_metadata = feature_md,
      provenance = list(
        dataset = "GSE281205", gsm = gsm, dpi = stage,
        original_116_provenance = input$provenance,
        age_assignment_status = ageflag,
        method = paste0(DECONTX_PACKAGE, "::decontX"),
        seed = seed, inferred_clusters = TRUE,
        z_input = NULL, background_input = NULL,
        background_note = "No independently measured empty-droplet matrix; inferred from other fitted cell clusters",
        modeling_rows_nonzero = length(expressed),
        n_genes_in_output = nrow(raw),
        raw_UMI_source_path = fin,
        raw_UMI_source_filename = basename(fin),
        raw_UMI_unchanged = TRUE,
        corrected_count_type = "fractional estimated native counts; DO NOT label raw UMI",
        next_step = "scDblFinder with explicit choice of original vs rounded corrected counts; separate original UMI for pseudobulk"
      )
    )
    saveRDS(saved, outfile, compress = "gzip")
    assert(file.exists(outfile), paste0("Failed to save output: ", outfile))
    output_rds[[i]] <- outfile

    if (SAVE_PLOTS) {
      tryCatch(plot_gsm(gsm, updated_md, raw_umi, corrected_umi),
               error = function(e) {
                 log_msg("WARNING: diagnostic figure failed for ", gsm,
                         ": ", conditionMessage(e))
               })
    }
    manifest[[i]] <- data.frame(
      gsm = gsm, dpi = stage, n_cells = ncol(raw),
      n_original_gene_features = nrow(raw),
      n_expressed_gene_features_modeled = length(expressed),
      source_116_rds = fin, output_118_rds = outfile,
      age_assignment_status = ageflag,
      background_available = FALSE,
      input_QC_passed = TRUE, stringsAsFactors = FALSE
    )
    per_cell[[i]] <- updated_md[, c(
      "cell_id", "gsm", "dpi", "replicate", "intended_group",
      "age_assignment_status", "nCount_RNA", "nFeature_RNA", "percent.mt",
      "decontX_contamination_fraction", "decontX_contamination_percent",
      "decontX_heuristic_cluster", "decontX_raw_UMI",
      "decontX_corrected_estimated_UMI", "decontX_retained_UMI_percent"
    ), drop = FALSE]
    sample_stats[[i]] <- data.frame(
      gsm = gsm, dpi = stage, n_cells = ncol(raw),
      age_assignment_status = ageflag,
      original_total_UMI = sum(raw_umi),
      corrected_total_estimated_UMI = sum(corrected_umi),
      corrected_over_raw_UMI_percent =
        100 * sum(corrected_umi) / sum(raw_umi),
      median_contamination_percent = med(100 * contamination),
      p25_contamination_percent = qtl(100 * contamination, 0.25),
      p75_contamination_percent = qtl(100 * contamination, 0.75),
      p90_contamination_percent = qtl(100 * contamination, 0.90),
      p95_contamination_percent = qtl(100 * contamination, 0.95),
      n_cells_contamination_ge_20pct = sum(contamination >= 0.2),
      n_cells_contamination_ge_50pct = sum(contamination >= 0.5),
      n_heuristic_clusters = length(unique(provisional_clusters[!is.na(provisional_clusters)])),
      n_gene_features = nrow(raw),
      n_nonzero_features_modeled = length(expressed),
      stringsAsFactors = FALSE
    )
    markers[[i]] <- marker_audit(gsm, stage, raw, corrected, feature_md)
    log_msg("Completed ", gsm,
            ": median contamination=", round(med(100 * contamination), 2),
            "%; total estimated UMI retained=",
            round(100 * sum(corrected_umi) / sum(raw_umi), 2), "%")

    rm(input, raw, model_input, corrected_fit, corrected, result, saved,
       md, feature_md, updated_md)
    invisible(gc(verbose = FALSE))
  }

  man <- do.call(rbind, manifest)
  cells <- do.call(rbind, per_cell)
  sm <- do.call(rbind, sample_stats)
  mark <- do.call(rbind, markers)
  write_table(man, "01_118_input_sample_audit_v01.csv")
  if (SAVE_PER_CELL) write_table(cells, "02_118_per_cell_contamination_v01.csv")
  write_table(sm, "03_118_per_sample_DecontX_summary_v01.csv")
  write_table(mark, "05_118_marker_detection_before_after_v01.csv")

  st <- do.call(rbind, lapply(split(sm, sm$dpi), function(x) {
    xx <- cells[cells$dpi == x$dpi[[1L]], , drop = FALSE]
    data.frame(
      dpi = x$dpi[[1L]], n_gsm_libraries = nrow(x),
      n_cells = sum(x$n_cells),
      median_cell_contamination_percent = med(xx$decontX_contamination_percent),
      original_total_UMI = sum(x$original_total_UMI),
      corrected_total_estimated_UMI = sum(x$corrected_total_estimated_UMI),
      corrected_over_raw_UMI_percent =
        100 * sum(x$corrected_total_estimated_UMI) / sum(x$original_total_UMI),
      metadata_conflicted_gsm = sum(grepl("PROVISIONAL", x$age_assignment_status)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(st) <- NULL
  write_table(st, "04_118_stage_DecontX_summary_v01.csv")
  if (SAVE_PLOTS) {
    tryCatch(plot_overall(sm), error = function(e) {
      log_msg("WARNING: overall figure failed: ", conditionMessage(e))
    })
  }
writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT_DIR, "06_118_REPORT_v01.txt"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()),
             file.path(OUT_DIR, "07_118_sessionInfo_v01.txt"), useBytes = TRUE)
  files <- c(
    "00_118_run_log_v01.log", "01_118_input_sample_audit_v01.csv",
    if (SAVE_PER_CELL) "02_118_per_cell_contamination_v01.csv",
    "03_118_per_sample_DecontX_summary_v01.csv",
    "04_118_stage_DecontX_summary_v01.csv",
    "05_118_marker_detection_before_after_v01.csv",
    "06_118_REPORT_v01.txt", "07_118_sessionInfo_v01.txt",
    file.path("RDS", basename(output_rds)),
    file.path("logs", list.files(MODEL_LOG_DIR, full.names = FALSE,
                                 pattern = "^118_GSM.*_DecontX_model_v01\\.log$")),
    if (SAVE_PLOTS) file.path("figures", list.files(FIG_DIR,
           full.names = FALSE, pattern = "^Fig_118_.*_v01\\.png$"))
  )
  idx <- data.frame(relative_path = files, file_exists = file.exists(file.path(OUT_DIR, files)))
  write_table(idx, "08_118_OUTPUT_INDEX_v01.csv")
  invisible(sm)
}

tryCatch(main(), error = function(e) {
  log_msg("\nERROR: ", conditionMessage(e))
  stop(e)
})

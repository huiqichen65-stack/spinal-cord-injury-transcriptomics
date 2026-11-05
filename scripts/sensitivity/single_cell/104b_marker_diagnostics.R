source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)
set.seed(1234)
PROJECT_DIR <- PROJECT_ROOT
VERSION <- "v04"
INPUT_RDS <- file.path(PROJECT_DIR, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
PREV_DIR <- file.path(PROJECT_DIR, "results", "104_Harmony_source_sensitivity_v03")
OUTPUT_DIR <- file.path(PROJECT_DIR, "results", "104_marker_diagnostics_v04")
FOCUS_LABEL <- "Clec7a+ lipid-associated Microglia"
SOURCES <- c("17702045", "GSE162610", "GSE196928")
MAIN_STAGES <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")
CORE <- c("Clec7a", "Lpl", "Apoe", "Trem2")
PANEL <- c(CORE, "Cst7", "Apoc4")
WRITE_CELL_LEVEL_GZ <- TRUE
MAKE_ZIP <- TRUE
ZIP_MAX_MB <- 30

packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})
if (!file.exists(INPUT_RDS)) stop("Missing input: ", INPUT_RDS)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUTPUT_DIR, "104_marker_diagnostics_v04.log")
if (file.exists(LOG_FILE)) unlink(LOG_FILE)
logit <- function(...) {
  message <- paste0(..., collapse = "")
  cat(message, "\n")
  cat(message, "\n", file = LOG_FILE, append = TRUE)
}
save_csv <- function(x, prefix) {
  path <- file.path(OUTPUT_DIR, paste0(prefix, "_v04.csv"))
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  logit("Saved ", basename(path), " | rows=", nrow(x))
  invisible(path)
}
pct <- function(n, d) ifelse(d > 0, 100 * n/d, NA_real_)
finite_mean <- function(v) if (any(is.finite(v))) mean(v[is.finite(v)]) else NA_real_
finite_median <- function(v) if (any(is.finite(v))) stats::median(v[is.finite(v)]) else NA_real_

logit("104 marker diagnostics v04 | ", as.character(Sys.time()))
logit("Input: ", INPUT_RDS)
logit("Output: ", OUTPUT_DIR)
logit("R ", as.character(getRversion()), "; Seurat ", as.character(utils::packageVersion("Seurat")))

obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) stop("Step 18 file is not a Seurat object")
if (!"RNA" %in% SeuratObject::Assays(obj)) stop("RNA assay missing")
SeuratObject::DefaultAssay(obj) <- "RNA"
md <- obj[[]]
required <- c("source_dataset", "timepoint_main", "myeloid_state")
if (length(setdiff(required, names(md)))) {
  stop("Missing required metadata fields: ", paste(setdiff(required, names(md)), collapse = ", "))
}
md$cell_id <- rownames(md)
md$source_dataset <- as.character(md$source_dataset)
md$timepoint_main <- as.character(md$timepoint_main)
md$myeloid_state <- as.character(md$myeloid_state)
if (anyNA(md[, required, drop = FALSE])) stop("Missing essential metadata values")
if (!all(md$source_dataset %in% SOURCES)) stop("Unexpected source labels")
if (!FOCUS_LABEL %in% md$myeloid_state) stop("Missing label: ", FOCUS_LABEL)
md$frozen_focus <- md$myeloid_state == FOCUS_LABEL
md$frozen_state <- ifelse(md$frozen_focus, "Frozen_Clec7a+", "Other_myeloid")
md$analysis_stage <- ifelse(md$source_dataset == "GSE196928" & md$timepoint_main == "Control", "Sham",
  md$timepoint_main)
md$is_main_stage <- md$timepoint_main %in% MAIN_STAGES
if ("treatment" %in% names(md)) {
  tx <- unique(as.character(md$treatment[md$source_dataset == "GSE196928"]))
  logit("GSE196928 recorded treatment: ", paste(tx, collapse = "; "))
  if (any(grepl("PLX5622", tx, ignore.case = TRUE), na.rm = TRUE)) {
    stop("PLX5622 rows detected; stop and check Step 18 atlas")
  }
}

obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) stop("JoinLayers failed: ",
  conditionMessage(e), "; inspect source-specific RNA assay layers"))
counts <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) stop("Cannot extract joined RNA counts: ",
  conditionMessage(e)))
if (!setequal(colnames(counts), md$cell_id) || ncol(counts) != nrow(md)) {
  stop("RNA counts columns do not exactly match metadata cell IDs")
}
counts <- counts[, md$cell_id, drop = FALSE]
if (!inherits(counts, "sparseMatrix") && !is.matrix(counts)) {
  stop("Unexpected joined RNA count class: ", paste(class(counts), collapse = "/"))
}
if ("x" %in% methods::slotNames(counts)) {
  bad_entries <- sum(!is.finite(methods::slot(counts, "x")))
} else {
  bad_entries <- sum(!is.finite(counts))
}
if (bad_entries != 0) stop("Non-finite matrix entries: ", bad_entries)
missing_core <- setdiff(CORE, rownames(counts))
if (length(missing_core)) stop("Core gene features absent from RNA: ", paste(missing_core, collapse = ", "))
available <- PANEL[PANEL %in% rownames(counts)]
missing_optional <- setdiff(PANEL, available)
logit("Detected marker features: ", paste(available, collapse = "; "))
if (length(missing_optional)) logit("Optional features not present: ", paste(missing_optional, collapse = "; "))

md$n_umi_joined <- as.numeric(Matrix::colSums(counts))
md$count_evaluable <- is.finite(md$n_umi_joined) & md$n_umi_joined > 0
if ("nCount_RNA" %in% names(md)) {
  nc <- suppressWarnings(as.numeric(as.character(md$nCount_RNA)))
  count_mismatch <- is.finite(nc) & abs(nc - md$n_umi_joined) > pmax(1, abs(nc) * 0.001)
} else {
  count_mismatch <- rep(FALSE, nrow(md))
}
source_count_audit <- do.call(rbind, lapply(split(seq_len(nrow(md)), md$source_dataset), function(ii) {
  data.frame(source_dataset = md$source_dataset[ii[1]], n_frozen_cells = length(ii), n_expression_evaluable = sum(md$count_evaluable[ii]),
    n_zero_count = sum(!md$count_evaluable[ii]), pct_zero_count = pct(sum(!md$count_evaluable[ii]),
      length(ii)), n_count_metadata_mismatch = sum(count_mismatch[ii]), median_umi_evaluable = finite_median(md$n_umi_joined[ii][md$count_evaluable[ii]]))
}))
rownames(source_count_audit) <- NULL
save_csv(source_count_audit, "00_104_source_RNA_count_audit")
if (any(md$count_evaluable & count_mismatch) && mean(count_mismatch) > 0.01) {
  stop("More than 1% RNA counts mismatch nCount_RNA; see count audit")
}
if (any(!md$count_evaluable & "nCount_RNA" %in% names(md))) {
  if ("nCount_RNA" %in% names(md)) {
    nc <- suppressWarnings(as.numeric(as.character(md$nCount_RNA)))
    if (any(!md$count_evaluable & is.finite(nc) & nc > 0)) {
      stop("Zero joined counts but positive nCount_RNA: possibly incomplete layers")
    }
  }
}
if (any(source_count_audit$pct_zero_count > 2)) {
  stop("More than 2% zero-count cells within one dataset; inspect joined RNA layer")
}

raw <- as.matrix(counts[available, , drop = FALSE])
rm(counts, obj)
gc(verbose = FALSE)
colnames(raw) <- md$cell_id
logexpr <- matrix(NA_real_, nrow = length(available), ncol = nrow(md), dimnames = list(available, md$cell_id))
valid <- which(md$count_evaluable)
logexpr[, valid] <- log1p(sweep(raw[, valid, drop = FALSE], 2, md$n_umi_joined[valid]/10000, "/"))
positives <- raw > 0
colnames(positives) <- md$cell_id

md$det_Clec7a <- positives["Clec7a", ] & md$count_evaluable
md$det_Lpl <- positives["Lpl", ] & md$count_evaluable
md$det_Apoe <- positives["Apoe", ] & md$count_evaluable
md$det_Trem2 <- positives["Trem2", ] & md$count_evaluable
md$det_Clec7a_Lpl <- md$det_Clec7a & md$det_Lpl
md$det_strict <- md$det_Clec7a_Lpl & (md$det_Apoe | md$det_Trem2)
md$n_core_detected <- colSums(positives[CORE, , drop = FALSE])
md$det_two_of_four <- md$n_core_detected >= 2 & md$count_evaluable

md$score_four <- rep(NA_real_, nrow(md))
md$score_four[valid] <- colMeans(logexpr[CORE, valid, drop = FALSE])
md$score_six_available <- rep(NA_real_, nrow(md))
md$score_six_available[valid] <- colMeans(logexpr[available, valid, drop = FALSE])
source_pct_rank <- function(score, source, eligible) {
  z <- rep(NA_real_, length(score))
  for (s in unique(source)) {
    ix <- which(source == s & eligible & is.finite(score))
    if (length(ix))
      z[ix] <- (rank(score[ix], ties.method = "average") - 0.5)/length(ix)
  }
  z
}
md$source_percentile_four <- source_pct_rank(md$score_four, md$source_dataset, md$count_evaluable)
md$source_top10_four <- is.finite(md$source_percentile_four) & md$source_percentile_four >= 0.9

marker_rows <- function(group_vars, id) {
  ix <- split(seq_len(nrow(md)), interaction(md[, group_vars, drop = FALSE], drop = TRUE, lex.order = TRUE))
  out <- vector("list", length(ix) * length(available))
  k <- 0L
  for (ii in ix) {
    base <- as.list(md[ii[1], group_vars, drop = FALSE])
    usable <- ii[md$count_evaluable[ii]]
    for (gene in available) {
      k <- k + 1L
      det <- if (length(usable))
        sum(positives[gene, usable]) else 0L
      ev <- logexpr[gene, usable]
      ans <- c(base, list(gene = gene, n_all_cells = length(ii), n_expression_evaluable = length(usable),
        n_detected = det, pct_detected = pct(det, length(usable)), mean_log1p_10k_all = finite_mean(ev),
        median_log1p_10k_all = finite_median(ev), mean_log1p_10k_detected = if (length(usable) &&
          any(positives[gene, usable])) {
          finite_mean(ev[positives[gene, usable]])
        } else NA_real_, median_umi = finite_median(md$n_umi_joined[usable])))
      out[[k]] <- as.data.frame(ans, stringsAsFactors = FALSE)
    }
  }
  if (!k)
    return(data.frame())
  z <- do.call(rbind, out[seq_len(k)])
  z$table_id <- id
  rownames(z) <- NULL
  z
}
by_source_gene <- marker_rows("source_dataset", "source")
by_source_stage_gene <- marker_rows(c("source_dataset", "timepoint_main"), "source_stage")
by_source_state_gene <- marker_rows(c("source_dataset", "frozen_state"), "source_frozen_state")
by_source_stage_state_gene <- marker_rows(c("source_dataset", "timepoint_main", "frozen_state"), "source_stage_frozen_state")
save_csv(by_source_gene, "01_104_marker_by_source")
save_csv(by_source_stage_gene, "02_104_marker_by_source_stage")
save_csv(by_source_state_gene, "03_104_marker_by_source_state")
save_csv(by_source_stage_state_gene, "04_104_marker_by_source_stage_state")

cascade <- function(keys) {
  ix <- split(seq_len(nrow(md)), interaction(md[, keys, drop = FALSE], drop = TRUE, lex.order = TRUE))
  out <- lapply(ix, function(ii) {
    usable <- ii[md$count_evaluable[ii]]
    n <- length(usable)
    data.frame(md[ii[1], keys, drop = FALSE], n_frozen_cells = length(ii), n_expression_evaluable = n,
      n_frozen_focus = sum(md$frozen_focus[ii]), n_Clec7a = sum(md$det_Clec7a[usable]), n_Lpl = sum(md$det_Lpl[usable]),
      n_Clec7a_Lpl = sum(md$det_Clec7a_Lpl[usable]), n_Clec7a_Lpl_Apoe_or_Trem2 = sum(md$det_strict[usable]),
      n_any_two_of_four = sum(md$det_two_of_four[usable]), pct_Clec7a = pct(sum(md$det_Clec7a[usable]),
        n), pct_Lpl = pct(sum(md$det_Lpl[usable]), n), pct_Clec7a_Lpl = pct(sum(md$det_Clec7a_Lpl[usable]),
        n), pct_strict = pct(sum(md$det_strict[usable]), n), pct_two_of_four = pct(sum(md$det_two_of_four[usable]),
        n), median_umi = finite_median(md$n_umi_joined[usable]), median_four_score = finite_median(md$score_four[usable]),
      median_source_percentile = finite_median(md$source_percentile_four[usable]), pct_source_top_decile = pct(sum(md$source_top10_four[usable]),
        n), stringsAsFactors = FALSE)
  })
  z <- do.call(rbind, out)
  rownames(z) <- NULL
  z
}
cascade_source <- cascade("source_dataset")
cascade_source_stage <- cascade(c("source_dataset", "timepoint_main"))
cascade_source_stage_state <- cascade(c("source_dataset", "timepoint_main", "frozen_state"))
save_csv(cascade_source, "05_104_codetection_cascade_source")
save_csv(cascade_source_stage, "06_104_codetection_cascade_source_stage")
save_csv(cascade_source_stage_state, "07_104_codetection_cascade_source_stage_state")

cutoffs <- as.numeric(stats::quantile(md$n_umi_joined[md$count_evaluable], probs = c(0, 0.25, 0.5, 0.75,
  1)))
if (length(unique(cutoffs)) < 5) {

  md$global_depth_bin <- NA_character_
  idx <- which(md$count_evaluable)
  md$global_depth_bin[idx] <- paste0("Q", pmin(4L, ceiling(4 * rank(md$n_umi_joined[idx], ties.method = "average")/length(idx))))
} else {
  md$global_depth_bin <- as.character(cut(md$n_umi_joined, breaks = cutoffs, include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4")))
}
depth_groups <- split(which(md$count_evaluable), interaction(md$source_dataset[md$count_evaluable], md$global_depth_bin[md$count_evaluable],
  drop = TRUE, lex.order = TRUE))
depth_table <- do.call(rbind, lapply(depth_groups, function(ii) {
  data.frame(source_dataset = md$source_dataset[ii[1]], global_depth_bin = md$global_depth_bin[ii[1]],
    n_cells = length(ii), median_umi = finite_median(md$n_umi_joined[ii]), pct_Clec7a = pct(sum(md$det_Clec7a[ii]),
      length(ii)), pct_Lpl = pct(sum(md$det_Lpl[ii]), length(ii)), pct_Apoe = pct(sum(md$det_Apoe[ii]),
      length(ii)), pct_Trem2 = pct(sum(md$det_Trem2[ii]), length(ii)), pct_strict = pct(sum(md$det_strict[ii]),
      length(ii)))
}))
rownames(depth_table) <- NULL
save_csv(depth_table, "08_104_depth_stratified_detection")

stage_summary <- cascade_source_stage
stage_summary$is_main_stage <- stage_summary$timepoint_main %in% MAIN_STAGES
stage_summary$timepoint_label <- ifelse(stage_summary$source_dataset == "GSE196928" & stage_summary$timepoint_main ==
  "Control", "Sham (Vehicle)", stage_summary$timepoint_main)
stage_summary$unit_of_inference <- "Cell-level descriptive, not independent animals"
save_csv(stage_summary, "09_104_source_stage_score_summary")

loo_cases <- list(All_sources = SOURCES, Drop_17702045 = setdiff(SOURCES, "17702045"), Drop_GSE162610 = setdiff(SOURCES,
  "GSE162610"), Drop_GSE196928 = setdiff(SOURCES, "GSE196928"))
loo_rows <- list()
for (nm in names(loo_cases)) {
  for (st in MAIN_STAGES) {
    ii <- which(md$source_dataset %in% loo_cases[[nm]] & md$timepoint_main == st)
    usable <- ii[md$count_evaluable[ii]]
    n <- length(ii)
    m <- length(usable)
    loo_rows[[length(loo_rows) + 1L]] <- data.frame(case = nm, source_included = paste(loo_cases[[nm]],
      collapse = ";"), timepoint_main = st, n_myeloid_cells = n, n_expression_evaluable = m, n_frozen_focus = sum(md$frozen_focus[ii]),
      pct_frozen_focus = pct(sum(md$frozen_focus[ii]), n), pct_strict_codetection = pct(sum(md$det_strict[usable]),
        m), pct_two_of_four = pct(sum(md$det_two_of_four[usable]), m), median_four_gene_score = finite_median(md$score_four[usable]),
      median_source_percentile = finite_median(md$source_percentile_four[usable]), pct_source_top_decile = pct(sum(md$source_top10_four[usable]),
        m), interpretation = if (n == 0)
        "NOT ESTIMABLE: missing stage after exclusion" else "-label composition / fixed RNA score; no reclustering", stringsAsFactors = FALSE)
  }
}
loo_table <- do.call(rbind, loo_rows)
save_csv(loo_table, "10_104_marker_score_leave_one_out")

why_rows <- list()
for (source in SOURCES) {
  a <- by_source_gene[by_source_gene$source_dataset == source, , drop = FALSE]
  a <- a[match(CORE, a$gene), , drop = FALSE]
  whole <- cascade_source[cascade_source$source_dataset == source, , drop = FALSE]
  zero_gene <- a$gene[a$n_detected == 0]
  weak_gene <- a$gene[a$pct_detected < 1 & a$n_detected > 0]
  cause <- if (length(zero_gene)) {
    paste0("At least one core marker has zero detected cells: ", paste(zero_gene, collapse = ", "),
      ". This prevents some co-detection combinations; may reflect upstream feature coverage.")
  } else if (whole$n_Clec7a_Lpl_Apoe_or_Trem2 == 0) {
    "All core genes detected separately, but strict co-detection is zero; inspect pairwise co-detection and depth."
  } else {
    "Strict co-detection is nonzero; inspect stage/state differences and source-specific rates."
  }
  why_rows[[length(why_rows) + 1L]] <- data.frame(source_dataset = source, n_evaluable_cells = whole$n_expression_evaluable,
    pct_strict_co_detection = whole$pct_strict, zero_detected_core_genes = paste(zero_gene, collapse = "; "),
    very_low_detection_core_genes_lt1pct = paste(weak_gene, collapse = "; "), description = cause,
    caveat = "RNA count zeros cannot distinguish biological absence from dropout, assay coverage, or upstream processing",
    stringsAsFactors = FALSE)
}
why <- do.call(rbind, why_rows)
save_csv(why, "11_104_strict_marker_failure_diagnosis")

v03_cluster_file <- file.path(PREV_DIR, "04_104_recluster_cluster_marker_audit_v03.csv")
if (file.exists(v03_cluster_file)) {
  old <- utils::read.csv(v03_cluster_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_v03 <- c("case", "arm", "marker_rich_proxy", "cluster")
  if (all(required_v03 %in% names(old))) {
    old$marker_rich_bool <- as.character(old$marker_rich_proxy) %in% c("TRUE", "True", "1")
    k <- split(seq_len(nrow(old)), interaction(old$case, old$arm, drop = TRUE))
    bridge <- do.call(rbind, lapply(k, function(ii) {
      data.frame(case = old$case[ii[1]], arm = old$arm[ii[1]], n_clusters_v03 = length(ii), n_old_proxy_positive_clusters = sum(old$marker_rich_bool[ii]),
        old_rule = "Mean 6-gene score AND strict raw-count co-detection; global cutoffs", what_v04_can_do = "Diagnose gene-specific source coverage; NOT reassign v03 clusters",
        required_for_de_novo_cluster_validation = "Per-cell cluster ID per v03 branch, or rerun v03 clustering with saved membership")
    }))
    rownames(bridge) <- NULL
    save_csv(bridge, "12_104_v03_cluster_proxy_caveat")
  } else {
    logit("v03 cluster audit exists but expected columns differ; skipping crosswalk")
  }
} else {
  logit("v03 cluster table not found; source-specific diagnostic remains complete")
}

if (WRITE_CELL_LEVEL_GZ) {
  cols <- c("cell_id", "source_dataset", "timepoint_main", "frozen_state", "n_umi_joined", "count_evaluable",
    "det_Clec7a", "det_Lpl", "det_Apoe", "det_Trem2", "det_Clec7a_Lpl", "det_strict", "n_core_detected",
    "det_two_of_four", "score_four", "score_six_available", "source_percentile_four", "source_top10_four")
  path <- file.path(OUTPUT_DIR, "13_104_cell_marker_scores_v04.csv.gz")
  cn <- gzfile(path, "wt", compression = 6)
  utils::write.csv(md[, cols, drop = FALSE], cn, row.names = FALSE, na = "")
  close(cn)
  logit("Saved ", basename(path), " | cells=", nrow(md))
}

pdf_path <- file.path(OUTPUT_DIR, "14_104_marker_diagnostic_figures_v04.pdf")
tryCatch({
    grDevices::pdf(pdf_path, width = 10.5, height = 7.4, onefile = TRUE)
    plot_dat <- by_source_state_gene
    plot_dat$source_dataset <- factor(plot_dat$source_dataset, levels = SOURCES)
    plot_dat$gene <- factor(plot_dat$gene, levels = rev(available))
    p1 <- ggplot(plot_dat, aes(x = source_dataset, y = gene, fill = pct_detected)) + geom_tile(color = "white", 
        linewidth = 0.5) + geom_text(aes(label = sprintf("%.1f%%", pct_detected)), size = 3) + facet_wrap(~frozen_state, 
        nrow = 1) + scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) + labs(title = "Marker detection by source and original state", 
        subtitle = "Raw RNA counts >0; denominators exclude audited zero-UMI cells", x = "Dataset", y = "Marker", 
        fill = "Detection %") + theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 20, 
        hjust = 1))
    print(p1)
    plot_dat2 <- by_source_stage_gene[by_source_stage_gene$gene %in% CORE & by_source_stage_gene$timepoint_main %in% 
        MAIN_STAGES, ]
    plot_dat2$timepoint_main <- factor(plot_dat2$timepoint_main, levels = MAIN_STAGES)
    p2 <- ggplot(plot_dat2, aes(x = timepoint_main, y = pct_detected, group = gene, color = gene)) + 
        geom_line(linewidth = 0.6, na.rm = TRUE) + geom_point(size = 1.8, na.rm = TRUE) + facet_wrap(~source_dataset, 
        nrow = 3, scales = "free_y") + labs(title = "Individual core-marker detection by source and stage", 
        subtitle = "Frequencies are cell-level descriptive (not animal-level replication)", x = "Stage (GSE196928 Control = Sham)", 
        y = "Detected cells (%)") + theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 40, 
        hjust = 1))
    print(p2)
    ss <- cascade_source_stage
    ss$timepoint_main <- factor(ss$timepoint_main, levels = c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", 
        "14 dpi", "28 dpi", "38 dpi"))
    p3 <- ggplot(ss, aes(x = timepoint_main, y = pct_strict, group = 1)) + geom_line(color = "#C25843", 
        linewidth = 0.7, na.rm = TRUE) + geom_point(color = "#C25843", size = 2, na.rm = TRUE) + facet_wrap(~source_dataset, 
        ncol = 1, scales = "free_x") + coord_cartesian(ylim = c(0, 100)) + labs(title = "Strict co-detection: Clec7a & Lpl & (Apoe or Trem2)", 
        subtitle = "Absence of strict co-detection is NOT evidence of state absence", x = "Stage (GSE196928 Control = Sham)", 
        y = "Cells satisfying rule (%)") + theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 40, 
        hjust = 1))
    print(p3)
    s4 <- cascade_source_stage_state
    s4$timepoint_main <- factor(s4$timepoint_main, levels = c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", 
        "14 dpi", "28 dpi", "38 dpi"))
    p4 <- ggplot(s4, aes(x = timepoint_main, y = median_source_percentile, group = frozen_state, color = frozen_state)) + 
        geom_line(linewidth = 0.7, na.rm = TRUE) + geom_point(size = 1.9, na.rm = TRUE) + facet_wrap(~source_dataset, 
        ncol = 1) + coord_cartesian(ylim = c(0, 1)) + labs(title = "Within-source percentile of fixed four-marker expression score", 
        subtitle = "Score is not a re-clustered identity call or biological replicate statistic", x = "Stage", 
        y = "Median source-relative percentile", color = "Original label") + theme_bw(base_size = 10) + 
        theme(axis.text.x = element_text(angle = 38, hjust = 1))
    print(p4)
    grDevices::dev.off()
    logit("Saved ", basename(pdf_path))
}, error = function(e) {
    if (grDevices::dev.cur() > 1) 
        try(grDevices::dev.off(), silent = TRUE)
    logit("WARNING: PDF figure creation failed: ", conditionMessage(e))
})

report_path <- file.path(OUTPUT_DIR, "16_104_REPORT_v04.txt")
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_path, useBytes = TRUE)
session_path <- file.path(OUTPUT_DIR, "15_104_sessionInfo_v04.txt")
writeLines(capture.output(utils::sessionInfo()), session_path, useBytes = TRUE)
logit("Saved ", basename(report_path), " and ", basename(session_path))

if (MAKE_ZIP) {
  if (requireNamespace("zip", quietly = TRUE)) {
    zip_name <- file.path(OUTPUT_DIR, "104_marker_diagnostics_results_v04.zip")
    candidates <- list.files(OUTPUT_DIR, full.names = TRUE)
    candidates <- candidates[!grepl("\\.zip$", basename(candidates), ignore.case = TRUE)]
    zip::zipr(zipfile = zip_name, files = basename(candidates), root = OUTPUT_DIR, include_directories = FALSE)
    mb <- file.info(zip_name)$size/1024^2
    if (!is.finite(mb) || mb > ZIP_MAX_MB) {
      unlink(zip_name)
      logit("ZIP exceeded ", ZIP_MAX_MB, " MB and was removed. All original results remain in: ",
        OUTPUT_DIR)
    } else {
      logit("ZIP ready: ", zip_name, " | ", round(mb, 2), " MB")
    }
  } else {
    logit("Optional zip package missing. Install with install.packages('zip')")
    logit("Results are saved uncompressed at ", OUTPUT_DIR)
  }
}
logit("104 marker diagnostics v04 finished successfully")

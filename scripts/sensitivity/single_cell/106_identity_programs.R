options(stringsAsFactors = FALSE, warn = 1)
ROOT <- "D:/zz"
INPUT <- file.path(ROOT, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
OUT <- file.path(ROOT, "results", "106_identity_programs")
STAGES <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")
SOURCES <- c("17702045", "GSE162610", "GSE196928")
FOCUS <- "Clec7a+ lipid-associated Microglia"
MIN_POSITIVE_CELLS_PER_SOURCE <- 10L
Z_CAP <- 3
MIN_PROGRAM_GENES <- 3L

PROGRAMS <- list(Microglia = c("P2ry12", "Tmem119", "Sall1", "Fcrls", "Siglech", "Gpr34", "Olfml3", "P2ry13"),
  Monocyte_derived_macrophage = c("Ccr2", "Ly6c2", "Plac8", "S100a8", "S100a9", "Cd14", "Mmp9"), Border_associated_macrophage = c("Mrc1",
    "Lyve1", "Cd163", "Folr2", "Pf4", "F13a1", "Stab1", "Timd4"))

pkgs <- c("Seurat", "SeuratObject", "Matrix", "ggplot2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Required R packages not installed: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})
if (!file.exists(INPUT)) stop("Input file not found: ", INPUT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "106_identity_programs_v01.log")
if (file.exists(LOG)) unlink(LOG)
log_msg <- function(...) {
  x <- paste0(..., collapse = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", x, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", x, "\n", sep = "", file = LOG, append = TRUE)
}
write_csv <- function(x, name) utils::write.csv(x, file.path(OUT, name), row.names = FALSE, na = "",
  fileEncoding = "UTF-8")
write_csv_gz <- function(x, name) {
  con <- gzfile(file.path(OUT, name), open = "wt", compression = 8)
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "")
}

o <- readRDS(INPUT)
if (!inherits(o, "Seurat")) stop("Input is not a Seurat object")
if (!"RNA" %in% Assays(o)) stop("RNA assay missing from input")
md <- o[[]]
needed <- c("source_dataset", "timepoint_main", "myeloid_state")
if (!all(needed %in% colnames(md))) {
  stop("Missing metadata fields: ", paste(setdiff(needed, colnames(md)), collapse = ", "))
}
md$cell_id <- rownames(md)
if (anyDuplicated(md$cell_id)) stop("Duplicated Seurat cell IDs")
md$source_dataset <- as.character(md$source_dataset)
md$timepoint_main <- as.character(md$timepoint_main)
md$myeloid_state <- as.character(md$myeloid_state)
if (!all(md$source_dataset %in% SOURCES)) {
  stop("Unexpected source_dataset labels: ", paste(setdiff(unique(md$source_dataset), SOURCES), collapse = ", "))
}
if (anyNA(md$myeloid_state) || anyNA(md$timepoint_main)) {
  stop("Missing stage or state labels: audit original Seurat metadata first")
}
if (!"sample_id" %in% names(md)) md$sample_id <- NA_character_
if (!"sample_id_source" %in% names(md)) md$sample_id_source <- NA_character_
if (!"celltype" %in% names(md)) md$celltype <- NA_character_
if (!"major_celltype" %in% names(md)) md$major_celltype <- NA_character_
md$sample_id <- as.character(md$sample_id)
md$sample_id_source <- as.character(md$sample_id_source)
md$celltype <- as.character(md$celltype)
md$major_celltype <- as.character(md$major_celltype)
md$sample_id[is.na(md$sample_id) | !nzchar(md$sample_id)] <- "Unrecoverable"
md$sample_id_source[is.na(md$sample_id_source) | !nzchar(md$sample_id_source)] <- "not_recoverable"
md$focus_group <- ifelse(md$myeloid_state == FOCUS, "Clec7a+ annotated focal state", ifelse(grepl("microglia",
  md$myeloid_state, ignore.case = TRUE), "Other microglia-labeled states", "Other myeloid states"))
md$sample_unit_interpretation <- ifelse(md$source_dataset == "17702045", "stage_group_not_biological_replicate",
  ifelse(md$source_dataset == "GSE196928", "pooled_library_not_individual_mouse", "sampling_unit_may_include_tissue_pool"))

DefaultAssay(o) <- "RNA"
o <- tryCatch(JoinLayers(o, assay = "RNA"), error = function(e) stop("Seurat v5 RNA JoinLayers failed: ",
  conditionMessage(e)))
cts <- tryCatch(GetAssayData(o, assay = "RNA", layer = "counts"), error = function(e) stop("RNA counts unavailable: ",
  conditionMessage(e)))
if (!inherits(cts, "sparseMatrix")) stop("RNA counts are not sparse; stop to avoid memory exhaustion")
if (!setequal(colnames(cts), md$cell_id)) stop("RNA count cell IDs do not match metadata")
cts <- cts[, md$cell_id, drop = FALSE]
if (ncol(cts) != nrow(md)) stop("Count matrix does not align to metadata")
rm(o)
invisible(gc(verbose = FALSE))

lib <- as.numeric(Matrix::colSums(cts))
md$rna_total_count_at_scoring <- lib
md$eligible_identity_scoring <- is.finite(lib) & lib > 0

all_genes <- unique(unlist(PROGRAMS, use.names = FALSE))
source_indices <- lapply(SOURCES, function(src) which(md$source_dataset == src & md$eligible_identity_scoring))
names(source_indices) <- SOURCES
manifest <- do.call(rbind, lapply(names(PROGRAMS), function(pr) {
  do.call(rbind, lapply(PROGRAMS[[pr]], function(g) {
    exists <- g %in% rownames(cts)
    do.call(rbind, lapply(SOURCES, function(src) {
      ii <- source_indices[[src]]
      detected <- if (exists)
        sum(cts[g, ii, drop = TRUE] > 0) else 0L
      data.frame(program = pr, gene = g, source_dataset = src, feature_in_merged_RNA = exists,
        n_source_cells = length(ii), n_cells_detected = detected, pct_source_cells_detected = if (length(ii))
          100 * detected/length(ii) else NA_real_, meets_detection_rule = exists && detected >= MIN_POSITIVE_CELLS_PER_SOURCE,
        stringsAsFactors = FALSE)
    }))
  }))
}))
shared_genes <- vapply(all_genes, function(g) {
  all(manifest$meets_detection_rule[manifest$gene == g])
}, logical(1))
shared_genes <- names(shared_genes)[shared_genes]
manifest$included_in_all_source_panel <- manifest$gene %in% shared_genes
manifest$interpretation_when_not_detected <- ifelse(manifest$n_cells_detected == 0L, "No detections in this source; technical absence cannot be excluded",
  ifelse(manifest$included_in_all_source_panel, "Included in fixed common panel", "Excluded from common panel: below detection threshold in >=1 source"))
write_csv(manifest, "01_106_identity_program_gene_manifest_v01.csv")

active <- lapply(PROGRAMS, function(gs) intersect(gs, shared_genes))
if (any(lengths(active) < MIN_PROGRAM_GENES)) {
  stop("Fewer than ", MIN_PROGRAM_GENES, " common detectable markers for a program. Review 01 gene manifest; no scores issued")
}

qc_source <- do.call(rbind, lapply(SOURCES, function(src) {
  ii <- which(md$source_dataset == src)
  data.frame(source_dataset = src, n_all_myeloid_cells = length(ii), n_scored_cells = sum(md$eligible_identity_scoring[ii]),
    n_zero_count_cells = sum(!md$eligible_identity_scoring[ii]), n_distinct_sample_id_strings = length(unique(md$sample_id[ii])),
    sample_id_provenance = paste(unique(md$sample_id_source[ii]), collapse = "; "), sample_unit_interpretation = unique(md$sample_unit_interpretation[ii])[1],
    stringsAsFactors = FALSE)
}))
write_csv(qc_source, "00_106_input_QC_v01.csv")

gene_mat <- as.matrix(cts[shared_genes, , drop = FALSE])
rm(cts)
invisible(gc(verbose = FALSE))
positive <- md$eligible_identity_scoring
log_cpk <- matrix(NA_real_, nrow = nrow(gene_mat), ncol = ncol(gene_mat), dimnames = dimnames(gene_mat))
log_cpk[, positive] <- log1p(sweep(gene_mat[, positive, drop = FALSE], 2, lib[positive], "/") * 10000)
rm(gene_mat)
invisible(gc(verbose = FALSE))
zexpr <- matrix(NA_real_, nrow = nrow(log_cpk), ncol = ncol(log_cpk), dimnames = dimnames(log_cpk))
for (src in SOURCES) {
  ii <- source_indices[[src]]
  for (g in shared_genes) {
    v <- log_cpk[g, ii]
    s <- stats::sd(v)
    if (is.finite(s) && s > 0) {
      zz <- (v - mean(v))/s
      zexpr[g, ii] <- pmax(-Z_CAP, pmin(Z_CAP, zz))
    }
  }
}
score_names <- c(Microglia = "MG_score", Monocyte_derived_macrophage = "MDM_score", Border_associated_macrophage = "BAM_score")
for (pr in names(active)) {
  genes <- active[[pr]]
  sn <- unname(score_names[[pr]])
  md[[sn]] <- colMeans(zexpr[genes, , drop = FALSE], na.rm = FALSE)
  md[[sn]][!positive] <- NA_real_
  if (any(!is.finite(md[[sn]][positive]))) {
    stop("Non-finite score for ", pr, "; inspect RNA input / gene audit")
  }
}

cell_export <- md[, c("cell_id", "source_dataset", "timepoint_main", "sample_id", "sample_id_source",
  "sample_unit_interpretation", "myeloid_state", "focus_group", "major_celltype", "celltype", "rna_total_count_at_scoring",
  "eligible_identity_scoring", "MG_score", "MDM_score", "BAM_score"), drop = FALSE]
write_csv_gz(cell_export, "02_106_per_cell_identity_scores_v01.csv.gz")

summarize_programs <- function(dat, group_vars) {
  if (!nrow(dat))
    return(data.frame())
  key <- do.call(paste, c(lapply(dat[group_vars], as.character), sep = "\037"))
  groups <- split(seq_len(nrow(dat)), key, drop = TRUE)
  out <- vector("list", length(groups) * length(score_names))
  j <- 0L
  for (ii in groups) {
    first <- dat[ii[1], group_vars, drop = FALSE]
    for (pr in names(score_names)) {
      sn <- unname(score_names[[pr]])
      values <- as.numeric(dat[[sn]][ii])
      good <- is.finite(values)
      x <- values[good]
      j <- j + 1L
      out[[j]] <- cbind(first, data.frame(program = pr, score_column = sn, n_all_cells = length(ii),
        n_scored_cells = length(x), mean_score = if (length(x))
          mean(x) else NA_real_, median_score = if (length(x))
          stats::median(x) else NA_real_, q25_score = if (length(x))
          unname(stats::quantile(x, 0.25)) else NA_real_, q75_score = if (length(x))
          unname(stats::quantile(x, 0.75)) else NA_real_, sd_score = if (length(x) > 1L)
          stats::sd(x) else NA_real_, n_program_genes = length(active[[pr]]), stringsAsFactors = FALSE))
    }
  }
  ans <- do.call(rbind, out[seq_len(j)])
  rownames(ans) <- NULL
  ans
}

by_source <- summarize_programs(cell_export, c("source_dataset", "myeloid_state", "focus_group"))
by_stage <- summarize_programs(cell_export, c("source_dataset", "timepoint_main", "myeloid_state", "focus_group"))
by_sample <- summarize_programs(cell_export, c("source_dataset", "sample_id", "sample_id_source", "sample_unit_interpretation",
  "timepoint_main", "myeloid_state", "focus_group"))
write_csv(by_source, "03_106_programs_by_source_state_v01.csv")
write_csv(by_stage, "04_106_programs_by_source_stage_state_v01.csv")
write_csv(by_sample, "05_106_programs_by_source_sample_stage_state_v01.csv")

marker_groups <- split(seq_len(nrow(cell_export)), do.call(paste, c(cell_export[c("source_dataset", "timepoint_main",
  "focus_group")], sep = "\037")), drop = TRUE)
marker_out <- vector("list", length(marker_groups) * length(shared_genes))
j <- 0L
for (ii in marker_groups) {
  d <- cell_export[ii[1], c("source_dataset", "timepoint_main", "focus_group"), drop = FALSE]
  ii <- ii[positive[ii]]
  for (g in shared_genes) {
    v <- log_cpk[g, ii]
    j <- j + 1L
    marker_out[[j]] <- cbind(d, data.frame(gene = g, program = names(active)[vapply(active, function(v) g %in%
      v, logical(1))][1], n_eligible_cells = length(ii), pct_cells_detected = if (length(v))
      mean(v > 0) * 100 else NA_real_, mean_log1p_CP10K = if (length(v))
      mean(v) else NA_real_, stringsAsFactors = FALSE))
  }
}
marker_stats <- do.call(rbind, marker_out[seq_len(j)])
rownames(marker_stats) <- NULL
write_csv(marker_stats, "06_106_marker_by_source_stage_focus_group_v01.csv")

program_levels <- names(active)
by_source$program <- factor(by_source$program, levels = program_levels)
by_stage$program <- factor(by_stage$program, levels = program_levels)
by_sample$program <- factor(by_sample$program, levels = program_levels)
state_order <- c(FOCUS, setdiff(sort(unique(md$myeloid_state)), FOCUS))
by_source$myeloid_state <- factor(by_source$myeloid_state, levels = rev(state_order))
by_stage$timepoint_main <- factor(by_stage$timepoint_main, levels = STAGES)
by_sample$timepoint_main <- factor(by_sample$timepoint_main, levels = STAGES)

pdf_path <- file.path(OUT, "07_106_identity_program_figures_v01.pdf")
graphics_ok <- TRUE
tryCatch({
  grDevices::pdf(pdf_path, width = 13.5, height = 7.5, onefile = TRUE, useDingbats = FALSE)
  p1 <- ggplot(by_source, aes(x = program, y = myeloid_state, fill = median_score)) + geom_tile(color = "white",
    linewidth = 0.35) + geom_text(aes(label = ifelse(n_scored_cells >= 1, sprintf("%.2f", median_score),
    "NA")), size = 2.6) + facet_wrap(~source_dataset, nrow = 1) + scale_fill_gradient2(low = "#4379AD",
    mid = "white", high = "#BB594D", midpoint = 0, name = "Median within-\nsource Z score") + labs(title = "Myeloid identity programs by original annotated state and source",
    subtitle = "Original Step18 labels retained; scores are standardized within each source, not lineage assignments",
    x = NULL, y = NULL) + theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30,
    hjust = 1), panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  print(p1)

  focal <- by_stage[as.character(by_stage$myeloid_state) == FOCUS & is.finite(by_stage$median_score),
    , drop = FALSE]
  if (nrow(focal)) {
    p2 <- ggplot(focal, aes(x = timepoint_main, y = median_score)) + geom_linerange(aes(ymin = q25_score,
      ymax = q75_score), linewidth = 0.65) + geom_point(aes(size = n_scored_cells), alpha = 0.9) +
      facet_grid(program ~ source_dataset, scales = "free_x") + labs(title = "Identity programs within the originally annotated Clec7a+ state",
      subtitle = "Point = cell-level median; vertical line = cell-level IQR (NOT biological-replicate uncertainty)",
      x = "Injury stage", y = "Source-relative identity score", size = "Cells") + theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
    print(p2)
  }

  smp <- by_sample[as.character(by_sample$myeloid_state) == FOCUS & by_sample$source_dataset != "17702045" &
    is.finite(by_sample$median_score), , drop = FALSE]
  if (nrow(smp)) {
    smp$source_dataset <- factor(smp$source_dataset, levels = c("GSE162610", "GSE196928"))
    smp$sample_id <- factor(smp$sample_id, levels = unique(smp$sample_id))
    p3 <- ggplot(smp, aes(x = sample_id, y = median_score, color = timepoint_main)) + geom_linerange(aes(ymin = q25_score,
      ymax = q75_score), linewidth = 0.5) + geom_point(aes(size = n_scored_cells), alpha = 0.9) +
      facet_grid(program ~ source_dataset, scales = "free_x", space = "free_x") + labs(title = "Focal-state identity programs by recoverable sample/library unit",
      subtitle = "GSE162610 units can include pools; GSE196928 B1/B3/B7 are pooled libraries; NOT animal replicates",
      x = "Reported sample or pooled library", y = "Source-relative identity score", color = "Stage",
      size = "Cells") + theme_bw(base_size = 8.5) + theme(axis.text.x = element_text(angle = 55,
      hjust = 1), panel.grid.minor = element_blank())
    print(p3)
  }
}, error = function(e) {
  graphics_ok <<- FALSE
  log_msg("Plotting failed (CSV summaries remain available): ", conditionMessage(e))
}, finally = {
  if (grDevices::dev.cur() > 1L)
    grDevices::dev.off()
})

notes <- c("106 | Source/stage/sample identity programs", paste("Input:", INPUT), paste("Output:", OUT),
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste("Cells:", nrow(md), "; score-eligible:",
    sum(positive), "; zero-RNA counts:", sum(!positive)), "", "PRESPECIFIED PROGRAMS (shared panel after source-level detection QC)",
  vapply(names(active), function(pr) paste0(pr, ": ", paste(active[[pr]], collapse = ", ")), character(1)),
  "", "METHOD:", "RNA counts from Step18 final myeloid Seurat object, joined across source layers.",
  "Normalize per cell using all RNA genes: log(1 + gene UMI / total RNA UMI * 10000).", "For each gene and source, Z-standardize over all nonzero-total-UMI myeloid cells; cap [-3,3].",
  "Identity-program score = per-cell mean of fixed, source-shared gene Z scores.", "A gene is included only when detected in >=10 cells in EVERY source; exclude all-source-unmatched genes.",
  "Report median, cell-level IQR, mean, SD and counts by source, stage, myeloid state, sample/library.",
  "Marker-level source/stage/focal-group detection/mean are additionally supplied.", "", "CRITICAL INTERPRETATION:",
  "Identity scores are source-relative within the fixed myeloid population, not cross-source calibrated.",
  "No identity reclassification, threshold-based exclusion or lineage inference is performed in Step106.",
  "Reduced P2ry12/Tmem119 after injury does not by itself establish monocyte origin.", "Shared injury/activation genes (Clec7a/Lpl/Apoe/Trem2) are NOT identity-defining here.",
  "Gene absence in original source data may masquerade as merged zero; detection screening is conservative.",
  "Within-cell summaries and score quantiles are NOT independent biological confidence intervals.",
  "The Step18 cell cohort and cell-state labels predate this audit; no independent reference mapping is claimed.",
  "Figshare 17702045 sample_id is stage-derived: do not count each stage label as an animal/library.",
  "GSE196928 B1/B3/B7 are pooled libraries, each corresponding to a condition/stage; sham != uninjured.",
  "GSE162610 sample IDs denote sampling units/libraries; some may contain pooled animals.", "GSE196928 PLX5622 arms were excluded previously; identity findings here concern retained Vehicle arms.",
  "The 14 dpi and 28 dpi groups are each represented by a single data source; temporal effects confounded.",
  "No statistical tests using cells as independent animal replicates are conducted.", "", "OUTPUT INDEX:",
  "00 input cell source, sample-provenance and zero-count audit", "01 gene identity program list, detection in each source, common panel inclusion",
  "02 compressed per-cell original identity and three identity-program scores", "03 per-source/per-original-state identity-program summaries",
  "04 per-source/per-stage/per-original-state identity-program summaries", "05 per-source/per-sample-or-pool/per-stage/per-original-state summaries",
  "06 source/stage/focal-vs-other gene-level detection and logRNA expression", "07 multi-page PDF figures (source/state, focus/stage, focus/sample)",
  "08 this report; 09 sessionInfo; run log", paste("Figure generation succeeded:", graphics_ok))
writeLines(enc2utf8(notes), file.path(OUT, "08_106_methods_and_limitations_v01.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUT, "09_106_sessionInfo_v01.txt"))

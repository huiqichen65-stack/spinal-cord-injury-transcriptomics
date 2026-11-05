source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)
SEED <- 1234L
N_HVG <- 3000L
N_PCS_COMPUTE <- 50L
N_PCS_COMPARE <- 30L
RESOLUTION <- 0.1
HARMONY_MAX_ITER <- 20L
MARKERS <- c("Clec7a", "Lpl", "Apoe", "Trem2", "Cst7", "Apoc4")
FOCUS_LABEL <- "Clec7a+ lipid-associated Microglia"
SOURCES <- c("17702045", "GSE162610", "GSE196928")
MAIN_STAGES <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")
FULL_STAGES <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")

RUN_SINGLE_SOURCE_RECLUSTER <- TRUE
RUN_UMAP_FOR_ALL_BRANCHES <- FALSE
RUN_UMAP_FOR_FULL <- TRUE
RUN_KNN_SOURCE_MIXING <- TRUE
MAX_KNN_CELLS_PER_SOURCE_STAGE <- 300L
MAX_PLOT_CELLS <- 12000L
SAVE_FULL_RDS <- FALSE
MAKE_SMALL_ZIP <- TRUE
ZIP_LIMIT_MB <- 100

MAX_ZERO_COUNT_FRACTION_PER_SOURCE <- 0.02

MAX_MISMATCH_FRACTION <- 0.01

set.seed(SEED)
req <- c("Seurat", "SeuratObject", "harmony", "Matrix", "ggplot2")
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

PROJECT_DIR <- if (dir.exists(PROJECT_ROOT)) PROJECT_ROOT else getwd()
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
INPUT <- file.path(PROJECT_DIR, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
OUTPUT <- file.path(PROJECT_DIR, "results", "104_Harmony_source_sensitivity_v03")
if (!file.exists(INPUT)) {
  stop("Cannot find input RDS: ", INPUT, "\nSet PROJECT_DIR above to the project root or the project root.")
}
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUTPUT, "104_Harmony_source_sensitivity_v03.log")
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)
logit <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}
write_csv <- function(x, nm) {
  fn <- file.path(OUTPUT, nm)
  utils::write.csv(x, fn, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  logit("Wrote: ", basename(fn), " (", nrow(x), " rows)")
  invisible(fn)
}
numfmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
logit("104 sensitivity analysis v03 | ", as.character(Sys.time()))
logit("PROJECT_DIR = ", PROJECT_DIR)
logit("INPUT = ", INPUT)
logit("R ", getRversion(), "; Seurat ", as.character(packageVersion("Seurat")), "; harmony ", as.character(packageVersion("harmony")))

obj0 <- readRDS(INPUT)
if (!inherits(obj0, "Seurat")) stop("Expected a Seurat object from Step 18")
if (!"RNA" %in% Assays(obj0)) stop("RNA assay missing from input")
DefaultAssay(obj0) <- "RNA"
md0 <- obj0[[]]
needed <- c("source_dataset", "timepoint_main", "myeloid_state")
if (length(setdiff(needed, colnames(md0)))) {
  stop("Missing metadata columns: ", paste(setdiff(needed, colnames(md0)), collapse = ", "))
}
md0$cell <- rownames(md0)
md0$source_dataset <- as.character(md0$source_dataset)
md0$timepoint_main <- as.character(md0$timepoint_main)
md0$myeloid_state <- as.character(md0$myeloid_state)
if (anyNA(md0$source_dataset) || any(!md0$source_dataset %in% SOURCES)) stop("Source dataset labels differ from 17702045 / GSE162610 / GSE196928")
if (anyNA(md0$timepoint_main) || any(!md0$timepoint_main %in% FULL_STAGES)) stop("Unexpected or missing timepoint_main values: inspect Step 18 metadata")
if (anyNA(md0$myeloid_state)) stop("Missing original myeloid_state labels")
if (!FOCUS_LABEL %in% md0$myeloid_state) stop("Expected focus label absent: ", FOCUS_LABEL)
if ("treatment" %in% colnames(md0)) {
  tx <- unique(trimws(as.character(md0$treatment[md0$source_dataset == "GSE196928"])))
  logit("GSE196928 treatment metadata: ", paste(tx, collapse = " / "))
  if (any(grepl("PLX5622", tx, ignore.case = TRUE)))
    stop("PLX5622 detected in final atlas! Verify upstream filtering.")
} else {
  logit("No treatment column in input; GSE196928 Vehicle origin is based on Step 05 selection")
}
if (nrow(md0) != 58533L) warning("Expected 58,533 final myeloid cells, observed ", nrow(md0), "; continuing but output must be reviewed before reporting")

md0$treatment_104 <- ifelse(md0$source_dataset == "GSE196928", "Vehicle (retained)", "Not separately reported in source metadata")
md0$stage_kind_104 <- ifelse(md0$source_dataset == "GSE196928" & md0$timepoint_main == "Control", "Sham; originally mapped to Control",
  as.character(md0$timepoint_main))
md0$is_main_stage <- md0$timepoint_main %in% MAIN_STAGES
md0$is_focus_frozen <- md0$myeloid_state == FOCUS_LABEL

obj0 <- tryCatch(JoinLayers(obj0, assay = "RNA"), error = function(e) {
  logit("JoinLayers returned error: ", conditionMessage(e))
  obj0
})
counts <- tryCatch(GetAssayData(obj0, assay = "RNA", layer = "counts"), error = function(e) stop("Joined RNA counts inaccessible: ",
  conditionMessage(e)))
if (ncol(counts) != nrow(md0)) stop("Counts/cell metadata dimensions disagree")
if (!setequal(colnames(counts), md0$cell)) stop("Counts cell barcodes mismatch")
counts <- counts[, md0$cell, drop = FALSE]
n_genes_original <- nrow(counts)
if (!all(c("Clec7a", "Lpl", "Apoe", "Trem2") %in% rownames(counts))) stop("Mandatory Clec7a/Lpl/Apoe/Trem2 markers missing from counts")
actual_markers <- intersect(MARKERS, rownames(counts))
logit("Available marker panel: ", paste(actual_markers, collapse = ", "))
logit("RNA counts class: ", paste(class(counts), collapse = "/"), "; genes=", nrow(counts), "; cells=",
  ncol(counts))

total_UMI <- Matrix::colSums(counts)
names(total_UMI) <- colnames(counts)
md0$joined_count_sum_104 <- as.numeric(total_UMI[md0$cell])
md0$invalid_total_count_104 <- !is.finite(md0$joined_count_sum_104) | md0$joined_count_sum_104 <= 0
md0$count_pass_104 <- !md0$invalid_total_count_104

nonfinite_entries <- if ("x" %in% methods::slotNames(counts)) {
  sum(!is.finite(methods::slot(counts, "x")))
} else if (is.matrix(counts)) {
  sum(!is.finite(counts))
} else {
  stop("Unsupported RNA count-matrix representation: ", paste(class(counts), collapse = "/"), ". Count integrity cannot be audited safely.")
}

if ("nCount_RNA" %in% names(md0)) {
  md0$stored_nCount_RNA_104 <- suppressWarnings(as.numeric(as.character(md0$nCount_RNA)))
  md0$nCount_discrepancy_104 <- is.finite(md0$stored_nCount_RNA_104) & is.finite(md0$joined_count_sum_104) &
    abs(md0$stored_nCount_RNA_104 - md0$joined_count_sum_104) > pmax(1, abs(md0$stored_nCount_RNA_104) *
      0.001)
} else {
  md0$stored_nCount_RNA_104 <- NA_real_
  md0$nCount_discrepancy_104 <- FALSE
}

qc_ix <- split(seq_len(nrow(md0)), interaction(md0$source_dataset, md0$timepoint_main, drop = TRUE, lex.order = TRUE))
qc_audit <- do.call(rbind, lapply(qc_ix, function(ii) data.frame(source_dataset = md0$source_dataset[ii[1]],
  timepoint_main = md0$timepoint_main[ii[1]], total_frozen_cells = length(ii), positive_count_cells = sum(md0$count_pass_104[ii]),
  zero_or_nonfinite_count_cells = sum(md0$invalid_total_count_104[ii]), zero_or_nonfinite_fraction = mean(md0$invalid_total_count_104[ii]),
  nCount_RNA_discrepancy_cells = sum(md0$nCount_discrepancy_104[ii]), median_joined_UMI_of_valid_cells = if (any(md0$count_pass_104[ii])) stats::median(md0$joined_count_sum_104[ii][md0$count_pass_104[ii]]) else NA_real_,
  stringsAsFactors = FALSE)))
rownames(qc_audit) <- NULL
write_csv(qc_audit, "00_104_RNA_count_QC_by_source_stage_v03.csv")

bad <- which(md0$invalid_total_count_104 | md0$nCount_discrepancy_104)
issue_details <- data.frame(cell = md0$cell[bad], source_dataset = md0$source_dataset[bad], timepoint_main = md0$timepoint_main[bad],
  myeloid_state = md0$myeloid_state[bad], stored_nCount_RNA = md0$stored_nCount_RNA_104[bad], joined_RNA_count_sum = md0$joined_count_sum_104[bad],
  zero_or_nonfinite = md0$invalid_total_count_104[bad], count_metadata_discrepancy = md0$nCount_discrepancy_104[bad],
  stringsAsFactors = FALSE)
write_csv(issue_details, "00_104_RNA_count_problem_cells_v03.csv")
source_ix <- split(seq_len(nrow(md0)), md0$source_dataset)
source_qc <- do.call(rbind, lapply(names(source_ix), function(src) {
  ii <- source_ix[[src]]
  data.frame(source_dataset = src, total_cells = length(ii), invalid_cells = sum(md0$invalid_total_count_104[ii]),
    invalid_fraction = mean(md0$invalid_total_count_104[ii]), metadata_discrepancy_cells = sum(md0$nCount_discrepancy_104[ii]),
    stringsAsFactors = FALSE)
}))
write_csv(source_qc, "00_104_RNA_count_QC_by_source_v03.csv")
logit("RNA count audit: ", sum(md0$invalid_total_count_104), " zero/nonfinite cells (", signif(mean(md0$invalid_total_count_104) *
  100, 3), "%), ", sum(md0$nCount_discrepancy_104), " count/metadata discrepancies")
for (i in seq_len(nrow(source_qc))) {
  logit("  ", source_qc$source_dataset[i], ": ", source_qc$invalid_cells[i], "/", source_qc$total_cells[i],
    " invalid, ", signif(source_qc$invalid_fraction[i] * 100, 3), "%")
}
if (is.finite(nonfinite_entries) && nonfinite_entries > 0) {
  stop("Nonfinite entries found inside joined RNA counts: ", nonfinite_entries, ". QC tables were saved in ",
    OUTPUT, ". Do not filter or impute these cells without investigating the layer source.")
}
if (any(source_qc$invalid_fraction > MAX_ZERO_COUNT_FRACTION_PER_SOURCE)) {
  stop("More than ", 100 * MAX_ZERO_COUNT_FRACTION_PER_SOURCE, "% zero-total-count cells within one source. Inspect 00_104_RNA_count_QC_by_source_v03.csv ",
    "and 00_104_RNA_count_problem_cells_v03.csv. ", "This may indicate incomplete count layers after JoinLayers; no cells were dropped.")
}
if (mean(md0$nCount_discrepancy_104) > MAX_MISMATCH_FRACTION || any(md0$invalid_total_count_104 & is.finite(md0$stored_nCount_RNA_104) &
  md0$stored_nCount_RNA_104 > 0)) {
  stop("Joined counts disagree with stored nCount_RNA for too many cells or ", "for a zero-count cell with originally positive UMI. ",
    "Check the saved 00_104_RNA_count_*.csv files before re-analysis.")
}
if (!any(md0$count_pass_104)) stop("No usable positive-count cells")

actual_valid <- which(md0$count_pass_104)
md0$marker_score_104 <- rep(NA_real_, nrow(md0))
md0$marker_codetect_104 <- rep(NA, nrow(md0))
valid_marker_counts <- as.matrix(counts[actual_markers, md0$cell[actual_valid], drop = FALSE])
valid_totals <- md0$joined_count_sum_104[actual_valid]
valid_expr <- log1p(sweep(valid_marker_counts, 2, valid_totals/10000, "/"))
md0$marker_score_104[actual_valid] <- colMeans(valid_expr)
md0$marker_codetect_104[actual_valid] <- (valid_marker_counts["Clec7a", ] > 0) & (valid_marker_counts["Lpl",
  ] > 0) & ((valid_marker_counts["Apoe", ] > 0) | (valid_marker_counts["Trem2", ] > 0))
rm(valid_expr, valid_marker_counts, valid_totals, counts, total_UMI)
gc(verbose = FALSE)

input_audit <- data.frame(category = c("input_file", "n_final_cells", "n_genes", "n_source_datasets",
  "n_frozen_focus_cells", "n_zero_count_cells_omitted_from_recluster", "available_markers", "reference_stage_label",
  "batch_correction", "same_final_cell_cohort", "independent_validation"), value = c(INPUT, as.character(nrow(md0)),
  as.character(n_genes_original), as.character(length(unique(md0$source_dataset))), as.character(sum(md0$is_focus_frozen)),
  as.character(sum(!md0$count_pass_104)), paste(actual_markers, collapse = ";"), "Control merges sham and uninjured origins; see Table S1/S2",
  "source_dataset only; never timepoint", "TRUE: Step 18 myeloid cohort", "FALSE"), stringsAsFactors = FALSE)
write_csv(input_audit, "00_104_INPUT_AUDIT_v03.csv")

summary_by <- function(d, fields) {
  if (!nrow(d))
    return(data.frame())
  ix <- split(seq_len(nrow(d)), interaction(d[, fields, drop = FALSE], drop = TRUE, lex.order = TRUE))
  res <- lapply(ix, function(i) {
    ans <- as.list(d[i[1], fields, drop = FALSE])
    ans$n_myeloid_cells <- length(i)
    ans$n_frozen_focus <- sum(d$is_focus_frozen[i])
    ans$pct_frozen_focus <- 100 * ans$n_frozen_focus/length(i)
    ans$n_marker_evaluable <- sum(!is.na(d$marker_codetect_104[i]))
    ans$n_marker_codetect <- sum(d$marker_codetect_104[i], na.rm = TRUE)
    ans$pct_marker_codetect <- if (ans$n_marker_evaluable > 0)
      100 * ans$n_marker_codetect/ans$n_marker_evaluable else NA_real_
    vals <- d$marker_score_104[i]
    ans$mean_marker_score <- if (any(is.finite(vals)))
      mean(vals, na.rm = TRUE) else NA_real_
    ans$median_marker_score <- if (any(is.finite(vals)))
      median(vals, na.rm = TRUE) else NA_real_
    ans$median_nCount_RNA <- if ("nCount_RNA" %in% names(d))
      median(d$nCount_RNA[i], na.rm = TRUE) else NA_real_
    as.data.frame(ans, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}
coverage <- summary_by(md0, c("source_dataset", "timepoint_main", "stage_kind_104", "treatment_104",
  "is_main_stage"))

coverage$n_library_units_confirmed <- NA_integer_
coverage$n_study_defined_sampling_units <- NA_integer_
for (k in seq_len(nrow(coverage))) {
  src <- coverage$source_dataset[k]
  st <- coverage$timepoint_main[k]
  if (src == "GSE162610") {
    n <- c(Control = 3L, `1 dpi` = 3L, `3 dpi` = 2L, `7 dpi` = 2L)[st]
    if (length(n) == 1 && !is.na(n)) {
      coverage$n_library_units_confirmed[k] <- unname(n)
      coverage$n_study_defined_sampling_units[k] <- unname(n)
    }
  } else if (src == "GSE196928") {
    coverage$n_library_units_confirmed[k] <- 1L
    coverage$n_study_defined_sampling_units[k] <- 1L
  }
}
coverage$unit_note <- ifelse(coverage$source_dataset == "17702045", "Stage/sex strata; specimens and library IDs unavailable",
  ifelse(coverage$source_dataset == "GSE162610", "GSM libraries may represent a single animal or an independently prepared 2-animal pool",
    "One pooled Vehicle library per included group; not animal-level replication"))
write_csv(coverage, "01_104_source_stage_treatment_coverage_v03.csv")
within_source <- summary_by(md0[md0$is_main_stage, ], c("source_dataset", "timepoint_main"))
within_source$unit_of_analysis <- "Cells (descriptive only; no specimen-level inference)"
write_csv(within_source, "02_104_frozen_state_within_source_v03.csv")

cases <- list(All_sources = SOURCES, Drop_17702045 = setdiff(SOURCES, "17702045"), Drop_GSE162610 = setdiff(SOURCES,
  "GSE162610"), Drop_GSE196928 = setdiff(SOURCES, "GSE196928"))
loo_tables <- lapply(names(cases), function(z) {
  d <- md0[md0$is_main_stage & md0$source_dataset %in% cases[[z]], , drop = FALSE]
  a <- summary_by(d, "timepoint_main")
  a$case <- z
  a$included_datasets <- paste(cases[[z]], collapse = "+")
  a
})
loo_tab <- do.call(rbind, loo_tables)
rownames(loo_tab) <- NULL

for (z in names(cases)) {
  absent <- setdiff(MAIN_STAGES, as.character(loo_tab$timepoint_main[loo_tab$case == z]))
  for (st in absent) {
    a <- as.list(loo_tab[1, , drop = FALSE])
    a[] <- NA
    a$timepoint_main <- st
    a$case <- z
    a$included_datasets <- paste(cases[[z]], collapse = "+")
    a$n_myeloid_cells <- 0L
    loo_tab <- rbind(loo_tab, as.data.frame(a, stringsAsFactors = FALSE))
  }
}
loo_tab$interpretation <- ifelse(loo_tab$n_myeloid_cells == 0, "NOT ESTIMABLE (no stage coverage)", "label; composition-only reweighting")
write_csv(loo_tab, "03_104_frozen_state_leave_one_out_v03.csv")

cluster_audits <- list()
branch_summaries <- list()
mix_audits <- list()
agreements <- list()
statuses <- list()
plots <- list()
logit("Starting independent re-clustering of the cell cohort")

ari <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]
  b <- b[ok]
  n <- length(a)
  if (n < 2)
    return(NA_real_)
  t <- table(a, b)
  choose2 <- function(x) x * (x - 1)/2
  observed <- sum(choose2(t))
  row_sum <- sum(choose2(rowSums(t)))
  col_sum <- sum(choose2(colSums(t)))
  expected <- row_sum * col_sum/choose2(n)
  maxindex <- (row_sum + col_sum)/2
  if (abs(maxindex - expected) < 1e-10)
    return(NA_real_)
  (observed - expected)/(maxindex - expected)
}

MARKER_HIGH_THRESHOLD <- mean(md0$marker_score_104, na.rm = TRUE) + 0.5 * stats::sd(md0$marker_score_104,
  na.rm = TRUE)
MARKER_POS_THRESHOLD <- mean(md0$marker_codetect_104, na.rm = TRUE)
cluster_table <- function(d, case_id, arm) {
  z <- summary_by(d, "cluster104")
  names(z)[names(z) == "cluster104"] <- "cluster"
  z$case <- case_id
  z$arm <- arm
  ord <- order(-z$mean_marker_score, -z$pct_marker_codetect)
  z$marker_rank <- integer(nrow(z))
  z$marker_rank[ord] <- seq_along(ord)
  z$frozen_focus_fraction <- z$pct_frozen_focus/100
  z$marker_rich_proxy <- z$n_myeloid_cells >= 30 & z$mean_marker_score >= MARKER_HIGH_THRESHOLD & z$pct_marker_codetect/100 >=
    MARKER_POS_THRESHOLD
  z$marker_high_threshold <- MARKER_HIGH_THRESHOLD
  z$co_detection_threshold <- MARKER_POS_THRESHOLD
  z
}

branch_summary <- function(d, case_id, arm, genes, nused, cluster_info) {

  high <- cluster_info$cluster[cluster_info$marker_rich_proxy]
  d$high_cluster <- d$cluster104 %in% high
  b <- summary_by(d, c("source_dataset", "timepoint_main"))
  key <- paste(d$source_dataset, d$timepoint_main, sep = "\r")
  sums <- tapply(as.integer(d$high_cluster), key, sum)
  bkey <- paste(b$source_dataset, b$timepoint_main, sep = "\r")
  b$n_marker_rich_reclustered <- as.integer(sums[bkey])
  b$pct_marker_rich_reclustered <- 100 * b$n_marker_rich_reclustered/b$n_myeloid_cells
  b$case <- case_id
  b$arm <- arm
  b$n_HVG <- length(genes)
  b$n_PCs <- nused
  b$marker_rich_definition <- paste0("cluster n>=30; mean score >=", signif(MARKER_HIGH_THRESHOLD,
    4), "; co-detect fraction >=", signif(MARKER_POS_THRESHOLD, 4))
  b
}

matched_stage_mixing <- function(embedding, d, case_id, arm) {
  if (!RUN_KNN_SOURCE_MIXING || !requireNamespace("RANN", quietly = TRUE))
    return(data.frame())
  if (length(unique(d$source_dataset)) < 2)
    return(data.frame())
  set.seed(SEED)
  d <- d[rownames(embedding), , drop = FALSE]
  combo <- interaction(d$timepoint_main, d$source_dataset, drop = TRUE)
  sampled <- unlist(lapply(split(seq_len(nrow(d)), combo), function(ii) {
    if (length(ii) <= MAX_KNN_CELLS_PER_SOURCE_STAGE)
      ii else sample(ii, MAX_KNN_CELLS_PER_SOURCE_STAGE)
  }), use.names = FALSE)
  dm <- d[sampled, , drop = FALSE]
  emb <- embedding[dm$cell, , drop = FALSE]
  out <- list()
  for (stage in unique(dm$timepoint_main)) {
    ii <- which(dm$timepoint_main == stage)
    s <- dm$source_dataset[ii]
    if (length(unique(s)) < 2 || length(ii) < 12)
      next
    k <- min(21L, length(ii))
    nn <- RANN::nn2(data = emb[ii, , drop = FALSE], k = k)$nn.idx
    nbr <- nn[, -1, drop = FALSE]
    same <- mean(matrix(s[nbr], nrow = nrow(nbr)) == matrix(s, nrow = length(s), ncol = ncol(nbr)))
    p <- prop.table(table(s))
    expected_same <- sum(p^2)
    out[[length(out) + 1L]] <- data.frame(case = case_id, arm = arm, stage = stage, sources = paste(sort(unique(s)),
      collapse = ";"), n_cells_sampled = length(ii), k = k - 1L, observed_same_source_neighbors = same,
      expected_same_source_random = expected_same, delta_observed_minus_random = same - expected_same,
      note = if (stage == "Control")
        "Sham vs uninjured sources are biologically non-equivalent" else "Stage-matched neighborhood; source can encode protocol differences")
  }
  if (!length(out))
    data.frame() else do.call(rbind, out)
}

make_umap_plot <- function(obj, d, case_id, arm, reduction_name) {
  if (!reduction_name %in% Reductions(obj))
    return(NULL)
  set.seed(SEED)
  keep <- sample(Cells(obj), min(ncol(obj), MAX_PLOT_CELLS))
  coords <- Embeddings(obj, reduction = reduction_name)[keep, 1:2, drop = FALSE]
  pdat <- data.frame(UMAP1 = coords[, 1], UMAP2 = coords[, 2], source = d[keep, "source_dataset"],
    frozen_state = ifelse(d[keep, "is_focus_frozen"], "Clec7a+", "Other"), de_novo_cluster = d[keep,
      "cluster104"], stringsAsFactors = FALSE)
  p <- ggplot(pdat, aes(x = UMAP1, y = UMAP2, color = source)) + geom_point(size = 0.12, alpha = 0.5) +
    coord_equal() + labs(title = paste(case_id, arm, "(recomputed)"), subtitle = "Colored by source; UMAP geometry cannot prove stage/source identifiability") +
    theme_classic(base_size = 11)
  p
}

run_one_case <- function(case_id, included_sources, single_source = FALSE) {
  start_time <- Sys.time()
  logit("\n================== ", case_id, " ==================")
  use_cells <- md0$cell[md0$source_dataset %in% included_sources & md0$count_pass_104]
  if (length(use_cells) < 200)
    stop("Too few cells for case ", case_id)
  o <- subset(obj0, cells = use_cells)
  DefaultAssay(o) <- "RNA"
  o@reductions <- list()
  o@graphs <- list()
  o@neighbors <- list()

  o <- tryCatch(JoinLayers(o, assay = "RNA"), error = function(e) {
    logit("JoinLayers warning: ", conditionMessage(e))
    o
  })
  test_count <- GetAssayData(o, assay = "RNA", layer = "counts")
  if (ncol(test_count) != ncol(o))
    stop("Joined RNA counts incomplete in case ", case_id)
  rm(test_count)
  o <- NormalizeData(o, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = FALSE)
  o <- FindVariableFeatures(o, assay = "RNA", selection.method = "vst", nfeatures = min(N_HVG, nrow(o)),
    verbose = FALSE)
  genes <- VariableFeatures(o)
  if (length(genes) < 100)
    stop("Unexpectedly few HVGs in ", case_id)
  o <- ScaleData(o, assay = "RNA", features = genes, verbose = FALSE)
  nc <- min(N_PCS_COMPUTE, length(genes) - 1L, ncol(o) - 1L)
  nused <- min(N_PCS_COMPARE, nc)
  o <- RunPCA(o, assay = "RNA", features = genes, npcs = nc, reduction.name = "pca_104", reduction.key = "r104PC_",
    seed.use = SEED, verbose = FALSE)
  d <- md0[Cells(o), , drop = FALSE]
  stopifnot(identical(d$cell, Cells(o)))

  o <- FindNeighbors(o, reduction = "pca_104", dims = seq_len(nused), graph.name = c("r104_pre_nn",
    "r104_pre_snn"), verbose = FALSE)
  o <- FindClusters(o, graph.name = "r104_pre_snn", resolution = RESOLUTION, cluster.name = "cluster_104_pre",
    random.seed = SEED, verbose = FALSE)
  d$cluster104 <- as.character(o$cluster_104_pre)
  audit_pre <- cluster_table(d, case_id, "PCA_only")
  cluster_audits[[paste0(case_id, "_PCA")]] <<- audit_pre
  branch_summaries[[paste0(case_id, "_PCA")]] <<- branch_summary(d, case_id, "PCA_only", genes, nused,
    audit_pre)

  pca_embeddings <- Embeddings(o, reduction = "pca_104")[, seq_len(nused), drop = FALSE]
  mix_audits[[paste0(case_id, "_PCA")]] <<- matched_stage_mixing(pca_embeddings, d, case_id, "PCA_only")
  cluster_pre <- setNames(d$cluster104, d$cell)
  if (RUN_UMAP_FOR_ALL_BRANCHES || (RUN_UMAP_FOR_FULL && case_id == "All_sources")) {
    o <- RunUMAP(o, reduction = "pca_104", dims = seq_len(nused), reduction.name = "umap_104_pre",
      reduction.key = "r104PRE_", seed.use = SEED, verbose = FALSE)
    plots[[paste0(case_id, "_PCA")]] <<- make_umap_plot(o, d, case_id, "PCA_only", "umap_104_pre")
  }

  if (!single_source && length(unique(d$source_dataset)) >= 2) {
    o <- harmony::RunHarmony(object = o, group.by.vars = "source_dataset", reduction.use = "pca_104",
      dims.use = seq_len(nused), reduction.save = "harmony_104", max_iter = HARMONY_MAX_ITER, plot_convergence = FALSE,
      verbose = FALSE)
    if (!"harmony_104" %in% Reductions(o))
      stop("Missing harmony_104 reduction")
    o <- FindNeighbors(o, reduction = "harmony_104", dims = seq_len(nused), graph.name = c("r104_harmony_nn",
      "r104_harmony_snn"), verbose = FALSE)
    o <- FindClusters(o, graph.name = "r104_harmony_snn", resolution = RESOLUTION, cluster.name = "cluster_104_harmony",
      random.seed = SEED, verbose = FALSE)
    d$cluster104 <- as.character(o$cluster_104_harmony)
    audit_h <- cluster_table(d, case_id, "Harmony")
    cluster_audits[[paste0(case_id, "_Harmony")]] <<- audit_h
    branch_summaries[[paste0(case_id, "_Harmony")]] <<- branch_summary(d, case_id, "Harmony", genes,
      nused, audit_h)
    harmony_emb <- Embeddings(o, reduction = "harmony_104")[, seq_len(nused), drop = FALSE]
    mix_audits[[paste0(case_id, "_Harmony")]] <<- matched_stage_mixing(harmony_emb, d, case_id, "Harmony")
    cluster_harmony <- setNames(d$cluster104, d$cell)
    matches <- intersect(names(cluster_pre), names(cluster_harmony))
    agreements[[case_id]] <<- data.frame(case = case_id, n_cells = length(matches), n_pre_clusters = length(unique(cluster_pre)),
      n_harmony_clusters = length(unique(cluster_harmony)), adjusted_rand_index = ari(cluster_pre[matches],
        cluster_harmony[matches]), comparison = "PCA vs Harmony de novo clusters; cohort", stringsAsFactors = FALSE)
    if (RUN_UMAP_FOR_ALL_BRANCHES || (RUN_UMAP_FOR_FULL && case_id == "All_sources")) {
      o <- RunUMAP(o, reduction = "harmony_104", dims = seq_len(nused), reduction.name = "umap_104_harmony",
        reduction.key = "r104HAR_", seed.use = SEED, verbose = FALSE)
      plots[[paste0(case_id, "_Harmony")]] <<- make_umap_plot(o, d, case_id, "Harmony", "umap_104_harmony")
    }
  }
  if (SAVE_FULL_RDS) {
    saveRDS(o, file.path(OUTPUT, paste0("RDS_104_", case_id, "_v03.rds")), compress = TRUE)
  }
  statuses[[case_id]] <<- data.frame(case = case_id, status = "success", included_sources = paste(included_sources,
    collapse = ";"), input_cells = length(use_cells), n_zero_count_excluded = sum(md0$source_dataset %in%
    included_sources & !md0$count_pass_104), n_HVG = length(genes), n_PCs = nused, independent_biological_validation = FALSE,
    elapsed_minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), note = if (single_source)
      "PCA only: Harmony not meaningful with one dataset" else "Recomputed HVG/PCA, PCA clusters, and source-adjusted Harmony clusters")
  logit("Completed ", case_id, ": ", numfmt(length(use_cells)), " cells")
  rm(o, d)
  gc(verbose = FALSE)
}

run_safely <- function(case_id, sources, single_source = FALSE) {
  tryCatch(run_one_case(case_id, sources, single_source), error = function(e) {
    logit("ERROR for ", case_id, ": ", conditionMessage(e))
    statuses[[case_id]] <<- data.frame(case = case_id, status = "failed", included_sources = paste(sources,
      collapse = ";"), input_cells = sum(md0$source_dataset %in% sources & md0$count_pass_104),
      n_zero_count_excluded = sum(md0$source_dataset %in% sources & !md0$count_pass_104), n_HVG = NA_integer_,
      n_PCs = NA_integer_, independent_biological_validation = FALSE, elapsed_minutes = NA_real_,
      note = conditionMessage(e), stringsAsFactors = FALSE)
    invisible(gc(verbose = FALSE))
  })
}
for (case_id in names(cases)) run_safely(case_id, cases[[case_id]])
if (RUN_SINGLE_SOURCE_RECLUSTER) {
  for (s in SOURCES) run_safely(paste0("Only_", s), s, single_source = TRUE)
}

bind_list <- function(x) {
  x <- Filter(function(v) is.data.frame(v) && nrow(v) > 0, x)
  if (!length(x))
    return(data.frame(note = "Analysis not available"))
  rownames_out <- do.call(rbind, x)
  rownames(rownames_out) <- NULL
  rownames_out
}
write_csv(bind_list(cluster_audits), "04_104_recluster_cluster_marker_audit_v03.csv")
write_csv(bind_list(branch_summaries), "05_104_recluster_source_stage_markers_v03.csv")
write_csv(bind_list(agreements), "06_104_pre_vs_harmony_cluster_agreement_v03.csv")
write_csv(bind_list(mix_audits), "07_104_stage_matched_source_mixing_v03.csv")
write_csv(bind_list(statuses), "08_104_branch_status_v03.csv")

not_identifiable <- data.frame(comparison = c("14 dpi vs 28 dpi (independent time effect)", "14 dpi after removing Figshare 17702045",
  "28 dpi after removing GSE196928", "Animal-level reproducibility in Figshare 17702045", "Animal-level variance for GSE196928 7 vs 28 dpi",
  "PLX5622 treatment effect in the analyzed atlas", "Cross-source Control vs Sham biological equivalence",
  "Pure source/batch effect separated from stage", "Independent external validation of Clec7a state",
  "Causal microglial lineage or lipid clearance"), verdict = c("Not identifiable", "No coverage", "No coverage",
  "Not testable", "Not testable", "Not estimable", "Not established", "Only partially identifiable",
  "Not provided by this analysis", "Not testable"), reason = c("14 dpi occurs only in Figshare; 28 dpi only in GSE196928; stage and source are confounded",
  "No 14 dpi cells remain in the retained sources", "No 28 dpi cells remain in the retained sources",
  "Individual animal/library IDs absent from public Figshare cell metadata", "Each retained GSE196928 condition corresponds to one pooled library; individual mice not separately sequenced",
  "Only Vehicle B1/B3/B7 were retained; no PLX5622 group in atlas", "GSE196928 baseline is surgically sham-treated; others contain uninjured samples",
  "Only partially overlapping stages; protocols, tissue selection and source vary together", "selection and source annotations originated in the same integrated project",
  "RNA expression/proxy scores do not establish origin or cellular function"), stringsAsFactors = FALSE)
write_csv(not_identifiable, "09_104_non_identifiable_contrasts_v03.csv")

pdf_file <- file.path(OUTPUT, "10_104_analysis_figures_v03.pdf")
pdf(pdf_file, width = 12, height = 8.5, onefile = TRUE)
tryCatch({

  ps <- ggplot(within_source, aes(x = factor(timepoint_main, levels = MAIN_STAGES), y = pct_frozen_focus,
    group = source_dataset, color = source_dataset)) + geom_line(linewidth = 0.75) + geom_point(size = 2) +
    labs(title = "Clec7a+ state: within-dataset descriptive proportions", subtitle = "Cannot infer animal-level reproducibility from cell counts",
      x = "Injury stage (categorical)", y = "% of retained myeloid cells") + theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  print(ps)
  pcd <- ggplot(within_source, aes(x = factor(timepoint_main, levels = MAIN_STAGES), y = pct_marker_codetect,
    group = source_dataset, color = source_dataset)) + geom_line(linewidth = 0.75) + geom_point(size = 2) +
    labs(title = "Fixed Clec7a/Lpl/(Apoe or Trem2) co-detection", subtitle = "Independent of Harmony coordinates; influenced by transcript capture",
      x = "Injury stage (categorical)", y = "% marker-co-detected cells") + theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  print(pcd)
  lp <- loo_tab[loo_tab$n_myeloid_cells > 0, , drop = FALSE]
  ploo <- ggplot(lp, aes(x = factor(timepoint_main, levels = MAIN_STAGES), y = pct_frozen_focus, group = case,
    color = case)) + geom_point(size = 1.8) + geom_line(linewidth = 0.7) + facet_wrap(~case, ncol = 2) +
    labs(title = "Leave-one-dataset-out: -label composition ONLY", subtitle = "Panels omit stages without coverage; NOT a new Harmony fit",
      x = "Stage", y = "Clec7a+ state (%)") + theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 40,
    hjust = 1), legend.position = "none")
  print(ploo)
  for (nm in names(plots)) if (!is.null(plots[[nm]]))
    print(plots[[nm]])
}, error = function(e) logit("Figure creation error: ", conditionMessage(e)), finally = dev.off())

status_df <- bind_list(statuses)
for (j in seq_len(nrow(status_df))) {
}
writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUTPUT, "11_104_REPORT_v03.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUTPUT, "12_104_sessionInfo_v03.txt"))

if (MAKE_SMALL_ZIP && requireNamespace("zip", quietly = TRUE)) {
  files <- list.files(OUTPUT, pattern = "\\.(csv|txt|log|pdf)$", full.names = TRUE)
  total_bytes <- sum(file.info(files)$size, na.rm = TRUE)
  zip_path <- file.path(OUTPUT, "104_Major1_sensitivity_results_v03.zip")
  if (length(files) && is.finite(total_bytes) && total_bytes < ZIP_LIMIT_MB * 1024^2) {
    oldwd <- getwd()
    tryCatch({
      setwd(OUTPUT)
      zip::zipr(zipfile = zip_path, files = basename(files), recurse = FALSE)
    }, error = function(e) logit("Optional zip failed: ", conditionMessage(e)), finally = setwd(oldwd))
    if (file.exists(zip_path)) {
      if (file.info(zip_path)$size > ZIP_LIMIT_MB * 1024^2) {
        unlink(zip_path)
        logit("ZIP exceeded max size; deleted ZIP (original results retained)")
      } else {
        logit("Created results ZIP: ", zip_path)
      }
    }
  } else logit("Skipping ZIP because uncompressed outputs exceed ", ZIP_LIMIT_MB, " MB")
} else if (MAKE_SMALL_ZIP) {
  logit("zip package not installed: files were saved but ZIP was not created")
}
logit("\n104 DONE. Please upload the small ZIP or CSV/PDF/log files.")
logit("IMPORTANT: Never claim successfully completed branches that 08 marks failed.")

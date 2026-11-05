options(stringsAsFactors = FALSE, warn = 1)
ROOT <- "D:/zz"
SEED <- 1234L
DBR_VALUES <- c(0.07, 0.05, 0.1)
EXPECTED_INPUT_N <- 19956L
EXPECTED_ORIGINAL_SINGLETS <- 18087L
VERSION <- "v01"

if (!dir.exists(ROOT)) stop("Project root does not exist: ", ROOT)
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
required <- c("Seurat", "SeuratObject", "Matrix", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scDblFinder)
})

OUT <- file.path(ROOT, "results", "140_Figshare17702045_doublet_rate_sensitivity")
TBL <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "140_Figshare17702045_doublet_rate_sensitivity_v01.log")
if (file.exists(LOG)) stop("Output log already exists. To avoid overwriting a prior run, move/rename output folder first: ",
  LOG)

logline <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n")
  cat(s, "\n", file = LOG, append = TRUE)
}
write_tab <- function(x, basename) {
  path <- file.path(TBL, basename)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  logline("WROTE ", path, " (", nrow(x), " rows)")
  invisible(path)
}
rate_name <- function(x) sprintf("dbr_%02d", as.integer(round(x * 100)))
as_chr <- function(x) as.character(x)
percent <- function(num, den) ifelse(den > 0, 100 * num/den, NA_real_)

norm_stage <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s <- gsub("[[:space:]_.-]+", "", s)
  out <- rep(NA_character_, length(s))
  out[s %in% c("uninjured", "uninj", "control", "ctrl", "0d", "0", "naive")] <- "Control"
  out[s %in% c("4h", "4hr", "4hours")] <- "4h"
  for (d in c(1L, 3L, 7L, 14L, 28L, 38L)) {
    out[s %in% c(paste0(d, "d"), paste0(d, "dpi"), paste0("d", d), paste0("day", d), paste0(d, "days"))] <- paste0(d,
      "dpi")
  }
  out
}
STAGES <- c("Control", "4h", "1dpi", "3dpi", "7dpi", "14dpi", "28dpi", "38dpi")
MAIN <- c("Control", "1dpi", "3dpi", "7dpi", "14dpi", "28dpi")
CLEC_STATE <- "Clec7a+ lipid-associated Microglia"

IN_DECONTX <- file.path(ROOT, "results", "03_DecontX", "03_DecontX_17702045_v01.rds")
IN_ORIG <- file.path(ROOT, "results", "04_doublet", "04_doublet_17702045_singlet_v03.rds")
IN_SUMMARY <- file.path(ROOT, "results", "04_doublet", "04_doublet_3scRNA_summary_v03.csv")
IN_FROZEN <- file.path(ROOT, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
IN_SCRIPT <- file.path(ROOT, "script", "04_doublet_scDblFinder_3scRNA_v03.R")
FILES <- c(decontx_pre_doublet = IN_DECONTX, original_7pct_singlets = IN_ORIG, original_summary = IN_SUMMARY,
  frozen_myeloid_atlas = IN_FROZEN, original_script = IN_SCRIPT)
if (!all(file.exists(FILES))) {
  stop("Missing required local files:\n", paste(paste(names(FILES), FILES)[!file.exists(FILES)], collapse = "\n"))
}
info <- file.info(FILES)
audit <- data.frame(role = names(FILES), path = as.character(FILES), size_bytes = as.numeric(info$size),
  modified = as.character(info$mtime), stringsAsFactors = FALSE)
write_tab(audit, "01_140_INPUT_AUDIT_v01.csv")
logline("Start: ", Sys.time(), "; project root=", ROOT)
logline("Only Figshare 17702045 is processed. Original files are READ ONLY.")
logline("scDblFinder version: ", as.character(packageVersion("scDblFinder")))
logline("Seurat version: ", as.character(packageVersion("Seurat")))
logline("Fixed values 0.05 / 0.07 / 0.10 are sensitivity assumptions, not measured capture rates.")
logline("No valid independent capture/library identifier for Figshare 17702045 was established.")
logline("Therefore stage/timepoint labels are NOT supplied as scDblFinder samples.")

obj <- readRDS(IN_DECONTX)
if (!inherits(obj, "Seurat")) stop("Expected Seurat object in DecontX input.")
if (!"RNA" %in% names(obj@assays)) stop("RNA assay missing.")
if (any(grepl("PLX5622", colnames(obj[[]]), ignore.case = TRUE))) {
  logline("Warning: unexpected metadata column containing PLX5622. Inspect input.")
}
meta <- obj[[]]
if (!identical(rownames(meta), colnames(obj))) stop("Metadata cell IDs are not aligned.")
if (anyDuplicated(rownames(meta))) stop("Non-unique input cell IDs.")

stage_candidates <- c("time", "timepoint", "timepoint_main", "time_raw")
stage_candidates <- intersect(stage_candidates, colnames(meta))
if (!length(stage_candidates)) stop("No stage metadata field on 17702045 DecontX object.")
best_stage_field <- stage_candidates[which.max(vapply(stage_candidates, function(k) sum(!is.na(norm_stage(meta[[k]]))),
  integer(1)))]
stage <- norm_stage(meta[[best_stage_field]])
if (anyNA(stage)) {
  bad <- unique(as.character(meta[[best_stage_field]])[is.na(stage)])
  stop("Cannot safely map time points: ", paste(bad, collapse = ", "))
}
ids <- colnames(obj)
if (length(ids) != EXPECTED_INPUT_N) {
  stop("17702045 DecontX cell count changed: expected ", EXPECTED_INPUT_N, ", got ", length(ids), ". Verify upstream object before testing sensitivity.")
}
DefaultAssay(obj) <- "RNA"
get_counts <- function(x) {
  ct <- tryCatch(SeuratObject::GetAssayData(x, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(ct)) {
    ct <- tryCatch(SeuratObject::GetAssayData(x, assay = "RNA", slot = "counts"), error = function(e) NULL)
  }
  if (is.null(ct)) {
    joined <- tryCatch(SeuratObject::JoinLayers(x, assay = "RNA"), error = function(e) NULL)
    if (!is.null(joined))
      ct <- SeuratObject::GetAssayData(joined, assay = "RNA", layer = "counts")
  }
  if (is.null(ct))
    stop("Cannot retrieve RNA counts from pre-doublet RDS.")
  as(ct, "dgCMatrix")
}
counts <- get_counts(obj)
if (!identical(colnames(counts), ids)) {
  if (!setequal(colnames(counts), ids))
    stop("Count matrix cell IDs do not match.")
  counts <- counts[, ids, drop = FALSE]
}
logline("Input: ", nrow(counts), " genes x ", ncol(counts), " cells; stage_col=", best_stage_field)

rm(obj, meta)
invisible(gc())

calls <- data.frame(cell_id = ids, stage = stage, stringsAsFactors = FALSE)
scenario <- list()
for (dbr in DBR_VALUES) {
  nm <- rate_name(dbr)
  logline("=== Running ", nm, " ===")
  set.seed(SEED)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce <- scDblFinder::scDblFinder(sce, dbr = dbr)
  out_id <- colnames(sce)
  if (!identical(out_id, ids))
    stop("Cell order/IDs changed during ", nm)
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!all(c("scDblFinder.class", "scDblFinder.score") %in% names(cd))) {
    stop("Expected class/score missing: ", nm)
  }
  cc <- as.character(cd[["scDblFinder.class"]])
  ss <- as.numeric(cd[["scDblFinder.score"]])
  if (anyNA(cc) || anyNA(ss) || !all(cc %in% c("singlet", "doublet"))) {
    stop("Missing/invalid class or score in ", nm)
  }
  calls[[paste0(nm, "_class")]] <- cc
  calls[[paste0(nm, "_score")]] <- ss
  n_d <- sum(cc == "doublet")
  scenario[[nm]] <- data.frame(dataset = "17702045", unit = "whole_dataset", dbr_expected = dbr, expected_dbr_pct = dbr *
    100, n_before = length(ids), n_called_doublets = n_d, observed_doublet_pct = percent(n_d, length(ids)),
    n_retained_singlets = sum(cc == "singlet"), score_median = stats::median(ss), score_p95 = as.numeric(stats::quantile(ss,
      0.95)), stringsAsFactors = FALSE)
  logline(nm, ": doublets=", n_d, "; singlets=", sum(cc == "singlet"), "; observed=", sprintf("%.3f%%",
    percent(n_d, length(ids))))
  rm(sce, cd, cc, ss)
  invisible(gc())
}
rm(counts)
invisible(gc())
write_tab(do.call(rbind, scenario), "02_140_DBR_SCENARIO_SUMMARY_v01.csv")

orig <- readRDS(IN_ORIG)
if (!inherits(orig, "Seurat")) stop("Original v03 singlet file is not a Seurat object.")
orig_ids <- colnames(orig)
rm(orig)
invisible(gc())
if (length(orig_ids) != EXPECTED_ORIGINAL_SINGLETS) {
  stop("Historical 7% file contains ", length(orig_ids), " singlets (expected ", EXPECTED_ORIGINAL_SINGLETS,
    "). Audit upstream provenance.")
}
if (!all(orig_ids %in% ids)) stop("Original 7% singlet IDs do not all match DecontX input.")
calls$historical_7pct_class <- ifelse(calls$cell_id %in% orig_ids, "singlet", "doublet")
original <- calls$historical_7pct_class
base <- calls$dbr_07_class
baseline_disagreement <- sum(original != base)
logline("Historic-vs-rerun 7% disagree on ", baseline_disagreement, " cells.")
if (baseline_disagreement) logline("CAUTION: 7% rerun is not exactly identical to original v03; interpretation must be qualified.")

pairwise <- list()
k <- 0L
for (refnm in c("historical_7pct", "dbr_07")) {
  ca <- if (refnm == "historical_7pct")
    original else calls$dbr_07_class
  for (dnm in c("dbr_05", "dbr_07", "dbr_10")) {
    cb <- calls[[paste0(dnm, "_class")]]
    n_stable <- sum(ca == cb)
    n_s_s <- sum(ca == "singlet" & cb == "singlet")
    n_s_d <- sum(ca == "singlet" & cb == "doublet")
    n_d_s <- sum(ca == "doublet" & cb == "singlet")
    n_d_d <- sum(ca == "doublet" & cb == "doublet")
    k <- k + 1L
    pairwise[[k]] <- data.frame(reference = refnm, scenario = dnm, n_cells = length(ids), same_class = n_stable,
      agreement_pct = percent(n_stable, length(ids)), singlet_to_singlet = n_s_s, singlet_to_doublet = n_s_d,
      doublet_to_singlet = n_d_s, doublet_to_doublet = n_d_d, singlet_jaccard = if ((n_s_s + n_s_d +
        n_d_s) > 0)
        n_s_s/(n_s_s + n_s_d + n_d_s) else NA_real_, stringsAsFactors = FALSE)
  }
}
write_tab(do.call(rbind, pairwise), "03_140_CELL_CALL_AGREEMENT_v01.csv")
con <- gzfile(file.path(TBL, "04_140_CELL_LEVEL_CALLS_v01.csv.gz"), open = "wt")
utils::write.csv(calls, con, row.names = FALSE)
close(con)
logline("WROTE 04_140_CELL_LEVEL_CALLS_v01.csv.gz")

stage_results <- list()
j <- 0L
for (nm in c("historical_7pct", "dbr_05", "dbr_07", "dbr_10")) {
  cls <- if (nm == "historical_7pct")
    original else calls[[paste0(nm, "_class")]]
  for (st in STAGES) {
    idx <- calls$stage == st
    if (!any(idx))
      next
    j <- j + 1L
    nd <- sum(cls[idx] == "doublet")
    stage_results[[j]] <- data.frame(scenario = nm, stage = st, unit_of_scDblFinder = "whole_dataset",
      stage_is_captured_library = FALSE, n_before = sum(idx), n_called_doublets = nd, n_singlets = sum(idx) -
        nd, observed_doublet_pct = percent(nd, sum(idx)), stringsAsFactors = FALSE)
  }
}
stage_df <- do.call(rbind, stage_results)
write_tab(stage_df, "05_140_STAGE_DOUBLETS_v01.csv")

logline("Loading Step18 myeloid metadata to evaluate conditional robustness...")
ref <- readRDS(IN_FROZEN)
if (!inherits(ref, "Seurat")) stop("Step 18 object is not Seurat.")
md <- ref[[]]
rm(ref)
invisible(gc())
need <- c("source_dataset", "myeloid_state", "timepoint_main")
if (!all(need %in% names(md))) stop("Step18 lacks metadata: ", paste(setdiff(need, names(md)), collapse = ", "))
if ("treatment" %in% names(md) && any(grepl("PLX5622", as_chr(md$treatment), ignore.case = TRUE), na.rm = TRUE)) {
  stop("PLX5622 found in atlas metadata; do not use this reference until audited.")
}
md$original_atlas_cell_id <- rownames(md)
md$source_dataset <- as_chr(md$source_dataset)
md$myeloid_state <- as_chr(md$myeloid_state)
md$stage <- norm_stage(md$timepoint_main)
md <- md[!is.na(md$stage) & nzchar(md$myeloid_state), , drop = FALSE]
f177 <- md[md$source_dataset == "17702045", , drop = FALSE]
if (!nrow(f177)) stop("atlas contains no 17702045 myeloid cells.")

possible <- list(prefixed_id = sub("^F177_", "", f177$original_atlas_cell_id), raw_atlas_id = f177$original_atlas_cell_id)
if ("cell_id_original" %in% names(f177)) {
  possible$explicit_metadata_id <- as_chr(f177$cell_id_original)
}
matched <- vapply(possible, function(v) sum(!is.na(v) & v %in% ids), integer(1))
chosen <- names(matched)[which.max(matched)]
f177$source_cell_id <- possible[[chosen]]
logline("177 cell key strategy=", chosen, "; matches=", max(matched), "/", nrow(f177))
if (max(matched) < 100L || anyDuplicated(f177$source_cell_id)) {
  stop("cell ID join is insufficient or ambiguous; inspect source IDs. No label inference is allowed.")
}
if (any(!f177$source_cell_id %in% ids)) {
  stop("source 177 cell IDs not all found in DecontX source after key normalization.")
}
if (any(!f177$source_cell_id %in% orig_ids)) {
  stop("Some 177 cells are not in historical v03 singlets: invalid reference provenance.")
}

ref_stats <- list()
m <- 0L
for (nm in c("historical_7pct", "dbr_05", "dbr_07", "dbr_10")) {
  cls <- if (nm == "historical_7pct")
    original else calls[[paste0(nm, "_class")]]
  keep_ids <- calls$cell_id[cls == "singlet"]
  keep_ref <- f177$source_cell_id %in% keep_ids
  ref_sub <- rbind(md[md$source_dataset != "17702045", c("stage", "myeloid_state"), drop = FALSE],
    f177[keep_ref, c("stage", "myeloid_state"), drop = FALSE])
  for (st in MAIN) {
    sub <- ref_sub[ref_sub$stage == st, , drop = FALSE]
    m <- m + 1L
    n_clec <- sum(sub$myeloid_state == CLEC_STATE)
    ref_stats[[m]] <- data.frame(scenario = nm, stage = st, n_frozen_annotated_myeloid = nrow(sub),
      n_frozen_annotated_clec7a = n_clec, clec7a_fraction_pct = percent(n_clec, nrow(sub)), interpretation = "-label subset only; new rescued singlets not annotated",
      stringsAsFactors = FALSE)
  }
}
props <- do.call(rbind, ref_stats)
write_tab(props, "06_140_FROZEN_ATLAS_STAGE_PROPORTIONS_v01.csv")

source_rows <- list()
src_i <- 0L
for (nm in c("historical_7pct", "dbr_05", "dbr_07", "dbr_10")) {
  cls <- if (nm == "historical_7pct")
    original else calls[[paste0(nm, "_class")]]
  kept_177 <- calls$cell_id[cls == "singlet"]
  sources <- unique(md$source_dataset)
  for (src in sources) {
    dat <- if (src == "17702045") {
      f177[f177$source_cell_id %in% kept_177, c("stage", "myeloid_state"), drop = FALSE]
    } else {
      md[md$source_dataset == src, c("stage", "myeloid_state"), drop = FALSE]
    }
    for (st in MAIN) {
      one <- dat[dat$stage == st, , drop = FALSE]
      n1 <- sum(one$myeloid_state == CLEC_STATE)
      src_i <- src_i + 1L
      source_rows[[src_i]] <- data.frame(scenario = nm, source_dataset = src, stage = st, n_frozen_myeloid = nrow(one),
        n_frozen_clec7a = n1, pct_frozen_clec7a = percent(n1, nrow(one)), stringsAsFactors = FALSE)
    }
  }
}
source_props <- do.call(rbind, source_rows)
write_tab(source_props, "06b_140_SOURCE_SPECIFIC_STAGE_PROPORTIONS_v01.csv")

phase_df <- list()
jj <- 0L
for (nm in unique(props$scenario)) {
  x <- props[props$scenario == nm, , drop = FALSE]
  for (phase in c("early_1_3_7dpi", "late_14_28dpi")) {
    stages <- if (phase == "early_1_3_7dpi")
      c("1dpi", "3dpi", "7dpi") else c("14dpi", "28dpi")
    sub <- x[x$stage %in% stages, , drop = FALSE]
    n0 <- sum(sub$n_frozen_annotated_myeloid)
    n1 <- sum(sub$n_frozen_annotated_clec7a)
    jj <- jj + 1L
    phase_df[[jj]] <- data.frame(scenario = nm, phase = phase, n_frozen_myeloid = n0, n_frozen_clec7a = n1,
      clec7a_pct = percent(n1, n0), stringsAsFactors = FALSE)
  }
}
phase_tab <- do.call(rbind, phase_df)
write_tab(phase_tab, "07_140_FROZEN_ATLAS_ENRICHMENT_v01.csv")

unknown_rows <- list()
t <- 0L
for (nm in c("dbr_05", "dbr_07", "dbr_10")) {
  alt <- calls[[paste0(nm, "_class")]]
  for (st in MAIN) {
    idx <- calls$stage == st
    new <- sum(idx & original == "doublet" & alt == "singlet")
    prev <- props[props$scenario == nm & props$stage == st, , drop = FALSE]
    den <- prev$n_frozen_annotated_myeloid
    numer <- prev$n_frozen_annotated_clec7a
    t <- t + 1L
    unknown_rows[[t]] <- data.frame(scenario = nm, stage = st, newly_rescued_unannotated_177_singlets = new,
      frozen_annotated_myeloid = den, frozen_annotated_clec7a = numer, known_only_pct = percent(numer,
        den), conservative_lower_pct = percent(numer, den + new), conservative_upper_pct = percent(numer +
        new, den + new), note = "Bounds assume all rescued cells are myeloid, none/all Clec7a; not full reannotation",
      stringsAsFactors = FALSE)
  }
}
write_tab(do.call(rbind, unknown_rows), "08_140_UNANNOTATED_NEW_SINGLET_BOUNDS_v01.csv")

qc <- data.frame(check = c("input_cells_match_original", "original_singlet_count_match", "rerun_7pct_exactly_reproduces_original",
  "source_177_frozen_cell_IDs_valid", "no_PLX5622_processing", "frozen_labels_reassigned_to_new_singlets",
  "specimen_level_replication_test", "rerun_all_downstream_integration"), passed = c(length(ids) ==
  EXPECTED_INPUT_N, length(orig_ids) == EXPECTED_ORIGINAL_SINGLETS, baseline_disagreement == 0, max(matched) ==
  nrow(f177), TRUE, FALSE, FALSE, FALSE), meaning = c("Expected DecontX input size", "Expected original 7% singlet total",
  "Must be TRUE before claiming exact historical reproducibility", "All source 177 cells have traceable unique IDs",
  "Only 17702045 inputs were accessed", "Intentional: new cells have no old annotation", "Not possible without independent specimen labels",
  "Outside scope: this is -annotation sensitivity"), stringsAsFactors = FALSE)
write_tab(qc, "09_140_FINAL_QC_v01.csv")

pdf(file.path(FIG, "10_140_DBR_SENSITIVITY_v01.pdf"), width = 11, height = 7.5)
tryCatch({
  ss <- do.call(rbind, scenario)
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 1))
  barplot(ss$observed_doublet_pct, names.arg = sprintf("%d%%", round(ss$expected_dbr_pct)), ylab = "Called doublets / QC cells (%)",
    xlab = "Expected dbr", main = "Figshare 17702045: doublet calls")
  barplot(ss$n_retained_singlets, names.arg = sprintf("%d%%", round(ss$expected_dbr_pct)), ylab = "Retained singlets",
    xlab = "Expected dbr", main = "Figshare 17702045: singlet retention")
  par(old_par)
  par(mfrow = c(1, 1), mar = c(6, 4, 4, 1))
  r <- props[props$scenario %in% c("dbr_05", "dbr_07", "dbr_10"), , drop = FALSE]
  mat <- sapply(c("dbr_05", "dbr_07", "dbr_10"), function(nm) r$clec7a_fraction_pct[match(paste(nm,
    MAIN), paste(r$scenario, r$stage))])
  matplot(seq_along(MAIN), mat, type = "b", pch = c(15, 16, 17), lty = 1, xlab = "Stage (six-stage subset)",
    ylab = "Clec7a+ / annotated myeloid (%)", xaxt = "n", main = "Conditional sensitivity: myeloid atlas")
  axis(1, at = seq_along(MAIN), labels = MAIN)
  legend("topleft", legend = c("5%", "7%", "10%"), pch = c(15, 16, 17), lty = 1, bty = "n")
  mtext("Newly rescued singlets are NOT assigned labels.", side = 1, line = 4, cex = 0.75)
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, paste("Plot error:", conditionMessage(e)))
  logline("Figure generation error: ", conditionMessage(e))
}, finally = grDevices::dev.off())

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "11_140_REPORT_v01.txt"), useBytes = TRUE)
utils::capture.output(sessionInfo(), file = file.path(OUT, "12_140_sessionInfo_v01.txt"))
all_files <- list.files(OUT, pattern = "\\.(csv|gz|txt|pdf|log)$", recursive = TRUE, full.names = TRUE)
file_idx <- data.frame(rel_path = substring(all_files, nchar(OUT) + 2), size_bytes = as.numeric(file.info(all_files)$size),
  stringsAsFactors = FALSE)
write_tab(file_idx, "00_140_OUTPUT_INDEX_v01.csv")
logline("FINISHED: ", Sys.time())
logline("Read tables/09_140_FINAL_QC_v01.csv and 11_140_REPORT_v01.txt first.")

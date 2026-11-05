source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
set.seed(20260930)

ROOT <- PROJECT_ROOT

STEP45_DIR <- file.path(ROOT, "results", "45_Bulk_batch_structure_and_analysis_object")

INPUT_RDS <- file.path(STEP45_DIR, "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "49_Bulk_WGCNA_network_modules")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "49_Bulk_WGCNA_network_modules_v01"

EXPECTED_SAMPLES <- 39L
EXPECTED_GENES <- 12977L

TIME_LEVELS <- c("Uninjured", "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi")

PHASE_LEVELS <- c("Early", "Transition", "Established late", "Extended late")

CORE_GENES <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
COR_TYPE <- "bicor"
MAX_P_OUTLIERS <- 0.1

POWER_VECTOR <- c(1:20, 22, 24, 26, 28, 30)
SFT_R2_TARGET <- 0.8
SFT_R2_FALLBACK <- 0.7

MIN_MODULE_SIZE <- 30L
DEEP_SPLIT <- 2L
MERGE_CUT_HEIGHT <- 0.25
MAX_BLOCK_SIZE <- 5000L

SAMPLE_CONNECTIVITY_Z_REVIEW <- -2.5

install_cran_if_missing <- require_cran_package

install_bioc_if_missing <- require_bioc_package

for (p in c("impute", "preprocessCore", "GO.db", "AnnotationDbi")) {
  install_bioc_if_missing(p)
}
install_cran_if_missing("WGCNA")

suppressPackageStartupMessages({
  library(WGCNA)
})

n_cores <- parallel::detectCores(logical = TRUE)
if (is.na(n_cores) || n_cores < 2L) {
  N_THREADS <- 1L
} else {
  N_THREADS <- min(8L, max(1L, n_cores - 1L))
}

try(WGCNA::enableWGCNAThreads(nThreads = N_THREADS), silent = TRUE)

write_csv <- function(x, filename) {
  utils::write.csv(x, file = file.path(OUT_DIR, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

safe_gene_symbol <- clean_gene_symbol

as_chr <- function(x) {
  as.character(x)
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), NA_character_, format(signif(x, digits), scientific = TRUE, trim = TRUE))
}

if (!file.exists(INPUT_RDS)) {
  stop("Step45 object not found:\n", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

required_fields <- c("expression_original", "expression_batch_adjusted_visualization", "metadata", "sample_outlier_diagnostics")

missing_fields <- setdiff(required_fields, names(obj))

if (length(missing_fields) > 0L) {
  stop("Step45 analysis object is missing required field(s): ", paste(missing_fields, collapse = ", "))
}

expr <- obj$expression_original
meta <- obj$metadata
step45_outlier <- obj$sample_outlier_diagnostics

if (!is.matrix(expr)) {
  expr <- as.matrix(expr)
}

storage.mode(expr) <- "double"

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("expression_original must have gene row names and sample column names.")
}

if (nrow(expr) != EXPECTED_GENES) {
  stop("Unexpected gene count in Step45 expression_original: ", nrow(expr), " ; expected ", EXPECTED_GENES)
}

if (ncol(expr) != EXPECTED_SAMPLES) {
  stop("Unexpected sample count in Step45 expression_original: ", ncol(expr), " ; expected ", EXPECTED_SAMPLES)
}

if (any(!is.finite(expr))) {
  stop("expression_original contains non-finite values. ", "Step49 will not silently impute or remove samples.")
}

if (anyDuplicated(rownames(expr))) {
  stop("Duplicated gene row names found in expression_original.")
}

if (anyDuplicated(colnames(expr))) {
  stop("Duplicated sample IDs found in expression_original.")
}

if (!"sample_id" %in% names(meta)) {
  stop("Step45 metadata lacks sample_id.")
}

if (!all(colnames(expr) %in% meta$sample_id)) {
  stop("Not all expression samples are present in Step45 metadata.")
}

meta <- meta[match(colnames(expr), meta$sample_id), , drop = FALSE]

rownames(meta) <- meta$sample_id

if (!identical(meta$sample_id, colnames(expr))) {
  stop("Metadata/sample order alignment failed.")
}

for (nm in c("batch", "timepoint", "phase")) {
  if (!nm %in% names(meta)) {
    stop("Step45 metadata lacks required column: ", nm)
  }
}

observed_times <- unique(as_chr(meta$timepoint))
if (!all(TIME_LEVELS %in% observed_times)) {
  stop("Not all time points are present. Observed: ", paste(observed_times, collapse = ", "))
}

review_flag_by_sample <- rep(FALSE, nrow(meta))
names(review_flag_by_sample) <- meta$sample_id

if (is.data.frame(step45_outlier) && all(c("sample_id", "review_flag") %in% names(step45_outlier))) {
    idx_flag <- match(meta$sample_id, step45_outlier$sample_id)
    found_flag <- !is.na(idx_flag)
    review_flag_by_sample[found_flag] <- as.logical(step45_outlier$review_flag[idx_flag[found_flag]])
}

datExpr0 <- t(expr)

gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)

if (!all(gsg$goodSamples)) {
    bad_samples <- rownames(datExpr0)[!gsg$goodSamples]
    write_csv(data.frame(sample_id = bad_samples, stringsAsFactors = FALSE), paste0("01A_", PREFIX, "_BAD_samples_STOP.csv"))
    stop("WGCNA goodSamplesGenes flagged sample(s): ", paste(bad_samples, collapse = ", "), "\nNo sample was removed automatically. Resolve upstream before rerunning.")
}

bad_genes <- colnames(datExpr0)[!gsg$goodGenes]

if (length(bad_genes) > 0L) {
  write_csv(data.frame(gene_id = bad_genes, gene_symbol_restored = safe_gene_symbol(bad_genes), reason = "WGCNA_goodSamplesGenes_failed",
    stringsAsFactors = FALSE), paste0("01B_", PREFIX, "_technical_bad_genes_removed.csv"))
}

datExpr <- datExpr0[, gsg$goodGenes, drop = FALSE]

if (nrow(datExpr) != EXPECTED_SAMPLES) {
  stop("Internal error: sample count changed during gene QC.")
}

if (ncol(datExpr) < 5000L) {
  stop("Unexpectedly few genes remain after WGCNA technical QC: ", ncol(datExpr))
}

sample_cor <- WGCNA::bicor(t(datExpr), use = "pairwise.complete.obs", maxPOutliers = MAX_P_OUTLIERS)

diag(sample_cor) <- NA_real_

mean_bicor_to_others <- rowMeans(sample_cor, na.rm = TRUE)

sample_connectivity_z <- as.numeric(scale(mean_bicor_to_others))

sample_tree <- stats::hclust(stats::dist(datExpr), method = "average")

sample_qc <- data.frame(sample_id = rownames(datExpr), batch = as_chr(meta[rownames(datExpr), "batch"]),
  timepoint = as_chr(meta[rownames(datExpr), "timepoint"]), phase = as_chr(meta[rownames(datExpr),
    "phase"]), step45_review_flag = as.logical(review_flag_by_sample[rownames(datExpr)]), mean_bicor_to_other_samples = mean_bicor_to_others,
  sample_connectivity_z = sample_connectivity_z, WGCNA_connectivity_review_flag = sample_connectivity_z <
    SAMPLE_CONNECTIVITY_Z_REVIEW, retained_in_network = TRUE, stringsAsFactors = FALSE)

sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = POWER_VECTOR, RsquaredCut = SFT_R2_TARGET, networkType = NETWORK_TYPE,
  corFnc = WGCNA::bicor, corOptions = list(use = "pairwise.complete.obs", maxPOutliers = MAX_P_OUTLIERS),
  blockSize = min(3000L, ncol(datExpr)), verbose = 5)

fit <- as.data.frame(sft$fitIndices, stringsAsFactors = FALSE)

required_sft_cols <- c("Power", "SFT.R.sq", "slope", "mean.k.")

missing_sft_cols <- setdiff(required_sft_cols, names(fit))

if (length(missing_sft_cols) > 0L) {
  stop("Unexpected pickSoftThreshold output; missing: ", paste(missing_sft_cols, collapse = ", "))
}

fit$signed_SFT_R2 <- -sign(fit$slope) * fit$SFT.R.sq
fit$negative_slope <- is.finite(fit$slope) & fit$slope < 0

cand_target <- which(fit$negative_slope & is.finite(fit$SFT.R.sq) & fit$SFT.R.sq >= SFT_R2_TARGET)

if (length(cand_target) > 0L) {
    chosen_row <- cand_target[which.min(fit$Power[cand_target])]
    selection_rule <- paste0("lowest power with negative slope and SFT.R.sq >= ", SFT_R2_TARGET)
} else {
    cand_fallback <- which(fit$negative_slope & is.finite(fit$SFT.R.sq) & fit$SFT.R.sq >= SFT_R2_FALLBACK)
    if (length(cand_fallback) > 0L) {
        chosen_row <- cand_fallback[which.min(fit$Power[cand_fallback])]
        selection_rule <- paste0("fallback: lowest power with negative slope and SFT.R.sq >= ", SFT_R2_FALLBACK)
    }
    else {
        cand_negative <- which(fit$negative_slope & is.finite(fit$SFT.R.sq))
        if (length(cand_negative) == 0L) {
            stop("No soft-threshold candidate has a negative scale-free slope. ", "Review the Step49 soft-threshold table.")
        }
        chosen_row <- cand_negative[which.max(fit$SFT.R.sq[cand_negative])]
        if (fit$SFT.R.sq[chosen_row] < 0.6) {
            write_csv(fit, paste0("02_", PREFIX, "_soft_threshold_scan.csv"))
            stop("No acceptable soft-threshold power was found. ", "Best negative-slope SFT.R.sq = ", 
                signif(fit$SFT.R.sq[chosen_row], 4), " at power ", fit$Power[chosen_row], ". Review file 02; network construction was NOT forced.")
        }
        selection_rule <- paste0("fallback: maximum SFT.R.sq among negative-slope powers; ", "accepted because best R2 >= 0.60")
    }
}

SOFT_POWER <- as.integer(fit$Power[chosen_row])

soft_selection <- data.frame(selected_power = SOFT_POWER, SFT_R_sq = fit$SFT.R.sq[chosen_row], signed_SFT_R2 = fit$signed_SFT_R2[chosen_row],
  slope = fit$slope[chosen_row], mean_connectivity = fit$mean.k.[chosen_row], selection_rule = selection_rule,
  network_type = NETWORK_TYPE, correlation = COR_TYPE, maxPOutliers = MAX_P_OUTLIERS, stringsAsFactors = FALSE)

tom_base <- file.path(OUT_DIR, paste0(PREFIX, "_TOM"))

net <- WGCNA::blockwiseModules(datExpr, power = SOFT_POWER, maxBlockSize = MAX_BLOCK_SIZE, networkType = NETWORK_TYPE,
  TOMType = TOM_TYPE, corType = COR_TYPE, maxPOutliers = MAX_P_OUTLIERS, deepSplit = DEEP_SPLIT, minModuleSize = MIN_MODULE_SIZE,
  mergeCutHeight = MERGE_CUT_HEIGHT, reassignThreshold = 0, pamRespectsDendro = FALSE, numericLabels = TRUE,
  saveTOMs = TRUE, saveTOMFileBase = tom_base, randomSeed = 20260930, verbose = 5)

module_labels <- net$colors
module_colors <- WGCNA::labels2colors(module_labels)

names(module_labels) <- colnames(datExpr)
names(module_colors) <- colnames(datExpr)

ME_list <- WGCNA::moduleEigengenes(datExpr, colors = module_colors, excludeGrey = FALSE)

MEs <- WGCNA::orderMEs(ME_list$eigengenes)

if (!identical(rownames(MEs), rownames(datExpr))) {
  MEs <- MEs[rownames(datExpr), , drop = FALSE]
}

kME_mat <- WGCNA::bicor(datExpr, MEs, use = "pairwise.complete.obs", maxPOutliers = MAX_P_OUTLIERS)

if (is.null(dim(kME_mat))) {
  stop("Unexpected scalar/vector result while computing kME.")
}

rownames(kME_mat) <- colnames(datExpr)

gene_ids <- colnames(datExpr)
gene_symbols <- safe_gene_symbol(gene_ids)

own_kME <- rep(NA_real_, length(gene_ids))
names(own_kME) <- gene_ids

for (i in seq_along(gene_ids)) {
    me_nm <- paste0("ME", module_colors[i])
    if (me_nm %in% colnames(kME_mat)) {
        own_kME[i] <- kME_mat[gene_ids[i], me_nm]
    }
}

gene_module <- data.frame(gene_id = gene_ids, gene_symbol_restored = gene_symbols, WGCNA_block = if (!is.null(net$blocks)) {
  as.integer(net$blocks)
} else {
  NA_integer_
}, module_label = as.integer(module_labels), module_color = module_colors, assigned_module = module_colors !=
  "grey", kME_own_module = as.numeric(own_kME), stringsAsFactors = FALSE)

module_size_tab <- as.data.frame(table(module_color = module_colors), stringsAsFactors = FALSE)

names(module_size_tab)[2] <- "n_genes"

module_size_tab$assigned_module <- module_size_tab$module_color != "grey"

label_map <- unique(gene_module[, c("module_color", "module_label"), drop = FALSE])

module_size_tab <- merge(module_size_tab, label_map, by = "module_color", all.x = TRUE, sort = FALSE)

module_size_tab <- module_size_tab[order(module_size_tab$assigned_module, module_size_tab$n_genes, decreasing = TRUE),
  , drop = FALSE]

ME_table <- data.frame(sample_id = rownames(MEs), batch = as_chr(meta[rownames(MEs), "batch"]), timepoint = as_chr(meta[rownames(MEs),
  "timepoint"]), phase = as_chr(meta[rownames(MEs), "phase"]), MEs, check.names = FALSE, stringsAsFactors = FALSE)

core_idx <- match(CORE_GENES, gene_module$gene_symbol_restored)

core_module <- data.frame(requested_gene = CORE_GENES, found_in_WGCNA_gene_set = !is.na(core_idx), gene_id = NA_character_,
  gene_symbol_restored = NA_character_, module_label = NA_integer_, module_color = NA_character_, kME_own_module = NA_real_,
  stringsAsFactors = FALSE)

found <- which(!is.na(core_idx))

if (length(found) > 0L) {
    core_module$gene_id[found] <- gene_module$gene_id[core_idx[found]]
    core_module$gene_symbol_restored[found] <- gene_module$gene_symbol_restored[core_idx[found]]
    core_module$module_label[found] <- gene_module$module_label[core_idx[found]]
    core_module$module_color[found] <- gene_module$module_color[core_idx[found]]
    core_module$kME_own_module[found] <- gene_module$kME_own_module[core_idx[found]]
}

n_modules_non_grey <- sum(module_size_tab$module_color != "grey")

n_genes_non_grey <- sum(gene_module$module_color != "grey")

n_grey <- sum(gene_module$module_color == "grey")

input_qc <- data.frame(metric = c("input_rds", "formal_expression_field", "batch_adjusted_field_used_for_primary_network",
  "expected_samples", "observed_samples", "expected_input_genes", "observed_input_genes", "technical_bad_genes_removed_n",
  "genes_used_for_WGCNA_n", "Step45_review_flagged_samples_retained_n", "WGCNA_sample_connectivity_review_flag_n",
  "network_type", "TOM_type", "correlation", "maxPOutliers", "soft_power", "minModuleSize", "deepSplit",
  "mergeCutHeight", "maxBlockSize", "non_grey_module_n", "genes_assigned_to_non_grey_modules_n", "grey_unassigned_gene_n",
  "automatic_sample_removal"), value = c(INPUT_RDS, "expression_original", "NO", EXPECTED_SAMPLES,
  nrow(datExpr), EXPECTED_GENES, ncol(datExpr0), length(bad_genes), ncol(datExpr), sum(review_flag_by_sample,
    na.rm = TRUE), sum(sample_qc$WGCNA_connectivity_review_flag, na.rm = TRUE), NETWORK_TYPE, TOM_TYPE,
  COR_TYPE, MAX_P_OUTLIERS, SOFT_POWER, MIN_MODULE_SIZE, DEEP_SPLIT, MERGE_CUT_HEIGHT, MAX_BLOCK_SIZE,
  n_modules_non_grey, n_genes_non_grey, n_grey, "NONE"), stringsAsFactors = FALSE)

files <- list(input_qc = paste0("01_", PREFIX, "_input_and_network_QC.csv"), soft_scan = paste0("02_",
  PREFIX, "_soft_threshold_scan.csv"), soft_selection = paste0("03_", PREFIX, "_soft_threshold_selection.csv"),
  gene_module = paste0("04_", PREFIX, "_gene_module_assignment.csv"), module_size = paste0("05_", PREFIX,
    "_module_size_summary.csv"), eigengenes = paste0("06_", PREFIX, "_module_eigengenes.csv"), sample_qc = paste0("07_",
    PREFIX, "_sample_network_QC.csv"), core_module = paste0("08_", PREFIX, "_core_gene_module_audit.csv"),
  analysis_object = paste0("09_", PREFIX, "_analysis_object.rds"), report = paste0("10_", PREFIX, "_REPORT.txt"),
  diagnostics = paste0("11_", PREFIX, "_diagnostics.pdf"), session = paste0("12_", PREFIX, "_sessionInfo.txt"))

write_csv(input_qc, files$input_qc)

write_csv(fit, files$soft_scan)

write_csv(soft_selection, files$soft_selection)

write_csv(gene_module, files$gene_module)

write_csv(module_size_tab, files$module_size)

write_csv(ME_table, files$eigengenes)

write_csv(sample_qc, files$sample_qc)

write_csv(core_module, files$core_module)

analysis_object <- list(step = 49, version = "v01", analysis = "Bulk WGCNA network construction and module detection",
  upstream_step45_file = INPUT_RDS, frozen_design = list(all_39_samples_retained = TRUE, all_8_timepoints_retained = TRUE,
    time_levels = TIME_LEVELS, phase_levels = PHASE_LEVELS, primary_network_expression = "Step45 expression_original",
    batch_adjusted_visualization_used = FALSE, no_DEG_preselection = TRUE, no_biological_module_selection_in_step49 = TRUE),
  settings = list(networkType = NETWORK_TYPE, TOMType = TOM_TYPE, corType = COR_TYPE, maxPOutliers = MAX_P_OUTLIERS,
    powerVector = POWER_VECTOR, softPower = SOFT_POWER, softThresholdSelectionRule = selection_rule,
    SFT_R2_target = SFT_R2_TARGET, SFT_R2_fallback = SFT_R2_FALLBACK, minModuleSize = MIN_MODULE_SIZE,
    deepSplit = DEEP_SPLIT, mergeCutHeight = MERGE_CUT_HEIGHT, maxBlockSize = MAX_BLOCK_SIZE, nThreads_requested = N_THREADS),
  metadata = meta, sample_qc = sample_qc, step45_review_flag_by_sample = review_flag_by_sample, technical_bad_genes_removed = bad_genes,
  genes_used = colnames(datExpr), gene_symbol_restored = safe_gene_symbol(colnames(datExpr)), soft_threshold_scan = fit,
  soft_threshold_selection = soft_selection, net = net, module_labels = module_labels, module_colors = module_colors,
  module_eigengenes = MEs, module_membership_bicor = kME_mat, gene_module_assignment = gene_module,
  module_size_summary = module_size_tab, core_gene_module_audit = core_module, TOM_file_base = tom_base,
  notes = c("Step49 constructs the co-expression network only.", "No module is selected as biologically key in Step49.",
    "Module-trait association is deferred to Step50.", "All samples are retained; review flags are audit information only.",
    "The primary network uses Step45 expression_original.", "The Step45 batch-adjusted visualization matrix is not used for the primary network.",
    "No DEG-based gene preselection is applied.", "Grey genes are unassigned, not a biological module."))

saveRDS(analysis_object, file = file.path(OUT_DIR, files$analysis_object), compress = "xz")

grDevices::pdf(file.path(OUT_DIR, files$diagnostics), width = 10, height = 7.5, onefile = TRUE, family = "sans")

plot(fit$Power, fit$signed_SFT_R2, xlab = "Soft-threshold power", ylab = "Signed scale-free topology fit (R^2)",
  main = "Step49 WGCNA soft-threshold scan", type = "n")

text(fit$Power, fit$signed_SFT_R2, labels = fit$Power, cex = 0.8)

abline(h = SFT_R2_TARGET, lty = 2)

abline(v = SOFT_POWER, lty = 3)

plot(fit$Power, fit$mean.k., xlab = "Soft-threshold power", ylab = "Mean connectivity", main = "Step49 mean connectivity",
  type = "b", pch = 16)

abline(v = SOFT_POWER, lty = 3)

plot(sample_tree, main = "Step49 sample clustering — all samples retained", xlab = "", sub = "", cex = 0.7)

if (!is.null(net$dendrograms) && !is.null(net$blockGenes)) {
    for (b in seq_along(net$dendrograms)) {
        block_gene_idx <- net$blockGenes[[b]]
        WGCNA::plotDendroAndColors(net$dendrograms[[b]], module_colors[block_gene_idx], groupLabels = "Final module", 
            main = paste0("Step49 gene dendrogram + module colors — block ", b), dendroLabels = FALSE, 
            hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    }
}

grDevices::dev.off()

module_summary_text <- paste0(module_size_tab$module_color, "=", module_size_tab$n_genes, collapse = "; ")

core_text <- paste0(core_module$requested_gene, ":", ifelse(core_module$found_in_WGCNA_gene_set, paste0(core_module$module_color,
  " (kME=", fmt_num(core_module$kME_own_module, 3), ")"), "not_found"), collapse = "; ")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = file.path(OUT_DIR, files$report), useBytes = TRUE)

sink(file.path(OUT_DIR, files$session))
print(sessionInfo())
sink()

tom_files <- list.files(OUT_DIR, pattern = paste0("^", PREFIX, "_TOM.*\\.RData$"), full.names = FALSE)

output_index <- data.frame(file = c(paste0("00_", PREFIX, "_OUTPUT_INDEX.csv"), unlist(files, use.names = FALSE),
  tom_files), description = c("Index of Step49 outputs", "Input/network QC and analysis settings",
  "Complete soft-threshold scan", "Selected soft-threshold power and explicit selection rule", "Complete gene-to-module assignment with own-module kME",
  "Module sizes including grey/unassigned genes", "Module eigengenes with sample metadata; no trait correlation performed here",
  "Sample-level network QC; all samples retained", "Clec7a/Trem2/Apoe/Lpl/Apoc4/Abca1 module audit only",
  "Complete Step49 WGCNA analysis object for Step50", "Plain-text Step49 report", "Diagnostic plots only; not a manuscript figure",
  "R sessionInfo", rep("Saved WGCNA topological overlap matrix block for reproducibility", length(tom_files))),
  stringsAsFactors = FALSE)

write_csv(output_index, paste0("00_", PREFIX, "_OUTPUT_INDEX.csv"))

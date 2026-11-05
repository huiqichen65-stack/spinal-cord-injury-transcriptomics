rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
set.seed(20261030)

ROOT <- PROJECT_ROOT

STEP45_RDS <- file.path(ROOT, "results", "45_Bulk_batch_structure_and_analysis_object", "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds")

STEP49_RDS <- file.path(ROOT, "results", "49_Bulk_WGCNA_network_modules", "09_49_Bulk_WGCNA_network_modules_v01_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "102_WGCNA_preprocessing_sensitivity")

PREFIX <- "102_WGCNA_preprocessing_sensitivity_v01"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(STEP45_RDS)) {
  stop("Required Step45 RDS not found:\n", STEP45_RDS)
}
if (!file.exists(STEP49_RDS)) {
  stop("Required Step49 RDS not found:\n", STEP49_RDS)
}

install_cran_if_missing <- require_cran_package

install_bioc_if_missing <- require_bioc_package

install_cran_if_missing("WGCNA")
install_bioc_if_missing("limma")

suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
})

n_cores <- parallel::detectCores(logical = TRUE)
if (is.na(n_cores) || n_cores < 2L) {
  N_THREADS <- 1L
} else {
  N_THREADS <- min(8L, max(1L, n_cores - 1L))
}

try(WGCNA::enableWGCNAThreads(nThreads = N_THREADS), silent = TRUE)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0(..., collapse = ""), "\n",
    sep = "")
}

write_csv <- function(x, filename) {
  utils::write.csv(x, file = file.path(OUT_DIR, filename), row.names = FALSE, na = "")
}

as_chr <- function(x) {
  as.character(x)
}

get_setting <- function(x, name, default) {
  if (is.null(x) || is.null(x[[name]]) || length(x[[name]]) == 0L) {
    return(default)
  }
  x[[name]][1]
}

fmt_num <- function(x, digits = 4L) {
  if (length(x) == 0L || is.na(x[1]))
    return("NA")
  formatC(as.numeric(x[1]), digits = digits, format = "fg", flag = "#")
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L)
    return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

obj45 <- readRDS(STEP45_RDS)
obj49 <- readRDS(STEP49_RDS)

required45 <- c("expression_original", "expression_batch_adjusted_visualization", "metadata", "sample_outlier_diagnostics",
  "design_time_only")

missing45 <- setdiff(required45, names(obj45))
if (length(missing45) > 0L) {
  stop("Step45 object lacks required field(s): ", paste(missing45, collapse = ", "))
}

required49 <- c("metadata", "module_eigengenes", "gene_module_assignment", "settings", "soft_threshold_selection",
  "sample_qc", "step45_review_flag_by_sample")

missing49 <- setdiff(required49, names(obj49))
if (length(missing49) > 0L) {
  stop("Step49 object lacks required field(s): ", paste(missing49, collapse = ", "))
}

expr_original <- as.matrix(obj45$expression_original)
expr_batch_adjusted <- as.matrix(obj45$expression_batch_adjusted_visualization)

storage.mode(expr_original) <- "double"
storage.mode(expr_batch_adjusted) <- "double"

if (is.null(rownames(expr_original)) || is.null(colnames(expr_original)) || is.null(rownames(expr_batch_adjusted)) ||
  is.null(colnames(expr_batch_adjusted))) {
  stop("Step45 expression matrices must have gene and sample dimnames.")
}

if (!identical(dim(expr_original), dim(expr_batch_adjusted))) {
  stop("Original and batch-adjusted matrices have different dimensions.")
}
if (!identical(rownames(expr_original), rownames(expr_batch_adjusted)) || !identical(colnames(expr_original),
  colnames(expr_batch_adjusted))) {
  stop("Original and batch-adjusted matrix dimnames are not identical.")
}
if (any(!is.finite(expr_original))) {
  stop("Step45 expression_original contains non-finite values.")
}
if (any(!is.finite(expr_batch_adjusted))) {
  stop("Step45 batch-adjusted expression contains non-finite values.")
}
if (anyDuplicated(rownames(expr_original))) {
  stop("Duplicated gene IDs in Step45 expression matrices.")
}
if (anyDuplicated(colnames(expr_original))) {
  stop("Duplicated sample IDs in Step45 expression matrices.")
}

meta <- as.data.frame(obj45$metadata, check.names = FALSE)
if (!"sample_id" %in% names(meta)) {
  stop("Step45 metadata lacks sample_id.")
}
if (!all(c("batch", "timepoint", "phase_primary") %in% names(meta))) {
  stop("Step45 metadata lacks batch, timepoint, or phase_primary.")
}

sample_ids <- colnames(expr_original)
if (!all(sample_ids %in% meta$sample_id)) {
  stop("Some expression sample IDs are absent from Step45 metadata.")
}

meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]
rownames(meta) <- meta$sample_id

if (!identical(meta$sample_id, sample_ids)) {
  stop("Step45 metadata could not be aligned exactly to expression columns.")
}
if (anyNA(meta$batch) || anyNA(meta$timepoint)) {
  stop("Missing batch/timepoint values in Step45 metadata.")
}

meta$batch <- factor(as_chr(meta$batch), levels = c("Batch_A", "Batch_B"))
meta$timepoint <- factor(as_chr(meta$timepoint), levels = c("Uninjured", "15 min", "1 dpi", "3 dpi",
  "7 dpi", "14 dpi", "28 dpi", "42 dpi"))
meta$phase_primary <- as_chr(meta$phase_primary)

if (anyNA(meta$batch) || anyNA(meta$timepoint)) {
  stop("Unexpected batch or timepoint level in Step45 metadata.")
}

primary_assignment <- as.data.frame(obj49$gene_module_assignment, stringsAsFactors = FALSE)
if (!all(c("gene_id", "module_color") %in% names(primary_assignment))) {
  stop("Step49 gene_module_assignment lacks gene_id/module_color.")
}

primary_MEs <- as.data.frame(obj49$module_eigengenes, check.names = FALSE)

if (is.null(rownames(primary_MEs))) {
  stop("Step49 module_eigengenes must have sample IDs as row names.")
}
if (!all(sample_ids %in% rownames(primary_MEs))) {
  stop("Step49 module eigengenes do not contain all 39 Step45 samples.")
}
primary_MEs <- primary_MEs[sample_ids, , drop = FALSE]

if (!all(sample_ids %in% obj49$metadata$sample_id)) {
  stop("Step49 metadata sample IDs do not match Step45 expression samples.")
}

if (!setequal(primary_assignment$gene_id, rownames(expr_original))) {
  stop("Step49 gene assignment does not match Step45 expression gene universe.")
}

primary_assignment <- primary_assignment[match(rownames(expr_original), primary_assignment$gene_id),
  , drop = FALSE]

if (anyNA(primary_assignment$module_color)) {
  stop("Missing module_color values in Step49 gene assignment.")
}

settings <- obj49$settings
SOFT_POWER <- as.integer(get_setting(obj49$soft_threshold_selection, "selected_power", get_setting(settings,
  "softPower", NA_integer_)))

if (is.na(SOFT_POWER)) {
  SOFT_POWER <- as.integer(get_setting(settings, "softPower", NA_integer_))
}
if (is.na(SOFT_POWER)) {
  stop("Could not recover the Step49 soft-threshold power.")
}

NETWORK_TYPE <- as.character(get_setting(settings, "networkType", "signed"))
TOM_TYPE <- as.character(get_setting(settings, "TOMType", "signed"))
COR_TYPE <- as.character(get_setting(settings, "corType", "bicor"))
MAX_P_OUTLIERS <- as.numeric(get_setting(settings, "maxPOutliers", 0.1))
MIN_MODULE_SIZE <- as.integer(get_setting(settings, "minModuleSize", 30L))
DEEP_SPLIT <- as.integer(get_setting(settings, "deepSplit", 2L))
MERGE_CUT_HEIGHT <- as.numeric(get_setting(settings, "mergeCutHeight", 0.25))
MAX_BLOCK_SIZE <- as.integer(get_setting(settings, "maxBlockSize", 5000L))

if (NETWORK_TYPE != "signed" || TOM_TYPE != "signed" || COR_TYPE != "bicor") {
  stop("Unexpected Step49 network settings. ", "This sensitivity script expects signed/signed/bicor.")
}

step44_cols <- grep("^step44_", names(meta), value = TRUE)

possible_source_batch_fields <- step44_cols[grepl("batch|run|lane|library|plate|seq", step44_cols, ignore.case = TRUE)]

official_batch_field_found <- length(possible_source_batch_fields) > 0L

batch_provenance <- data.frame(item = c("Step45 batch field", "Operational Batch_A rule", "Operational Batch_B rule",
  "Independent Step44 batch field found", "Step44 candidate batch-related fields", "Interpretation"),
  value = c("metadata$batch", "Sample IDs beginning jerry. or normal.", "Sample IDs beginning D", ifelse(official_batch_field_found,
    "YES", "NO"), ifelse(official_batch_field_found, paste(possible_source_batch_fields, collapse = "; "),
    "None detected among Step44 metadata columns"), paste("Batch_A/Batch_B are operational labels inherited from Step45;",
    "experimental batch identity requires confirmation from the source", "records or laboratory metadata.")),
  stringsAsFactors = FALSE)

write_csv(batch_provenance, paste0("01_", PREFIX, "_batch_provenance_audit.csv"))

batch_time_table <- as.data.frame(as.table(table(batch = meta$batch, timepoint = meta$timepoint)), stringsAsFactors = FALSE)
names(batch_time_table)[3] <- "n_samples"
write_csv(batch_time_table, paste0("02_", PREFIX, "_batch_by_timepoint_counts.csv"))

timepoint_batch_presence <- do.call(rbind, lapply(levels(meta$timepoint), function(tp) {
  z <- table(factor(meta$batch, levels = c("Batch_A", "Batch_B"))[meta$timepoint == tp])
  data.frame(timepoint = tp, n_Batch_A = as.integer(z["Batch_A"]), n_Batch_B = as.integer(z["Batch_B"]),
    both_batches_present = all(z > 0L), stringsAsFactors = FALSE)
}))

write_csv(timepoint_batch_presence, paste0("03_", PREFIX, "_shared_batch_timepoints.csv"))

batch_design <- model.matrix(~timepoint + batch, data = meta)
batch_design_rank <- qr(batch_design)$rank

batch_design_audit <- data.frame(model = "~ timepoint + batch", n_samples = nrow(meta), n_parameters = ncol(batch_design),
  rank = batch_design_rank, full_rank = batch_design_rank == ncol(batch_design), shared_timepoints_with_both_batches = sum(timepoint_batch_presence$both_batches_present),
  total_timepoints = nrow(timepoint_batch_presence), stringsAsFactors = FALSE)

write_csv(batch_design_audit, paste0("04_", PREFIX, "_batch_design_rank_audit.csv"))

if (batch_design_rank < ncol(batch_design)) {
  stop("The model ~ timepoint + batch is not full rank. ", "Batch-adjusted module association tests cannot proceed safely.")
}

step45_outlier <- as.data.frame(obj45$sample_outlier_diagnostics, stringsAsFactors = FALSE)
if (!all(c("sample_id", "review_flag") %in% names(step45_outlier))) {
  stop("Step45 sample_outlier_diagnostics lacks sample_id/review_flag.")
}

step49_sample_qc <- as.data.frame(obj49$sample_qc, stringsAsFactors = FALSE)
if (!"sample_id" %in% names(step49_sample_qc)) {
  stop("Step49 sample_qc lacks sample_id.")
}

sample_outlier_audit <- merge(step45_outlier, step49_sample_qc, by = "sample_id", all = TRUE, suffixes = c("_Step45",
  "_Step49"), sort = FALSE)

sample_outlier_audit <- sample_outlier_audit[match(sample_ids, sample_outlier_audit$sample_id), , drop = FALSE]

write_csv(sample_outlier_audit, paste0("05_", PREFIX, "_sample_outlier_audit.csv"))

flag_idx <- match(sample_ids, step45_outlier$sample_id)
review_flag <- rep(FALSE, length(sample_ids))
review_flag[!is.na(flag_idx)] <- as.logical(step45_outlier$review_flag[flag_idx[!is.na(flag_idx)]])
review_flag[is.na(review_flag)] <- FALSE
names(review_flag) <- sample_ids

flagged_samples <- names(review_flag)[review_flag]
retained_samples_after_flag_exclusion <- names(review_flag)[!review_flag]

flag49 <- obj49$step45_review_flag_by_sample
if (!is.null(names(flag49))) {
  flag49_aligned <- as.logical(flag49[sample_ids])
  flag49_aligned[is.na(flag49_aligned)] <- FALSE
  if (!identical(unname(flag49_aligned), unname(review_flag))) {
    stop("Step45 and Step49 review-flag definitions do not match. ", "Resolve the upstream sample audit before sensitivity analysis.")
  }
}

sample_set_audit <- data.frame(sample_id = sample_ids, batch = as_chr(meta$batch), timepoint = as_chr(meta$timepoint),
  Step45_review_flag = as.logical(review_flag), used_primary_network = TRUE, used_batch_adjusted_network = TRUE,
  used_outlier_excluded_network = !review_flag, stringsAsFactors = FALSE)

write_csv(sample_set_audit, paste0("06_", PREFIX, "_sample_set_by_network.csv"))

msg("Step45 review-flagged samples: ", length(flagged_samples), if (length(flagged_samples) > 0L) {
  paste0(" [", paste(flagged_samples, collapse = ", "), "]")
} else {
  ""
})

if (length(retained_samples_after_flag_exclusion) < 30L) {
  stop("Fewer than 30 samples remain after excluding Step45 review flags. ", "Outlier-excluded network will not be constructed.")
}

datExpr_adjusted_for_scan <- t(expr_batch_adjusted)

power_vector <- c(1:20, 22, 24, 26, 28, 30)

sft_adjusted <- WGCNA::pickSoftThreshold(datExpr_adjusted_for_scan, powerVector = power_vector, RsquaredCut = 0.8,
  networkType = NETWORK_TYPE, corFnc = WGCNA::bicor, corOptions = list(use = "pairwise.complete.obs",
    maxPOutliers = MAX_P_OUTLIERS), blockSize = min(3000L, ncol(datExpr_adjusted_for_scan)), verbose = 3)

sft_adjusted_table <- as.data.frame(sft_adjusted$fitIndices, stringsAsFactors = FALSE)

if (!"Power" %in% names(sft_adjusted_table)) {
  stop("Unexpected pickSoftThreshold output for batch-adjusted matrix.")
}

if (!"SFT.R.sq" %in% names(sft_adjusted_table)) {
  sft_adjusted_table$SFT.R.sq <- NA_real_
}
if (!"slope" %in% names(sft_adjusted_table)) {
  sft_adjusted_table$slope <- NA_real_
}
if (!"mean.k." %in% names(sft_adjusted_table)) {
  sft_adjusted_table$mean.k. <- NA_real_
}

sft_adjusted_table$negative_slope <- is.finite(sft_adjusted_table$slope) & sft_adjusted_table$slope <
  0

select_power_same_rule <- function(fit) {
    ix <- which(fit$negative_slope & is.finite(fit$SFT.R.sq) & fit$SFT.R.sq >= 0.8)
    if (length(ix) > 0L) {
        j <- ix[which.min(fit$Power[ix])]
        return(list(power = as.integer(fit$Power[j]), rule = "lowest power with negative slope and SFT.R.sq >= 0.80"))
    }
    ix <- which(fit$negative_slope & is.finite(fit$SFT.R.sq) & fit$SFT.R.sq >= 0.7)
    if (length(ix) > 0L) {
        j <- ix[which.min(fit$Power[ix])]
        return(list(power = as.integer(fit$Power[j]), rule = "fallback: lowest power with negative slope and SFT.R.sq >= 0.70"))
    }
    ix <- which(fit$negative_slope & is.finite(fit$SFT.R.sq))
    if (length(ix) > 0L) {
        j <- ix[which.max(fit$SFT.R.sq[ix])]
        if (fit$SFT.R.sq[j] >= 0.6) {
            return(list(power = as.integer(fit$Power[j]), rule = "fallback: maximum SFT.R.sq among negative-slope powers; R2 >= 0.60"))
        }
    }
    list(power = NA_integer_, rule = "No acceptable power by the Step49 selection rule")
}

adjusted_power_selection <- select_power_same_rule(sft_adjusted_table)

power14_idx <- which(sft_adjusted_table$Power == SOFT_POWER)

power14_r2 <- if (length(power14_idx) == 1L) {
  sft_adjusted_table$SFT.R.sq[power14_idx]
} else {
  NA_real_
}

sft_adjusted_table$selected_by_same_rule <- (sft_adjusted_table$Power == adjusted_power_selection$power)

sft_adjusted_table$fixed_primary_power <- (sft_adjusted_table$Power == SOFT_POWER)

write_csv(sft_adjusted_table, paste0("07_", PREFIX, "_batch_adjusted_soft_threshold_scan.csv"))

soft_threshold_summary <- data.frame(frozen_primary_power = SOFT_POWER, batch_adjusted_matrix_power_selected_by_same_rule = adjusted_power_selection$power,
  batch_adjusted_selection_rule = adjusted_power_selection$rule, batch_adjusted_SFT_R2_at_frozen_primary_power = power14_r2,
  fixed_power_used_for_batch_adjusted_network = SOFT_POWER, interpretation = "Network sensitivity uses the primary power to isolate preprocessing effects; scan is diagnostic only.",
  stringsAsFactors = FALSE)

write_csv(soft_threshold_summary, paste0("08_", PREFIX, "_batch_adjusted_soft_threshold_summary.csv"))

run_wgcna_sensitivity <- function(expression_gene_by_sample, metadata, network_name, seed) {
    expression_gene_by_sample <- as.matrix(expression_gene_by_sample)
    storage.mode(expression_gene_by_sample) <- "double"
    datExpr0 <- t(expression_gene_by_sample)
    if (!identical(rownames(datExpr0), metadata$sample_id)) {
        stop("Sample order mismatch for network ", network_name)
    }
    gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)
    if (!all(gsg$goodSamples)) {
        bad_samples <- rownames(datExpr0)[!gsg$goodSamples]
        stop("WGCNA goodSamplesGenes flagged sample(s) in ", network_name, ": ", paste(bad_samples, collapse = ", "), 
            ". No sample will be removed automatically.")
    }
    bad_genes <- colnames(datExpr0)[!gsg$goodGenes]
    datExpr <- datExpr0[, gsg$goodGenes, drop = FALSE]
    if (ncol(datExpr) < 5000L) {
        stop("Fewer than 5,000 genes remain after WGCNA technical QC in ", network_name)
    }
    net <- WGCNA::blockwiseModules(datExpr, power = SOFT_POWER, maxBlockSize = MAX_BLOCK_SIZE, networkType = NETWORK_TYPE, 
        TOMType = TOM_TYPE, corType = COR_TYPE, maxPOutliers = MAX_P_OUTLIERS, deepSplit = DEEP_SPLIT, 
        minModuleSize = MIN_MODULE_SIZE, mergeCutHeight = MERGE_CUT_HEIGHT, reassignThreshold = 0, pamRespectsDendro = FALSE, 
        numericLabels = TRUE, saveTOMs = FALSE, randomSeed = seed, verbose = 4)
    module_colors <- WGCNA::labels2colors(net$colors)
    names(module_colors) <- colnames(datExpr)
    ME_list <- WGCNA::moduleEigengenes(datExpr, colors = module_colors, excludeGrey = FALSE)
    MEs <- WGCNA::orderMEs(ME_list$eigengenes)
    if (!identical(rownames(MEs), rownames(datExpr))) {
        MEs <- MEs[rownames(datExpr), , drop = FALSE]
    }
    gene_assignment <- data.frame(gene_id = colnames(datExpr), module_color = as.character(module_colors), 
        stringsAsFactors = FALSE)
    module_size <- as.data.frame(table(module_color = as.character(module_colors)), stringsAsFactors = FALSE)
    names(module_size)[2] <- "n_genes"
    list(network_name = network_name, sample_ids = rownames(datExpr), genes_used = colnames(datExpr), 
        technical_bad_genes = bad_genes, net = net, gene_module_assignment = gene_assignment, module_eigengenes = MEs, 
        module_size_summary = module_size, settings = list(softPower = SOFT_POWER, networkType = NETWORK_TYPE, 
            TOMType = TOM_TYPE, corType = COR_TYPE, maxPOutliers = MAX_P_OUTLIERS, minModuleSize = MIN_MODULE_SIZE, 
            deepSplit = DEEP_SPLIT, mergeCutHeight = MERGE_CUT_HEIGHT, maxBlockSize = MAX_BLOCK_SIZE, 
            saveTOMs = FALSE, randomSeed = seed))
}

network_batch_adjusted <- run_wgcna_sensitivity(expression_gene_by_sample = expr_batch_adjusted, metadata = meta,
  network_name = "batch_adjusted_all_samples", seed = 20261031L)

network_outlier_excluded <- NULL

if (length(flagged_samples) > 0L) {
    keep_samples <- retained_samples_after_flag_exclusion
    network_outlier_excluded <- run_wgcna_sensitivity(expression_gene_by_sample = expr_original[, keep_samples, 
        drop = FALSE], metadata = meta[match(keep_samples, meta$sample_id), , drop = FALSE], network_name = "original_expression_excluding_review_flags", 
        seed = 20261032L)
} else {
    msg("No Step45 review-flagged samples; outlier-excluded network skipped.")
}

compare_modules <- function(primary_assignment, alternative_assignment, network_name) {
    common_genes <- intersect(primary_assignment$gene_id, alternative_assignment$gene_id)
    if (length(common_genes) < 5000L) {
        stop("Too few common genes for module comparison in ", network_name)
    }
    pmap <- setNames(as.character(primary_assignment$module_color), primary_assignment$gene_id)
    amap <- setNames(as.character(alternative_assignment$module_color), alternative_assignment$gene_id)
    pmods <- setdiff(unique(pmap[common_genes]), "grey")
    amods <- setdiff(unique(amap[common_genes]), "grey")
    rows <- list()
    for (pm in pmods) {
        pg <- common_genes[pmap[common_genes] == pm]
        candidates <- lapply(amods, function(am) {
            ag <- common_genes[amap[common_genes] == am]
            ov <- length(intersect(pg, ag))
            un <- length(union(pg, ag))
            den <- length(common_genes)
            a <- ov
            b <- length(pg) - ov
            c <- length(ag) - ov
            d <- den - a - b - c
            fisher_p <- tryCatch(stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")$p.value, 
                error = function(e) NA_real_)
            data.frame(network_name = network_name, primary_module = pm, alternative_module = am, primary_module_n_common_universe = length(pg), 
                alternative_module_n_common_universe = length(ag), overlap_n = ov, union_n = un, jaccard = if (un > 
                  0L) 
                  ov/un
                else NA_real_, primary_module_retention = if (length(pg) > 0L) {
                  ov/length(pg)
                }
                else {
                  NA_real_
                }, alternative_module_purity = if (length(ag) > 0L) {
                  ov/length(ag)
                }
                else {
                  NA_real_
                }, overlap_coefficient = if (min(length(pg), length(ag)) > 0L) {
                  ov/min(length(pg), length(ag))
                }
                else {
                  NA_real_
                }, fisher_P = fisher_p, stringsAsFactors = FALSE)
        })
        if (length(candidates) > 0L) {
            candidates <- do.call(rbind, candidates)
            candidates <- candidates[order(candidates$jaccard, candidates$overlap_n, decreasing = TRUE), 
                , drop = FALSE]
            rows[[length(rows) + 1L]] <- candidates[1, , drop = FALSE]
        }
    }
    out <- if (length(rows) > 0L) {
        do.call(rbind, rows)
    }
    else {
        data.frame()
    }
    if (nrow(out) > 0L) {
        out$FDR_across_primary_modules <- stats::p.adjust(out$fisher_P, method = "BH")
    }
    out
}

network_list <- list(batch_adjusted_all_samples = network_batch_adjusted)

if (!is.null(network_outlier_excluded)) {
  network_list$original_expression_excluding_review_flags <- network_outlier_excluded
}

module_overlap_all <- list()
eigengene_cor_all <- list()
brown_gene_membership_all <- list()
core_gene_membership_all <- list()

core_genes <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

primary_symbols <- if ("gene_symbol_restored" %in% names(primary_assignment)) {
  as.character(primary_assignment$gene_symbol_restored)
} else {
  primary_assignment$gene_id
}

primary_assignment$gene_symbol_restored <- primary_symbols

for (nm in names(network_list)) {
  alt <- network_list[[nm]]

  overlap <- compare_modules(primary_assignment, alt$gene_module_assignment, nm)
  module_overlap_all[[nm]] <- overlap

  primary_ME_names <- grep("^ME", colnames(primary_MEs), value = TRUE)
  alt_ME_names <- grep("^ME", colnames(alt$module_eigengenes), value = TRUE)

  common_samples <- intersect(rownames(primary_MEs), rownames(alt$module_eigengenes))

  cor_rows <- list()
  for (pme in primary_ME_names) {
    for (ame in alt_ME_names) {
      cor_rows[[length(cor_rows) + 1L]] <- data.frame(network_name = nm, primary_module = sub("^ME",
        "", pme), alternative_module = sub("^ME", "", ame), n_samples = length(common_samples),
        pearson_r = safe_cor(primary_MEs[common_samples, pme], alt$module_eigengenes[common_samples,
          ame], method = "pearson"), spearman_rho = safe_cor(primary_MEs[common_samples, pme],
          alt$module_eigengenes[common_samples, ame], method = "spearman"), stringsAsFactors = FALSE)
    }
  }
  eigengene_cor_all[[nm]] <- do.call(rbind, cor_rows)

  alt_map <- setNames(as.character(alt$gene_module_assignment$module_color), alt$gene_module_assignment$gene_id)
  brown_ids <- primary_assignment$gene_id[primary_assignment$module_color == "brown"]
  brown_idx <- match(brown_ids, primary_assignment$gene_id)

  brown_gene_membership_all[[nm]] <- data.frame(network_name = nm, gene_id = brown_ids, gene_symbol = primary_assignment$gene_symbol_restored[brown_idx],
    primary_module = "brown", alternative_module = unname(alt_map[brown_ids]), retained_in_alternative_gene_universe = brown_ids %in%
      names(alt_map), stringsAsFactors = FALSE)

  alt_symbol_map <- setNames(as.character(alt$gene_module_assignment$module_color), alt$gene_module_assignment$gene_id)

  for (cg in core_genes) {
    p_hit <- which(primary_assignment$gene_symbol_restored == cg)
    if (length(p_hit) == 0L) {
      primary_id <- NA_character_
      primary_mod <- NA_character_
      alt_mod <- NA_character_
    } else {
      primary_id <- primary_assignment$gene_id[p_hit[1]]
      primary_mod <- primary_assignment$module_color[p_hit[1]]
      alt_mod <- unname(alt_symbol_map[primary_id])
      if (length(alt_mod) == 0L)
        alt_mod <- NA_character_
    }

    core_gene_membership_all[[length(core_gene_membership_all) + 1L]] <- data.frame(network_name = nm,
      core_gene = cg, gene_id = primary_id, primary_module = primary_mod, alternative_module = alt_mod,
      stringsAsFactors = FALSE)
  }
}

module_overlap <- do.call(rbind, module_overlap_all)
eigengene_correlations <- do.call(rbind, eigengene_cor_all)
brown_gene_membership <- do.call(rbind, brown_gene_membership_all)
core_gene_membership <- do.call(rbind, core_gene_membership_all)

primary_gene_membership <- data.frame(network_name = "primary_original_expression_all_samples", gene_id = primary_assignment$gene_id,
  gene_symbol = primary_assignment$gene_symbol_restored, module_color = as.character(primary_assignment$module_color),
  stringsAsFactors = FALSE)

alternative_gene_membership <- list()
symbol_by_gene <- setNames(primary_assignment$gene_symbol_restored, primary_assignment$gene_id)

for (nm in names(network_list)) {
  ga <- network_list[[nm]]$gene_module_assignment
  alternative_gene_membership[[nm]] <- data.frame(network_name = nm, gene_id = ga$gene_id, gene_symbol = unname(symbol_by_gene[ga$gene_id]),
    module_color = as.character(ga$module_color), stringsAsFactors = FALSE)
}

gene_membership_all_networks <- do.call(rbind, c(list(primary_gene_membership), alternative_gene_membership))

network_input_summary <- data.frame(network_name = c("primary_original_expression_all_samples", names(network_list)),
  expression_source = c("Step45 expression_original", ifelse(names(network_list) == "batch_adjusted_all_samples",
    "Step45 expression_batch_adjusted_visualization", "Step45 expression_original")), n_samples = c(length(sample_ids),
    vapply(network_list, function(x) length(x$sample_ids), integer(1))), n_genes_used = c(nrow(primary_assignment),
    vapply(network_list, function(x) length(x$genes_used), integer(1))), stringsAsFactors = FALSE)

write_csv(network_input_summary, paste0("17_", PREFIX, "_network_input_sample_gene_summary.csv"))
write_csv(gene_membership_all_networks, paste0("18_", PREFIX, "_exact_gene_module_membership_all_networks.csv"))

write_csv(module_overlap, paste0("09_", PREFIX, "_module_overlap_best_matches.csv"))
write_csv(eigengene_correlations, paste0("10_", PREFIX, "_module_eigengene_correlations_all_pairs.csv"))
write_csv(brown_gene_membership, paste0("11_", PREFIX, "_primary_brown_gene_membership_sensitivity.csv"))
write_csv(core_gene_membership, paste0("12_", PREFIX, "_six_core_gene_module_membership.csv"))

run_batch_association <- function(MEs, metadata, network_name) {
    m <- metadata[match(rownames(MEs), metadata$sample_id), , drop = FALSE]
    if (anyNA(m$sample_id)) {
        stop("Could not align metadata for batch association: ", network_name)
    }
    b <- factor(as_chr(m$batch), levels = c("Batch_A", "Batch_B"))
    tp <- factor(as_chr(m$timepoint), levels = levels(meta$timepoint))
    design <- model.matrix(~tp + b)
    design_rank <- qr(design)$rank
    full_rank <- design_rank == ncol(design)
    audit <- data.frame(network_name = network_name, model = "~ timepoint + batch", n_samples = nrow(m), 
        n_parameters = ncol(design), rank = design_rank, full_rank = full_rank, stringsAsFactors = FALSE)
    if (!full_rank) {
        stop("Batch association design is not full rank for ", network_name)
    }
    coef_name <- "bBatch_B"
    if (!coef_name %in% colnames(design)) {
        stop("Could not identify bBatch_B coefficient for ", network_name)
    }
    me_cols <- grep("^ME", colnames(MEs), value = TRUE)
    me_cols <- setdiff(me_cols, "MEgrey")
    fit <- limma::lmFit(t(as.matrix(MEs[, me_cols, drop = FALSE])), design)
    fit <- limma::eBayes(fit, trend = FALSE, robust = TRUE)
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "none")
    tt$ME <- rownames(tt)
    tt$module_color <- sub("^ME", "", tt$ME)
    out <- data.frame(network_name = network_name, comparison = "Batch_B_vs_Batch_A_adjusted_for_timepoint", 
        n_samples = nrow(m), module_color = tt$module_color, estimate_BatchB_minus_BatchA = tt$logFC, 
        t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
    list(result = out, audit = audit)
}

run_phase_association <- function(MEs, metadata, network_name) {
    m_all <- metadata[match(rownames(MEs), metadata$sample_id), , drop = FALSE]
    keep <- !is.na(m_all$phase_primary) & as_chr(m_all$phase_primary) %in% c("Early", "Established late")
    m <- m_all[keep, , drop = FALSE]
    me <- MEs[rownames(MEs) %in% m$sample_id, , drop = FALSE]
    me <- me[m$sample_id, , drop = FALSE]
    b <- factor(as_chr(m$batch), levels = c("Batch_A", "Batch_B"))
    g <- factor(ifelse(as_chr(m$phase_primary) == "Established late", "EstablishedLate", "Early"), levels = c("Early", 
        "EstablishedLate"))
    design <- model.matrix(~b + g)
    design_rank <- qr(design)$rank
    full_rank <- design_rank == ncol(design)
    audit <- data.frame(network_name = network_name, model = "~ batch + phase_primary", n_samples = nrow(m), 
        n_parameters = ncol(design), rank = design_rank, full_rank = full_rank, stringsAsFactors = FALSE)
    if (!full_rank) {
        stop("Primary-phase association design is not full rank for ", network_name)
    }
    coef_name <- grep("^g", colnames(design), value = TRUE)
    if (length(coef_name) != 1L) {
        stop("Could not uniquely identify phase coefficient for ", network_name)
    }
    me_cols <- grep("^ME", colnames(me), value = TRUE)
    me_cols <- setdiff(me_cols, "MEgrey")
    fit <- limma::lmFit(t(as.matrix(me[, me_cols, drop = FALSE])), design)
    fit <- limma::eBayes(fit, trend = FALSE, robust = TRUE)
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "none")
    tt$ME <- rownames(tt)
    tt$module_color <- sub("^ME", "", tt$ME)
    out <- data.frame(network_name = network_name, comparison = "EstablishedLate_vs_Early_adjusted_for_batch", 
        n_samples = nrow(m), module_color = tt$module_color, estimate_EstablishedLate_minus_Early = tt$logFC, 
        t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
    list(result = out, audit = audit)
}

batch_results <- list()
phase_results <- list()
model_audits <- list()

primary_batch <- run_batch_association(primary_MEs, meta, "primary_original_expression_all_samples")
primary_phase <- run_phase_association(primary_MEs, meta, "primary_original_expression_all_samples")

batch_results[[length(batch_results) + 1L]] <- primary_batch$result
phase_results[[length(phase_results) + 1L]] <- primary_phase$result
model_audits[[length(model_audits) + 1L]] <- primary_batch$audit
model_audits[[length(model_audits) + 1L]] <- primary_phase$audit

adj_batch <- run_batch_association(network_batch_adjusted$module_eigengenes, meta, "batch_adjusted_expression_all_samples")
adj_phase <- run_phase_association(network_batch_adjusted$module_eigengenes, meta, "batch_adjusted_expression_all_samples")

batch_results[[length(batch_results) + 1L]] <- adj_batch$result
phase_results[[length(phase_results) + 1L]] <- adj_phase$result
model_audits[[length(model_audits) + 1L]] <- adj_batch$audit
model_audits[[length(model_audits) + 1L]] <- adj_phase$audit

if (!is.null(network_outlier_excluded)) {
    outlier_meta <- meta[match(network_outlier_excluded$sample_ids, meta$sample_id), , drop = FALSE]
    outlier_batch <- run_batch_association(network_outlier_excluded$module_eigengenes, outlier_meta, 
        "original_expression_excluding_review_flags")
    outlier_phase <- run_phase_association(network_outlier_excluded$module_eigengenes, outlier_meta, 
        "original_expression_excluding_review_flags")
    batch_results[[length(batch_results) + 1L]] <- outlier_batch$result
    phase_results[[length(phase_results) + 1L]] <- outlier_phase$result
    model_audits[[length(model_audits) + 1L]] <- outlier_batch$audit
    model_audits[[length(model_audits) + 1L]] <- outlier_phase$audit
}

batch_association_all <- do.call(rbind, batch_results)
phase_association_all <- do.call(rbind, phase_results)
model_audit_all <- do.call(rbind, model_audits)

write_csv(batch_association_all, paste0("13_", PREFIX, "_module_batch_association_all_networks.csv"))
write_csv(phase_association_all, paste0("14_", PREFIX, "_module_phase_association_all_networks.csv"))
write_csv(model_audit_all, paste0("15_", PREFIX, "_association_model_rank_audit.csv"))

get_association_row <- function(x, network_name, module_color) {
  z <- x[x$network_name == network_name & x$module_color == module_color, , drop = FALSE]
  if (nrow(z) == 0L)
    return(NULL)
  z[1, , drop = FALSE]
}

brown_summary_rows <- list()

primary_brown_n <- sum(primary_assignment$module_color == "brown")

primary_batch_brown <- get_association_row(batch_association_all, "primary_original_expression_all_samples",
  "brown")
primary_phase_brown <- get_association_row(phase_association_all, "primary_original_expression_all_samples",
  "brown")

brown_summary_rows[[length(brown_summary_rows) + 1L]] <- data.frame(network_name = "primary_original_expression_all_samples",
  primary_brown_n = primary_brown_n, matched_alternative_module = "brown", alternative_module_n = primary_brown_n,
  overlap_n = primary_brown_n, jaccard = 1, primary_brown_retention = 1, alternative_module_purity = 1,
  matched_ME_pearson_r = 1, batch_estimate_B_minus_A = if (!is.null(primary_batch_brown)) {
    primary_batch_brown$estimate_BatchB_minus_BatchA
  } else NA_real_, batch_FDR = if (!is.null(primary_batch_brown)) {
    primary_batch_brown$adj.P.Val
  } else NA_real_, phase_estimate_late_minus_early = if (!is.null(primary_phase_brown)) {
    primary_phase_brown$estimate_EstablishedLate_minus_Early
  } else NA_real_, phase_FDR = if (!is.null(primary_phase_brown)) {
    primary_phase_brown$adj.P.Val
  } else NA_real_, stringsAsFactors = FALSE)

for (nm in names(network_list)) {
  ov <- module_overlap[module_overlap$network_name == nm & module_overlap$primary_module == "brown",
    , drop = FALSE]

  if (nrow(ov) == 0L)
    next
  ov <- ov[1, , drop = FALSE]
  matched_alt <- as.character(ov$alternative_module)

  ec <- eigengene_correlations[eigengene_correlations$network_name == nm & eigengene_correlations$primary_module ==
    "brown" & eigengene_correlations$alternative_module == matched_alt, , drop = FALSE]

  ba <- get_association_row(batch_association_all, if (nm == "batch_adjusted_all_samples") {
    "batch_adjusted_expression_all_samples"
  } else {
    "original_expression_excluding_review_flags"
  }, matched_alt)

  pa <- get_association_row(phase_association_all, if (nm == "batch_adjusted_all_samples") {
    "batch_adjusted_expression_all_samples"
  } else {
    "original_expression_excluding_review_flags"
  }, matched_alt)

  brown_summary_rows[[length(brown_summary_rows) + 1L]] <- data.frame(network_name = nm, primary_brown_n = primary_brown_n,
    matched_alternative_module = matched_alt, alternative_module_n = ov$alternative_module_n_common_universe,
    overlap_n = ov$overlap_n, jaccard = ov$jaccard, primary_brown_retention = ov$primary_module_retention,
    alternative_module_purity = ov$alternative_module_purity, matched_ME_pearson_r = if (nrow(ec) >
      0L)
      ec$pearson_r[1] else NA_real_, batch_estimate_B_minus_A = if (!is.null(ba)) {
      ba$estimate_BatchB_minus_BatchA
    } else NA_real_, batch_FDR = if (!is.null(ba))
      ba$adj.P.Val else NA_real_, phase_estimate_late_minus_early = if (!is.null(pa)) {
      pa$estimate_EstablishedLate_minus_Early
    } else NA_real_, phase_FDR = if (!is.null(pa))
      pa$adj.P.Val else NA_real_, stringsAsFactors = FALSE)
}

brown_stability_summary <- do.call(rbind, brown_summary_rows)

write_csv(brown_stability_summary, paste0("16_", PREFIX, "_brown_module_stability_summary.csv"))

sensitivity_object <- list(step = 102, version = "v01", purpose = "WGCNA preprocessing and outlier sensitivity",
  input_step45_rds = STEP45_RDS, input_step49_rds = STEP49_RDS, frozen_primary_settings = list(softPower = SOFT_POWER,
    networkType = NETWORK_TYPE, TOMType = TOM_TYPE, corType = COR_TYPE, maxPOutliers = MAX_P_OUTLIERS,
    minModuleSize = MIN_MODULE_SIZE, deepSplit = DEEP_SPLIT, mergeCutHeight = MERGE_CUT_HEIGHT, maxBlockSize = MAX_BLOCK_SIZE),
  batch_provenance = batch_provenance, batch_timepoint_counts = batch_time_table, timepoint_batch_presence = timepoint_batch_presence,
  sample_outlier_audit = sample_outlier_audit, flagged_samples = flagged_samples, batch_adjusted_soft_threshold_scan = sft_adjusted_table,
  batch_adjusted_soft_threshold_selection = adjusted_power_selection, batch_adjusted_network = network_batch_adjusted,
  outlier_excluded_network = network_outlier_excluded, module_overlap = module_overlap, eigengene_correlations = eigengene_correlations,
  module_batch_association = batch_association_all, module_phase_association = phase_association_all,
  brown_stability_summary = brown_stability_summary)

saveRDS(sensitivity_object, file = file.path(OUT_DIR, paste0("20_", PREFIX, "_analysis_object.rds")),
  compress = "xz")

pdf(file = file.path(OUT_DIR, paste0("21_", PREFIX, "_diagnostics.pdf")), width = 10, height = 7, onefile = TRUE)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))

barplot(t(as.matrix(table(meta$timepoint, meta$batch))), beside = TRUE, las = 2, ylab = "Number of samples",
  main = "Operational batch by timepoint", legend.text = TRUE, args.legend = list(x = "topright", bty = "n",
    cex = 0.8))

brown_plot <- brown_stability_summary[brown_stability_summary$network_name != "primary_original_expression_all_samples",
  , drop = FALSE]

if (nrow(brown_plot) > 0L) {
  barplot(brown_plot$jaccard, names.arg = gsub("_", "\n", brown_plot$network_name), las = 2, ylim = c(0,
    1), ylab = "Jaccard index", main = "Primary Brown module: gene overlap")
} else {
  plot.new()
  title("No alternative network available for Brown comparison")
}

if (nrow(brown_plot) > 0L) {
  nm0 <- brown_plot$network_name[1]
  alt_color0 <- brown_plot$matched_alternative_module[1]
  alt_obj0 <- network_list[[nm0]]
  alt_me_name0 <- paste0("ME", alt_color0)

  if (!is.null(alt_obj0) && "MEbrown" %in% colnames(primary_MEs) && alt_me_name0 %in% colnames(alt_obj0$module_eigengenes)) {
    common0 <- intersect(rownames(primary_MEs), rownames(alt_obj0$module_eigengenes))
    plot(primary_MEs[common0, "MEbrown"], alt_obj0$module_eigengenes[common0, alt_me_name0], xlab = "Primary MEbrown",
      ylab = paste0("Alternative ", alt_me_name0), main = paste("Brown eigengene comparison:",
        nm0), pch = 19)
    abline(lm(alt_obj0$module_eigengenes[common0, alt_me_name0] ~ primary_MEs[common0, "MEbrown"]))
  } else {
    plot.new()
    title("Matched Brown eigengene unavailable")
  }
} else {
  plot.new()
  title("No alternative eigengene comparison")
}

plot(sft_adjusted_table$Power, sft_adjusted_table$SFT.R.sq, type = "b", pch = 19, xlab = "Soft-threshold power",
  ylab = "Scale-free topology fit (R-squared)", main = "Batch-adjusted matrix: soft-threshold scan")
abline(h = 0.8, lty = 2)
abline(v = SOFT_POWER, lty = 3)

par(old_par)
dev.off()

get_brown_row <- function(network_name) {
  brown_stability_summary[brown_stability_summary$network_name == network_name, , drop = FALSE]
}

brown_primary_row <- get_brown_row("primary_original_expression_all_samples")

brown_adj_row <- get_brown_row("batch_adjusted_all_samples")

brown_outlier_row <- get_brown_row("original_expression_excluding_review_flags")

primary_brown_batch <- get_association_row(batch_association_all, "primary_original_expression_all_samples",
  "brown")
adjusted_brown_batch <- if (nrow(brown_adj_row) > 0L) {
  get_association_row(batch_association_all, "batch_adjusted_expression_all_samples", as.character(brown_adj_row$matched_alternative_module))
} else NULL

primary_brown_phase <- get_association_row(phase_association_all, "primary_original_expression_all_samples",
  "brown")
adjusted_brown_phase <- if (nrow(brown_adj_row) > 0L) {
  get_association_row(phase_association_all, "batch_adjusted_expression_all_samples", as.character(brown_adj_row$matched_alternative_module))
} else NULL

writeLines("Numerical results are provided in the accompanying output tables.", 
    con = file.path(OUT_DIR, paste0("22_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

sink(file.path(OUT_DIR, paste0("23_", PREFIX, "_sessionInfo.txt")))
print(sessionInfo())
sink()

output_files <- list.files(OUT_DIR, full.names = TRUE, recursive = FALSE)

output_index <- data.frame(file_name = basename(output_files), size_bytes = as.numeric(file.info(output_files)$size),
  size_MB = round(as.numeric(file.info(output_files)$size)/1024^2, 3), stringsAsFactors = FALSE)

write_csv(output_index, paste0("00_", PREFIX, "_OUTPUT_INDEX.csv"))

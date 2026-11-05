source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)

set.seed(109)

find_project_root <- function() {
    env_root <- Sys.getenv("SCI_PROJECT_ROOT", unset = "")
    candidates <- unique(c(if (nzchar(env_root)) env_root else character(0), PROJECT_ROOT, getwd()))
    is_project <- function(x) {
        if (!dir.exists(x)) 
            return(FALSE)
        file.exists(file.path(x, "results", "45_Bulk_batch_structure_and_analysis_object", "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds"))
    }
    ok <- candidates[vapply(candidates, is_project, logical(1))]
    if (!length(ok)) {
        stop("Could not identify the zz project root.\n", "Expected the Step45 analysis object under results/.\n", 
            "Open the zz RStudio project or set ", "Sys.setenv(SCI_PROJECT_ROOT='/path/to/project').")
    }
    normalizePath(ok[1], winslash = "/", mustWork = TRUE)
}

PROJECT_ROOT <- find_project_root()
path_in <- function(...) file.path(PROJECT_ROOT, ...)

STEP45_RDS <- path_in("results", "45_Bulk_batch_structure_and_analysis_object", "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds")

OUT_DIR <- path_in("results", "109_Bulk_timepoint_sensitivity")
TABLES_DIR <- file.path(OUT_DIR, "tables")

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "109_Bulk_timepoint_sensitivity_v01"

num_name <- function(num, suffix) {
  paste0(sprintf("%02d", num), "_", PREFIX, "_", suffix)
}

LOG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".log"))
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)

log_msg <- function(...) {
  z <- paste0(...)
  cat(z, "\n")
  cat(z, "\n", file = LOG_FILE, append = TRUE)
}

write_csv <- function(x, filename) {
  utils::write.csv(x, file.path(TABLES_DIR, filename), row.names = FALSE, fileEncoding = "UTF-8")
}

EARLY_TIMEPOINTS <- c("Uninjured", "15 min", "1 dpi", "3 dpi")

ESTABLISHED_LATE_TIMEPOINTS <- c("14 dpi", "28 dpi")

CORE_GENES <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

ADJ_P_CUTOFF <- 0.05
ABS_LOGFC_CUTOFF <- 0.5
MIN_SAMPLES_PER_GROUP <- 2L
MIN_RESIDUAL_DF <- 2L

required_pkgs <- c("limma")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs)) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them before running Step109. No input files were modified.")
}

suppressPackageStartupMessages({
  library(limma)
})

safe_gene_symbol <- function(x) {
  ifelse(grepl("^X[0-9]", x), sub("^X", "", x), x)
}

annotate_deg <- function(tt) {
    tt$stat_sig <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < ADJ_P_CUTOFF
    tt$effect_sig <- tt$stat_sig & !is.na(tt$logFC) & abs(tt$logFC) >= ABS_LOGFC_CUTOFF
    tt$direction <- "NS"
    tt$direction[tt$stat_sig & !is.na(tt$logFC) & tt$logFC >= ABS_LOGFC_CUTOFF] <- "Up"
    tt$direction[tt$stat_sig & !is.na(tt$logFC) & tt$logFC <= -ABS_LOGFC_CUTOFF] <- "Down"
    tt
}

add_gene_columns <- function(tt) {
    if (is.null(rownames(tt))) {
        stop("limma result has no gene row names.")
    }
    out <- data.frame(gene_id = rownames(tt), gene_symbol_restored = safe_gene_symbol(rownames(tt)), 
        tt, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(out) <- NULL
    annotate_deg(out)
}

safe_rbind <- function(x) {

  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, x)
  if (!length(x))
    return(data.frame())

  all_cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing_cols <- setdiff(all_cols, names(z))
    for (nm in missing_cols) z[[nm]] <- NA
    z <- z[, all_cols, drop = FALSE]
    rownames(z) <- NULL
    z
  })

  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L)
    return(NA_real_)
  if (length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

count_degs <- function(tt, analysis_id, analysis_label) {
  data.frame(analysis_id = analysis_id, analysis = analysis_label, n_genes_tested = nrow(tt), n_adjP_lt_0.05 = sum(tt$adj.P.Val <
    ADJ_P_CUTOFF, na.rm = TRUE), n_up_adjP_lt_0.05_logFC_ge_0.5 = sum(tt$adj.P.Val < ADJ_P_CUTOFF &
    tt$logFC >= ABS_LOGFC_CUTOFF, na.rm = TRUE), n_down_adjP_lt_0.05_logFC_le_minus0.5 = sum(tt$adj.P.Val <
    ADJ_P_CUTOFF & tt$logFC <= -ABS_LOGFC_CUTOFF, na.rm = TRUE), stringsAsFactors = FALSE)
}

extract_core_genes <- function(tt, analysis_id, analysis_label) {
    out <- lapply(CORE_GENES, function(g) {
        idx <- match(g, tt$gene_symbol_restored)
        if (is.na(idx)) {
            return(data.frame(analysis_id = analysis_id, analysis = analysis_label, requested_gene = g, 
                found = FALSE, gene_id = NA_character_, logFC = NA_real_, AveExpr = NA_real_, t = NA_real_, 
                P.Value = NA_real_, adj.P.Val = NA_real_, direction = NA_character_, effect_sig = NA, 
                stringsAsFactors = FALSE))
        }
        z <- tt[idx, , drop = FALSE]
        data.frame(analysis_id = analysis_id, analysis = analysis_label, requested_gene = g, found = TRUE, 
            gene_id = as.character(z$gene_id[1]), logFC = as.numeric(z$logFC[1]), AveExpr = if ("AveExpr" %in% 
                names(z)) 
                as.numeric(z$AveExpr[1])
            else NA_real_, t = if ("t" %in% names(z)) 
                as.numeric(z$t[1])
            else NA_real_, P.Value = as.numeric(z$P.Value[1]), adj.P.Val = as.numeric(z$adj.P.Val[1]), 
            direction = as.character(z$direction[1]), effect_sig = as.logical(z$effect_sig[1]), stringsAsFactors = FALSE)
    })
    safe_rbind(out)
}

build_design <- function(md) {
    md$analysis_group <- factor(as.character(md$analysis_group), levels = c("Early", "Target"))
    batch_factor <- droplevels(factor(as.character(md$batch)))
    group_mm <- stats::model.matrix(~0 + analysis_group, data = md)
    colnames(group_mm) <- c("Group_Early", "Group_Target")
    if (nlevels(batch_factor) > 1L) {
        batch_mm <- stats::model.matrix(~batch_factor)
        batch_mm <- batch_mm[, -1, drop = FALSE]
        colnames(batch_mm) <- make.names(colnames(batch_mm))
        design <- cbind(group_mm, batch_mm)
        batch_adjusted <- TRUE
    }
    else {
        design <- group_mm
        batch_adjusted <- FALSE
    }
    list(design = design, batch_adjusted = batch_adjusted, batch_levels = levels(batch_factor))
}

fit_one_analysis <- function(analysis_id, analysis_label, target_timepoints, excluded_timepoints = character(0), 
    is_reference = FALSE) {
    contrast_timepoints <- unique(c(EARLY_TIMEPOINTS, target_timepoints))
    candidate <- meta[meta$timepoint %in% contrast_timepoints & !(meta$timepoint %in% excluded_timepoints), 
        , drop = FALSE]
    excluded_md <- meta[meta$timepoint %in% excluded_timepoints & meta$timepoint %in% contrast_timepoints, 
        , drop = FALSE]
    candidate$analysis_group <- ifelse(candidate$timepoint %in% EARLY_TIMEPOINTS, "Early", "Target")
    candidate$analysis_group <- factor(candidate$analysis_group, levels = c("Early", "Target"))
    sample_ids <- rownames(candidate)
    expr_sub <- expr[, sample_ids, drop = FALSE]
    n_early <- sum(candidate$analysis_group == "Early")
    n_target <- sum(candidate$analysis_group == "Target")
    n_samples <- nrow(candidate)
    sample_summary <- data.frame(analysis_id = analysis_id, analysis = analysis_label, is_reference_model = is_reference, 
        target_timepoints = paste(target_timepoints, collapse = "; "), early_timepoints_included = paste(intersect(EARLY_TIMEPOINTS, 
            unique(candidate$timepoint)), collapse = "; "), sensitivity_excluded_timepoints = if (length(excluded_timepoints)) {
            paste(excluded_timepoints, collapse = "; ")
        }
        else {
            "None"
        }, n_early = n_early, n_target = n_target, n_total = n_samples, n_sensitivity_excluded = nrow(excluded_md), 
        sensitivity_excluded_sample_ids = if (nrow(excluded_md)) {
            paste(rownames(excluded_md), collapse = "; ")
        }
        else {
            "None"
        }, included_sample_ids = if (n_samples) {
            paste(sample_ids, collapse = "; ")
        }
        else {
            "None"
        }, other_timepoints_not_in_this_contrast = paste(setdiff(unique(meta$timepoint), contrast_timepoints), 
            collapse = "; "), stringsAsFactors = FALSE)
    if (n_samples) {
        dist <- as.data.frame(table(analysis_group = candidate$analysis_group, timepoint = candidate$timepoint, 
            batch = candidate$batch), stringsAsFactors = FALSE)
        dist <- dist[dist$Freq > 0, , drop = FALSE]
        dist$analysis_id <- analysis_id
        dist$analysis <- analysis_label
        names(dist)[names(dist) == "Freq"] <- "n_samples"
        dist <- dist[, c("analysis_id", "analysis", "analysis_group", "timepoint", "batch", "n_samples"), 
            drop = FALSE]
    }
    else {
        dist <- data.frame()
    }
    base_diag <- data.frame(analysis_id = analysis_id, analysis = analysis_label, status = "NOT_ESTIMABLE", 
        n_samples = n_samples, n_early = n_early, n_target = n_target, n_genes = nrow(expr_sub), n_batches = length(unique(candidate$batch)), 
        batch_adjusted = NA, n_parameters = NA_integer_, design_rank = NA_integer_, residual_df = NA_integer_, 
        full_rank = FALSE, reason = NA_character_, stringsAsFactors = FALSE)
    if (n_early < MIN_SAMPLES_PER_GROUP || n_target < MIN_SAMPLES_PER_GROUP) {
        base_diag$reason <- paste0("At least ", MIN_SAMPLES_PER_GROUP, " samples per group are required; observed Early=", 
            n_early, ", Target=", n_target, ".")
        return(list(table = NULL, sample_summary = sample_summary, distribution = dist, diagnostics = base_diag, 
            weights = data.frame()))
    }
    if (anyNA(candidate$batch)) {
        base_diag$reason <- "Missing batch labels in the selected samples."
        return(list(table = NULL, sample_summary = sample_summary, distribution = dist, diagnostics = base_diag, 
            weights = data.frame()))
    }
    design_info <- build_design(candidate)
    design <- design_info$design
    q <- qr(design)
    residual_df <- nrow(design) - q$rank
    base_diag$batch_adjusted <- design_info$batch_adjusted
    base_diag$n_parameters <- ncol(design)
    base_diag$design_rank <- q$rank
    base_diag$residual_df <- residual_df
    base_diag$full_rank <- q$rank == ncol(design)
    if (q$rank != ncol(design)) {
        base_diag$reason <- paste("The group and batch terms are not independently estimable", "(design matrix is not full rank).")
        return(list(table = NULL, sample_summary = sample_summary, distribution = dist, diagnostics = base_diag, 
            weights = data.frame()))
    }
    if (residual_df < MIN_RESIDUAL_DF) {
        base_diag$reason <- paste0("Residual degrees of freedom < ", MIN_RESIDUAL_DF, "; model not fitted.")
        return(list(table = NULL, sample_summary = sample_summary, distribution = dist, diagnostics = base_diag, 
            weights = data.frame()))
    }
    contrast <- matrix(0, nrow = ncol(design), ncol = 1, dimnames = list(colnames(design), "Target_vs_Early"))
    contrast["Group_Target", 1] <- 1
    contrast["Group_Early", 1] <- -1
    fit_result <- tryCatch({
        aw <- limma::arrayWeights(expr_sub, design = design)
        fit0 <- limma::lmFit(expr_sub, design = design, weights = aw)
        fit1 <- limma::contrasts.fit(fit0, contrast)
        fit1 <- limma::eBayes(fit1, trend = TRUE, robust = TRUE)
        tt <- limma::topTable(fit1, coef = "Target_vs_Early", number = Inf, sort.by = "P", adjust.method = "BH")
        tt <- add_gene_columns(tt)
        aw_df <- data.frame(analysis_id = analysis_id, analysis = analysis_label, sample_id = sample_ids, 
            analysis_group = as.character(candidate$analysis_group), timepoint = as.character(candidate$timepoint), 
            batch = as.character(candidate$batch), array_weight = as.numeric(aw), stringsAsFactors = FALSE)
        list(table = tt, weights = aw_df, error = NULL)
    }, error = function(e) {
        list(table = NULL, weights = data.frame(), error = conditionMessage(e))
    })
    if (!is.null(fit_result$error)) {
        base_diag$status <- "ERROR"
        base_diag$reason <- fit_result$error
    }
    else {
        base_diag$status <- "FITTED"
        base_diag$reason <- if (design_info$batch_adjusted) {
            "Full-rank batch-adjusted model fitted."
        }
        else {
            "Full-rank model fitted; only one batch was represented, so no batch term was estimable."
        }
    }
    list(table = fit_result$table, sample_summary = sample_summary, distribution = dist, diagnostics = base_diag, 
        weights = fit_result$weights)
}

compare_to_reference <- function(test_table, reference_table, analysis_id, analysis_label) {
    if (is.null(test_table) || is.null(reference_table)) {
        return(data.frame(analysis_id = analysis_id, analysis = analysis_label, status = "NOT_AVAILABLE", 
            n_common_genes = NA_integer_, pearson_logFC = NA_real_, spearman_logFC = NA_real_, same_direction_fraction = NA_real_, 
            reference_effect_sig_n = NA_integer_, sensitivity_effect_sig_n = NA_integer_, shared_effect_sig_n = NA_integer_, 
            effect_sig_jaccard = NA_real_, stringsAsFactors = FALSE))
    }
    ref <- reference_table[, c("gene_id", "logFC", "effect_sig"), drop = FALSE]
    test <- test_table[, c("gene_id", "logFC", "effect_sig"), drop = FALSE]
    names(ref)[2:3] <- c("logFC_reference", "effect_sig_reference")
    names(test)[2:3] <- c("logFC_sensitivity", "effect_sig_sensitivity")
    m <- merge(ref, test, by = "gene_id", all = FALSE)
    ref_sig <- m$gene_id[m$effect_sig_reference %in% TRUE]
    test_sig <- m$gene_id[m$effect_sig_sensitivity %in% TRUE]
    union_sig <- union(ref_sig, test_sig)
    shared_sig <- intersect(ref_sig, test_sig)
    data.frame(analysis_id = analysis_id, analysis = analysis_label, status = "COMPARED_WITH_REFERENCE", 
        n_common_genes = nrow(m), pearson_logFC = safe_cor(m$logFC_reference, m$logFC_sensitivity, "pearson"), 
        spearman_logFC = safe_cor(m$logFC_reference, m$logFC_sensitivity, "spearman"), same_direction_fraction = mean(sign(m$logFC_reference) == 
            sign(m$logFC_sensitivity), na.rm = TRUE), reference_effect_sig_n = length(ref_sig), sensitivity_effect_sig_n = length(test_sig), 
        shared_effect_sig_n = length(shared_sig), effect_sig_jaccard = if (length(union_sig)) {
            length(shared_sig)/length(union_sig)
        }
        else {
            NA_real_
        }, stringsAsFactors = FALSE)
}

if (!file.exists(STEP45_RDS)) {
  stop("Step45 input not found:\n", STEP45_RDS)
}

obj <- readRDS(STEP45_RDS)

if (!all(c("expression_original", "metadata") %in% names(obj))) {
  stop("Step45 analysis object must contain 'expression_original' and 'metadata'.")
}

expr <- as.matrix(obj$expression_original)
storage.mode(expr) <- "numeric"
meta <- as.data.frame(obj$metadata, stringsAsFactors = FALSE)

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix must have gene row names and sample column names.")
}
if (anyDuplicated(rownames(expr))) {
  stop("Duplicated gene identifiers were detected in expression_original.")
}
if (is.null(rownames(meta))) {
  stop("Step45 metadata must have sample IDs as row names.")
}
if (!all(colnames(expr) %in% rownames(meta))) {
  stop("Expression sample IDs do not all match Step45 metadata row names.")
}

required_meta <- c("batch", "timepoint")
missing_meta <- setdiff(required_meta, colnames(meta))
if (length(missing_meta)) {
  stop("Step45 metadata is missing required field(s): ", paste(missing_meta, collapse = ", "))
}

meta <- meta[colnames(expr), , drop = FALSE]
meta$batch <- as.character(meta$batch)
meta$timepoint <- as.character(meta$timepoint)

if (anyNA(expr) || any(!is.finite(expr))) {
  stop("Non-finite values were detected in expression_original.")
}
if (anyNA(meta$timepoint) || any(!nzchar(meta$timepoint))) {
  stop("Missing or empty timepoint labels were detected.")
}
if (anyNA(meta$batch) || any(!nzchar(meta$batch))) {
  stop("Missing or empty batch labels were detected.")
}

required_timepoints <- unique(c(EARLY_TIMEPOINTS, ESTABLISHED_LATE_TIMEPOINTS))
missing_timepoints <- setdiff(required_timepoints, unique(meta$timepoint))
if (length(missing_timepoints)) {
  stop("Expected timepoint label(s) not found in Step45 metadata: ", paste(missing_timepoints, collapse = ", "),
    "\nObserved labels: ", paste(sort(unique(meta$timepoint)), collapse = "; "))
}

reference <- fit_one_analysis(analysis_id = "REFERENCE_primary_all_samples", analysis_label = "Reference: Established late (14+28 dpi) vs Early",
  target_timepoints = ESTABLISHED_LATE_TIMEPOINTS, excluded_timepoints = character(0), is_reference = TRUE)

if (reference$diagnostics$status != "FITTED") {
  stop("The reference Step46 primary contrast could not be reproduced. ", "Reason: ", reference$diagnostics$reason)
}

reference_file <- num_name(4, "reference_primary_all_samples.csv")
write_csv(reference$table, reference_file)

analysis_specs <- list(list(id = "S1_exclude_Uninjured", label = "Sensitivity 1: exclude Uninjured",
  target = ESTABLISHED_LATE_TIMEPOINTS, exclude = "Uninjured", output_no = 5L, output_suffix = "sensitivity_exclude_Uninjured.csv"),
  list(id = "S2_exclude_15min", label = "Sensitivity 2: exclude 15 min", target = ESTABLISHED_LATE_TIMEPOINTS,
    exclude = "15 min", output_no = 6L, output_suffix = "sensitivity_exclude_15min.csv"), list(id = "S3_exclude_Uninjured_and_15min",
    label = "Sensitivity 3: exclude Uninjured and 15 min", target = ESTABLISHED_LATE_TIMEPOINTS,
    exclude = c("Uninjured", "15 min"), output_no = 7L, output_suffix = "sensitivity_exclude_Uninjured_and_15min.csv"),
  list(id = "S4_14dpi_vs_Early", label = "Sensitivity 4: 14 dpi vs Early", target = "14 dpi", exclude = character(0),
    output_no = 8L, output_suffix = "sensitivity_14dpi_vs_Early.csv"), list(id = "S5_28dpi_vs_Early",
    label = "Sensitivity 5: 28 dpi vs Early", target = "28 dpi", exclude = character(0), output_no = 9L,
    output_suffix = "sensitivity_28dpi_vs_Early.csv"))

results <- list()
results[[reference$diagnostics$analysis_id]] <- reference

for (i in seq_along(analysis_specs)) {
    sp <- analysis_specs[[i]]
    res <- fit_one_analysis(analysis_id = sp$id, analysis_label = sp$label, target_timepoints = sp$target, 
        excluded_timepoints = sp$exclude, is_reference = FALSE)
    results[[sp$id]] <- res
    if (!is.null(res$table)) {
        write_csv(res$table, num_name(sp$output_no, sp$output_suffix))
    }
    else {
    }
}

sample_summary <- safe_rbind(lapply(results, `[[`, "sample_summary"))
sample_distribution <- safe_rbind(lapply(results, `[[`, "distribution"))
design_diagnostics <- safe_rbind(lapply(results, `[[`, "diagnostics"))
array_weights <- safe_rbind(lapply(results, `[[`, "weights"))

write_csv(sample_summary, num_name(1, "sample_inclusion_summary.csv"))
write_csv(sample_distribution, num_name(2, "sample_distribution_by_group_batch_timepoint.csv"))
write_csv(design_diagnostics, num_name(3, "design_diagnostics_and_estimability.csv"))

if (nrow(array_weights)) {
  write_csv(array_weights, num_name(13, "arrayWeights_by_analysis.csv"))
}

deg_count_rows <- list()
core_gene_rows <- list()
concordance_rows <- list()

for (nm in names(results)) {
    z <- results[[nm]]
    if (!is.null(z$table)) {
        deg_count_rows[[length(deg_count_rows) + 1L]] <- count_degs(z$table, nm, z$diagnostics$analysis)
        core_gene_rows[[length(core_gene_rows) + 1L]] <- extract_core_genes(z$table, nm, z$diagnostics$analysis)
    }
    if (nm != "REFERENCE_primary_all_samples") {
        concordance_rows[[length(concordance_rows) + 1L]] <- compare_to_reference(z$table, reference$table, 
            nm, z$diagnostics$analysis)
    }
}

deg_summary <- safe_rbind(deg_count_rows)
core_summary <- safe_rbind(core_gene_rows)
concordance_summary <- safe_rbind(concordance_rows)

write_csv(deg_summary, num_name(10, "DEG_count_summary.csv"))
write_csv(concordance_summary, num_name(11, "gene_effect_concordance_vs_reference.csv"))
write_csv(core_summary, num_name(12, "core_gene_sensitivity_results.csv"))

audit_checks <- data.frame(check = c("Step45_expression_original_present", "Step45_metadata_present",
  "Expression_sample_ids_match_metadata", "Required_timepoint_labels_present", "Reference_primary_model_fitted",
  "Exactly_five_sensitivity_models_defined", "All_fitted_models_full_rank", "No_source_objects_modified"),
  result = c(TRUE, TRUE, all(colnames(expr) %in% rownames(meta)), all(required_timepoints %in% unique(meta$timepoint)),
    reference$diagnostics$status == "FITTED", length(analysis_specs) == 5L, all(design_diagnostics$full_rank[design_diagnostics$status ==
      "FITTED"]), TRUE), detail = c("Read from Step45 analysis object.", "Read from Step45 analysis object.",
    "All expression columns were matched to metadata rows.", paste(required_timepoints, collapse = "; "),
    reference$diagnostics$reason, paste(vapply(analysis_specs, `[[`, character(1), "id"), collapse = "; "),
    "Any non-estimable sensitivity model is reported and not forced.", "The script only reads Step45 and writes new Step109 outputs."),
  stringsAsFactors = FALSE)

write_csv(audit_checks, num_name(14, "audit_checks.csv"))

output_files <- data.frame(No = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
  "13", "14"), file = c(num_name(1, "sample_inclusion_summary.csv"), num_name(2, "sample_distribution_by_group_batch_timepoint.csv"),
  num_name(3, "design_diagnostics_and_estimability.csv"), num_name(4, "reference_primary_all_samples.csv"),
  num_name(5, "sensitivity_exclude_Uninjured.csv"), num_name(6, "sensitivity_exclude_15min.csv"), num_name(7,
    "sensitivity_exclude_Uninjured_and_15min.csv"), num_name(8, "sensitivity_14dpi_vs_Early.csv"),
  num_name(9, "sensitivity_28dpi_vs_Early.csv"), num_name(10, "DEG_count_summary.csv"), num_name(11,
    "gene_effect_concordance_vs_reference.csv"), num_name(12, "core_gene_sensitivity_results.csv"),
  num_name(13, "arrayWeights_by_analysis.csv"), num_name(14, "audit_checks.csv")), purpose = c("Included/excluded samples and group sizes for each model.",
  "Sample counts by analysis, group, time point, and batch.", "Model rank, residual degrees of freedom, and estimability.",
  "Full-sample Step46 primary contrast; reference only.", "Full gene-level results after excluding Uninjured.",
  "Full gene-level results after excluding 15 min.", "Full gene-level results after excluding both early time points.",
  "Full gene-level results for 14 dpi vs Early.", "Full gene-level results for 28 dpi vs Early.", "DEG counts under the Step46 thresholds.",
  "Effect-size and DEG-set concordance against the reference.", "Clec7a/lipid-handling core genes across fitted models.",
  "Sample-level array weights for each fitted model.", "Input, design, and output integrity checks."),
  status = c(rep("WRITTEN", 4), vapply(analysis_specs, function(sp) {
    z <- results[[sp$id]]
    if (is.null(z$table)) z$diagnostics$status else "WRITTEN"
  }, character(1)), rep("WRITTEN", 5)), stringsAsFactors = FALSE)

utils::write.csv(output_files, file.path(TABLES_DIR, paste0("00_", PREFIX, "_OUTPUT_INDEX.csv")), row.names = FALSE,
  fileEncoding = "UTF-8")

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT_DIR, paste0(PREFIX, "_REPORT.txt")), useBytes = TRUE)

writeLines(capture.output(utils::sessionInfo()), file.path(OUT_DIR, paste0(PREFIX, "_sessionInfo.txt")),
  useBytes = TRUE)

log_msg("Sensitivity models not estimable or failed: ", sum(design_diagnostics$status != "FITTED" & design_diagnostics$analysis_id !=
  "REFERENCE_primary_all_samples"), "/5")

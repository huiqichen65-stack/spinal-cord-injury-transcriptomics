source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

STEP45_DIR <- file.path(ROOT, "results", "45_Bulk_batch_structure_and_analysis_object")

INPUT_RDS <- file.path(STEP45_DIR, "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "46_Bulk_limma_temporal_DEG_analysis")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "46_Bulk_limma_temporal_DEG_analysis_v06"

num_name <- function(num, suffix) {
  paste0(sprintf("%02d", num), "_", PREFIX, "_", suffix)
}

ADJ_P_CUTOFF <- 0.05
ABS_LOGFC_CUTOFF <- 0.5

CORE_GENES <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

REVIEW_FLAGGED_SAMPLES <- c("jerry.1.8", "jerry.1.10", "D28.2.19")

required_pkgs <- c("limma")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall before rerunning.")
}

suppressPackageStartupMessages({
  library(limma)
})

safe_gene_symbol <- clean_gene_symbol

add_gene_columns <- function(tt) {
    if (is.null(rownames(tt))) {
        stop("Result table has no row names.")
    }
    out <- data.frame(gene_id = rownames(tt), gene_symbol_restored = safe_gene_symbol(rownames(tt)), 
        tt, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(out) <- NULL
    out
}

annotate_deg <- function(tt) {
    if (!"adj.P.Val" %in% colnames(tt)) {
        return(tt)
    }
    if (!"logFC" %in% colnames(tt)) {
        tt$stat_sig <- tt$adj.P.Val < ADJ_P_CUTOFF
        return(tt)
    }
    tt$stat_sig <- (tt$adj.P.Val < ADJ_P_CUTOFF)
    tt$effect_sig <- (tt$adj.P.Val < ADJ_P_CUTOFF & abs(tt$logFC) >= ABS_LOGFC_CUTOFF)
    tt$direction <- "NS"
    tt$direction[tt$adj.P.Val < ADJ_P_CUTOFF & tt$logFC >= ABS_LOGFC_CUTOFF] <- "Up"
    tt$direction[tt$adj.P.Val < ADJ_P_CUTOFF & tt$logFC <= -ABS_LOGFC_CUTOFF] <- "Down"
    tt
}

write_csv <- function(x, filename) {
  write.csv(x, file.path(OUT_DIR, filename), row.names = FALSE, fileEncoding = "UTF-8")
}

safe_rbind <- function(lst) {
    lst <- Filter(function(x) {
        !is.null(x) && is.data.frame(x)
    }, lst)
    if (length(lst) == 0) {
        return(data.frame())
    }
    all_names <- unique(unlist(lapply(lst, colnames), use.names = FALSE))
    lst2 <- lapply(lst, function(x) {
        miss <- setdiff(all_names, colnames(x))
        if (length(miss) > 0) {
            for (nm in miss) {
                x[[nm]] <- NA
            }
        }
        x <- x[, all_names, drop = FALSE]
        rownames(x) <- NULL
        x
    })
    out <- do.call(rbind, lst2)
    rownames(out) <- NULL
    out
}

design_rank <- function(design, model_name) {
    q <- qr(design)
    data.frame(model = model_name, n_samples = nrow(design), n_parameters = ncol(design), rank = q$rank, 
        full_rank = q$rank == ncol(design), stringsAsFactors = FALSE)
}

fit_limma_weighted <- function(expr_mat, design, model_label) {
    if (ncol(expr_mat) != nrow(design)) {
        stop("Sample-number mismatch in ", model_label)
    }
    if (qr(design)$rank != ncol(design)) {
        stop("Design is not full rank for ", model_label)
    }
    aw <- limma::arrayWeights(expr_mat, design = design)
    fit <- limma::lmFit(expr_mat, design = design, weights = aw)
    list(fit = fit, array_weights = aw)
}

fit_contrast_ebayes <- function(fit, contrast_matrix) {
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    limma::eBayes(fit2, trend = TRUE, robust = TRUE)
}

get_single_coef_table <- function(fit, coef_name) {
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = "BH")
    tt <- add_gene_columns(tt)
    annotate_deg(tt)
}

count_degs <- function(tt, analysis_name) {
    data.frame(analysis = analysis_name, n_tested = nrow(tt), n_adjP_lt_0.05 = sum(tt$adj.P.Val < ADJ_P_CUTOFF, 
        na.rm = TRUE), n_up_adjP_lt_0.05_logFC_ge_0.5 = sum(tt$adj.P.Val < ADJ_P_CUTOFF & tt$logFC >= 
        ABS_LOGFC_CUTOFF, na.rm = TRUE), n_down_adjP_lt_0.05_logFC_le_minus0.5 = sum(tt$adj.P.Val < ADJ_P_CUTOFF & 
        tt$logFC <= -ABS_LOGFC_CUTOFF, na.rm = TRUE), stringsAsFactors = FALSE)
}

extract_core <- function(tt, analysis_name) {
    idx <- match(CORE_GENES, tt$gene_symbol_restored)
    out <- data.frame(analysis = analysis_name, requested_gene = CORE_GENES, found = !is.na(idx), stringsAsFactors = FALSE)
    cols_keep <- intersect(c("gene_id", "gene_symbol_restored", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", 
        "B", "direction"), colnames(tt))
    tmp <- as.data.frame(matrix(NA, nrow = length(CORE_GENES), ncol = length(cols_keep)), stringsAsFactors = FALSE)
    colnames(tmp) <- cols_keep
    for (i in seq_along(idx)) {
        if (!is.na(idx[i])) {
            tmp[i, ] <- tt[idx[i], cols_keep, drop = FALSE]
        }
    }
    cbind(out, tmp)
}

build_time_design <- function(meta_sub) {
    meta_sub <- droplevels(meta_sub)
    design <- model.matrix(~0 + timepoint + batch, data = meta_sub)
    cn <- colnames(design)
    cn <- gsub("^timepointUninjured$", "T_Uninjured", cn)
    cn <- gsub("^timepoint15 min$", "T_15min", cn)
    cn <- gsub("^timepoint1 dpi$", "T_1dpi", cn)
    cn <- gsub("^timepoint3 dpi$", "T_3dpi", cn)
    cn <- gsub("^timepoint7 dpi$", "T_7dpi", cn)
    cn <- gsub("^timepoint14 dpi$", "T_14dpi", cn)
    cn <- gsub("^timepoint28 dpi$", "T_28dpi", cn)
    cn <- gsub("^timepoint42 dpi$", "T_42dpi", cn)
    cn <- gsub("^batchBatch_B$", "Batch_B", cn)
    colnames(design) <- cn
    design
}

if (!file.exists(INPUT_RDS)) {
  stop("Cannot find Step45 analysis object:\n", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

required_elements <- c("expression_original", "metadata")

if (!all(required_elements %in% names(obj))) {
  stop("Step45 RDS does not contain required elements.")
}

expr <- as.matrix(obj$expression_original)

storage.mode(expr) <- "numeric"

meta <- obj$metadata

if (!all(colnames(expr) %in% rownames(meta))) {
  stop("Expression sample IDs do not match metadata.")
}

meta <- meta[colnames(expr), , drop = FALSE]

if (!all(c("sample_id", "batch", "timepoint", "phase", "phase_primary") %in% colnames(meta))) {
  stop("Step45 metadata is missing required fields.")
}

meta$sample_id <- rownames(meta)

ids_all <- rownames(meta)

ids_primary <- rownames(meta)[as.character(meta$phase) %in% c("Early", "Established late")]

ids_transition <- rownames(meta)[as.character(meta$phase) %in% c("Early", "Transition")]

ids_extended <- rownames(meta)[as.character(meta$phase) %in% c("Early", "Extended late")]

ids_14_28 <- rownames(meta)[as.character(meta$timepoint) %in% c("14 dpi", "28 dpi")]

sample_set_summary <- data.frame(analysis = c("Full temporal model", "Primary: Established late vs Early",
  "Transition: 7 dpi vs Early", "Extended late: 42 dpi vs Early", "Sensitivity: 28 dpi vs 14 dpi"),
  n_samples = c(length(ids_all), length(ids_primary), length(ids_transition), length(ids_extended),
    length(ids_14_28)), stringsAsFactors = FALSE)

write_csv(sample_set_summary, num_name(13, "sample_set_summary.csv"))

primary_balance <- as.data.frame(table(Phase = meta[ids_primary, "phase"], Batch = meta[ids_primary,
  "batch"]))

write_csv(primary_balance, num_name(14, "primary_phase_by_batch.csv"))

meta_all <- droplevels(meta[ids_all, , drop = FALSE])

expr_all <- expr[, ids_all, drop = FALSE]

design_time <- build_time_design(meta_all)

design_time_diag <- design_rank(design_time, "Full 8-time-point: ~ 0 + timepoint + batch")

if (!design_time_diag$full_rank) {
  stop("Full temporal design is not full rank.")
}

fit_time_base <- fit_limma_weighted(expr_all, design_time, "Full temporal model")

aw_time <- data.frame(sample_id = colnames(expr_all), batch = as.character(meta_all$batch), timepoint = as.character(meta_all$timepoint),
  array_weight = as.numeric(fit_time_base$array_weights), stringsAsFactors = FALSE)

write_csv(aw_time, num_name(15, "arrayWeights_full_time_model.csv"))

time_contrast_matrix <- limma::makeContrasts(`15min_vs_Uninjured` = T_15min - T_Uninjured, `1dpi_vs_Uninjured` = T_1dpi -
  T_Uninjured, `3dpi_vs_Uninjured` = T_3dpi - T_Uninjured, `7dpi_vs_Uninjured` = T_7dpi - T_Uninjured,
  `14dpi_vs_Uninjured` = T_14dpi - T_Uninjured, `28dpi_vs_Uninjured` = T_28dpi - T_Uninjured, `42dpi_vs_Uninjured` = T_42dpi -
    T_Uninjured, `28dpi_vs_14dpi` = T_28dpi - T_14dpi, levels = design_time)

fit_time <- fit_contrast_ebayes(fit_time_base$fit, time_contrast_matrix)

time_contrast_names <- colnames(time_contrast_matrix)

time_tables <- list()

for (cn in time_contrast_names) {
    tt <- get_single_coef_table(fit_time, cn)
    time_tables[[cn]] <- tt
    write_csv(tt, num_name(20 + match(cn, time_contrast_names), paste0(cn, ".csv")))
}

overall_coef_names <- c("15min_vs_Uninjured", "1dpi_vs_Uninjured", "3dpi_vs_Uninjured", "7dpi_vs_Uninjured",
  "14dpi_vs_Uninjured", "28dpi_vs_Uninjured", "42dpi_vs_Uninjured")

tt_F <- limma::topTable(fit_time, coef = overall_coef_names, number = Inf, sort.by = "F", adjust.method = "BH")

tt_F <- add_gene_columns(tt_F)

tt_F$stat_sig <- (tt_F$adj.P.Val < ADJ_P_CUTOFF)

write_csv(tt_F, num_name(16, "overall_temporal_Ftest.csv"))

meta_primary <- droplevels(meta[ids_primary, , drop = FALSE])

expr_primary <- expr[, ids_primary, drop = FALSE]

primary_group <- ifelse(as.character(meta_primary$phase) == "Early", "Early", "EstablishedLate")

meta_primary$primary_group <- factor(primary_group, levels = c("Early", "EstablishedLate"))

design_primary <- model.matrix(~0 + primary_group + batch, data = meta_primary)

colnames(design_primary) <- c("Early", "EstablishedLate", "Batch_B")

primary_diag <- design_rank(design_primary, "Primary: Established late vs Early")

if (!primary_diag$full_rank) {
  stop("Primary design is not full rank.")
}

fit_primary_base <- fit_limma_weighted(expr_primary, design_primary, "Primary model")

aw_primary <- data.frame(sample_id = colnames(expr_primary), batch = as.character(meta_primary$batch),
  phase = as.character(meta_primary$phase), array_weight = as.numeric(fit_primary_base$array_weights),
  stringsAsFactors = FALSE)

write_csv(aw_primary, num_name(17, "arrayWeights_primary_model.csv"))

primary_contrast <- limma::makeContrasts(EstablishedLate_vs_Early = EstablishedLate - Early, levels = design_primary)

fit_primary <- fit_contrast_ebayes(fit_primary_base$fit, primary_contrast)

tt_primary <- get_single_coef_table(fit_primary, "EstablishedLate_vs_Early")

write_csv(tt_primary, num_name(4, "PRIMARY_EstablishedLate_vs_Early.csv"))

meta_transition <- droplevels(meta[ids_transition, , drop = FALSE])

expr_transition <- expr[, ids_transition, drop = FALSE]

meta_transition$two_group <- factor(ifelse(as.character(meta_transition$phase) == "Transition", "Target",
  "Early"), levels = c("Early", "Target"))

design_transition <- model.matrix(~0 + two_group + batch, data = meta_transition)

colnames(design_transition) <- c("Early", "Target", "Batch_B")

transition_diag <- design_rank(design_transition, "Transition: 7 dpi vs Early")

if (!transition_diag$full_rank) {
  stop("Transition design is not full rank.")
}

fit_transition_base <- fit_limma_weighted(expr_transition, design_transition, "Transition model")

contrast_transition <- limma::makeContrasts(Transition7dpi_vs_Early = Target - Early, levels = design_transition)

fit_transition <- fit_contrast_ebayes(fit_transition_base$fit, contrast_transition)

tt_transition <- get_single_coef_table(fit_transition, "Transition7dpi_vs_Early")

write_csv(tt_transition, num_name(5, "Transition7dpi_vs_Early.csv"))

meta_extended <- droplevels(meta[ids_extended, , drop = FALSE])

expr_extended <- expr[, ids_extended, drop = FALSE]

meta_extended$two_group <- factor(ifelse(as.character(meta_extended$phase) == "Extended late", "Target",
  "Early"), levels = c("Early", "Target"))

design_extended <- model.matrix(~0 + two_group + batch, data = meta_extended)

colnames(design_extended) <- c("Early", "Target", "Batch_B")

extended_diag <- design_rank(design_extended, "Extended late: 42 dpi vs Early")

if (!extended_diag$full_rank) {
  stop("Extended-late design is not full rank.")
}

fit_extended_base <- fit_limma_weighted(expr_extended, design_extended, "Extended-late model")

contrast_extended <- limma::makeContrasts(Extended42dpi_vs_Early = Target - Early, levels = design_extended)

fit_extended <- fit_contrast_ebayes(fit_extended_base$fit, contrast_extended)

tt_extended <- get_single_coef_table(fit_extended, "Extended42dpi_vs_Early")

write_csv(tt_extended, num_name(6, "Extended42dpi_vs_Early.csv"))

meta_14_28 <- droplevels(meta[ids_14_28, , drop = FALSE])

expr_14_28 <- expr[, ids_14_28, drop = FALSE]

meta_14_28$late_time <- factor(as.character(meta_14_28$timepoint), levels = c("14 dpi", "28 dpi"))

design_14_28 <- model.matrix(~0 + late_time + batch, data = meta_14_28)

colnames(design_14_28) <- c("T14", "T28", "Batch_B")

late_diag <- design_rank(design_14_28, "Sensitivity: 28 dpi vs 14 dpi")

if (!late_diag$full_rank) {
  stop("14-vs-28 design is not full rank.")
}

fit_14_28_base <- fit_limma_weighted(expr_14_28, design_14_28, "14 vs 28 model")

contrast_14_28 <- limma::makeContrasts(`28dpi_vs_14dpi` = T28 - T14, levels = design_14_28)

fit_14_28 <- fit_contrast_ebayes(fit_14_28_base$fit, contrast_14_28)

tt_14_28 <- get_single_coef_table(fit_14_28, "28dpi_vs_14dpi")

write_csv(tt_14_28, num_name(7, "Sensitivity_28dpi_vs_14dpi.csv"))

ids_primary_sensitivity <- setdiff(ids_primary, REVIEW_FLAGGED_SAMPLES)

actually_excluded_primary <- intersect(ids_primary, REVIEW_FLAGGED_SAMPLES)

meta_primary_sens <- droplevels(meta[ids_primary_sensitivity, , drop = FALSE])

expr_primary_sens <- expr[, ids_primary_sensitivity, drop = FALSE]

meta_primary_sens$primary_group <- factor(ifelse(as.character(meta_primary_sens$phase) == "Early", "Early",
  "EstablishedLate"), levels = c("Early", "EstablishedLate"))

design_primary_sens <- model.matrix(~0 + primary_group + batch, data = meta_primary_sens)

colnames(design_primary_sens) <- c("Early", "EstablishedLate", "Batch_B")

primary_sens_diag <- design_rank(design_primary_sens, "Primary sensitivity excluding Step45 review flags")

if (!primary_sens_diag$full_rank) {
  stop("Primary sensitivity design is not full rank.")
}

fit_primary_sens_base <- fit_limma_weighted(expr_primary_sens, design_primary_sens, "Primary sensitivity model")

contrast_primary_sens <- limma::makeContrasts(EstablishedLate_vs_Early = EstablishedLate - Early, levels = design_primary_sens)

fit_primary_sens <- fit_contrast_ebayes(fit_primary_sens_base$fit, contrast_primary_sens)

tt_primary_sens <- get_single_coef_table(fit_primary_sens, "EstablishedLate_vs_Early")

write_csv(tt_primary_sens, num_name(8, "PRIMARY_sensitivity_excluding_review_flags.csv"))

merge_cols <- c("gene_id", "gene_symbol_restored", "logFC", "t", "P.Value", "adj.P.Val", "direction")

primary_merge <- tt_primary[, merge_cols, drop = FALSE]

sens_merge <- tt_primary_sens[, merge_cols, drop = FALSE]

colnames(primary_merge)[-(1:2)] <- paste0(colnames(primary_merge)[-(1:2)], "_primary")

colnames(sens_merge)[-(1:2)] <- paste0(colnames(sens_merge)[-(1:2)], "_sensitivity")

concordance_df <- merge(primary_merge, sens_merge, by = c("gene_id", "gene_symbol_restored"), all = FALSE)

write_csv(concordance_df, num_name(9, "PRIMARY_sensitivity_gene_concordance.csv"))

concordance_summary <- data.frame(metric = c("n_common_genes", "Pearson_logFC", "Spearman_logFC", "Pearson_t",
  "Spearman_t", "same_logFC_direction_fraction", "primary_effect_sig_n", "sensitivity_effect_sig_n",
  "effect_sig_overlap_n"), value = c(nrow(concordance_df), cor(concordance_df$logFC_primary, concordance_df$logFC_sensitivity,
  method = "pearson", use = "complete.obs"), cor(concordance_df$logFC_primary, concordance_df$logFC_sensitivity,
  method = "spearman", use = "complete.obs"), cor(concordance_df$t_primary, concordance_df$t_sensitivity,
  method = "pearson", use = "complete.obs"), cor(concordance_df$t_primary, concordance_df$t_sensitivity,
  method = "spearman", use = "complete.obs"), mean(sign(concordance_df$logFC_primary) == sign(concordance_df$logFC_sensitivity),
  na.rm = TRUE), sum(tt_primary$effect_sig, na.rm = TRUE), sum(tt_primary_sens$effect_sig, na.rm = TRUE),
  sum(tt_primary$gene_id[tt_primary$effect_sig] %in% tt_primary_sens$gene_id[tt_primary_sens$effect_sig])),
  stringsAsFactors = FALSE)

write_csv(concordance_summary, num_name(3, "PRIMARY_sensitivity_concordance_summary.csv"))

deg_count_list <- list(count_degs(tt_primary, "PRIMARY_EstablishedLate_vs_Early"), count_degs(tt_transition,
  "Transition7dpi_vs_Early"), count_degs(tt_extended, "Extended42dpi_vs_Early"), count_degs(tt_14_28,
  "Sensitivity_28dpi_vs_14dpi"))

for (nm in overall_coef_names) {
  deg_count_list[[length(deg_count_list) + 1]] <- count_degs(time_tables[[nm]], nm)
}

deg_count_summary <- safe_rbind(deg_count_list)

expected_deg_rows <- 4 + length(overall_coef_names)

if (nrow(deg_count_summary) != expected_deg_rows) {
  stop("DEG summary integrity check failed: expected ", expected_deg_rows, " rows but obtained ", nrow(deg_count_summary),
    ".")
}

write_csv(deg_count_summary, num_name(1, "DEG_count_summary.csv"))

core_tables <- list(extract_core(tt_primary, "PRIMARY_EstablishedLate_vs_Early"), extract_core(tt_primary_sens,
  "PRIMARY_sensitivity_excluding_review_flags"), extract_core(tt_transition, "Transition7dpi_vs_Early"),
  extract_core(tt_extended, "Extended42dpi_vs_Early"), extract_core(tt_14_28, "Sensitivity_28dpi_vs_14dpi"))

for (nm in overall_coef_names) {
  core_tables[[length(core_tables) + 1]] <- extract_core(time_tables[[nm]], nm)
}

core_stats <- safe_rbind(core_tables)

expected_core_rows <- length(CORE_GENES) * (5 + length(overall_coef_names))

if (nrow(core_stats) != expected_core_rows) {
  stop("Core-gene summary integrity check failed: expected ", expected_core_rows, " rows but obtained ",
    nrow(core_stats), ".")
}

write_csv(core_stats, num_name(2, "core_gene_all_contrasts.csv"))

core_idx <- match(CORE_GENES, safe_gene_symbol(rownames(expr)))

core_time_means <- list()

for (g_i in seq_along(CORE_GENES)) {
    if (is.na(core_idx[g_i])) {
        next
    }
    vals <- expr[core_idx[g_i], , drop = TRUE]
    tmp <- data.frame(sample_id = colnames(expr), value = vals, timepoint = as.character(meta[colnames(expr), 
        "timepoint"]), batch = as.character(meta[colnames(expr), "batch"]), stringsAsFactors = FALSE)
    agg_mean <- aggregate(value ~ timepoint, data = tmp, FUN = mean)
    agg_sd <- aggregate(value ~ timepoint, data = tmp, FUN = sd)
    agg_n <- aggregate(value ~ timepoint, data = tmp, FUN = length)
    colnames(agg_mean)[2] <- "mean_expression"
    colnames(agg_sd)[2] <- "sd_expression"
    colnames(agg_n)[2] <- "n"
    one <- Reduce(function(x, y) {
        merge(x, y, by = "timepoint", all = TRUE)
    }, list(agg_mean, agg_sd, agg_n))
    one$gene <- CORE_GENES[g_i]
    one$time_order <- match(one$timepoint, c("Uninjured", "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", 
        "28 dpi", "42 dpi"))
    one <- one[order(one$time_order), , drop = FALSE]
    core_time_means[[length(core_time_means) + 1]] <- one
}

core_time_means <- safe_rbind(core_time_means)

write_csv(core_time_means, num_name(10, "core_gene_raw_timepoint_means.csv"))

design_diagnostics <- safe_rbind(list(design_time_diag, primary_diag, transition_diag, extended_diag,
  late_diag, primary_sens_diag))

write_csv(design_diagnostics, num_name(11, "design_diagnostics.csv"))

contrast_notes <- data.frame(analysis = c("Overall temporal F-test", "15min_vs_Uninjured", "1dpi_vs_Uninjured",
  "3dpi_vs_Uninjured", "7dpi_vs_Uninjured", "14dpi_vs_Uninjured", "28dpi_vs_Uninjured", "42dpi_vs_Uninjured",
  "PRIMARY_EstablishedLate_vs_Early", "Transition7dpi_vs_Early", "Extended42dpi_vs_Early", "Sensitivity_28dpi_vs_14dpi"),
  role = c("Global temporal test", "Exploratory temporal contrast", "Exploratory temporal contrast",
    "Temporal contrast", "Temporal/transition contrast", "Temporal contrast", "Temporal contrast",
    "Temporal/persistence contrast", "PRIMARY biological contrast", "Secondary transition analysis",
    "Secondary persistence analysis", "Sensitivity comparison"), caution = c("Batch-adjusted additive model",
    "Uninjured is Batch_A only; 15 min is Batch_B only. Interpret cautiously.", "Uninjured is Batch_A only; 1 dpi is Batch_B only. Interpret cautiously.",
    "Both batches present at 3 dpi, but Uninjured is Batch_A only.", "Both batches present at 7 dpi, but Uninjured is Batch_A only.",
    "Both batches present at 14 dpi, but Uninjured is Batch_A only.", "Both batches present at 28 dpi, but Uninjured is Batch_A only.",
    "Both batches present at 42 dpi, but Uninjured is Batch_A only.", "Primary phase comparison is batch-adjusted and phase groups contain both batches.",
    "Early group contains both batches; 7 dpi contains both batches.", "Early group contains both batches; 42 dpi contains both batches.",
    "14 and 28 dpi both contain Batch_A and Batch_B."), stringsAsFactors = FALSE)

write_csv(contrast_notes, num_name(12, "contrast_interpretation_notes.csv"))

step46_object <- list(input_step45 = INPUT_RDS, expression_original = expr, metadata = meta, full_time = list(design = design_time,
  array_weights = fit_time_base$array_weights, contrast_matrix = time_contrast_matrix, fit = fit_time,
  overall_F_table = tt_F, contrast_tables = time_tables), primary = list(sample_ids = ids_primary,
  design = design_primary, array_weights = fit_primary_base$array_weights, fit = fit_primary, table = tt_primary),
  transition = list(sample_ids = ids_transition, design = design_transition, fit = fit_transition,
    table = tt_transition), extended_late = list(sample_ids = ids_extended, design = design_extended,
    fit = fit_extended, table = tt_extended), late_14_vs_28 = list(sample_ids = ids_14_28, design = design_14_28,
    fit = fit_14_28, table = tt_14_28), primary_sensitivity = list(excluded_review_flags = actually_excluded_primary,
    sample_ids = ids_primary_sensitivity, design = design_primary_sens, fit = fit_primary_sens, table = tt_primary_sens,
    concordance_summary = concordance_summary), thresholds = list(adj_P_Val = ADJ_P_CUTOFF, abs_logFC = ABS_LOGFC_CUTOFF),
  core_genes = CORE_GENES, notes = c("Formal testing used the original transformed expression matrix.",
    "DESeq2 was not used because the matrix is not raw counts.", "All formal models include batch.",
    "arrayWeights were used to down-weight lower-quality samples without deleting them.", "eBayes was run with trend=TRUE and robust=TRUE.",
    "Primary contrast: Established late (14+28 dpi) vs Early (Uninjured+15 min+1+3 dpi).", "7 dpi is treated as Transition.",
    "42 dpi is treated as Extended late/persistence.", "Primary sensitivity analysis excludes Step45 review-flagged samples only for robustness assessment.",
    "No manuscript figure was generated in Step46."))

saveRDS(step46_object, file.path(OUT_DIR, num_name(30, "analysis_object.rds")), compress = FALSE)

output_index <- data.frame(No = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
  "13", "14", "15", "16", "17", sprintf("%02d", 21:(20 + length(time_contrast_names))), "30", "31",
  "32"), file = c(num_name(1, "DEG_count_summary.csv"), num_name(2, "core_gene_all_contrasts.csv"),
  num_name(3, "PRIMARY_sensitivity_concordance_summary.csv"), num_name(4, "PRIMARY_EstablishedLate_vs_Early.csv"),
  num_name(5, "Transition7dpi_vs_Early.csv"), num_name(6, "Extended42dpi_vs_Early.csv"), num_name(7,
    "Sensitivity_28dpi_vs_14dpi.csv"), num_name(8, "PRIMARY_sensitivity_excluding_review_flags.csv"),
  num_name(9, "PRIMARY_sensitivity_gene_concordance.csv"), num_name(10, "core_gene_raw_timepoint_means.csv"),
  num_name(11, "design_diagnostics.csv"), num_name(12, "contrast_interpretation_notes.csv"), num_name(13,
    "sample_set_summary.csv"), num_name(14, "primary_phase_by_batch.csv"), num_name(15, "arrayWeights_full_time_model.csv"),
  num_name(16, "overall_temporal_Ftest.csv"), num_name(17, "arrayWeights_primary_model.csv"), vapply(seq_along(time_contrast_names),
    function(i) num_name(20 + i, paste0(time_contrast_names[i], ".csv")), character(1)), num_name(30,
    "analysis_object.rds"), num_name(31, "REPORT.txt"), num_name(32, "sessionInfo.txt")), purpose = c("DEG count summary",
  "Core genes across contrasts", "Primary sensitivity concordance summary", "Primary Established late vs Early full table",
  "7 dpi vs Early", "42 dpi vs Early", "28 dpi vs 14 dpi", "Primary sensitivity full table", "Gene-level primary/sensitivity concordance",
  "Core-gene time-point means", "Design diagnostics", "Contrast interpretation notes", "Sample-set summary",
  "Primary phase by batch", "Full-time array weights", "Overall temporal F-test", "Primary-model array weights",
  rep("Time-point contrast table", length(time_contrast_names)), "Step46 analysis object", "Step46 report",
  "R session info"), stringsAsFactors = FALSE)

write.csv(output_index, file.path(OUT_DIR, paste0("00_", PREFIX, "_OUTPUT_INDEX.csv")), row.names = FALSE,
  fileEncoding = "UTF-8")

primary_counts <- count_degs(tt_primary, "PRIMARY")

transition_counts <- count_degs(tt_transition, "Transition")

extended_counts <- count_degs(tt_extended, "Extended")

late_counts <- count_degs(tt_14_28, "28vs14")

writeLines("Numerical results are provided in the accompanying output tables.", 
    con = file.path(OUT_DIR, num_name(31, "REPORT.txt")))

capture.output(sessionInfo(), file = file.path(OUT_DIR, num_name(32, "sessionInfo.txt")))

print(table(meta_primary$primary_group, meta_primary$batch))

print(primary_counts, row.names = FALSE)

print(concordance_summary, row.names = FALSE)

print(extract_core(tt_primary, "PRIMARY")[, intersect(c("requested_gene", "found", "logFC", "t", "P.Value",
  "adj.P.Val", "direction"), colnames(extract_core(tt_primary, "PRIMARY"))), drop = FALSE], row.names = FALSE)

cat("1) ", file.path(OUT_DIR, num_name(1, "DEG_count_summary.csv")), "\n", sep = "")

cat("2) ", file.path(OUT_DIR, num_name(2, "core_gene_all_contrasts.csv")), "\n", sep = "")

cat("3) ", file.path(OUT_DIR, num_name(3, "PRIMARY_sensitivity_concordance_summary.csv")), "\n", sep = "")

cat("4) ", file.path(OUT_DIR, num_name(4, "PRIMARY_EstablishedLate_vs_Early.csv")), "\n", sep = "")

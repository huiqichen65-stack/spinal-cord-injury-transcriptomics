rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
set.seed(20260930)

ROOT <- PROJECT_ROOT

STEP45_RDS <- file.path(ROOT, "results", "45_Bulk_batch_structure_and_analysis_object", "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds")

STEP47_RDS <- file.path(ROOT, "results", "47_Bulk_scRNA_Clec7a_signature_GSEA", "30_47_Bulk_scRNA_Clec7a_signature_GSEA_v02_analysis_object.rds")

STEP48_RDS <- file.path(ROOT, "results", "48_Bulk_functional_enrichment", "10_48_Bulk_functional_enrichment_v04_analysis_object.rds")

STEP49_RDS <- file.path(ROOT, "results", "49_Bulk_WGCNA_network_modules", "09_49_Bulk_WGCNA_network_modules_v01_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "50_Bulk_WGCNA_module_trait_integration")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "50_Bulk_WGCNA_module_trait_integration_v01"

EXPECTED_SAMPLES <- 39L
EXPECTED_GENES <- 12977L
EXPECTED_SIGNATURE_N <- 1320L
EXPECTED_NON_GREY_MODULES <- 11L
EXPECTED_SOFT_POWER <- 14L

TIME_LEVELS <- c("Uninjured", "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi")

PHASE_LEVELS <- c("Early", "Transition", "Established late", "Extended late")

CORE_GENES <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

FDR_CUTOFF <- 0.05
MIN_GS_SIZE <- 10L
MAX_GS_SIZE <- 500L

install_cran_if_missing <- require_cran_package

install_bioc_if_missing <- require_bioc_package

install_cran_if_missing("statmod")

for (p in c("limma", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi")) {
  install_bioc_if_missing(p)
}

suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

write_csv <- function(x, filename) {
  utils::write.csv(x, file = file.path(OUT_DIR, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

split_symbol_string <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
        return(character(0))
    }
    y <- unlist(strsplit(x, "/", fixed = TRUE), use.names = FALSE)
    unique(y[!is.na(y) & nzchar(y)])
}

empty_enrichment_table <- function() {
  data.frame(ID = character(0), Description = character(0), GeneRatio = character(0), BgRatio = character(0),
    pvalue = numeric(0), p.adjust = numeric(0), qvalue = numeric(0), geneID = character(0), Count = integer(0),
    geneID_SYMBOL = character(0), stringsAsFactors = FALSE)
}

for (f in c(STEP45_RDS, STEP47_RDS, STEP48_RDS, STEP49_RDS)) {
  if (!file.exists(f)) {
    stop("Required upstream object not found:\n", f)
  }
}

obj45 <- readRDS(STEP45_RDS)
obj47 <- readRDS(STEP47_RDS)
obj48 <- readRDS(STEP48_RDS)
obj49 <- readRDS(STEP49_RDS)

required49 <- c("metadata", "module_eigengenes", "gene_module_assignment", "module_size_summary", "soft_threshold_selection",
  "step45_review_flag_by_sample")

missing49 <- setdiff(required49, names(obj49))

if (length(missing49) > 0L) {
  stop("Step49 analysis object lacks required field(s): ", paste(missing49, collapse = ", "))
}

meta <- obj49$metadata
MEs_all <- as.data.frame(obj49$module_eigengenes, check.names = FALSE)
gene_module <- obj49$gene_module_assignment
module_size <- obj49$module_size_summary

if (nrow(meta) != EXPECTED_SAMPLES) {
  stop("Step49 metadata sample count mismatch: ", nrow(meta), " vs expected ", EXPECTED_SAMPLES)
}

if (nrow(gene_module) != EXPECTED_GENES) {
  stop("Step49 gene-module table count mismatch: ", nrow(gene_module), " vs expected ", EXPECTED_GENES)
}

non_grey_modules <- as.character(module_size$module_color[as.logical(module_size$assigned_module)])

non_grey_modules <- unique(non_grey_modules[non_grey_modules != "grey"])

if (length(non_grey_modules) != EXPECTED_NON_GREY_MODULES) {
  stop("Unexpected non-grey module count: ", length(non_grey_modules), " vs ", EXPECTED_NON_GREY_MODULES)
}

soft_power <- as.integer(obj49$soft_threshold_selection$selected_power[1])

if (length(soft_power) != 1L || is.na(soft_power) || soft_power != EXPECTED_SOFT_POWER) {
  stop("Step49 soft power mismatch. Observed ", soft_power, "; expected power ", EXPECTED_SOFT_POWER)
}

core_rows <- gene_module[gene_module$gene_symbol_restored %in% CORE_GENES, , drop = FALSE]

if (nrow(core_rows) != length(CORE_GENES)) {
  stop("Step49 core-gene audit mismatch; expected all six core genes.")
}

core_modules <- unique(core_rows$module_color)

if (length(core_modules) != 1L || core_modules[1] != "brown") {
  stop("Step49 core-gene module audit changed unexpectedly. ", "Expected all six core genes in brown.")
}

if (is.null(rownames(MEs_all))) {
  stop("Step49 module eigengenes lack sample row names.")
}

if (!"sample_id" %in% names(meta)) {
  stop("Step49 metadata lacks sample_id.")
}

meta <- meta[match(rownames(MEs_all), meta$sample_id), , drop = FALSE]

if (anyNA(meta$sample_id) || !identical(meta$sample_id, rownames(MEs_all))) {
  stop("Step49 module eigengene / metadata alignment failed.")
}

me_cols <- paste0("ME", non_grey_modules)

missing_me_cols <- setdiff(me_cols, colnames(MEs_all))

if (length(missing_me_cols) > 0L) {
  stop("Missing expected module eigengene column(s): ", paste(missing_me_cols, collapse = ", "))
}

MEs <- MEs_all[, me_cols, drop = FALSE]

if (!"signature_genes" %in% names(obj47)) {
  stop("Step47 object lacks signature_genes.")
}

signature_genes <- unique(as.character(obj47$signature_genes))

signature_genes <- signature_genes[!is.na(signature_genes) & nzchar(signature_genes)]

if (length(signature_genes) != EXPECTED_SIGNATURE_N) {
  stop("Frozen Step47 signature size mismatch: ", length(signature_genes), " vs expected ", EXPECTED_SIGNATURE_N)
}

if (!"result_tables" %in% names(obj48) || !"key_pathway_summary" %in% names(obj48$result_tables)) {
  stop("Step48 object lacks result_tables$key_pathway_summary.")
}

key48 <- obj48$result_tables$key_pathway_summary

if (!"expression_original" %in% names(obj45)) {
  stop("Step45 object lacks expression_original.")
}

expr45 <- obj45$expression_original

if (nrow(expr45) != EXPECTED_GENES || ncol(expr45) != EXPECTED_SAMPLES) {
  stop("Step45 expression dimensions changed unexpectedly: ", nrow(expr45), " x ", ncol(expr45))
}

review_flag <- obj49$step45_review_flag_by_sample

if (is.null(names(review_flag))) {
  stop("Step49 review-flag vector lacks sample names.")
}

review_flag <- as.logical(review_flag[meta$sample_id])

names(review_flag) <- meta$sample_id

if (anyNA(review_flag)) {
  stop("Step45 review flags could not be aligned to all Step50 samples.")
}

run_binary_ME_comparison <- function(reference_label, target_label, group_vector, sample_keep, comparison_name) {
    sample_ids <- meta$sample_id[sample_keep]
    if (length(sample_ids) < 8L) {
        stop("Too few samples for comparison ", comparison_name)
    }
    m <- meta[sample_keep, , drop = FALSE]
    g <- factor(as.character(group_vector[sample_keep]), levels = c(reference_label, target_label))
    if (anyNA(g)) {
        stop("Unexpected group value in comparison ", comparison_name)
    }
    b <- factor(as.character(m$batch), levels = c("Batch_A", "Batch_B"))
    design <- model.matrix(~b + g)
    coef_candidates <- grep("^g", colnames(design), value = TRUE)
    if (length(coef_candidates) != 1L) {
        stop("Could not uniquely identify target coefficient for ", comparison_name, ". Design columns: ", 
            paste(colnames(design), collapse = ", "))
    }
    coef_name <- coef_candidates[1]
    me_matrix <- t(as.matrix(MEs[sample_ids, , drop = FALSE]))
    fit <- limma::lmFit(me_matrix, design)
    fit <- limma::eBayes(fit, trend = FALSE, robust = TRUE)
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "none")
    tt$ME <- rownames(tt)
    tt$module_color <- sub("^ME", "", tt$ME)
    out <- data.frame(comparison = comparison_name, reference = reference_label, target = target_label, 
        n_samples = length(sample_ids), n_reference = sum(g == reference_label), n_target = sum(g == 
            target_label), module_color = tt$module_color, estimate_target_minus_reference = tt$logFC, 
        AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, B = tt$B, stringsAsFactors = FALSE)
    out
}

run_global_time_test <- function() {
    m <- meta
    b <- factor(as.character(m$batch), levels = c("Batch_A", "Batch_B"))
    tp <- factor(as.character(m$timepoint), levels = TIME_LEVELS)
    design <- model.matrix(~b + tp)
    time_coef <- grep("^tp", colnames(design), value = TRUE)
    if (length(time_coef) != 7L) {
        stop("Unexpected number of time coefficients in global model: ", length(time_coef))
    }
    me_matrix <- t(as.matrix(MEs))
    fit <- limma::lmFit(me_matrix, design)
    fit <- limma::eBayes(fit, trend = FALSE, robust = TRUE)
    tt <- limma::topTable(fit, coef = time_coef, number = Inf, adjust.method = "BH", sort.by = "none")
    tt$ME <- rownames(tt)
    tt$module_color <- sub("^ME", "", tt$ME)
    data.frame(comparison = "Overall_timepoint_effect_adjusted_for_batch", n_samples = nrow(meta), module_color = tt$module_color, 
        F = tt$F, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
}

run_batch_audit <- function() {
    m <- meta
    b <- factor(as.character(m$batch), levels = c("Batch_A", "Batch_B"))
    tp <- factor(as.character(m$timepoint), levels = TIME_LEVELS)
    design <- model.matrix(~tp + b)
    coef_name <- "bBatch_B"
    if (!coef_name %in% colnames(design)) {
        stop("Could not identify batch coefficient in all-time batch audit.")
    }
    me_matrix <- t(as.matrix(MEs))
    fit <- limma::lmFit(me_matrix, design)
    fit <- limma::eBayes(fit, trend = FALSE, robust = TRUE)
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "none")
    tt$ME <- rownames(tt)
    tt$module_color <- sub("^ME", "", tt$ME)
    data.frame(comparison = "Batch_B_vs_Batch_A_adjusted_for_timepoint", n_samples = nrow(meta), module_color = tt$module_color, 
        estimate_BatchB_minus_BatchA = tt$logFC, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, 
        stringsAsFactors = FALSE)
}

phase_chr <- as.character(meta$phase)

time_chr <- as.character(meta$timepoint)

keep_primary <- phase_chr %in% c("Early", "Established late")

primary_full <- run_binary_ME_comparison(reference_label = "Early", target_label = "Established late",
  group_vector = phase_chr, sample_keep = keep_primary, comparison_name = "EstablishedLate_vs_Early")

keep_primary_sens <- keep_primary & !review_flag

primary_sens <- run_binary_ME_comparison(reference_label = "Early", target_label = "Established late",
  group_vector = phase_chr, sample_keep = keep_primary_sens, comparison_name = "EstablishedLate_vs_Early_excluding_Step45_review_samples")

keep_transition <- phase_chr %in% c("Early", "Transition")

transition_vs_early <- run_binary_ME_comparison(reference_label = "Early", target_label = "Transition",
  group_vector = phase_chr, sample_keep = keep_transition, comparison_name = "Transition_7dpi_vs_Early")

keep_extended <- phase_chr %in% c("Early", "Extended late")

extended_vs_early <- run_binary_ME_comparison(reference_label = "Early", target_label = "Extended late",
  group_vector = phase_chr, sample_keep = keep_extended, comparison_name = "ExtendedLate_42dpi_vs_Early")

keep_28_14 <- time_chr %in% c("14 dpi", "28 dpi")

d28_vs_d14 <- run_binary_ME_comparison(reference_label = "14 dpi", target_label = "28 dpi", group_vector = time_chr,
  sample_keep = keep_28_14, comparison_name = "28dpi_vs_14dpi")

global_time <- run_global_time_test()
batch_audit <- run_batch_audit()

secondary_module_trait <- rbind(transition_vs_early, extended_vs_early, d28_vs_d14)

primary_merge <- merge(primary_full, primary_sens[, c("module_color", "n_samples", "estimate_target_minus_reference",
  "t", "P.Value", "adj.P.Val"), drop = FALSE], by = "module_color", all = TRUE, suffixes = c("_primary",
  "_sensitivity"), sort = FALSE)

names(primary_merge)[names(primary_merge) == "n_samples_sensitivity"] <- "sensitivity_n_samples"

primary_merge$module_order <- match(primary_merge$module_color, non_grey_modules)

primary_merge <- primary_merge[order(primary_merge$module_order), , drop = FALSE]

primary_merge$module_order <- NULL

primary_sensitivity_pearson <- suppressWarnings(cor(primary_merge$estimate_target_minus_reference_primary,
  primary_merge$estimate_target_minus_reference_sensitivity, method = "pearson", use = "complete.obs"))

primary_sensitivity_spearman <- suppressWarnings(cor(primary_merge$estimate_target_minus_reference_primary,
  primary_merge$estimate_target_minus_reference_sensitivity, method = "spearman", use = "complete.obs"))

primary_sensitivity_same_direction <- mean(sign(primary_merge$estimate_target_minus_reference_primary) ==
  sign(primary_merge$estimate_target_minus_reference_sensitivity), na.rm = TRUE)

ME_long_list <- list()

for (module in non_grey_modules) {
    me_col <- paste0("ME", module)
    for (tp in TIME_LEVELS) {
        idx <- time_chr == tp
        v <- MEs[idx, me_col]
        ME_long_list[[length(ME_long_list) + 1L]] <- data.frame(module_color = module, timepoint = tp, 
            n = sum(idx), ME_mean = mean(v, na.rm = TRUE), ME_sd = stats::sd(v, na.rm = TRUE), ME_sem = stats::sd(v, 
                na.rm = TRUE)/sqrt(sum(idx)), stringsAsFactors = FALSE)
    }
}

module_time_summary <- do.call(rbind, ME_long_list)

symbol_module_split <- split(as.character(gene_module$module_color), as.character(gene_module$gene_symbol_restored))

ambiguous_symbols <- names(Filter(function(z) {
  length(unique(z[!is.na(z) & nzchar(z)])) > 1L
}, symbol_module_split))

if (length(ambiguous_symbols) > 0L) {
  stop("Restored gene symbols map to multiple WGCNA modules; ", "module gene-set enrichment would be ambiguous. First symbols: ",
    paste(head(ambiguous_symbols, 10), collapse = ", "))
}

universe_symbols <- unique(as.character(gene_module$gene_symbol_restored))

universe_symbols <- universe_symbols[!is.na(universe_symbols) & nzchar(universe_symbols)]

signature_in_universe <- intersect(signature_genes, universe_symbols)

hyper_module_set <- function(module, gene_set, set_name) {
    module_genes <- unique(as.character(gene_module$gene_symbol_restored[gene_module$module_color == 
        module]))
    module_genes <- intersect(module_genes, universe_symbols)
    gene_set <- intersect(unique(gene_set), universe_symbols)
    a <- length(intersect(module_genes, gene_set))
    module_n <- length(module_genes)
    set_n <- length(gene_set)
    universe_n <- length(universe_symbols)
    b <- module_n - a
    c <- set_n - a
    d <- universe_n - a - b - c
    if (any(c(a, b, c, d) < 0)) {
        stop("Invalid contingency table for module ", module, " / set ", set_name)
    }
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
    data.frame(module_color = module, gene_set = set_name, universe_n = universe_n, module_n = module_n, 
        gene_set_n_in_universe = set_n, overlap_n = a, overlap_fraction_of_module = ifelse(module_n > 
            0, a/module_n, NA_real_), overlap_fraction_of_gene_set = ifelse(set_n > 0, a/set_n, NA_real_), 
        odds_ratio = unname(ft$estimate), P.Value = ft$p.value, overlap_genes = paste(sort(intersect(module_genes, 
            gene_set)), collapse = "/"), stringsAsFactors = FALSE)
}

sig_overlap_list <- lapply(non_grey_modules, function(module) {
  hyper_module_set(module = module, gene_set = signature_in_universe, set_name = "Step47_Clec7a_UP_signature")
})

signature_overlap <- do.call(rbind, sig_overlap_list)

signature_overlap$adj.P.Val <- p.adjust(signature_overlap$P.Value, method = "BH")

signature_overlap <- signature_overlap[order(signature_overlap$adj.P.Val, -signature_overlap$odds_ratio),
  , drop = FALSE]

required_key48_cols <- c("key_group", "direction", "significant_BH_0.05", "core_enrichment_SYMBOL")

missing_key48 <- setdiff(required_key48_cols, names(key48))

if (length(missing_key48) > 0L) {
  stop("Step48 key pathway summary lacks required column(s): ", paste(missing_key48, collapse = ", "))
}

is_gsea <- rep(FALSE, nrow(key48))

if ("source_table" %in% names(key48)) {
  is_gsea <- is_gsea | grepl("^GSEA", key48$source_table)
}

if ("analysis" %in% names(key48)) {
  is_gsea <- is_gsea | grepl("^GSEA", key48$analysis)
}

theme_specs <- data.frame(key_group = c("lipid_handling_transport", "lysosomal_processing", "myeloid_activation",
  "sterol_biosynthesis"), required_direction = c("positive", "positive", "positive", "negative"), stringsAsFactors = FALSE)

theme_gene_sets <- list()
theme_qc_list <- list()

for (i in seq_len(nrow(theme_specs))) {
    grp <- theme_specs$key_group[i]
    dr <- theme_specs$required_direction[i]
    keep <- is_gsea & key48$key_group == grp & key48$direction == dr & as.logical(key48$significant_BH_0.05)
    terms_n <- sum(keep, na.rm = TRUE)
    genes <- split_symbol_string(key48$core_enrichment_SYMBOL[keep])
    genes_in_universe <- intersect(genes, universe_symbols)
    theme_gene_sets[[grp]] <- genes_in_universe
    theme_qc_list[[length(theme_qc_list) + 1L]] <- data.frame(key_group = grp, required_direction = dr, 
        significant_GSEA_terms_n = terms_n, unique_core_SYMBOL_n = length(genes), unique_core_SYMBOL_in_WGCNA_universe_n = length(genes_in_universe), 
        stringsAsFactors = FALSE)
}

theme_qc <- do.call(rbind, theme_qc_list)

if (any(theme_qc$unique_core_SYMBOL_in_WGCNA_universe_n < 5L)) {
  stop("One or more Step48 mechanism theme gene sets ", "has fewer than 5 WGCNA-universe genes.")
}

mechanism_overlap_list <- list()

for (grp in names(theme_gene_sets)) {

  tmp <- lapply(non_grey_modules, function(module) {
    hyper_module_set(module = module, gene_set = theme_gene_sets[[grp]], set_name = paste0("Step48_",
      grp))
  })

  tmp <- do.call(rbind, tmp)

  tmp$adj.P.Val <- p.adjust(tmp$P.Value, method = "BH")

  mechanism_overlap_list[[grp]] <- tmp
}

mechanism_overlap <- do.call(rbind, mechanism_overlap_list)

rownames(mechanism_overlap) <- NULL

mechanism_overlap <- mechanism_overlap[order(mechanism_overlap$gene_set, mechanism_overlap$adj.P.Val,
  -mechanism_overlap$odds_ratio), , drop = FALSE]

signature_for_merge <- signature_overlap[, c("module_color", "gene_set_n_in_universe", "overlap_n", "odds_ratio",
  "P.Value", "adj.P.Val"), drop = FALSE]

names(signature_for_merge) <- c("module_color", "Step47_signature_n_in_universe", "Step47_signature_overlap_n",
  "Step47_signature_odds_ratio", "Step47_signature_P.Value", "Step47_signature_adj.P.Val")

evidence <- merge(primary_merge, signature_for_merge, by = "module_color", all.x = TRUE, sort = FALSE)

batch_for_merge <- batch_audit[, c("module_color", "estimate_BatchB_minus_BatchA", "P.Value", "adj.P.Val"),
  drop = FALSE]

names(batch_for_merge) <- c("module_color", "estimate_BatchB_minus_BatchA", "batch_audit_P.Value", "batch_audit_adj.P.Val")

evidence <- merge(evidence, batch_for_merge, by = "module_color", all.x = TRUE, sort = FALSE)

core_count <- table(factor(core_rows$module_color, levels = non_grey_modules))

core_count_df <- data.frame(module_color = non_grey_modules, frozen_core_gene_n = as.integer(core_count[non_grey_modules]),
  stringsAsFactors = FALSE)

evidence <- merge(evidence, core_count_df, by = "module_color", all.x = TRUE, sort = FALSE)

for (grp in names(theme_gene_sets)) {
    set_nm <- paste0("Step48_", grp)
    tmp <- mechanism_overlap[mechanism_overlap$gene_set == set_nm, c("module_color", "overlap_n", "odds_ratio", 
        "adj.P.Val"), drop = FALSE]
    names(tmp)[-1] <- paste0(grp, c("_overlap_n", "_odds_ratio", "_adj.P.Val"))
    evidence <- merge(evidence, tmp, by = "module_color", all.x = TRUE, sort = FALSE)
}

evidence$Tier1_eligible <- evidence$estimate_target_minus_reference_primary > 0 & evidence$adj.P.Val_primary <
  FDR_CUTOFF & evidence$estimate_target_minus_reference_sensitivity > 0 & evidence$adj.P.Val_sensitivity <
  FDR_CUTOFF & evidence$Step47_signature_odds_ratio > 1 & evidence$Step47_signature_adj.P.Val < FDR_CUTOFF

evidence$Tier2_eligible <- evidence$estimate_target_minus_reference_primary > 0 & evidence$adj.P.Val_primary <
  FDR_CUTOFF & evidence$estimate_target_minus_reference_sensitivity > 0 & evidence$Step47_signature_odds_ratio >
  1 & evidence$Step47_signature_adj.P.Val < FDR_CUTOFF

selected_module <- NA_character_
selection_tier <- "NONE"
selection_note <- paste("No Tier1/Tier2 module satisfied the pre-specified rule;", "no key module was forced.")

tier1 <- evidence[evidence$Tier1_eligible, , drop = FALSE]

tier2 <- evidence[evidence$Tier2_eligible, , drop = FALSE]

choose_by_rule <- function(df) {
    if (nrow(df) == 0L) {
        return(NA_character_)
    }
    df <- df[order(df$Step47_signature_adj.P.Val, -df$Step47_signature_odds_ratio, -df$t_primary), , 
        drop = FALSE]
    df$module_color[1]
}

if (nrow(tier1) > 0L) {
    selected_module <- choose_by_rule(tier1)
    selection_tier <- "Tier1"
    selection_note <- paste("Selected from Tier1 modules by smallest Step47 signature BH FDR,", "then larger signature odds ratio, then larger primary t statistic.")
} else if (nrow(tier2) > 0L) {
    selected_module <- choose_by_rule(tier2)
    selection_tier <- "Tier2_fallback"
    selection_note <- paste("No Tier1 module; selected from Tier2 fallback by smallest Step47", "signature BH FDR, then larger odds ratio, then larger primary t statistic.")
}

evidence$selected_key_module <- !is.na(selected_module) & evidence$module_color == selected_module

evidence$selection_tier <- ifelse(evidence$selected_key_module, selection_tier, "")

evidence$module_order <- match(evidence$module_color, non_grey_modules)

evidence <- evidence[order(evidence$selected_key_module, decreasing = TRUE), , drop = FALSE]

if (sum(!evidence$selected_key_module) > 1L) {
    sel_part <- evidence[evidence$selected_key_module, , drop = FALSE]
    nonsel_part <- evidence[!evidence$selected_key_module, , drop = FALSE]
    nonsel_part <- nonsel_part[order(nonsel_part$module_order), , drop = FALSE]
    evidence <- rbind(sel_part, nonsel_part)
}

evidence$module_order <- NULL

key_module_gene_evidence <- data.frame(gene_id = character(0), gene_symbol_restored = character(0), module_color = character(0),
  kME_own_module = numeric(0), in_Step47_Clec7a_UP_signature = logical(0), in_Step48_lipid_handling_transport_core = logical(0),
  in_Step48_lysosomal_processing_core = logical(0), in_Step48_myeloid_activation_core = logical(0),
  in_Step48_sterol_biosynthesis_negative_core = logical(0), Step46_primary_logFC = numeric(0), Step46_primary_t = numeric(0),
  Step46_primary_adj.P.Val = numeric(0), stringsAsFactors = FALSE)

if (!is.na(selected_module)) {
    key_module_gene_evidence <- gene_module[gene_module$module_color == selected_module, c("gene_id", 
        "gene_symbol_restored", "module_color", "kME_own_module"), drop = FALSE]
    key_module_gene_evidence$in_Step47_Clec7a_UP_signature <- key_module_gene_evidence$gene_symbol_restored %in% 
        signature_in_universe
    key_module_gene_evidence$in_Step48_lipid_handling_transport_core <- key_module_gene_evidence$gene_symbol_restored %in% 
        theme_gene_sets$lipid_handling_transport
    key_module_gene_evidence$in_Step48_lysosomal_processing_core <- key_module_gene_evidence$gene_symbol_restored %in% 
        theme_gene_sets$lysosomal_processing
    key_module_gene_evidence$in_Step48_myeloid_activation_core <- key_module_gene_evidence$gene_symbol_restored %in% 
        theme_gene_sets$myeloid_activation
    key_module_gene_evidence$in_Step48_sterol_biosynthesis_negative_core <- key_module_gene_evidence$gene_symbol_restored %in% 
        theme_gene_sets$sterol_biosynthesis
    if ("input_primary" %in% names(obj48) && all(c("gene_symbol", "logFC", "t", "adj.P.Val") %in% names(obj48$input_primary))) {
        p48 <- obj48$input_primary[, c("gene_symbol", "logFC", "t", "adj.P.Val"), drop = FALSE]
        p48 <- p48[!duplicated(p48$gene_symbol), , drop = FALSE]
        idx <- match(key_module_gene_evidence$gene_symbol_restored, p48$gene_symbol)
        key_module_gene_evidence$Step46_primary_logFC <- p48$logFC[idx]
        key_module_gene_evidence$Step46_primary_t <- p48$t[idx]
        key_module_gene_evidence$Step46_primary_adj.P.Val <- p48$adj.P.Val[idx]
    }
    else {
        key_module_gene_evidence$Step46_primary_logFC <- NA_real_
        key_module_gene_evidence$Step46_primary_t <- NA_real_
        key_module_gene_evidence$Step46_primary_adj.P.Val <- NA_real_
    }
    key_module_gene_evidence <- key_module_gene_evidence[order(-abs(key_module_gene_evidence$kME_own_module)), 
        , drop = FALSE]
}

mapping_raw <- AnnotationDbi::select(org.Mm.eg.db, keys = universe_symbols, keytype = "SYMBOL", columns = "ENTREZID")

mapping_raw <- mapping_raw[!is.na(mapping_raw$SYMBOL) & nzchar(mapping_raw$SYMBOL) & !is.na(mapping_raw$ENTREZID) &
  nzchar(mapping_raw$ENTREZID), , drop = FALSE]

mapping_raw <- unique(mapping_raw[, c("SYMBOL", "ENTREZID"), drop = FALSE])

background_entrez <- unique(mapping_raw$ENTREZID)

key_entrez <- character(0)

if (!is.na(selected_module)) {
    key_symbols <- unique(key_module_gene_evidence$gene_symbol_restored)
    key_entrez <- unique(mapping_raw$ENTREZID[mapping_raw$SYMBOL %in% key_symbols])
}

enrichment_mapping_qc <- data.frame(metric = c("WGCNA_universe_unique_SYMBOL_n", "WGCNA_universe_mapped_SYMBOL_n",
  "WGCNA_universe_unique_ENTREZ_n", "selected_key_module", "selected_key_module_SYMBOL_n", "selected_key_module_mapped_SYMBOL_n",
  "selected_key_module_unique_ENTREZ_n"), value = c(length(universe_symbols), length(unique(mapping_raw$SYMBOL)),
  length(background_entrez), ifelse(is.na(selected_module), "NONE", selected_module), ifelse(is.na(selected_module),
    0, nrow(key_module_gene_evidence)), ifelse(is.na(selected_module), 0, length(unique(mapping_raw$SYMBOL[mapping_raw$SYMBOL %in%
    key_module_gene_evidence$gene_symbol_restored]))), length(key_entrez)), stringsAsFactors = FALSE)

key_go_obj <- NULL
key_kegg_obj <- NULL
key_go <- empty_enrichment_table()
key_kegg <- empty_enrichment_table()
kegg_error <- NA_character_

if (!is.na(selected_module)) {

  if (length(key_entrez) < MIN_GS_SIZE) {
    stop("Selected key module has too few mapped Entrez genes for enrichment: ", length(key_entrez))
  }

  key_go_obj <- clusterProfiler::enrichGO(gene = key_entrez, universe = background_entrez, OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS_SIZE,
    maxGSSize = MAX_GS_SIZE, readable = FALSE)

  key_go <- as.data.frame(key_go_obj)

  key_kegg_obj <- tryCatch({
    clusterProfiler::enrichKEGG(gene = key_entrez, universe = background_entrez, organism = "mmu",
      keyType = "ncbi-geneid", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE)
  }, error = function(e) {
    kegg_error <<- conditionMessage(e)
    NULL
  })

  if (!is.null(key_kegg_obj)) {
    key_kegg <- as.data.frame(key_kegg_obj)
  }

  entrez_to_symbol <- split(mapping_raw$SYMBOL, mapping_raw$ENTREZID)

  entrez_string_to_symbol <- function(x) {

    if (is.na(x) || !nzchar(x)) {
      return(NA_character_)
    }

    ids <- strsplit(x, "/", fixed = TRUE)[[1]]

    syms <- unique(unlist(entrez_to_symbol[ids], use.names = FALSE))

    syms <- syms[!is.na(syms) & nzchar(syms)]

    paste(syms, collapse = "/")
  }

  if (nrow(key_go) > 0L && "geneID" %in% names(key_go)) {
    key_go$geneID_SYMBOL <- vapply(key_go$geneID, entrez_string_to_symbol, character(1))
  }

  if (nrow(key_kegg) > 0L && "geneID" %in% names(key_kegg)) {
    key_kegg$geneID_SYMBOL <- vapply(key_kegg$geneID, entrez_string_to_symbol, character(1))
  }
}

key_enrichment_summary <- rbind(if (nrow(key_go) > 0L) {
  data.frame(database = "GO_BP", key_go[, intersect(c("ID", "Description", "pvalue", "p.adjust", "Count",
    "geneID_SYMBOL"), names(key_go)), drop = FALSE], stringsAsFactors = FALSE)
} else {
  NULL
}, if (nrow(key_kegg) > 0L) {
  data.frame(database = "KEGG", key_kegg[, intersect(c("ID", "Description", "pvalue", "p.adjust", "Count",
    "geneID_SYMBOL"), names(key_kegg)), drop = FALSE], stringsAsFactors = FALSE)
} else {
  NULL
})

if (!is.null(key_enrichment_summary) && nrow(key_enrichment_summary) > 0L) {
    key_enrichment_summary$significant_BH_0.05 <- !is.na(key_enrichment_summary$p.adjust) & key_enrichment_summary$p.adjust < 
        FDR_CUTOFF
    key_enrichment_summary <- key_enrichment_summary[order(!key_enrichment_summary$significant_BH_0.05, 
        key_enrichment_summary$p.adjust, key_enrichment_summary$pvalue), , drop = FALSE]
} else {
    key_enrichment_summary <- data.frame(database = character(0), ID = character(0), Description = character(0), 
        pvalue = numeric(0), p.adjust = numeric(0), Count = integer(0), geneID_SYMBOL = character(0), 
        significant_BH_0.05 = logical(0), stringsAsFactors = FALSE)
}

input_qc <- data.frame(metric = c("Step45_RDS", "Step47_RDS", "Step48_RDS", "Step49_RDS", "samples_n",
  "WGCNA_genes_n", "non_grey_modules_n", "Step49_soft_power", "Step47_signature_n", "Step47_signature_in_WGCNA_universe_n",
  "restored_SYMBOL_cross_module_ambiguity_n", "Step45_review_flagged_samples_n", "primary_samples_n",
  "primary_sensitivity_samples_n", "primary_sensitivity_estimate_Pearson", "primary_sensitivity_estimate_Spearman",
  "primary_sensitivity_same_direction_fraction", "key_module_selection_tier", "selected_key_module",
  "key_module_hard_coded", "Step48_mechanism_used_for_key_module_selection"), value = c(STEP45_RDS,
  STEP47_RDS, STEP48_RDS, STEP49_RDS, nrow(meta), nrow(gene_module), length(non_grey_modules), soft_power,
  length(signature_genes), length(signature_in_universe), length(ambiguous_symbols), sum(review_flag),
  sum(keep_primary), sum(keep_primary_sens), primary_sensitivity_pearson, primary_sensitivity_spearman,
  primary_sensitivity_same_direction, selection_tier, ifelse(is.na(selected_module), "NONE", selected_module),
  "NO", "NO"), stringsAsFactors = FALSE)

files <- list(input_qc = paste0("01_", PREFIX, "_input_and_integration_QC.csv"), primary_trait = paste0("02_",
  PREFIX, "_PRIMARY_module_trait_and_sensitivity.csv"), secondary_trait = paste0("03_", PREFIX, "_secondary_module_trait.csv"),
  global_batch = paste0("04_", PREFIX, "_global_time_and_batch_audit.csv"), time_summary = paste0("05_",
    PREFIX, "_module_eigengene_time_summary.csv"), signature_overlap = paste0("06_", PREFIX, "_Step47_Clec7a_signature_module_enrichment.csv"),
  theme_qc = paste0("07_", PREFIX, "_Step48_mechanism_theme_gene_set_QC.csv"), mechanism_overlap = paste0("08_",
    PREFIX, "_Step48_mechanism_module_enrichment.csv"), selection = paste0("09_", PREFIX, "_module_evidence_and_key_selection.csv"),
  key_genes = paste0("10_", PREFIX, "_key_module_gene_evidence.csv"), key_go = paste0("11_", PREFIX,
    "_key_module_ORA_GO_BP.csv"), key_kegg = paste0("12_", PREFIX, "_key_module_ORA_KEGG.csv"), enrichment_qc = paste0("13_",
    PREFIX, "_key_module_enrichment_mapping_QC.csv"), enrichment_summary = paste0("14_", PREFIX,
    "_key_module_enrichment_summary.csv"), analysis_object = paste0("30_", PREFIX, "_analysis_object.rds"),
  report = paste0("31_", PREFIX, "_REPORT.txt"), session = paste0("32_", PREFIX, "_sessionInfo.txt"))

write_csv(input_qc, files$input_qc)

write_csv(primary_merge, files$primary_trait)

write_csv(secondary_module_trait, files$secondary_trait)

all_cols <- unique(unlist(lapply(list(data.frame(analysis_type = "global_time", global_time, stringsAsFactors = FALSE),
  data.frame(analysis_type = "batch_audit", batch_audit, stringsAsFactors = FALSE)), names)))

pad_df <- function(x, cols) {
  miss <- setdiff(cols, names(x))
  for (nm in miss) {
    x[[nm]] <- NA
  }
  x[, cols, drop = FALSE]
}

global_batch_combined <- rbind(pad_df(data.frame(analysis_type = "global_time", global_time, stringsAsFactors = FALSE),
  all_cols), pad_df(data.frame(analysis_type = "batch_audit", batch_audit, stringsAsFactors = FALSE),
  all_cols))

write_csv(global_batch_combined, files$global_batch)

write_csv(module_time_summary, files$time_summary)

write_csv(signature_overlap, files$signature_overlap)

write_csv(theme_qc, files$theme_qc)

write_csv(mechanism_overlap, files$mechanism_overlap)

write_csv(evidence, files$selection)

write_csv(key_module_gene_evidence, files$key_genes)

write_csv(key_go, files$key_go)

write_csv(key_kegg, files$key_kegg)

write_csv(enrichment_mapping_qc, files$enrichment_qc)

write_csv(key_enrichment_summary, files$enrichment_summary)

analysis_object <- list(step = 50, version = "v01", analysis = "Bulk WGCNA module-trait and cross-omics integration",
  frozen_upstream = list(Step45 = STEP45_RDS, Step47 = STEP47_RDS, Step48 = STEP48_RDS, Step49 = STEP49_RDS),
  frozen_network = list(soft_power = soft_power, non_grey_modules = non_grey_modules, gene_module_assignment = gene_module,
    module_eigengenes = MEs), selection_rule = list(Tier1 = "Primary late>early FDR<0.05 positive; sensitivity FDR<0.05 positive; Step47 signature enrichment FDR<0.05 OR>1",
    Tier2 = "Primary late>early FDR<0.05 positive; sensitivity same positive direction; Step47 signature enrichment FDR<0.05 OR>1",
    tie_break = "smallest Step47 signature FDR, then larger odds ratio, then larger primary t", Step48_mechanism_used_for_selection = FALSE),
  module_trait = list(primary = primary_full, primary_sensitivity = primary_sens, primary_merged = primary_merge,
    transition_vs_early = transition_vs_early, extended_vs_early = extended_vs_early, d28_vs_d14 = d28_vs_d14,
    global_time = global_time, batch_audit = batch_audit, time_summary = module_time_summary), Step47_signature = list(signature_genes = signature_genes,
    signature_in_WGCNA_universe = signature_in_universe, module_enrichment = signature_overlap),
  Step48_mechanism = list(theme_specs = theme_specs, theme_gene_sets = theme_gene_sets, theme_qc = theme_qc,
    module_enrichment = mechanism_overlap), key_module = list(selected_module = selected_module,
    selection_tier = selection_tier, selection_note = selection_note, all_module_evidence = evidence,
    gene_evidence = key_module_gene_evidence, enrichment_mapping_qc = enrichment_mapping_qc, ORA_GO_BP = key_go_obj,
    ORA_GO_BP_table = key_go, ORA_KEGG = key_kegg_obj, ORA_KEGG_table = key_kegg, enrichment_summary = key_enrichment_summary,
    KEGG_error = kegg_error), notes = c("No Step49 network reconstruction was performed.", "Key module was not hard-coded by color.",
    "Step48 mechanism themes were not used to select the key module.", "The primary module trait is Established late vs Early with batch adjustment.",
    "A sensitivity analysis excludes the three Step45 review-flagged samples.", "Step47 signature enrichment is a module-level overrepresentation test.",
    "Step48 mechanism linkage uses significant GSEA core genes in the directions.", "Key-module ORA background is the WGCNA-tested universe, not the whole genome.",
    "No manuscript figure is generated in Step50."))

saveRDS(analysis_object, file.path(OUT_DIR, files$analysis_object), compress = "xz")

selected_row <- if (!is.na(selected_module)) {
  evidence[evidence$module_color == selected_module, , drop = FALSE]
} else {
  evidence[FALSE, , drop = FALSE]
}

go_sig_n <- if (nrow(key_go) > 0L && "p.adjust" %in% names(key_go)) {
  sum(!is.na(key_go$p.adjust) & key_go$p.adjust < FDR_CUTOFF)
} else {
  0L
}

kegg_sig_n <- if (nrow(key_kegg) > 0L && "p.adjust" %in% names(key_kegg)) {
  sum(!is.na(key_kegg$p.adjust) & key_kegg$p.adjust < FDR_CUTOFF)
} else {
  0L
}

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT_DIR, files$report), useBytes = TRUE)

sink(file.path(OUT_DIR, files$session))
print(sessionInfo())
sink()

output_index <- data.frame(file = c(paste0("00_", PREFIX, "_OUTPUT_INDEX.csv"), unlist(files, use.names = FALSE)),
  description = c("Index of Step50 outputs", "-input and integration QC", "Primary Established late vs Early batch-adjusted module associations plus review-sample sensitivity",
    "Transition, persistence, and 28-vs-14 secondary module associations", "Overall timepoint F-test and batch-after-time audit",
    "Descriptive module eigengene mean/SD/SEM across all 8 time points", "Step47 Clec7a-UP signature overrepresentation across all non-grey modules",
    "Construction QC for Step48 mechanism-theme GSEA core gene sets", "Step48 mechanism-theme overrepresentation across all non-grey modules",
    "All-module evidence table and objective key-module selection", "Selected key-module gene evidence with kME, Step47, Step48, and Step46 annotations",
    "Complete selected key-module GO BP ORA", "Complete selected key-module KEGG ORA", "Selected key-module SYMBOL/ENTREZ enrichment mapping QC",
    "Selected key-module significant enrichment summary", "Complete Step50 analysis object", "Plain-text Step50 report",
    "R sessionInfo"), stringsAsFactors = FALSE)

write_csv(output_index, paste0("00_", PREFIX, "_OUTPUT_INDEX.csv"))

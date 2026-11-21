rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)
set.seed(20261118)

ROOT <- PROJECT_ROOT
STEP45_RDS <- file.path( ROOT, "results", "45_Bulk_batch_structure_and_analysis_object", "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds" )

STEP46_PRIMARY_CSV <- file.path( ROOT, "results", "46_Bulk_limma_temporal_DEG_analysis", "04_46_Bulk_limma_temporal_DEG_analysis_v06_PRIMARY_EstablishedLate_vs_Early.csv" )

SIGNATURE_FILE <- file.path( ROOT, "results", "21_Clec7a_DEG_volcano_myeloid", "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv" )

OUT_DIR <- file.path( ROOT, "results", "148_Reviewer1_RemainingMajor4_proxy_sensitivity" )
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "148_Reviewer1_RemainingMajor4_proxy_sensitivity_v01"
FDR_CUTOFF <- 0.05
EXPECTED_GENES <- 12977L
EXPECTED_ALL_SAMPLES <- 39L
EXPECTED_PRIMARY_SAMPLES <- 29L
EXPECTED_SIGNATURE_N <- 1320L

CORE_GENES <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")
TIMEPOINT_ORDER <- c( "Uninjured", "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi" )

required_pkgs <- c("limma", "fgsea")
missing_pkgs <- required_pkgs[ !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE) ]
if (length(missing_pkgs) > 0L) {
  stop( "Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall before running Step148.",
    if ("fgsea" %in% missing_pkgs) "\nFor fgsea: BiocManager::install('fgsea')" else "" )
}
suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
})

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

num_name <- function(num, suffix) {
  paste0(sprintf("%02d", as.integer(num)), "_", PREFIX, "_", suffix)
}

out_path <- function(num, suffix) file.path(OUT_DIR, num_name(num, suffix))

write_csv <- function(x, num, suffix) {
  utils::write.csv( x, out_path(num, suffix), row.names = FALSE, na = "", fileEncoding = "UTF-8" )
}

safe_gene_symbol <- function(x) {
  ifelse(grepl("^X[0-9]", x), sub("^X", "", x), x)
}

design_rank_audit <- function(design, label) {
  q <- qr(design)
  data.frame( model = label, n_samples = nrow(design), n_columns = ncol(design), rank = q$rank, full_rank = q$rank == ncol(design), design_columns = paste(colnames(design), collapse = ";"),
    stringsAsFactors = FALSE )
}

fit_primary_model <- function(expr_mat, meta_df, include_proxy = FALSE) {
  if (include_proxy) {
    design <- model.matrix(~ 0 + primary_group + operational_group_proxy, data = meta_df)
    if (ncol(design) != 3L) stop("WITH-proxy design expected 3 columns.")
    colnames(design) <- c("Early", "EstablishedLate", "Proxy_B")
    model_label <- "WITH proxy"
  } else {
    design <- model.matrix(~ 0 + primary_group, data = meta_df)
    if (ncol(design) != 2L) stop("WITHOUT-proxy design expected 2 columns.")
    colnames(design) <- c("Early", "EstablishedLate")
    model_label <- "WITHOUT proxy"
  }

  rank_audit <- design_rank_audit(design, model_label)
  if (!rank_audit$full_rank) stop(model_label, " design is not full rank.")

  aw <- limma::arrayWeights(expr_mat, design = design)
  fit0 <- limma::lmFit(expr_mat, design = design, weights = aw)
  contrast_matrix <- limma::makeContrasts( EstablishedLate_vs_Early = EstablishedLate - Early, levels = design )
  fit1 <- limma::contrasts.fit(fit0, contrast_matrix)
  fit2 <- limma::eBayes(fit1, trend = TRUE, robust = TRUE)
  tt <- limma::topTable( fit2, coef = "EstablishedLate_vs_Early", number = Inf, sort.by = "P", adjust.method = "BH" )
  if (is.null(rownames(tt))) stop(model_label, " result table has no gene row names.")

  out <- data.frame( gene_id = rownames(tt), gene_symbol_restored = safe_gene_symbol(rownames(tt)), transformed_scale_difference = as.numeric(tt$logFC), AveExpr = as.numeric(tt$AveExpr),
    t = as.numeric(tt$t), P.Value = as.numeric(tt$P.Value), adj.P.Val = as.numeric(tt$adj.P.Val), B = as.numeric(tt$B), stringsAsFactors = FALSE )

  out$status <- "Not significant"
  out$status[ !is.na(out$adj.P.Val) & out$adj.P.Val < FDR_CUTOFF & out$transformed_scale_difference > 0 ] <- "Higher"
  out$status[ !is.na(out$adj.P.Val) & out$adj.P.Val < FDR_CUTOFF & out$transformed_scale_difference < 0 ] <- "Lower"
  out <- out[order(out$gene_id), , drop = FALSE]
  rownames(out) <- NULL

  weight_df <- data.frame( sample_id = colnames(expr_mat), model = model_label, operational_group_proxy = as.character(meta_df$operational_group_proxy),
    primary_group = as.character(meta_df$primary_group), array_weight = as.numeric(aw), stringsAsFactors = FALSE )

  list( label = model_label, design = design, rank_audit = rank_audit, array_weights = weight_df, table = out, fit = fit2 )
}

count_model <- function(tab, model_label) {
  n_higher <- sum(tab$status == "Higher", na.rm = TRUE)
  n_lower <- sum(tab$status == "Lower", na.rm = TRUE)
  data.frame( model = model_label, tested_genes_n = nrow(tab), FDR_cutoff = FDR_CUTOFF, absolute_effect_threshold = "None", higher_FDR_lt_0.05_n = n_higher, lower_FDR_lt_0.05_n = n_lower,
    total_FDR_lt_0.05_n = n_higher + n_lower, not_significant_n = sum(tab$status == "Not significant", na.rm = TRUE), stringsAsFactors = FALSE )
}

jaccard_index <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

make_cross_tab_df <- function(row_var, col_var, row_name, row_levels = NULL, col_levels = NULL) {
  if (!is.null(row_levels)) row_var <- factor(as.character(row_var), levels = row_levels)
  if (!is.null(col_levels)) col_var <- factor(as.character(col_var), levels = col_levels)
  z <- as.data.frame.matrix(table(row_var, col_var, useNA = "ifany"))
  z[[row_name]] <- rownames(z)
  rownames(z) <- NULL
  z[, c(row_name, setdiff(names(z), row_name)), drop = FALSE]
}

collapse_rank_table <- function(tt) {
  required_cols <- c( "gene_symbol_restored", "transformed_scale_difference", "t", "P.Value", "adj.P.Val" )
  if (!all(required_cols %in% colnames(tt))) {
    stop("Rank table missing: ", paste(setdiff(required_cols, colnames(tt)), collapse = ", "))
  }
  x <- data.frame( gene = safe_gene_symbol(as.character(tt$gene_symbol_restored)), transformed_scale_difference = as.numeric(tt$transformed_scale_difference), t = as.numeric(tt$t),
    P.Value = as.numeric(tt$P.Value), adj.P.Val = as.numeric(tt$adj.P.Val), stringsAsFactors = FALSE )
  x <- x[!is.na(x$gene) & nzchar(x$gene) & is.finite(x$t), , drop = FALSE]
  x <- x[order(abs(x$t), decreasing = TRUE), , drop = FALSE]
  x <- x[!duplicated(x$gene), , drop = FALSE]
  x <- x[order(x$t, decreasing = TRUE), , drop = FALSE]
  rownames(x) <- NULL
  x
}

run_signature_gsea <- function(tt, signature_genes, model_label) {
  rank_table <- collapse_rank_table(tt)
  ranks <- rank_table$t
  names(ranks) <- rank_table$gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  overlap <- intersect(signature_genes, names(ranks))
  if (length(overlap) < 15L) {
    stop("Too few frozen-signature genes overlap bulk ranks for ", model_label, ": ", length(overlap))
  }
  pathways <- list(Clec7a_UP_scRNA = overlap)
  set.seed(20261118)
  fg <- fgsea::fgseaMultilevel( pathways = pathways, stats = ranks, minSize = 15, maxSize = 5000, eps = 0, scoreType = "std" )
  if (nrow(fg) != 1L) stop("Unexpected fgsea row count for ", model_label)
  leading <- fg$leadingEdge[[1]]
  summary <- data.frame( model = model_label, pathway = as.character(fg$pathway), ES = as.numeric(fg$ES), NES = as.numeric(fg$NES), pval = as.numeric(fg$pval),
    fgsea_padj_single_pathway = as.numeric(fg$padj), signature_n_total = length(signature_genes), signature_n_in_bulk = length(overlap), leading_edge_n = length(leading),
    stringsAsFactors = FALSE )
  leading_df <- rank_table[match(leading, rank_table$gene), , drop = FALSE]
  leading_df$model <- model_label
  leading_df$leading_edge_order <- seq_len(nrow(leading_df))
  leading_df <- leading_df[, c( "model", "leading_edge_order", "gene", "transformed_scale_difference", "t", "P.Value", "adj.P.Val" ), drop = FALSE]
  list( summary = summary, leading_edge = leading_df, leading_genes = leading, rank_table = rank_table )
}

for (f in c(STEP45_RDS, STEP46_PRIMARY_CSV, SIGNATURE_FILE)) {
  if (!file.exists(f)) stop("Required input not found:\n", f)
}

msg("Reading frozen Step45 analysis object...")
obj45 <- readRDS(STEP45_RDS)
if (!all(c("expression_original", "metadata") %in% names(obj45))) {
  stop("Step45 object lacks expression_original and/or metadata.")
}

expr <- as.matrix(obj45$expression_original)
storage.mode(expr) <- "numeric"
meta <- obj45$metadata

if (nrow(expr) != EXPECTED_GENES || ncol(expr) != EXPECTED_ALL_SAMPLES) {
  stop( "Unexpected Step45 expression dimensions: ", nrow(expr), " x ", ncol(expr), "; expected ", EXPECTED_GENES, " x ", EXPECTED_ALL_SAMPLES, "." )
}
if (!all(colnames(expr) %in% rownames(meta))) {
  stop("Expression sample IDs do not match Step45 metadata row names.")
}
meta <- meta[colnames(expr), , drop = FALSE]
required_meta_cols <- c("batch", "timepoint", "phase")
if (!all(required_meta_cols %in% colnames(meta))) {
  stop("Step45 metadata lacks: ", paste(setdiff(required_meta_cols, colnames(meta)), collapse = ", "))
}
meta$sample_id <- rownames(meta)

batch_values <- unique(as.character(meta$batch))
if (!all(c("Batch_A", "Batch_B") %in% batch_values)) {
  stop("Expected Step45 batch proxy levels Batch_A and Batch_B; found: ", paste(sort(batch_values), collapse = ", "))
}
meta$operational_group_proxy <- factor( as.character(meta$batch), levels = c("Batch_A", "Batch_B") )
meta$timepoint <- factor(as.character(meta$timepoint), levels = TIMEPOINT_ORDER)

ids_primary <- rownames(meta)[ as.character(meta$phase) %in% c("Early", "Established late") ]
if (length(ids_primary) != EXPECTED_PRIMARY_SAMPLES) {
  stop("Primary sample count mismatch: found ", length(ids_primary), "; expected ", EXPECTED_PRIMARY_SAMPLES)
}
meta_primary <- droplevels(meta[ids_primary, , drop = FALSE])
expr_primary <- expr[, ids_primary, drop = FALSE]
meta_primary$primary_group <- factor( ifelse(as.character(meta_primary$phase) == "Early", "Early", "EstablishedLate"), levels = c("Early", "EstablishedLate") )
meta_primary$operational_group_proxy <- factor( as.character(meta_primary$operational_group_proxy), levels = c("Batch_A", "Batch_B") )

sample_audit <- data.frame( sample_id = rownames(meta), timepoint = as.character(meta$timepoint), phase = as.character(meta$phase),
  operational_group_proxy = as.character(meta$operational_group_proxy), in_primary_29_sample_comparison = rownames(meta) %in% ids_primary, primary_group = ifelse(
  rownames(meta) %in% ids_primary, ifelse(as.character(meta$phase) == "Early", "Early", "EstablishedLate"), NA_character_ ), stringsAsFactors = FALSE )
write_csv(sample_audit, 1, "sample_set_audit.csv")

proxy_by_timepoint <- make_cross_tab_df( meta$timepoint, meta$operational_group_proxy, "timepoint", row_levels = TIMEPOINT_ORDER, col_levels = c("Batch_A", "Batch_B") )
names(proxy_by_timepoint)[names(proxy_by_timepoint) == "Batch_A"] <- "Group_A_n"
names(proxy_by_timepoint)[names(proxy_by_timepoint) == "Batch_B"] <- "Group_B_n"
proxy_by_timepoint$both_proxy_groups_represented <-
  proxy_by_timepoint$Group_A_n > 0 & proxy_by_timepoint$Group_B_n > 0
proxy_by_timepoint$identifiability_note <- ifelse( proxy_by_timepoint$both_proxy_groups_represented, "Both proxy groups represented at this time point",
  "Only one proxy group represented; stage and proxy effects cannot be independently identified within this time point" )
write_csv(proxy_by_timepoint, 2, "proxy_by_timepoint_cross_tab.csv")

proxy_by_primary_phase <- make_cross_tab_df( meta_primary$primary_group, meta_primary$operational_group_proxy, "primary_group",
  row_levels = c("Early", "EstablishedLate"), col_levels = c("Batch_A", "Batch_B") )
names(proxy_by_primary_phase)[names(proxy_by_primary_phase) == "Batch_A"] <- "Group_A_n"
names(proxy_by_primary_phase)[names(proxy_by_primary_phase) == "Batch_B"] <- "Group_B_n"
proxy_by_primary_phase$total_n <- proxy_by_primary_phase$Group_A_n + proxy_by_primary_phase$Group_B_n
write_csv(proxy_by_primary_phase, 3, "proxy_by_primary_phase_cross_tab.csv")

msg("Fitting primary model WITHOUT proxy...")
fit_without <- fit_primary_model(expr_primary, meta_primary, include_proxy = FALSE)
msg("Fitting primary model WITH proxy...")
fit_with <- fit_primary_model(expr_primary, meta_primary, include_proxy = TRUE)

design_audit <- rbind(fit_without$rank_audit, fit_with$rank_audit)
both_timepoints_n <- sum(proxy_by_timepoint$both_proxy_groups_represented)
single_group_timepoints <- proxy_by_timepoint$timepoint[ !proxy_by_timepoint$both_proxy_groups_represented ]

identifiability_audit <- data.frame( item = c( "WITH-proxy primary design full rank", "WITHOUT-proxy primary design full rank", "Time points with both proxy groups represented",
  "Total time points", "Time points represented by only one proxy group", "Uninjured proxy distribution", "15 min proxy distribution", "1 dpi proxy distribution",
  "Primary Early proxy distribution", "Primary EstablishedLate proxy distribution", "Interpretation" ), value = c( as.character(fit_with$rank_audit$full_rank),
  as.character(fit_without$rank_audit$full_rank), as.character(both_timepoints_n), as.character(nrow(proxy_by_timepoint)), paste(single_group_timepoints, collapse = "; "),
  paste0("A=", proxy_by_timepoint$Group_A_n[proxy_by_timepoint$timepoint == "Uninjured"], "; B=", proxy_by_timepoint$Group_B_n[proxy_by_timepoint$timepoint == "Uninjured"]),
  paste0("A=", proxy_by_timepoint$Group_A_n[proxy_by_timepoint$timepoint == "15 min"], "; B=", proxy_by_timepoint$Group_B_n[proxy_by_timepoint$timepoint == "15 min"]),
  paste0("A=", proxy_by_timepoint$Group_A_n[proxy_by_timepoint$timepoint == "1 dpi"], "; B=", proxy_by_timepoint$Group_B_n[proxy_by_timepoint$timepoint == "1 dpi"]),
  paste0("A=", proxy_by_primary_phase$Group_A_n[proxy_by_primary_phase$primary_group == "Early"], "; B=", proxy_by_primary_phase$Group_B_n[proxy_by_primary_phase$primary_group == "Early"]),
  paste0("A=", proxy_by_primary_phase$Group_A_n[proxy_by_primary_phase$primary_group == "EstablishedLate"],
  "; B=", proxy_by_primary_phase$Group_B_n[proxy_by_primary_phase$primary_group == "EstablishedLate"]), paste( "Full-rank estimation does not establish biological separability.",
  "Uninjured occurs only in Group A, whereas 15 min and 1 dpi occur only in Group B;", "therefore those individual early-stage contrasts are confounded with the unverified proxy.",
  "The primary Early vs EstablishedLate comparison is more balanced but remains a sensitivity analysis to an unverified operational grouping." ) ), stringsAsFactors = FALSE )

design_and_identifiability <- rbind( data.frame( item = paste0("Design rank: ", design_audit$model), value = paste0( "rank ", design_audit$rank, "/", design_audit$n_columns,
  "; full_rank=", design_audit$full_rank, "; columns=", design_audit$design_columns ), stringsAsFactors = FALSE ), identifiability_audit )
write_csv(design_and_identifiability, 4, "design_rank_and_identifiability_audit.csv")

aw_cmp <- merge( fit_without$array_weights[, c("sample_id", "operational_group_proxy", "primary_group", "array_weight")], fit_with$array_weights[, c("sample_id", "array_weight")],
  by = "sample_id", suffixes = c("_without_proxy", "_with_proxy"), all = TRUE, sort = FALSE )
aw_cmp$array_weight_difference_with_minus_without <-
  aw_cmp$array_weight_with_proxy - aw_cmp$array_weight_without_proxy
write_csv(aw_cmp, 5, "arrayWeights_with_vs_without_proxy.csv")

write_csv(fit_with$table, 6, "PRIMARY_with_proxy_all_genes.csv")
write_csv(fit_without$table, 7, "PRIMARY_without_proxy_all_genes.csv")

with_tab <- fit_with$table
without_tab <- fit_without$table

rename_cols <- c( "transformed_scale_difference", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "status" )
idx <- names(with_tab) %in% rename_cols
names(with_tab)[idx] <- paste0(names(with_tab)[idx], "_with_proxy")
idx <- names(without_tab) %in% rename_cols
names(without_tab)[idx] <- paste0(names(without_tab)[idx], "_without_proxy")

paired <- merge( with_tab, without_tab, by = c("gene_id", "gene_symbol_restored"), all = TRUE, sort = FALSE )
if (nrow(paired) != EXPECTED_GENES) {
  stop("Paired gene table has ", nrow(paired), " rows; expected ", EXPECTED_GENES)
}
paired$coefficient_difference_with_minus_without <-
  paired$transformed_scale_difference_with_proxy - paired$transformed_scale_difference_without_proxy
paired$same_direction <-
  sign(paired$transformed_scale_difference_with_proxy) == sign(paired$transformed_scale_difference_without_proxy)
paired$FDR_significant_with_proxy <-
  !is.na(paired$adj.P.Val_with_proxy) & paired$adj.P.Val_with_proxy < FDR_CUTOFF
paired$FDR_significant_without_proxy <-
  !is.na(paired$adj.P.Val_without_proxy) & paired$adj.P.Val_without_proxy < FDR_CUTOFF
paired$FDR_significant_in_both <-
  paired$FDR_significant_with_proxy & paired$FDR_significant_without_proxy
paired$FDR_status_concordant <- paired$status_with_proxy == paired$status_without_proxy
paired <- paired[order(paired$gene_id), , drop = FALSE]
write_csv(paired, 9, "gene_level_with_vs_without_proxy.csv")

model_counts <- rbind( count_model(fit_with$table, "WITH proxy"), count_model(fit_without$table, "WITHOUT proxy") )
write_csv(model_counts, 8, "model_FDR_count_summary.csv")

sig_with <- paired$gene_id[paired$FDR_significant_with_proxy]
sig_without <- paired$gene_id[paired$FDR_significant_without_proxy]
higher_with <- paired$gene_id[paired$status_with_proxy == "Higher"]
higher_without <- paired$gene_id[paired$status_without_proxy == "Higher"]
lower_with <- paired$gene_id[paired$status_with_proxy == "Lower"]
lower_without <- paired$gene_id[paired$status_without_proxy == "Lower"]

coef_pearson <- cor( paired$transformed_scale_difference_with_proxy, paired$transformed_scale_difference_without_proxy, method = "pearson", use = "complete.obs" )
coef_spearman <- cor( paired$transformed_scale_difference_with_proxy, paired$transformed_scale_difference_without_proxy, method = "spearman", use = "complete.obs" )
t_pearson <- cor( paired$t_with_proxy, paired$t_without_proxy, method = "pearson", use = "complete.obs" )
same_direction_pct <- 100 * mean(paired$same_direction, na.rm = TRUE)
FDR_status_concordance_pct <- 100 * mean(paired$FDR_status_concordant, na.rm = TRUE)

comparison_summary <- data.frame( metric = c( "Pearson_r_transformed_scale_difference", "Spearman_rho_transformed_scale_difference", "Pearson_r_moderated_t", "Same_direction_genes_percent",
  "FDR_status_concordance_percent", "FDR_significant_WITH_proxy_n", "FDR_significant_WITHOUT_proxy_n", "FDR_significant_overlap_n", "FDR_significant_union_n",
  "FDR_significant_gene_set_Jaccard", "Higher_WITH_proxy_n", "Higher_WITHOUT_proxy_n", "Higher_overlap_n", "Higher_gene_set_Jaccard", "Lower_WITH_proxy_n", "Lower_WITHOUT_proxy_n",
  "Lower_overlap_n", "Lower_gene_set_Jaccard" ), value = c( coef_pearson, coef_spearman, t_pearson, same_direction_pct, FDR_status_concordance_pct, length(sig_with), length(sig_without),
  length(intersect(sig_with, sig_without)), length(union(sig_with, sig_without)), jaccard_index(sig_with, sig_without), length(higher_with), length(higher_without),
  length(intersect(higher_with, higher_without)), jaccard_index(higher_with, higher_without), length(lower_with), length(lower_without), length(intersect(lower_with, lower_without)),
  jaccard_index(lower_with, lower_without) ), stringsAsFactors = FALSE )
write_csv(comparison_summary, 10, "with_vs_without_proxy_concordance_summary.csv")

fdr_overlap_summary <- data.frame( set = c("All FDR-significant genes", "Higher genes", "Lower genes"), with_proxy_n = c(length(sig_with), length(higher_with), length(lower_with)),
  without_proxy_n = c(length(sig_without), length(higher_without), length(lower_without)), overlap_n = c( length(intersect(sig_with, sig_without)),
  length(intersect(higher_with, higher_without)), length(intersect(lower_with, lower_without)) ), union_n = c( length(union(sig_with, sig_without)),
  length(union(higher_with, higher_without)), length(union(lower_with, lower_without)) ), Jaccard = c( jaccard_index(sig_with, sig_without), jaccard_index(higher_with, higher_without),
  jaccard_index(lower_with, lower_without) ), stringsAsFactors = FALSE )
write_csv(fdr_overlap_summary, 11, "FDR_significant_gene_set_overlap.csv")

core_rows <- paired[paired$gene_symbol_restored %in% CORE_GENES, , drop = FALSE]
core_rows$core_order <- match(core_rows$gene_symbol_restored, CORE_GENES)
core_rows <- core_rows[order(core_rows$core_order), , drop = FALSE]
core_export <- core_rows[, c( "gene_symbol_restored", "gene_id", "transformed_scale_difference_with_proxy", "t_with_proxy", "P.Value_with_proxy", "adj.P.Val_with_proxy", "status_with_proxy",
  "transformed_scale_difference_without_proxy", "t_without_proxy", "P.Value_without_proxy", "adj.P.Val_without_proxy", "status_without_proxy",
  "coefficient_difference_with_minus_without", "same_direction" ), drop = FALSE]
write_csv(core_export, 12, "six_core_genes_with_vs_without_proxy.csv")

msg("Reading frozen Step21 Clec7a-UP signature...")
sig_raw <- read.csv(SIGNATURE_FILE, check.names = FALSE, stringsAsFactors = FALSE)
required_sig_cols <- c("gene", "avg_log2FC", "p_val_adj")
if (!all(required_sig_cols %in% colnames(sig_raw))) {
  stop("Signature file lacks: ", paste(setdiff(required_sig_cols, colnames(sig_raw)), collapse = ", "))
}
sig <- sig_raw[!is.na(sig_raw$gene) & nzchar(sig_raw$gene), , drop = FALSE]
sig$gene <- safe_gene_symbol(as.character(sig$gene))
sig <- sig[!duplicated(sig$gene), , drop = FALSE]
if (nrow(sig) != EXPECTED_SIGNATURE_N) {
  stop("Frozen signature size mismatch: found ", nrow(sig), "; expected ", EXPECTED_SIGNATURE_N)
}
if (any(sig$p_val_adj >= 0.05, na.rm = TRUE) || any(sig$avg_log2FC < 0.25, na.rm = TRUE)) {
  stop("Frozen Step21 signature does not satisfy its original scRNA definition.")
}
signature_genes <- unique(sig$gene)

msg("Running frozen-signature GSEA WITH proxy...")
gsea_with <- run_signature_gsea(fit_with$table, signature_genes, "WITH proxy")
msg("Running frozen-signature GSEA WITHOUT proxy...")
gsea_without <- run_signature_gsea(fit_without$table, signature_genes, "WITHOUT proxy")

gsea_summary <- rbind(gsea_with$summary, gsea_without$summary)
gsea_summary$BH_FDR_across_two_proxy_models <- p.adjust(gsea_summary$pval, method = "BH")
gsea_summary$enrichment_direction <- ifelse( gsea_summary$NES > 0, "Positive: signature enriched toward genes higher in EstablishedLate",
  "Negative: signature enriched toward genes lower in EstablishedLate" )
write_csv(gsea_summary, 13, "Clec7a_signature_GSEA_with_vs_without_proxy.csv")

leading_edges <- rbind(gsea_with$leading_edge, gsea_without$leading_edge)
write_csv(leading_edges, 14, "Clec7a_signature_GSEA_leading_edges.csv")

le_with <- gsea_with$leading_genes
le_without <- gsea_without$leading_genes
gsea_concordance <- data.frame( metric = c( "NES_WITH_proxy", "NES_WITHOUT_proxy", "Absolute_NES_difference", "Leading_edge_WITH_proxy_n", "Leading_edge_WITHOUT_proxy_n",
  "Leading_edge_overlap_n", "Leading_edge_union_n", "Leading_edge_Jaccard" ), value = c( gsea_with$summary$NES, gsea_without$summary$NES,
  abs(gsea_with$summary$NES - gsea_without$summary$NES), length(le_with), length(le_without), length(intersect(le_with, le_without)), length(union(le_with, le_without)),
  jaccard_index(le_with, le_without) ), stringsAsFactors = FALSE )
write_csv(gsea_concordance, 15, "Clec7a_signature_GSEA_concordance.csv")

msg("Auditing WITH-proxy model against frozen Step46 primary CSV...")
step46_primary <- read.csv(STEP46_PRIMARY_CSV, check.names = FALSE, stringsAsFactors = FALSE)
required_step46_cols <- c("gene_id", "gene_symbol_restored", "logFC", "t", "P.Value", "adj.P.Val")
if (!all(required_step46_cols %in% colnames(step46_primary))) {
  stop("Frozen Step46 primary CSV lacks: ", paste(setdiff(required_step46_cols, colnames(step46_primary)), collapse = ", "))
}
step46_cmp <- merge( step46_primary[, required_step46_cols, drop = FALSE], fit_with$table[, c("gene_id", "transformed_scale_difference", "t", "P.Value", "adj.P.Val"), drop = FALSE],
  by = "gene_id", suffixes = c("_Step46", "_Step148_with_proxy"), all = FALSE )
reproduction_audit <- data.frame( metric = c( "genes_compared_n", "Pearson_r_coefficient", "Pearson_r_moderated_t", "max_abs_coefficient_difference", "max_abs_t_difference",
  "max_abs_P_value_difference", "max_abs_adjusted_P_difference" ), value = c( nrow(step46_cmp), cor(step46_cmp$logFC, step46_cmp$transformed_scale_difference, use = "complete.obs"),
  cor(step46_cmp$t_Step46, step46_cmp$t_Step148_with_proxy, use = "complete.obs"), max(abs(step46_cmp$logFC - step46_cmp$transformed_scale_difference), na.rm = TRUE),
  max(abs(step46_cmp$t_Step46 - step46_cmp$t_Step148_with_proxy), na.rm = TRUE), max(abs(step46_cmp$P.Value_Step46 - step46_cmp$P.Value_Step148_with_proxy), na.rm = TRUE),
  max(abs(step46_cmp$adj.P.Val_Step46 - step46_cmp$adj.P.Val_Step148_with_proxy), na.rm = TRUE) ), stringsAsFactors = FALSE )
write_csv(reproduction_audit, 16, "Step46_with_proxy_reproduction_audit.csv")

interpretation_notes <- data.frame( topic = c( "Operational-group definition", "Proxy status", "Primary comparison", "WITH-proxy model", "WITHOUT-proxy model", "Coefficient interpretation",
  "Significance definition", "Absolute effect threshold", "Cross-tab limitation", "Identifiability limitation", "Primary-phase balance", "Permitted conclusion" ), note = c(
  "Group A = Step45 Batch_A (sample IDs beginning with jerry. or normal.); Group B = Step45 Batch_B (sample IDs beginning with D).",
  "The operational grouping is sample-ID-derived and is NOT a verified experimental batch assignment.",
  "Established late (14 + 28 dpi) versus Early (Uninjured + 15 min + 1 dpi + 3 dpi), same 29 samples in both models.",
  "~ 0 + primary_group + operational_group_proxy; arrayWeights; limma; eBayes(trend=TRUE, robust=TRUE).", "~ 0 + primary_group; arrayWeights; limma; eBayes(trend=TRUE, robust=TRUE).",
  "Difference on the supplied transformed expression scale; not interpreted as log2 fold change.", "BH-adjusted P < 0.05; Higher/Lower determined only by coefficient sign.", "None.",
  "Both proxy groups occur at 5 of 8 time points. Uninjured is represented only in Group A; 15 min and 1 dpi only in Group B.",
  "A full-rank model establishes mathematical estimability but does not prove that stage and technical effects are biologically separable. Individual Uninjured-vs-15-min and Uninjured-vs-1-dpi comparisons cannot independently identify stage and proxy effects.",
  "Within the 29-sample primary comparison, Early and EstablishedLate both contain both proxy groups, making the phase-level contrast more balanced than the fully confounded individual early time points.",
  "Compare WITH- and WITHOUT-proxy results as a sensitivity analysis. Do not describe the proxy as a validated batch correction or claim that the comparison proves removal of technical confounding."
  ), stringsAsFactors = FALSE )
write_csv(interpretation_notes, 17, "interpretation_notes.csv")

scatter_pdf <- out_path(20, "transformed_scale_difference_concordance.pdf")
scatter_png <- out_path(20, "transformed_scale_difference_concordance.png")
x <- paired$transformed_scale_difference_without_proxy
y <- paired$transformed_scale_difference_with_proxy
finite_xy <- is.finite(x) & is.finite(y)
xy_lim <- range(c(x[finite_xy], y[finite_xy]), finite = TRUE)

plot_scatter <- function() {
  plot( x[finite_xy], y[finite_xy], pch = 16, cex = 0.45, col = grDevices::adjustcolor("grey35", alpha.f = 0.30), xlab = "WITHOUT proxy: difference on supplied transformed scale",
    ylab = "WITH proxy: difference on supplied transformed scale", xlim = xy_lim, ylim = xy_lim, asp = 1, main = "Primary EstablishedLate vs Early" )
  abline(a = 0, b = 1, lty = 2, lwd = 1.2, col = "grey45")
  core_idx <- which(paired$gene_symbol_restored %in% CORE_GENES)
  points(x[core_idx], y[core_idx], pch = 19, cex = 0.8)
  text( x[core_idx], y[core_idx], labels = paired$gene_symbol_restored[core_idx], pos = 4, cex = 0.75 )
  legend( "topleft", legend = c( paste0("Pearson r = ", sprintf("%.4f", coef_pearson)), paste0("Spearman rho = ", sprintf("%.4f", coef_spearman)) ), bty = "n" )
}

grDevices::pdf(scatter_pdf, width = 6.5, height = 6.5)
plot_scatter()
grDevices::dev.off()

grDevices::png(scatter_png, width = 1950, height = 1950, res = 300)
plot_scatter()
grDevices::dev.off()

output_index <- data.frame( No = c( "00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "20a", "20b", "31", "32" ), file = c(
  paste0("00_", PREFIX, "_OUTPUT_INDEX.csv"), num_name(1, "sample_set_audit.csv"), num_name(2, "proxy_by_timepoint_cross_tab.csv"), num_name(3, "proxy_by_primary_phase_cross_tab.csv"),
  num_name(4, "design_rank_and_identifiability_audit.csv"), num_name(5, "arrayWeights_with_vs_without_proxy.csv"), num_name(6, "PRIMARY_with_proxy_all_genes.csv"),
  num_name(7, "PRIMARY_without_proxy_all_genes.csv"), num_name(8, "model_FDR_count_summary.csv"), num_name(9, "gene_level_with_vs_without_proxy.csv"),
  num_name(10, "with_vs_without_proxy_concordance_summary.csv"), num_name(11, "FDR_significant_gene_set_overlap.csv"), num_name(12, "six_core_genes_with_vs_without_proxy.csv"),
  num_name(13, "Clec7a_signature_GSEA_with_vs_without_proxy.csv"), num_name(14, "Clec7a_signature_GSEA_leading_edges.csv"), num_name(15, "Clec7a_signature_GSEA_concordance.csv"),
  num_name(16, "Step46_with_proxy_reproduction_audit.csv"), num_name(17, "interpretation_notes.csv"), num_name(20, "transformed_scale_difference_concordance.pdf"),
  num_name(20, "transformed_scale_difference_concordance.png"), num_name(31, "REPORT.txt"), num_name(32, "sessionInfo.txt") ), purpose = c( "Output index",
  "All 39 samples and primary-comparison membership", "Reviewer-requested proxy x timepoint cross-tabulation", "Proxy x primary-phase cross-tabulation",
  "Design rank and identifiability limitations", "Sample-specific arrayWeights in both primary models", "Primary model WITH proxy: all genes", "Primary model WITHOUT proxy: all genes",
  "FDR-only Higher/Lower counts in both models", "Paired all-gene comparison", "Coefficient/t-statistic/significance concordance summary", "FDR-significant gene-set overlap and Jaccard",
  "Six prespecified lipid-handling genes in both models", "Frozen Clec7a-UP signature GSEA in both models", "GSEA leading-edge genes", "GSEA NES and leading-edge concordance",
  "Audit that WITH-proxy model reproduces frozen Step46", "Reviewer-facing interpretation and identifiability notes", "Diagnostic coefficient concordance scatter PDF",
  "Diagnostic coefficient concordance scatter PNG", "Text report", "R session information" ), stringsAsFactors = FALSE )
write.csv( output_index, out_path(0, "OUTPUT_INDEX.csv"), row.names = FALSE, fileEncoding = "UTF-8" )

with_counts <- count_model(fit_with$table, "WITH proxy")
without_counts <- count_model(fit_without$table, "WITHOUT proxy")

report_lines <- c( "Reviewer 1 Remaining Major Comment 4 - operational-group proxy sensitivity", paste0("Script: ", PREFIX, ".R"), "", "SCOPE",
  "- Same 29-sample primary EstablishedLate vs Early contrast WITH and WITHOUT the unverified operational-group proxy.",
  "- The only intended design change is inclusion/exclusion of the proxy term.", "- Both models use the original supplied transformed expression matrix from frozen Step45.",
  "- Both models use model-specific arrayWeights and limma eBayes(trend=TRUE, robust=TRUE).", "- No ZIP archive is used.", "", "COEFFICIENT INTERPRETATION",
  "- Exported coefficient: transformed_scale_difference", "- Interpretation: difference on the supplied transformed expression scale", "- No log2 fold-change interpretation is assigned.",
  "- No absolute effect-size threshold is used.", "- Statistical significance: BH-adjusted P < 0.05.", "", "PROXY DEFINITION",
  "- Group A corresponds to Step45 Batch_A: sample IDs beginning with jerry. or normal.", "- Group B corresponds to Step45 Batch_B: sample IDs beginning with D.",
  "- This is an unverified sample-ID-derived operational grouping, not a validated experimental batch.", "", "PROXY x TIMEPOINT CROSS-TAB",
  paste(capture.output(print(proxy_by_timepoint, row.names = FALSE)), collapse = "\n"), "", "PRIMARY PHASE x PROXY CROSS-TAB",
  paste(capture.output(print(proxy_by_primary_phase, row.names = FALSE)), collapse = "\n"), "", "IDENTIFIABILITY",
  paste0("- Both proxy groups are represented at ", both_timepoints_n, " of ", nrow(proxy_by_timepoint), " time points."),
  paste0("- Time points represented by only one proxy group: ", paste(single_group_timepoints, collapse = ", "), "."),
  "- Uninjured occurs only in Group A; 15 min and 1 dpi occur only in Group B.",
  "- Therefore individual Uninjured-vs-15-min and Uninjured-vs-1-dpi comparisons cannot independently identify stage and proxy effects.",
  "- A full-rank design establishes mathematical estimability, not biological separation of stage and technical effects.",
  "- The primary 29-sample phase contrast contains both proxy groups in both Early and EstablishedLate and is more balanced, but proxy adjustment remains a sensitivity analysis.", "",
  "PRIMARY MODEL COUNTS", paste(capture.output(print(rbind(with_counts, without_counts), row.names = FALSE)), collapse = "\n"), "", "WITH vs WITHOUT PROXY CONCORDANCE",
  paste0("- Pearson r of transformed-scale differences: ", signif(coef_pearson, 6)), paste0("- Spearman rho of transformed-scale differences: ", signif(coef_spearman, 6)),
  paste0("- Pearson r of moderated t statistics: ", signif(t_pearson, 6)), paste0("- Same-direction genes: ", sprintf("%.2f%%", same_direction_pct)),
  paste0("- FDR-status concordance: ", sprintf("%.2f%%", FDR_status_concordance_pct)), paste0("- FDR-significant gene-set Jaccard: ", signif(jaccard_index(sig_with, sig_without), 6)), "",
  "FROZEN Clec7a-UP SIGNATURE GSEA", paste(capture.output(print( gsea_summary[, c("model", "NES", "pval", "BH_FDR_across_two_proxy_models", "leading_edge_n"), drop = FALSE],
  row.names = FALSE )), collapse = "\n"), "", "STEP46 REPRODUCTION AUDIT", paste(capture.output(print(reproduction_audit, row.names = FALSE)), collapse = "\n"), "", "INTERPRETIVE LIMIT",
  "These results assess sensitivity to an unverified identifier-derived proxy. They do not validate the proxy as an experimental batch, do not prove complete removal of technical confounding, and do not establish independent identifiability of stage and proxy effects at time points represented by only one proxy group."
  )
writeLines(report_lines, out_path(31, "REPORT.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), out_path(32, "sessionInfo.txt"), useBytes = TRUE)

msg("Step148 completed successfully.")
msg("Output directory: ", OUT_DIR)
msg("Coefficient Pearson r (WITH vs WITHOUT proxy): ", signif(coef_pearson, 6))
msg("Coefficient Spearman rho: ", signif(coef_spearman, 6))
msg("FDR-significant set Jaccard: ", signif(jaccard_index(sig_with, sig_without), 6))
msg("Time points with both proxy groups: ", both_timepoints_n, "/", nrow(proxy_by_timepoint))
msg("WITH proxy FDR-significant genes: ", length(sig_with))
msg("WITHOUT proxy FDR-significant genes: ", length(sig_without))
msg("No absolute effect-size cutoff was used.")
msg("Bulk coefficients are interpreted only on the supplied transformed scale.")

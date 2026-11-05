source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

STEP44_DIR <- file.path(ROOT, "results", "44_Bulk_17702045_input_QC")

EXPR_FILE <- file.path(STEP44_DIR, "44_Bulk_17702045_input_QC_v04_expression_matrix.rds")

META_FILE <- file.path(STEP44_DIR, "44_Bulk_17702045_input_QC_v04_metadata_clean.rds")

OUT_DIR <- file.path(ROOT, "results", "45_Bulk_batch_structure_and_analysis_object")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "45_Bulk_batch_structure_and_analysis_object_v02"

TIME_LEVELS <- c("Uninjured", "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi")

PHASE_LEVELS <- c("Early", "Transition", "Established late", "Extended late")

PRIMARY_PHASE_LEVELS <- c("Early", "Established late")

N_TOP_PCA_GENES <- 3000

required_pkgs <- c("limma", "ggplot2")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them before rerunning.")
}

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

make_pca <- function(mat, top_genes) {
    x <- mat[top_genes, , drop = FALSE]
    pca_input <- t(x)
    keep_cols <- apply(pca_input, 2, function(z) {
        s <- sd(z, na.rm = TRUE)
        is.finite(s) && s > 0
    })
    pca_input <- pca_input[, keep_cols, drop = FALSE]
    prcomp(pca_input, center = TRUE, scale. = TRUE)
}

pca_table <- function(pca, meta, label) {
    ve <- 100 * pca$sdev^2/sum(pca$sdev^2)
    out <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 
        3], PC4 = pca$x[, 4], PC5 = pca$x[, 5], matrix_type = label, stringsAsFactors = FALSE)
    out$timepoint <- as.character(meta[out$sample_id, "timepoint"])
    out$batch <- as.character(meta[out$sample_id, "batch"])
    out$phase <- as.character(meta[out$sample_id, "phase"])
    attr(out, "variance_explained") <- ve
    out
}

plot_pca <- function(df, variance_explained, title_text) {
    ggplot(df, aes(x = PC1, y = PC2, color = timepoint, shape = batch, label = sample_id)) + geom_point(size = 3) + 
        geom_text(size = 2.4, vjust = -0.75, check_overlap = TRUE) + labs(title = title_text, x = sprintf("PC1 (%.1f%%)", 
        variance_explained[1]), y = sprintf("PC2 (%.1f%%)", variance_explained[2]), color = "Time", shape = "Batch") + 
        theme_classic(base_family = "sans", base_size = 10) + theme(text = element_text(family = "sans", 
        color = "black"))
}

write_matrix_csv <- function(x, path) {
    write.csv(x, path, row.names = TRUE, fileEncoding = "UTF-8")
}

model_diagnostics <- function(design, model_name) {
    qr_obj <- qr(design)
    data.frame(model = model_name, n_samples = nrow(design), n_parameters = ncol(design), rank = qr_obj$rank, 
        full_rank = qr_obj$rank == ncol(design), stringsAsFactors = FALSE)
}

if (!file.exists(EXPR_FILE)) {
  stop("Cannot find Step44 expression matrix:\n", EXPR_FILE)
}

if (!file.exists(META_FILE)) {
  stop("Cannot find Step44 metadata:\n", META_FILE)
}

expr <- readRDS(EXPR_FILE)

meta0 <- readRDS(META_FILE)

expr <- as.matrix(expr)

storage.mode(expr) <- "numeric"

if (nrow(expr) == 0 || ncol(expr) == 0) {
  stop("Expression matrix is empty.")
}

if (!all(colnames(expr) %in% rownames(meta0))) {
  stop("Expression sample IDs do not fully match Step44 metadata.")
}

meta0 <- meta0[colnames(expr), , drop = FALSE]

sample_id <- colnames(expr)

batch <- rep(NA_character_, length(sample_id))

batch[grepl("^(jerry\\.|normal\\.)", sample_id)] <- "Batch_A"

batch[grepl("^D", sample_id)] <- "Batch_B"

if (anyNA(batch)) {
    bad <- sample_id[is.na(batch)]
    stop("The following sample IDs could not be assigned to Batch_A/Batch_B:\n", paste(bad, collapse = "\n"), 
        "\nEdit the batch definition before continuing.")
}

if (!"timepoint_std" %in% colnames(meta0)) {
  stop("Step44 metadata does not contain timepoint_std.")
}

timepoint <- as.character(meta0$timepoint_std)

if (anyNA(timepoint) || any(!timepoint %in% TIME_LEVELS)) {
  stop("Some timepoint_std values are missing or outside the TIME_LEVELS.")
}

phase <- rep(NA_character_, length(timepoint))

phase[timepoint %in% c("Uninjured", "15 min", "1 dpi", "3 dpi")] <- "Early"

phase[timepoint == "7 dpi"] <- "Transition"

phase[timepoint %in% c("14 dpi", "28 dpi")] <- "Established late"

phase[timepoint == "42 dpi"] <- "Extended late"

phase_primary <- ifelse(phase %in% PRIMARY_PHASE_LEVELS, phase, NA_character_)

meta <- data.frame(sample_id = sample_id, batch = factor(batch, levels = c("Batch_A", "Batch_B")), timepoint = factor(timepoint,
  levels = TIME_LEVELS), phase = factor(phase, levels = PHASE_LEVELS), phase_primary = factor(phase_primary,
  levels = PRIMARY_PHASE_LEVELS), stringsAsFactors = FALSE)

rownames(meta) <- meta$sample_id

for (nm in colnames(meta0)) {
    new_nm <- paste0("step44_", nm)
    meta[[new_nm]] <- meta0[rownames(meta), nm]
}

write.csv(meta, file.path(OUT_DIR, paste0(PREFIX, "_formal_metadata.csv")), row.names = FALSE, fileEncoding = "UTF-8")

batch_time <- as.data.frame.matrix(table(meta$batch, meta$timepoint))

batch_time$batch <- rownames(batch_time)

batch_time <- batch_time[, c("batch", TIME_LEVELS), drop = FALSE]

write.csv(batch_time, file.path(OUT_DIR, paste0(PREFIX, "_batch_by_timepoint.csv")), row.names = FALSE)

batch_phase <- as.data.frame.matrix(table(meta$batch, meta$phase))

batch_phase$batch <- rownames(batch_phase)

write.csv(batch_phase, file.path(OUT_DIR, paste0(PREFIX, "_batch_by_phase.csv")), row.names = FALSE)

meta$batch <- relevel(meta$batch, ref = "Batch_A")

meta$timepoint <- relevel(meta$timepoint, ref = "Uninjured")

design_full_time <- model.matrix(~batch + timepoint, data = meta)

primary_keep <- !is.na(meta$phase_primary)

meta_primary <- droplevels(meta[primary_keep, , drop = FALSE])

meta_primary$batch <- relevel(meta_primary$batch, ref = "Batch_A")

meta_primary$phase_primary <- relevel(meta_primary$phase_primary, ref = "Early")

design_primary <- model.matrix(~batch + phase_primary, data = meta_primary)

design_time_only <- model.matrix(~timepoint, data = meta)

design_diag <- rbind(model_diagnostics(design_full_time, "~ batch + timepoint"), model_diagnostics(design_primary,
  "~ batch + phase_primary [Early vs Established late]"), model_diagnostics(design_time_only, "~ timepoint"))

write.csv(design_diag, file.path(OUT_DIR, paste0(PREFIX, "_design_rank_diagnostics.csv")), row.names = FALSE)

write_matrix_csv(design_full_time, file.path(OUT_DIR, paste0(PREFIX, "_design_full_time.csv")))

write_matrix_csv(design_primary, file.path(OUT_DIR, paste0(PREFIX, "_design_primary_phase.csv")))

if (any(!design_diag$full_rank[design_diag$model %in% c("~ batch + timepoint", "~ batch + phase_primary [Early vs Established late]")])) {
    stop(paste0("At least one formal downstream design is not full rank.\n", "Inspect *_design_rank_diagnostics.csv before proceeding."))
}

gene_names <- rownames(expr)

x_numeric <- grepl("^X[0-9]", gene_names)

gene_name_audit <- data.frame(original_gene_name = gene_names, starts_X_numeric = x_numeric, proposed_restored_name = ifelse(x_numeric,
  sub("^X", "", gene_names), gene_names), stringsAsFactors = FALSE)

write.csv(gene_name_audit[gene_name_audit$starts_X_numeric, , drop = FALSE], file.path(OUT_DIR, paste0(PREFIX,
  "_Xnumeric_gene_name_audit.csv")), row.names = FALSE)

gene_var <- apply(expr, 1, var, na.rm = TRUE)

gene_var[!is.finite(gene_var)] <- 0

variable_genes <- names(sort(gene_var, decreasing = TRUE))

variable_genes <- variable_genes[gene_var[variable_genes] > 0]

top_genes <- head(variable_genes, min(N_TOP_PCA_GENES, length(variable_genes)))

if (length(top_genes) < 100) {
  stop("Too few variable genes for PCA.")
}

writeLines(top_genes, con = file.path(OUT_DIR, paste0(PREFIX, "_PCA_top_variable_genes.txt")))

pca_raw <- make_pca(expr, top_genes)

pca_raw_df <- pca_table(pca_raw, meta, "Original transformed expression")

ve_raw <- attr(pca_raw_df, "variance_explained")

write.csv(pca_raw_df, file.path(OUT_DIR, paste0(PREFIX, "_PCA_original.csv")), row.names = FALSE)

expr_batch_adjusted <- limma::removeBatchEffect(expr, batch = meta$batch, design = design_time_only)

pca_adj <- make_pca(expr_batch_adjusted, top_genes)

pca_adj_df <- pca_table(pca_adj, meta, "Batch-adjusted visualization matrix")

ve_adj <- attr(pca_adj_df, "variance_explained")

write.csv(pca_adj_df, file.path(OUT_DIR, paste0(PREFIX, "_PCA_batch_adjusted.csv")), row.names = FALSE)

cor_raw <- cor(expr, use = "pairwise.complete.obs", method = "pearson")

cor_adj <- cor(expr_batch_adjusted, use = "pairwise.complete.obs", method = "pearson")

write_matrix_csv(cor_raw, file.path(OUT_DIR, paste0(PREFIX, "_correlation_original.csv")))

write_matrix_csv(cor_adj, file.path(OUT_DIR, paste0(PREFIX, "_correlation_batch_adjusted.csv")))

median_cor_all <- vapply(seq_len(ncol(cor_raw)), function(i) {
  median(cor_raw[i, -i], na.rm = TRUE)
}, numeric(1))

names(median_cor_all) <- colnames(cor_raw)

median_cor_same_time <- vapply(rownames(meta), function(sid) {
    tp <- meta[sid, "timepoint"]
    same_ids <- rownames(meta)[meta$timepoint == tp]
    same_ids <- setdiff(same_ids, sid)
    if (length(same_ids) == 0) {
        return(NA_real_)
    }
    median(cor_raw[sid, same_ids], na.rm = TRUE)
}, numeric(1))

pc5 <- pca_raw_df[, c("PC1", "PC2", "PC3", "PC4", "PC5")]

rownames(pc5) <- pca_raw_df$sample_id

pca_distance_to_time_centroid <- vapply(rownames(meta), function(sid) {
    tp <- as.character(meta[sid, "timepoint"])
    ids <- rownames(meta)[as.character(meta$timepoint) == tp]
    centroid <- colMeans(pc5[ids, , drop = FALSE])
    sqrt(sum((as.numeric(pc5[sid, , drop = TRUE]) - centroid)^2))
}, numeric(1))

outlier_df <- data.frame(sample_id = rownames(meta), batch = as.character(meta$batch), timepoint = as.character(meta$timepoint),
  phase = as.character(meta$phase), median_correlation_all = median_cor_all[rownames(meta)], median_correlation_same_time = median_cor_same_time[rownames(meta)],
  pca5_distance_to_time_centroid = pca_distance_to_time_centroid[rownames(meta)], stringsAsFactors = FALSE)

within_cor_cut <- quantile(outlier_df$median_correlation_same_time, 0.25, na.rm = TRUE) - 1.5 * IQR(outlier_df$median_correlation_same_time,
  na.rm = TRUE)

pca_dist_cut <- quantile(outlier_df$pca5_distance_to_time_centroid, 0.75, na.rm = TRUE) + 1.5 * IQR(outlier_df$pca5_distance_to_time_centroid,
  na.rm = TRUE)

outlier_df$review_low_same_time_correlation <- outlier_df$median_correlation_same_time < within_cor_cut

outlier_df$review_high_PCA_distance <- outlier_df$pca5_distance_to_time_centroid > pca_dist_cut

outlier_df$review_flag <- outlier_df$review_low_same_time_correlation | outlier_df$review_high_PCA_distance

write.csv(outlier_df, file.path(OUT_DIR, paste0(PREFIX, "_sample_outlier_diagnostics.csv")), row.names = FALSE)

pdf(file.path(OUT_DIR, paste0(PREFIX, "_diagnostics.pdf")), width = 10, height = 7.5, onefile = TRUE,
  family = "sans")

bt_long <- as.data.frame(table(Batch = meta$batch, Time = meta$timepoint))

p_bt <- ggplot(bt_long, aes(x = Time, y = Batch, size = Freq)) + geom_point() + geom_text(aes(label = Freq),
  vjust = -1, size = 3) + labs(title = "Bulk sample structure: batch by time point", x = NULL, y = NULL,
  size = "n") + theme_classic(base_family = "sans", base_size = 10) + theme(axis.text.x = element_text(angle = 45,
  hjust = 1))

print(p_bt)

print(plot_pca(pca_raw_df, ve_raw, "PCA before batch adjustment"))

print(plot_pca(pca_adj_df, ve_adj, "PCA after batch adjustment (visualization only)"))

p_out1 <- ggplot(outlier_df, aes(x = reorder(sample_id, median_correlation_same_time), y = median_correlation_same_time,
  shape = batch)) + geom_point(size = 2.6) + geom_hline(yintercept = within_cor_cut, linetype = 2) +
  coord_flip() + labs(title = "Median correlation with samples from the same time point", x = NULL,
  y = "Median Pearson correlation", shape = "Batch") + theme_classic(base_family = "sans", base_size = 9)

print(p_out1)

p_out2 <- ggplot(outlier_df, aes(x = reorder(sample_id, pca5_distance_to_time_centroid), y = pca5_distance_to_time_centroid,
  shape = batch)) + geom_point(size = 2.6) + geom_hline(yintercept = pca_dist_cut, linetype = 2) +
  coord_flip() + labs(title = "Distance to time-point centroid in first five PCs", x = NULL, y = "Euclidean distance",
  shape = "Batch") + theme_classic(base_family = "sans", base_size = 9)

print(p_out2)

op <- par(mar = c(7, 7, 3, 2))

heatmap(cor_raw, Rowv = NA, Colv = NA, scale = "none", margins = c(8, 8), main = "Sample correlation before batch adjustment")

heatmap(cor_adj, Rowv = NA, Colv = NA, scale = "none", margins = c(8, 8), main = "Sample correlation after batch adjustment")

par(op)

dev.off()

analysis_object <- list(expression_original = expr, expression_batch_adjusted_visualization = expr_batch_adjusted,
  metadata = meta, metadata_primary = meta_primary, design_full_time = design_full_time, design_primary = design_primary,
  design_time_only = design_time_only, pca_top_genes = top_genes, pca_original = pca_raw, pca_batch_adjusted = pca_adj,
  sample_correlation_original = cor_raw, sample_correlation_batch_adjusted = cor_adj, sample_outlier_diagnostics = outlier_df,
  time_levels = TIME_LEVELS, phase_levels = PHASE_LEVELS, primary_phase_levels = PRIMARY_PHASE_LEVELS,
  notes = c("Input expression matrix is transformed continuous expression, not raw counts.", "Use limma, not DESeq2, for downstream differential analysis.",
    "Batch-adjusted expression is for visualization/QC only.", "Formal testing must use expression_original with batch in the design.",
    "Primary comparison is Established late (14+28 dpi) vs Early (Uninjured+15 min+1+3 dpi).", "7 dpi is Transition and 42 dpi is Extended late/persistence.",
    "No samples were automatically excluded in Step45."))

saveRDS(analysis_object, file.path(OUT_DIR, paste0(PREFIX, "_analysis_object.rds")), compress = FALSE)

flagged <- outlier_df[outlier_df$review_flag, , drop = FALSE]

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = file.path(OUT_DIR, paste0(PREFIX, "_REPORT.txt")))

capture.output(sessionInfo(), file = file.path(OUT_DIR, paste0(PREFIX, "_sessionInfo.txt")))

print(table(meta$batch, meta$timepoint))

print(design_diag, row.names = FALSE)

if (nrow(flagged) == 0) {
} else {
  print(flagged[, c("sample_id", "batch", "timepoint", "median_correlation_same_time", "pca5_distance_to_time_centroid",
    "review_flag"), drop = FALSE], row.names = FALSE)
}

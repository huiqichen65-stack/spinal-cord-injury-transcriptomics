source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

SIGNATURE_FILE <- file.path(ROOT, "results", "21_Clec7a_DEG_volcano_myeloid", "21_Clec7a_vs_all_myeloid_DEG_up_v01.csv")

STEP46_FILE <- file.path(ROOT, "results", "46_Bulk_limma_temporal_DEG_analysis", "30_46_Bulk_limma_temporal_DEG_analysis_v06_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "47_Bulk_scRNA_Clec7a_signature_GSEA")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "47_Bulk_scRNA_Clec7a_signature_GSEA_v02"

EXPECTED_SIGNATURE_N <- 1320L

MATCHED_TIME_CONTRASTS <- c("1dpi_vs_Uninjured", "3dpi_vs_Uninjured", "7dpi_vs_Uninjured", "14dpi_vs_Uninjured",
  "28dpi_vs_Uninjured")

MATCHED_TIME_LABELS <- c("1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

required_pkgs <- c("fgsea")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall fgsea first, e.g. with BiocManager::install('fgsea').")
}

suppressPackageStartupMessages({
  library(fgsea)
})

num_name <- function(num, suffix) {
  paste0(sprintf("%02d", num), "_", PREFIX, "_", suffix)
}

write_csv <- function(x, num, suffix) {
  write.csv(x, file.path(OUT_DIR, num_name(num, suffix)), row.names = FALSE, fileEncoding = "UTF-8")
}

safe_gene_symbol <- clean_gene_symbol

collapse_rank_table <- function(tt) {

  required_cols <- c("gene_symbol_restored", "logFC", "t", "P.Value", "adj.P.Val")

  if (!all(required_cols %in% colnames(tt))) {
    stop("Bulk contrast table is missing required columns: ", paste(setdiff(required_cols, colnames(tt)),
      collapse = ", "))
  }

  x <- data.frame(gene = safe_gene_symbol(as.character(tt$gene_symbol_restored)), logFC = as.numeric(tt$logFC),
    t = as.numeric(tt$t), P.Value = as.numeric(tt$P.Value), adj.P.Val = as.numeric(tt$adj.P.Val),
    stringsAsFactors = FALSE)

  x <- x[!is.na(x$gene) & nzchar(x$gene) & is.finite(x$t), , drop = FALSE]

  ord <- order(abs(x$t), decreasing = TRUE)

  x <- x[ord, , drop = FALSE]

  x <- x[!duplicated(x$gene), , drop = FALSE]

  x <- x[order(x$t, decreasing = TRUE), , drop = FALSE]

  rownames(x) <- NULL

  x
}

make_rank_vector <- function(rank_table) {
    ranks <- rank_table$t
    names(ranks) <- rank_table$gene
    ranks <- ranks[is.finite(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    ranks
}

run_one_gsea <- function(tt, contrast_id, contrast_label, role, caution, signature_genes) {
    rank_table <- collapse_rank_table(tt)
    ranks <- make_rank_vector(rank_table)
    overlap <- intersect(signature_genes, names(ranks))
    if (length(overlap) < 15) {
        stop("Too few signature genes overlap the ranked bulk list for ", contrast_id, ": ", length(overlap))
    }
    pathways <- list(Clec7a_UP_scRNA = overlap)
    set.seed(20260929)
    fg <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15, maxSize = 5000, eps = 0, 
        scoreType = "std")
    if (nrow(fg) != 1) {
        stop("Unexpected fgsea result row count for ", contrast_id, ": ", nrow(fg))
    }
    leading <- fg$leadingEdge[[1]]
    fg_out <- data.frame(contrast_id = contrast_id, contrast_label = contrast_label, role = role, caution = caution, 
        pathway = as.character(fg$pathway), NES = as.numeric(fg$NES), ES = as.numeric(fg$ES), pval = as.numeric(fg$pval), 
        fgsea_padj_within_call = as.numeric(fg$padj), signature_n_total = length(signature_genes), signature_n_in_bulk = length(overlap), 
        signature_overlap_fraction = length(overlap)/length(signature_genes), leading_edge_n = length(leading), 
        leading_edge_fraction_of_overlap = length(leading)/length(overlap), leading_edge_genes = paste(leading, 
            collapse = ";"), stringsAsFactors = FALSE)
    leading_df <- rank_table[match(leading, rank_table$gene), , drop = FALSE]
    leading_df$contrast_id <- contrast_id
    leading_df$contrast_label <- contrast_label
    leading_df$leading_edge_order <- seq_len(nrow(leading_df))
    leading_df <- leading_df[, c("contrast_id", "contrast_label", "leading_edge_order", "gene", "logFC", 
        "t", "P.Value", "adj.P.Val"), drop = FALSE]
    list(summary = fg_out, leading_edge = leading_df, rank_table = rank_table, overlap = overlap)
}

bind_rows_base <- function(lst) {
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
        x[, all_names, drop = FALSE]
    })
    out <- do.call(rbind, lst2)
    rownames(out) <- NULL
    out
}

if (!file.exists(SIGNATURE_FILE)) {
  stop("Cannot find frozen Step21 UP signature file:\n", SIGNATURE_FILE)
}

if (!file.exists(STEP46_FILE)) {
  stop("Cannot find Step46 analysis object:\n", STEP46_FILE)
}

sig_raw <- read.csv(SIGNATURE_FILE, check.names = FALSE, stringsAsFactors = FALSE)

required_sig_cols <- c("gene", "avg_log2FC", "p_val_adj")

if (!all(required_sig_cols %in% colnames(sig_raw))) {
  stop("Step21 UP signature file is missing required columns: ", paste(setdiff(required_sig_cols, colnames(sig_raw)),
    collapse = ", "))
}

sig <- sig_raw[!is.na(sig_raw$gene) & nzchar(sig_raw$gene), , drop = FALSE]

sig$gene <- safe_gene_symbol(as.character(sig$gene))

sig <- sig[!duplicated(sig$gene), , drop = FALSE]

if (nrow(sig) != EXPECTED_SIGNATURE_N) {
  stop("Frozen Step21 signature size mismatch: expected ", EXPECTED_SIGNATURE_N, " unique UP genes, but found ",
    nrow(sig), ". Do not continue until this is resolved.")
}

if (any(sig$p_val_adj >= 0.05, na.rm = TRUE)) {
  stop("The Step21 UP signature file contains p_val_adj >= 0.05.")
}

if (any(sig$avg_log2FC < 0.25, na.rm = TRUE)) {
  stop("The Step21 UP signature file contains avg_log2FC < 0.25.")
}

signature_genes <- unique(sig$gene)

bulk <- readRDS(STEP46_FILE)

required_bulk_slots <- c("primary", "primary_sensitivity", "transition", "extended_late", "full_time")

if (!all(required_bulk_slots %in% names(bulk))) {
  stop("Step46 object is missing required slots: ", paste(setdiff(required_bulk_slots, names(bulk)),
    collapse = ", "))
}

if (is.null(bulk$full_time$contrast_tables)) {
  stop("Step46 full_time$contrast_tables is missing.")
}

missing_time <- setdiff(MATCHED_TIME_CONTRASTS, names(bulk$full_time$contrast_tables))

if (length(missing_time) > 0) {
  stop("Step46 object is missing matched temporal contrast(s): ", paste(missing_time, collapse = ", "))
}

primary_rank_table <- collapse_rank_table(bulk$primary$table)

bulk_genes <- primary_rank_table$gene

sig$in_bulk <- sig$gene %in% bulk_genes

sig$bulk_primary_t <- primary_rank_table$t[match(sig$gene, primary_rank_table$gene)]

sig$bulk_primary_logFC <- primary_rank_table$logFC[match(sig$gene, primary_rank_table$gene)]

sig$bulk_primary_adjP <- primary_rank_table$adj.P.Val[match(sig$gene, primary_rank_table$gene)]

overlap_n <- sum(sig$in_bulk)

signature_audit <- data.frame(metric = c("Frozen Step21 Clec7a-UP signature size", "Genes found in Step46 bulk universe",
  "Genes absent from Step46 bulk universe", "Bulk overlap fraction", "Signature definition p_val_adj threshold",
  "Signature definition avg_log2FC threshold", "Bulk GSEA ranking metric"), value = c(as.character(length(signature_genes)),
  as.character(overlap_n), as.character(length(signature_genes) - overlap_n), sprintf("%.6f", overlap_n/length(signature_genes)),
  "< 0.05", ">= 0.25", "limma moderated t statistic"), stringsAsFactors = FALSE)

write_csv(signature_audit, 1, "signature_audit.csv")

write_csv(sig, 2, "signature_bulk_overlap.csv")

if (overlap_n < 500) {
  warning("Only ", overlap_n, " of ", length(signature_genes), " signature genes are present in bulk.")
}

contrast_specs <- list(list(id = "PRIMARY_EstablishedLate_vs_Early", label = "Established late (14+28 dpi) vs Early",
  role = "PRIMARY", caution = "Primary biological contrast; batch-adjusted limma model.", table = bulk$primary$table),
  list(id = "PRIMARY_sensitivity_excluding_review_flags", label = "Established late vs Early; review-flagged samples excluded",
    role = "ROBUSTNESS", caution = "Sensitivity analysis only; not the primary result.", table = bulk$primary_sensitivity$table),
  list(id = "Transition7dpi_vs_Early", label = "7 dpi vs Early", role = "TRANSITION", caution = "Secondary transition-phase analysis.",
    table = bulk$transition$table), list(id = "Extended42dpi_vs_Early", label = "42 dpi vs Early",
    role = "PERSISTENCE", caution = "Secondary persistence analysis; 42 dpi is not present in scRNA.",
    table = bulk$extended_late$table))

for (i in seq_along(MATCHED_TIME_CONTRASTS)) {
    cid <- MATCHED_TIME_CONTRASTS[i]
    tlab <- MATCHED_TIME_LABELS[i]
    caution_text <- if (cid == "1dpi_vs_Uninjured") {
        "Exploratory only because 1 dpi and Uninjured are strongly confounded with operational batch."
    }
    else if (cid == "3dpi_vs_Uninjured") {
        "Matched scRNA-bulk temporal profile; both batches present at 3 dpi, Uninjured is Batch_A only."
    }
    else {
        "Matched scRNA-bulk temporal profile."
    }
    contrast_specs[[length(contrast_specs) + 1]] <- list(id = cid, label = paste0(tlab, " vs Uninjured"), 
        role = "MATCHED_TEMPORAL_PROFILE", caution = caution_text, table = bulk$full_time$contrast_tables[[cid]])
}

gsea_results <- list()

for (i in seq_along(contrast_specs)) {
    sp <- contrast_specs[[i]]
    gsea_results[[sp$id]] <- run_one_gsea(tt = sp$table, contrast_id = sp$id, contrast_label = sp$label, 
        role = sp$role, caution = sp$caution, signature_genes = signature_genes)
}

gsea_summary <- bind_rows_base(lapply(gsea_results, function(x) {
  x$summary
}))

gsea_summary$BH_FDR_across_planned_contrasts <- p.adjust(gsea_summary$pval, method = "BH")

gsea_summary$enrichment_direction <- ifelse(gsea_summary$NES > 0, "Positive: signature enriched toward genes increased in target group",
  "Negative: signature enriched toward genes decreased in target group")

planned_order <- vapply(contrast_specs, function(x) {
  x$id
}, character(1))

gsea_summary$order_tmp <- match(gsea_summary$contrast_id, planned_order)

gsea_summary <- gsea_summary[order(gsea_summary$order_tmp), , drop = FALSE]

gsea_summary$order_tmp <- NULL
rownames(gsea_summary) <- NULL

write_csv(gsea_summary, 3, "GSEA_summary_all_contrasts.csv")

primary_id <- "PRIMARY_EstablishedLate_vs_Early"

primary_gsea <- gsea_summary[gsea_summary$contrast_id == primary_id, , drop = FALSE]

if (nrow(primary_gsea) != 1) {
  stop("Primary GSEA result is missing or duplicated.")
}

write_csv(primary_gsea, 4, "GSEA_PRIMARY_EstablishedLate_vs_Early.csv")

primary_le <- gsea_results[[primary_id]]$leading_edge

write_csv(primary_le, 5, "GSEA_PRIMARY_leading_edge.csv")

primary_rank <- gsea_results[[primary_id]]$rank_table

primary_rank$is_Clec7a_UP_signature <- primary_rank$gene %in% signature_genes

primary_rank$is_leading_edge <- primary_rank$gene %in% primary_le$gene

write_csv(primary_rank, 6, "GSEA_PRIMARY_ranked_genes.csv")

temporal_profile <- gsea_summary[gsea_summary$contrast_id %in% MATCHED_TIME_CONTRASTS, , drop = FALSE]

temporal_profile$timepoint <- MATCHED_TIME_LABELS[match(temporal_profile$contrast_id, MATCHED_TIME_CONTRASTS)]

temporal_profile$time_order <- match(temporal_profile$contrast_id, MATCHED_TIME_CONTRASTS)

temporal_profile <- temporal_profile[order(temporal_profile$time_order), , drop = FALSE]

rownames(temporal_profile) <- NULL

write_csv(temporal_profile, 7, "GSEA_matched_scRNA_bulk_temporal_profile.csv")

phase_ids <- c("PRIMARY_EstablishedLate_vs_Early", "PRIMARY_sensitivity_excluding_review_flags", "Transition7dpi_vs_Early",
  "Extended42dpi_vs_Early")

phase_summary <- gsea_summary[gsea_summary$contrast_id %in% phase_ids, , drop = FALSE]

write_csv(phase_summary, 8, "GSEA_phase_and_robustness_summary.csv")

all_leading <- bind_rows_base(lapply(gsea_results, function(x) {
  x$leading_edge
}))

write_csv(all_leading, 9, "GSEA_leading_edge_all_contrasts.csv")

primary_sens_id <- "PRIMARY_sensitivity_excluding_review_flags"

primary_row <- gsea_summary[gsea_summary$contrast_id == primary_id, , drop = FALSE]

sens_row <- gsea_summary[gsea_summary$contrast_id == primary_sens_id, , drop = FALSE]

primary_le_genes <- gsea_results[[primary_id]]$leading_edge$gene

sens_le_genes <- gsea_results[[primary_sens_id]]$leading_edge$gene

le_union <- union(primary_le_genes, sens_le_genes)

le_intersect <- intersect(primary_le_genes, sens_le_genes)

primary_sensitivity_concordance <- data.frame(metric = c("PRIMARY_NES", "Sensitivity_NES", "Absolute_NES_difference",
  "PRIMARY_pval", "Sensitivity_pval", "PRIMARY_BH_FDR_across_planned_contrasts", "Sensitivity_BH_FDR_across_planned_contrasts",
  "PRIMARY_leading_edge_n", "Sensitivity_leading_edge_n", "Leading_edge_overlap_n", "Leading_edge_Jaccard"),
  value = c(primary_row$NES, sens_row$NES, abs(primary_row$NES - sens_row$NES), primary_row$pval, sens_row$pval,
    primary_row$BH_FDR_across_planned_contrasts, sens_row$BH_FDR_across_planned_contrasts, length(primary_le_genes),
    length(sens_le_genes), length(le_intersect), if (length(le_union) == 0) {
      NA_real_
    } else {
      length(le_intersect)/length(le_union)
    }), stringsAsFactors = FALSE)

write_csv(primary_sensitivity_concordance, 10, "GSEA_PRIMARY_sensitivity_concordance.csv")

expected_gsea_rows <- length(contrast_specs)

if (nrow(gsea_summary) != expected_gsea_rows) {
  stop("GSEA summary integrity check failed: expected ", expected_gsea_rows, " rows, found ", nrow(gsea_summary))
}

if (nrow(temporal_profile) != length(MATCHED_TIME_CONTRASTS)) {
  stop("Temporal-profile integrity check failed.")
}

if (nrow(primary_le) == 0) {
  stop("Primary leading edge is empty.")
}

step47_object <- list(signature_file = SIGNATURE_FILE, step46_file = STEP46_FILE, signature_table = sig,
  signature_genes = signature_genes, signature_audit = signature_audit, planned_contrasts = contrast_specs,
  gsea_results = gsea_results, gsea_summary = gsea_summary, temporal_profile = temporal_profile, phase_summary = phase_summary,
  all_leading_edges = all_leading, primary_sensitivity_concordance = primary_sensitivity_concordance,
  notes = c("Frozen Step21 Clec7a-UP signature contains 1320 genes.", "Bulk ranking metric is the limma moderated t statistic.",
    "Positive NES means the signature is enriched toward genes increased in the target group.", "Primary comparison is Established late (14+28 dpi) vs Early.",
    "7 dpi vs Early is the transition analysis.", "42 dpi vs Early is the persistence analysis.",
    "Matched temporal profile includes 1, 3, 7, 14 and 28 dpi vs Uninjured.", "1 dpi vs Uninjured is exploratory due to operational-batch confounding.",
    "No manuscript figure is generated in Step47."))

saveRDS(step47_object, file.path(OUT_DIR, num_name(30, "analysis_object.rds")), compress = FALSE)

output_index <- data.frame(No = c("00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "30",
  "31", "32"), file = c(paste0("00_", PREFIX, "_OUTPUT_INDEX.csv"), num_name(1, "signature_audit.csv"),
  num_name(2, "signature_bulk_overlap.csv"), num_name(3, "GSEA_summary_all_contrasts.csv"), num_name(4,
    "GSEA_PRIMARY_EstablishedLate_vs_Early.csv"), num_name(5, "GSEA_PRIMARY_leading_edge.csv"), num_name(6,
    "GSEA_PRIMARY_ranked_genes.csv"), num_name(7, "GSEA_matched_scRNA_bulk_temporal_profile.csv"),
  num_name(8, "GSEA_phase_and_robustness_summary.csv"), num_name(9, "GSEA_leading_edge_all_contrasts.csv"),
  num_name(10, "GSEA_PRIMARY_sensitivity_concordance.csv"), num_name(30, "analysis_object.rds"), num_name(31,
    "REPORT.txt"), num_name(32, "sessionInfo.txt")), purpose = c("Output index", "Frozen signature audit",
  "1320-gene signature-to-bulk overlap", "All planned GSEA results", "Primary GSEA result", "Primary leading-edge genes",
  "Primary ranked bulk gene list", "Matched scRNA-bulk temporal GSEA profile", "Primary/transition/persistence/robustness summary",
  "Leading-edge genes for all contrasts", "Primary vs sensitivity GSEA concordance", "Step47 analysis object",
  "Step47 text report", "R session info"), stringsAsFactors = FALSE)

write.csv(output_index, file.path(OUT_DIR, paste0("00_", PREFIX, "_OUTPUT_INDEX.csv")), row.names = FALSE,
  fileEncoding = "UTF-8")

primary_NES <- primary_row$NES
primary_p <- primary_row$pval
primary_FDR <- primary_row$BH_FDR_across_planned_contrasts

writeLines("Numerical results are provided in the accompanying output tables.", 
    con = file.path(OUT_DIR, num_name(31, "REPORT.txt")))

capture.output(sessionInfo(), file = file.path(OUT_DIR, num_name(32, "sessionInfo.txt")))

print(primary_gsea[, c("contrast_label", "NES", "pval", "BH_FDR_across_planned_contrasts", "leading_edge_n"),
  drop = FALSE], row.names = FALSE)

print(temporal_profile[, c("timepoint", "NES", "pval", "BH_FDR_across_planned_contrasts"), drop = FALSE],
  row.names = FALSE)

send_files <- c(num_name(1, "signature_audit.csv"), num_name(3, "GSEA_summary_all_contrasts.csv"), num_name(4,
  "GSEA_PRIMARY_EstablishedLate_vs_Early.csv"), num_name(5, "GSEA_PRIMARY_leading_edge.csv"), num_name(7,
  "GSEA_matched_scRNA_bulk_temporal_profile.csv"))

for (i in seq_along(send_files)) {
}

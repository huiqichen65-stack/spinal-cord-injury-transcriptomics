rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

STEP47_RDS <- file.path(ROOT, "results", "47_Bulk_scRNA_Clec7a_signature_GSEA", "30_47_Bulk_scRNA_Clec7a_signature_GSEA_v02_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "54_Figure3D_scRNA_Clec7a_signature_GSEA")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "54_Figure3D_scRNA_Clec7a_signature_GSEA_v01"

PRIMARY_ID <- "PRIMARY_EstablishedLate_vs_Early"

EXPECTED_SIGNATURE_TOTAL <- 1320L
EXPECTED_SIGNATURE_IN_BULK <- 1028L
EXPECTED_LEADING_EDGE <- 471L

EXPECTED_NES <- 2.30285
EXPECTED_ES <- 0.52723

NES_TOL <- 0.005
ES_TOL <- 0.005

CURVE_COLOR <- "#C94C4C"
LEADING_EDGE_COLOR <- "#C94C4C"
HIT_COLOR <- "#222222"
RANK_POS_COLOR <- "#D9A0A0"
RANK_NEG_COLOR <- "#91ADD0"

install_if_missing <- require_package

for (p in c("ggplot2", "patchwork", "svglite")) {
  install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

msg <- timestamp_message

write_csv <- write_result_csv

sci_text <- function(x, digits = 2) {
    if (!is.finite(x)) {
        return("NA")
    }
    format(x, scientific = TRUE, digits = digits, trim = TRUE)
}

if (!file.exists(STEP47_RDS)) {
  stop("Step47 object not found:\n", STEP47_RDS)
}

obj47 <- readRDS(STEP47_RDS)

required_top <- c("signature_genes", "gsea_results", "gsea_summary")

missing_top <- setdiff(required_top, names(obj47))

if (length(missing_top) > 0L) {
  stop("Step47 object lacks required field(s): ", paste(missing_top, collapse = ", "))
}

if (!PRIMARY_ID %in% names(obj47$gsea_results)) {
  stop("Step47 gsea_results lacks Primary result: ", PRIMARY_ID)
}

primary <- obj47$gsea_results[[PRIMARY_ID]]

required_primary <- c("summary", "leading_edge", "rank_table", "overlap")

missing_primary <- setdiff(required_primary, names(primary))

if (length(missing_primary) > 0L) {
  stop("Step47 Primary result lacks field(s): ", paste(missing_primary, collapse = ", "))
}

signature_genes <- unique(as.character(obj47$signature_genes))

signature_genes <- signature_genes[!is.na(signature_genes) & nzchar(signature_genes)]

rank_table <- primary$rank_table
overlap_genes <- unique(as.character(primary$overlap))

leading_df <- primary$leading_edge
leading_genes <- unique(as.character(leading_df$gene))

summary_row <- obj47$gsea_summary[obj47$gsea_summary$contrast_id == PRIMARY_ID, , drop = FALSE]

if (length(signature_genes) != EXPECTED_SIGNATURE_TOTAL) {
  stop("Frozen signature size mismatch: ", length(signature_genes), " vs expected ", EXPECTED_SIGNATURE_TOTAL)
}

if (length(overlap_genes) != EXPECTED_SIGNATURE_IN_BULK) {
  stop("Frozen signature-to-Bulk overlap mismatch: ", length(overlap_genes), " vs expected ", EXPECTED_SIGNATURE_IN_BULK)
}

if (length(leading_genes) != EXPECTED_LEADING_EDGE) {
  stop("leading-edge count mismatch: ", length(leading_genes), " vs expected ", EXPECTED_LEADING_EDGE)
}

if (nrow(summary_row) != 1L) {
  stop("Primary GSEA summary row is missing or duplicated.")
}

required_rank_cols <- c("gene", "t", "logFC", "P.Value", "adj.P.Val")

missing_rank_cols <- setdiff(required_rank_cols, names(rank_table))

if (length(missing_rank_cols) > 0L) {
  stop("Primary rank table lacks column(s): ", paste(missing_rank_cols, collapse = ", "))
}

if (anyDuplicated(rank_table$gene)) {
  stop("Primary rank table contains duplicated gene symbols.")
}

if (any(!is.finite(rank_table$t)) || any(diff(rank_table$t) > 0)) {
  stop("rank table is not a finite descending moderated-t ranking.")
}

if (!all(overlap_genes %in% rank_table$gene)) {
  stop("Frozen signature overlap contains genes absent from ranked list.")
}

if (!all(leading_genes %in% overlap_genes)) {
  stop("Frozen leading-edge genes are not a subset of the signature overlap.")
}

frozen_NES <- as.numeric(summary_row$NES)

frozen_ES <- as.numeric(summary_row$ES)

frozen_p <- as.numeric(summary_row$pval)

if (!"BH_FDR_across_planned_contrasts" %in% names(summary_row)) {
  stop("Step47 summary lacks BH_FDR_across_planned_contrasts.")
}

frozen_FDR <- as.numeric(summary_row$BH_FDR_across_planned_contrasts)

if (abs(frozen_NES - EXPECTED_NES) > NES_TOL || abs(frozen_ES - EXPECTED_ES) > ES_TOL) {
  stop("Primary GSEA statistic mismatch.\n", "Observed NES=", frozen_NES, ", ES=", frozen_ES, "\nExpected approximately NES=",
    EXPECTED_NES, ", ES=", EXPECTED_ES)
}

rank_genes <- as.character(rank_table$gene)

rank_stats <- as.numeric(rank_table$t)

N <- length(rank_genes)

is_hit <- rank_genes %in% overlap_genes

Nh <- sum(is_hit)

if (Nh != EXPECTED_SIGNATURE_IN_BULK) {
  stop("Hit count mismatch while constructing enrichment curve.")
}

hit_weights <- abs(rank_stats[is_hit])

hit_weight_sum <- sum(hit_weights)

if (!is.finite(hit_weight_sum) || hit_weight_sum <= 0) {
  stop("Invalid hit-weight sum for GSEA curve.")
}

miss_n <- N - Nh

if (miss_n <= 0L) {
  stop("Invalid GSEA ranked-list size / hit count.")
}

step_value <- numeric(N)

step_value[is_hit] <- abs(rank_stats[is_hit])/hit_weight_sum

step_value[!is_hit] <- -1/miss_n

running_ES <- cumsum(step_value)

curve_ES <- if (abs(max(running_ES)) >= abs(min(running_ES))) {
  max(running_ES)
} else {
  min(running_ES)
}

if (abs(curve_ES - frozen_ES) > 0.01) {
  stop("Reconstructed plotting ES does not match Step47 ES closely enough.\n", "Reconstructed ES = ",
    signif(curve_ES, 6), "; ES = ", signif(frozen_ES, 6))
}

peak_index <- if (frozen_ES >= 0) {
  which.max(running_ES)
} else {
  which.min(running_ES)
}

curve_df <- data.frame(rank = seq_len(N), running_ES = running_ES, moderated_t = rank_stats, gene = rank_genes,
  is_signature = is_hit, is_leading_edge = rank_genes %in% leading_genes, stringsAsFactors = FALSE)

hit_df <- curve_df[curve_df$is_signature, , drop = FALSE]

zero_index <- which.min(abs(rank_stats))

positive_n <- sum(rank_stats > 0)

negative_n <- sum(rank_stats < 0)

stats_label <- paste0("NES = ", sprintf("%.2f", frozen_NES), "\nFDR = ", sci_text(frozen_FDR, digits = 3),
  "\nLeading edge = ", length(leading_genes), "/", length(overlap_genes))

p_top <- ggplot(curve_df, aes(x = rank, y = running_ES)) + geom_hline(yintercept = 0, linewidth = 0.35,
  color = "#777777") + geom_line(linewidth = 0.9, color = CURVE_COLOR) + geom_vline(xintercept = peak_index,
  linewidth = 0.35, linetype = "dashed", color = "#777777") + annotate("text", x = N * 0.62, y = max(running_ES) *
  0.76, label = stats_label, family = "Arial", size = 3.2, hjust = 0, vjust = 1, color = "black") +
  scale_x_continuous(limits = c(1, N), expand = c(0, 0)) + labs(x = NULL, y = "Enrichment\nscore") +
  theme_classic(base_family = "Arial", base_size = 10) + theme(plot.title = element_blank(), plot.tag = element_blank(),
  axis.title.y = element_text(family = "Arial", size = 9.5, color = "black"), axis.text.y = element_text(family = "Arial",
    size = 8.2, color = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.line.x = element_blank(), axis.line.y = element_line(linewidth = 0.5, color = "black"), panel.grid = element_blank(),
  plot.margin = margin(t = 5, r = 6, b = 0, l = 6, unit = "pt"))

p_hits <- ggplot(hit_df, aes(x = rank)) + geom_segment(data = hit_df[!hit_df$is_leading_edge, , drop = FALSE],
  aes(x = rank, xend = rank, y = 0, yend = 1), linewidth = 0.28, color = HIT_COLOR, alpha = 0.55) +
  geom_segment(data = hit_df[hit_df$is_leading_edge, , drop = FALSE], aes(x = rank, xend = rank, y = 0,
    yend = 1), linewidth = 0.42, color = LEADING_EDGE_COLOR, alpha = 0.95) + scale_x_continuous(limits = c(1,
  N), expand = c(0, 0)) + coord_cartesian(ylim = c(0, 1)) + labs(x = NULL, y = NULL) + theme_void(base_family = "Arial") +
  theme(plot.margin = margin(t = 0, r = 6, b = 0, l = 6, unit = "pt"))

rank_metric_df <- data.frame(rank = seq_len(N), moderated_t = rank_stats, sign = ifelse(rank_stats >=
  0, "Established late", "Early"), stringsAsFactors = FALSE)

p_rank <- ggplot(rank_metric_df, aes(x = rank, y = moderated_t)) + geom_hline(yintercept = 0, linewidth = 0.3,
  color = "#777777") + geom_area(data = rank_metric_df[rank_metric_df$moderated_t >= 0, , drop = FALSE],
  aes(y = moderated_t), fill = RANK_POS_COLOR, alpha = 0.85) + geom_area(data = rank_metric_df[rank_metric_df$moderated_t <
  0, , drop = FALSE], aes(y = moderated_t), fill = RANK_NEG_COLOR, alpha = 0.85) + geom_vline(xintercept = zero_index,
  linewidth = 0.3, color = "#777777") + annotate("text", x = N * 0.02, y = min(max(rank_stats) * 0.72,
  max(rank_stats) - 0.2), label = "Established late", family = "Arial", size = 2.9, hjust = 0, color = "#9D3B3B") +
  annotate("text", x = N * 0.98, y = max(min(rank_stats) * 0.72, min(rank_stats) + 0.2), label = "Early",
    family = "Arial", size = 2.9, hjust = 1, color = "#416B9B") + scale_x_continuous(limits = c(1,
  N), breaks = scales::pretty_breaks(n = 5), expand = c(0, 0)) + labs(x = "Rank in ordered Bulk gene list",
  y = "Moderated t") + theme_classic(base_family = "Arial", base_size = 10) + theme(axis.title = element_text(family = "Arial",
  size = 9.5, color = "black"), axis.text = element_text(family = "Arial", size = 8, color = "black"),
  axis.line = element_line(linewidth = 0.5, color = "black"), panel.grid = element_blank(), plot.margin = margin(t = 0,
    r = 6, b = 5, l = 6, unit = "pt"))

combined <- p_top/p_hits/p_rank + patchwork::plot_layout(heights = c(3.5, 0.55, 1.6))

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))

PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave(filename = SVG_FILE, plot = combined, device = svglite::svglite, width = 6, height = 5.2, units = "in",
  bg = "white")

ggsave(filename = PDF_FILE, plot = combined, device = grDevices::cairo_pdf, width = 6, height = 5.2,
  units = "in", family = "Arial", bg = "white")

ggsave(filename = PNG_FILE, plot = combined, width = 6, height = 5.2, units = "in", dpi = 600, bg = "white")

qc <- data.frame(metric = c("Step47_source", "primary_contrast", "ranking_metric", "ranked_genes_n",
  "signature_total_n", "signature_in_Bulk_n", "leading_edge_n", "frozen_ES", "reconstructed_plot_ES",
  "frozen_NES", "frozen_p", "frozen_BH_FDR_across_planned_contrasts", "positive_ranked_genes_n", "negative_ranked_genes_n",
  "panel_letter_inside_plot", "internal_plot_title", "GSEA_recomputed"), value = c(STEP47_RDS, "Established late (14+28 dpi) vs Early",
  "limma moderated t", N, length(signature_genes), Nh, length(leading_genes), frozen_ES, curve_ES,
  frozen_NES, frozen_p, frozen_FDR, positive_n, negative_n, "NO", "NO", "NO; running curve reconstructed for visualization only"),
  stringsAsFactors = FALSE)

write_csv(qc, paste0("01_", PREFIX, "_plot_QC.csv"))

write_csv(curve_df, paste0("02_", PREFIX, "_running_enrichment_curve.csv"))

hit_export <- hit_df[, c("rank", "gene", "moderated_t", "is_leading_edge"), drop = FALSE]

write_csv(hit_export, paste0("03_", PREFIX, "_signature_hit_positions.csv"))

write_csv(leading_df, paste0("04_", PREFIX, "_frozen_leading_edge_genes.csv"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, paste0("05_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

msg("NES = ", signif(frozen_NES, 6), "; FDR = ", sci_text(frozen_FDR, 3))

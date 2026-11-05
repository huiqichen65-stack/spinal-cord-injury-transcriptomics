rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

STATE_MARKER_FILE <- file.path(ROOT, "results", "18_myeloid_annotation_cleanup", "18_myeloid_state_markers_v01.csv")

STEP50_RDS <- file.path(ROOT, "results", "50_Bulk_WGCNA_module_trait_integration", "30_50_Bulk_WGCNA_module_trait_integration_v01_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "57_Figure3G_Bulk_WGCNA_module_trait_heatmap")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "57_Figure3G_Bulk_WGCNA_module_trait_heatmap_v05"

EXPECTED_NON_GREY_MODULES <- 11L
EXPECTED_STATES <- 9L
EXPECTED_SELECTED_MODULE <- "brown"

TARGET_STATE_RAW <- "Clec7a+ lipid-associated Microglia"

FIG_WIDTH <- 12.5
FIG_HEIGHT <- 6.6

STATE_LEVELS <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

STATE_LABELS <- c(`Homeostatic Microglia` = "Homeostatic\nMicroglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling\nMicroglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory\nMicroglia", `Ccl12+ acute Microglia` = "Ccl12+ acute\nMicroglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory\nMicroglia", `Proliferating Microglia` = "Proliferating\nMicroglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ AP\nMacrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ tissue-repair\nMacrophage", 
    `Acute-response Macrophage` = "Acute-response\nMacrophage")

install_if_missing <- require_package

for (p in c("ggplot2", "dplyr", "patchwork", "svglite", "scales")) {
  install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(svglite)
  library(scales)
})

write_csv <- write_result_csv

fdr_symbol <- function(x) {
  out <- rep("ns", length(x))
  out[!is.na(x) & x < 0.05] <- "*"
  out[!is.na(x) & x < 0.01] <- "**"
  out[!is.na(x) & x < 0.001] <- "***"
  out[!is.na(x) & x < 0.0001] <- "****"
  out
}

safe_log2_or <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- log2(x[ok])
  out[!is.na(x) & x == 0] <- -Inf
  out[is.infinite(x) & x > 0] <- Inf
  out
}

format_log2or <- function(x) {
  out <- rep("NA", length(x))
  out[is.finite(x)] <- sprintf("%.2f", x[is.finite(x)])
  out[is.infinite(x) & x > 0] <- "Inf"
  out[is.infinite(x) & x < 0] <- "-Inf"
  out
}

hyper_module_set <- function(module, gene_set, set_name, gene_module, universe_symbols) {
    module_genes <- unique(as.character(gene_module$gene_symbol_restored[gene_module$module_color == 
        module]))
    module_genes <- intersect(module_genes, universe_symbols)
    gene_set <- intersect(unique(as.character(gene_set)), universe_symbols)
    universe_n <- length(universe_symbols)
    module_n <- length(module_genes)
    set_n <- length(gene_set)
    a <- length(intersect(module_genes, gene_set))
    b <- module_n - a
    c <- set_n - a
    d <- universe_n - a - b - c
    if (any(c(a, b, c, d) < 0)) {
        stop("Invalid Fisher table for module=", module, ", state=", set_name)
    }
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "two.sided")
    data.frame(module_color = module, state = set_name, universe_n = universe_n, module_n = module_n, 
        state_marker_n_in_universe = set_n, overlap_n = a, odds_ratio = unname(ft$estimate), P.Value = ft$p.value, 
        overlap_genes = paste(sort(intersect(module_genes, gene_set)), collapse = "/"), stringsAsFactors = FALSE)
}

if (!file.exists(STATE_MARKER_FILE)) {
 stop("Step18 state-marker file not found:\n", STATE_MARKER_FILE)
}
if (!file.exists(STEP50_RDS)) {
 stop("Step50 analysis object not found:\n", STEP50_RDS)
}

state_markers <- read.csv(STATE_MARKER_FILE, check.names = FALSE, stringsAsFactors = FALSE)

obj50 <- readRDS(STEP50_RDS)

required_marker_cols <- c("cluster", "gene")
miss <- setdiff(required_marker_cols, names(state_markers))
if (length(miss) > 0L) {
  stop("Step18 state-marker file lacks: ", paste(miss, collapse = ", "))
}

state_markers$cluster <- as.character(state_markers$cluster)
state_markers$gene <- as.character(state_markers$gene)

missing_states <- setdiff(STATE_LEVELS, unique(state_markers$cluster))
if (length(missing_states) > 0L) {
    stop("Step18 marker file is missing state(s): ", paste(missing_states, collapse = "; "))
}

state_markers <- state_markers[state_markers$cluster %in% STATE_LEVELS & !is.na(state_markers$gene) & 
    nzchar(state_markers$gene), , drop = FALSE]

state_markers <- state_markers[!duplicated(state_markers[, c("cluster", "gene")]), , drop = FALSE]

required_top <- c("frozen_network", "module_trait", "key_module")
miss <- setdiff(required_top, names(obj50))
if (length(miss) > 0L) {
  stop("Step50 object lacks: ", paste(miss, collapse = ", "))
}

module_order <- as.character(obj50$frozen_network$non_grey_modules)
gene_module <- obj50$frozen_network$gene_module_assignment

if (length(module_order) != EXPECTED_NON_GREY_MODULES) {
  stop("Unexpected non-grey module count.")
}

required_gm <- c("gene_symbol_restored", "module_color")
miss <- setdiff(required_gm, names(gene_module))
if (length(miss) > 0L) {
 stop("WGCNA gene assignment lacks: ", paste(miss, collapse = ", "))
}

selected_module <- as.character(obj50$key_module$selected_module)

if (!identical(selected_module, EXPECTED_SELECTED_MODULE)) {
    stop("selected module mismatch: observed ", selected_module, "; expected ", EXPECTED_SELECTED_MODULE)
}

symbol_module_split <- split(as.character(gene_module$module_color), as.character(gene_module$gene_symbol_restored))

ambiguous_symbols <- names(Filter(function(z) {
  length(unique(z[!is.na(z)])) > 1L
}, symbol_module_split))

ambiguous_symbols <- ambiguous_symbols[!is.na(ambiguous_symbols) & nzchar(ambiguous_symbols)]

if (length(ambiguous_symbols) > 0L) {
 stop("Some restored SYMBOLs map to multiple WGCNA modules.")
}

universe_symbols <- unique(as.character(gene_module$gene_symbol_restored))
universe_symbols <- universe_symbols[!is.na(universe_symbols) & nzchar(universe_symbols)]

state_gene_sets <- lapply(STATE_LEVELS, function(st) {
  unique(state_markers$gene[state_markers$cluster == st])
})
names(state_gene_sets) <- STATE_LEVELS

state_gene_sets_in_universe <- lapply(state_gene_sets, function(g) intersect(g, universe_symbols))

state_set_qc <- data.frame(state = STATE_LEVELS, marker_genes_frozen = vapply(state_gene_sets, length, 
    integer(1)), marker_genes_in_WGCNA_universe = vapply(state_gene_sets_in_universe, length, integer(1)), 
    stringsAsFactors = FALSE)

state_set_qc$fraction_in_WGCNA_universe <- state_set_qc$marker_genes_in_WGCNA_universe/state_set_qc$marker_genes_frozen

if (any(state_set_qc$marker_genes_in_WGCNA_universe < 5L)) {
    bad <- state_set_qc$state[state_set_qc$marker_genes_in_WGCNA_universe < 5L]
    stop("Too few WGCNA-universe markers for: ", paste(bad, collapse = "; "))
}

write_csv(state_set_qc, paste0("01_", PREFIX, "_state_signature_QC.csv"))

state_enrichment_list <- list()

for (st in STATE_LEVELS) {

  tmp <- lapply(module_order, function(module) {
    hyper_module_set(
      module = module,
      gene_set = state_gene_sets_in_universe[[st]],
      set_name = st,
      gene_module = gene_module,
      universe_symbols = universe_symbols
    )
  })

  tmp <- do.call(rbind, tmp)

  tmp$adj.P.Val <- p.adjust(tmp$P.Value, method = "BH")

  state_enrichment_list[[st]] <- tmp
}

state_enrichment <- do.call(rbind, state_enrichment_list)
rownames(state_enrichment) <- NULL

if (nrow(state_enrichment) != EXPECTED_STATES * EXPECTED_NON_GREY_MODULES) {
  stop("Unexpected state-enrichment row count.")
}

state_enrichment$log2_OR <- safe_log2_or(state_enrichment$odds_ratio)
state_enrichment$significance <- fdr_symbol(state_enrichment$adj.P.Val)

state_enrichment$enrichment_direction <- ifelse(state_enrichment$adj.P.Val < 0.05 & state_enrichment$odds_ratio > 
    1, "Significant enrichment", ifelse(state_enrichment$adj.P.Val < 0.05 & state_enrichment$odds_ratio < 
    1, "Significant depletion", "Not significant"))

state_enrichment$cell_label <- paste0(format_log2or(state_enrichment$log2_OR), "\n(", state_enrichment$significance, 
    ")")

state_enrichment$state_order <- match(state_enrichment$state, STATE_LEVELS)
state_enrichment$module_factor <- factor(state_enrichment$module_color, levels = rev(module_order))
state_enrichment$x_numeric <- state_enrichment$state_order

required_trait_fields <- c("primary", "transition_vs_early", "extended_vs_early", "d28_vs_d14")

miss <- setdiff(required_trait_fields, names(obj50$module_trait))
if (length(miss) > 0L) {
  stop("Step50 module_trait lacks: ", paste(miss, collapse = ", "))
}

contrast_spec <- list(list(id = "EstablishedLate_vs_Early", label = "Established late\nvs Early", order = 1L, 
    data = obj50$module_trait$primary), list(id = "Transition_7dpi_vs_Early", label = "7 dpi\nvs Early", 
    order = 2L, data = obj50$module_trait$transition_vs_early), list(id = "ExtendedLate_42dpi_vs_Early", 
    label = "42 dpi\nvs Early", order = 3L, data = obj50$module_trait$extended_vs_early), list(id = "28dpi_vs_14dpi", 
    label = "28 dpi\nvs 14 dpi", order = 4L, data = obj50$module_trait$d28_vs_d14))

required_trait_cols <- c("comparison", "reference", "target", "module_color", "estimate_target_minus_reference", 
    "P.Value", "adj.P.Val")

temporal_list <- vector("list", length(contrast_spec))

for (i in seq_along(contrast_spec)) {
    sp <- contrast_spec[[i]]
    z <- sp$data
    miss <- setdiff(required_trait_cols, names(z))
    if (length(miss) > 0L) {
        stop("contrast ", sp$id, " lacks: ", paste(miss, collapse = ", "))
    }
    if (nrow(z) != EXPECTED_NON_GREY_MODULES) {
        stop("Unexpected module row count for ", sp$id)
    }
    z$contrast_id <- sp$id
    z$contrast_label <- sp$label
    z$contrast_order <- sp$order
    temporal_list[[i]] <- z
}

temporal_df <- do.call(rbind, temporal_list)
rownames(temporal_df) <- NULL

temporal_df$estimate <- temporal_df$estimate_target_minus_reference
temporal_df$significance <- fdr_symbol(temporal_df$adj.P.Val)
temporal_df$cell_label <- paste0(sprintf("%.2f", temporal_df$estimate), "\n(", temporal_df$significance, 
    ")")
temporal_df$x_numeric <- temporal_df$contrast_order
temporal_df$module_factor <- factor(temporal_df$module_color, levels = rev(module_order))

module_annotation <- data.frame(module_color = module_order, stringsAsFactors = FALSE)

module_annotation$module_factor <- factor(module_annotation$module_color, levels = rev(module_order))
module_annotation$module_label <- paste0("ME", module_annotation$module_color)
module_annotation$label_face <- ifelse(module_annotation$module_color == selected_module, "bold", "plain")
module_annotation$x <- 1

finite_or <- state_enrichment$log2_OR[is.finite(state_enrichment$log2_OR)]

if (length(finite_or) == 0L) {
  stop("No finite state-signature log2 odds ratios.")
}

left_limit <- max(1, min(6, ceiling(max(abs(finite_or), na.rm = TRUE) * 2)/2))

state_enrichment$left_fill_value <- NA_real_

left_sig <- (!is.na(state_enrichment$adj.P.Val) & state_enrichment$adj.P.Val < 0.05 & !is.na(state_enrichment$odds_ratio))

state_enrichment$left_fill_value[left_sig] <- state_enrichment$log2_OR[left_sig]

state_enrichment$left_fill_value[is.infinite(state_enrichment$left_fill_value) & state_enrichment$left_fill_value > 
    0] <- left_limit

state_enrichment$left_fill_value[is.infinite(state_enrichment$left_fill_value) & state_enrichment$left_fill_value < 
    0] <- -left_limit

state_enrichment$left_fill_value <- pmax(-left_limit, pmin(state_enrichment$left_fill_value, left_limit))

right_limit <- ceiling(max(abs(temporal_df$estimate), na.rm = TRUE) * 100)/100

if (!is.finite(right_limit) || right_limit <= 0) {
  stop("Invalid temporal estimate range.")
}

left_highlight <- state_enrichment[state_enrichment$module_color == selected_module & state_enrichment$state == 
    TARGET_STATE_RAW, , drop = FALSE]

right_highlight <- temporal_df[temporal_df$module_color == selected_module & temporal_df$contrast_id == 
    "EstablishedLate_vs_Early", , drop = FALSE]

if (nrow(left_highlight) != 1L) {
  stop("Could not uniquely identify brown × Clec7a state cell.")
}

if (nrow(right_highlight) != 1L) {
  stop("Could not uniquely identify brown × Primary temporal cell.")
}

if (abs(right_highlight$estimate - 0.25446) > 0.01) {
 stop("brown Primary estimate mismatch.")
}

left_labels <- unname(STATE_LABELS[STATE_LEVELS])

p_left <- ggplot(
  state_enrichment,
  aes(x = x_numeric, y = module_factor)
) +

  geom_tile(
    width = 0.94,
    height = 0.92,
    fill = "#F0F0F0",
    color = "white",
    linewidth = 0.70
  ) +

  geom_tile(
    aes(fill = left_fill_value),
    width = 0.94,
    height = 0.92,
    color = NA
  ) +

  geom_text(
    aes(label = cell_label),
    family = "Arial",
    size = 2.35,
    lineheight = 0.88,
    color = "black"
  ) +

  geom_tile(
    data = left_highlight,
    aes(x = x_numeric, y = module_factor),
    inherit.aes = FALSE,
    width = 0.94,
    height = 0.92,
    fill = NA,
    color = "black",
    linewidth = 1.10,
    linetype = "dashed"
  ) +

  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-left_limit, left_limit),
    oob = scales::squish,
    na.value = "transparent",
    name = expression(log[2] * "(OR)")
  ) +
  scale_x_continuous(
    breaks = seq_along(STATE_LEVELS),
    labels = left_labels,
    limits = c(0.5, EXPECTED_STATES + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "scRNA state-signature enrichment",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(
      family = "Arial", size = 10.5, face = "bold",
      hjust = 0.5, color = "black", margin = margin(b = 6)
    ),
    axis.text.x = element_text(
      family = "Arial", size = 7.3, color = "black",
      angle = 60, hjust = 1, vjust = 1, lineheight = 0.90
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(family = "Arial", size = 8.5),
    legend.text = element_text(family = "Arial", size = 7.7),
    plot.margin = margin(t = 7, r = 3, b = 3, l = 3, unit = "pt")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(35, "mm"),
      barheight = grid::unit(3.5, "mm"),
      frame.colour = NA
    )
  )

p_center <- ggplot(module_annotation, aes(x = x, y = module_factor)) + geom_hline(yintercept = seq(0.5, 
    EXPECTED_NON_GREY_MODULES + 0.5, by = 1), color = "#777777", linewidth = 0.35) + geom_point(aes(color = module_color), 
    size = 5.1, shape = 16, show.legend = FALSE) + scale_color_identity() + geom_text(aes(x = 1.22, label = module_label, 
    fontface = label_face), family = "Arial", size = 3.65, hjust = 0, color = "black") + scale_x_continuous(limits = c(0.2, 
    2.35), expand = c(0, 0)) + scale_y_discrete(drop = FALSE) + labs(title = "WGCNA modules", x = NULL, 
    y = NULL) + theme_void(base_family = "Arial") + theme(plot.title = element_text(family = "Arial", 
    size = 10.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 6)), plot.margin = margin(t = 7, 
    r = 4, b = 3, l = 4, unit = "pt"))

right_labels <- vapply(contrast_spec, function(x) x$label, character(1))

p_right <- ggplot(temporal_df, aes(x = x_numeric, y = module_factor, fill = estimate)) + geom_tile(width = 0.94, 
    height = 0.92, color = "white", linewidth = 0.65) + geom_text(aes(label = cell_label), family = "Arial", 
    size = 2.8, lineheight = 0.88, color = "black") + geom_tile(data = right_highlight, aes(x = x_numeric, 
    y = module_factor), inherit.aes = FALSE, width = 0.94, height = 0.92, fill = NA, color = "black", 
    linewidth = 1, linetype = "dashed") + scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", 
    midpoint = 0, limits = c(-right_limit, right_limit), oob = scales::squish, name = "ME estimate") + 
    scale_x_continuous(breaks = seq_along(right_labels), labels = right_labels, limits = c(0.5, length(right_labels) + 
        0.5), expand = c(0, 0)) + scale_y_discrete(drop = FALSE) + labs(title = "Temporal association", 
    x = NULL, y = NULL) + theme_minimal(base_family = "Arial", base_size = 10) + theme(plot.title = element_text(family = "Arial", 
    size = 10.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 6)), axis.text.x = element_text(family = "Arial", 
    size = 8, color = "black", lineheight = 0.92, margin = margin(t = 5)), axis.text.y = element_blank(), 
    axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.title = element_text(family = "Arial", 
        size = 8.5), legend.text = element_text(family = "Arial", size = 7.7), plot.margin = margin(t = 7, 
        r = 3, b = 3, l = 3, unit = "pt")) + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
    barwidth = grid::unit(28, "mm"), barheight = grid::unit(3.5, "mm"), frame.colour = NA))

p_combined <- p_left + p_center + p_right + patchwork::plot_layout(widths = c(5.4, 2.15, 2.7))

ggsave(file.path(OUT_DIR, paste0(PREFIX, ".svg")), p_combined, device = svglite::svglite, width = FIG_WIDTH, 
    height = FIG_HEIGHT, units = "in", bg = "white")

ggsave(file.path(OUT_DIR, paste0(PREFIX, ".pdf")), p_combined, device = grDevices::cairo_pdf, width = FIG_WIDTH, 
    height = FIG_HEIGHT, units = "in", family = "Arial", bg = "white")

ggsave(file.path(OUT_DIR, paste0(PREFIX, ".png")), p_combined, width = FIG_WIDTH, height = FIG_HEIGHT, 
    units = "in", dpi = 600, bg = "white")

ggsave(file.path(OUT_DIR, paste0(PREFIX, "_LEFT_scRNA_state_enrichment.svg")), p_left, device = svglite::svglite, 
    width = 7.2, height = 6.2, units = "in", bg = "white")

ggsave(file.path(OUT_DIR, paste0(PREFIX, "_RIGHT_temporal_association.svg")), p_right, device = svglite::svglite, 
    width = 4, height = 6.2, units = "in", bg = "white")

write_csv(state_enrichment, paste0("02_", PREFIX, "_scRNA_state_WGCNA_enrichment.csv"))

write_csv(temporal_df[, c("comparison", "reference", "target", "module_color", "estimate_target_minus_reference", 
    "P.Value", "adj.P.Val", "significance", "contrast_id", "contrast_label", "contrast_order"), drop = FALSE], 
    paste0("03_", PREFIX, "_frozen_temporal_association.csv"))

dual_target_summary <- data.frame(evidence = c("scRNA_state_signature_enrichment", "frozen_temporal_association"), 
    module_color = selected_module, target = c("Clec7a+ lipid-handling Microglia", "Established late vs Early"), 
    effect_metric = c("log2 odds ratio", "batch-adjusted ME estimate"), effect_value = c(left_highlight$log2_OR, 
        right_highlight$estimate), raw_P = c(left_highlight$P.Value, right_highlight$P.Value), FDR = c(left_highlight$adj.P.Val, 
        right_highlight$adj.P.Val), stringsAsFactors = FALSE)

write_csv(dual_target_summary, paste0("04_", PREFIX, "_brown_dual_association_summary.csv"))

qc <- data.frame(metric = c("Step18_state_marker_source", "Step50_source", "frozen_scRNA_states_n", "frozen_WGCNA_modules_n", 
    "selected_key_module", "left_analysis", "left_universe", "left_FDR_adjustment", "right_analysis", 
    "v03_left_color_rule", "v03_right_palette", "scRNA_marker_analysis_rerun", "WGCNA_rerun", "temporal_limma_rerun", 
    "key_module_reselected", "v05_change"), value = c(STATE_MARKER_FILE, STEP50_RDS, length(STATE_LEVELS), 
    length(module_order), selected_module, "Cross-omics two-sided Fisher test of state markers within WGCNA modules", 
    "WGCNA-tested restored-SYMBOL universe", "BH across 11 modules within each scRNA state", "Step50 batch-adjusted module eigengene contrasts", 
    "FDR<0.05 & OR>1 red; FDR<0.05 & OR<1 blue; all non-significant cells neutral gray", "Stronger blue-white-red diverging palette", 
    "NO", "NO", "NO", "NO", "layout only: center WGCNA-module panel widened; all statistics unchanged"), 
    stringsAsFactors = FALSE)

write_csv(qc, paste0("05_", PREFIX, "_plot_QC.csv"))

n_left_enriched <- sum(state_enrichment$enrichment_direction == "Significant enrichment")
n_left_depleted <- sum(state_enrichment$enrichment_direction == "Significant depletion")
n_left_ns <- sum(state_enrichment$enrichment_direction == "Not significant")

writeLines("Numerical results are provided in the accompanying output tables.", 
    file.path(OUT_DIR, paste0("06_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

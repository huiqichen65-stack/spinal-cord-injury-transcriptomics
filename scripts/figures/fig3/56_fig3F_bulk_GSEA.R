source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

GSEA_FILE <- file.path(ROOT, "results", "48_Bulk_functional_enrichment", "06_48_Bulk_functional_enrichment_v04_GSEA_GO_BP.csv")

STEP48_RDS <- file.path(ROOT, "results", "48_Bulk_functional_enrichment", "10_48_Bulk_functional_enrichment_v04_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "56_Figure3F_Bulk_GO_BP_GSEA_ridgeline")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "56_Figure3F_Bulk_GO_BP_GSEA_ridgeline_v02"

FDR_CUTOFF <- 0.05

FIG_WIDTH <- 8.2
FIG_HEIGHT <- 7.6

MANUAL_POSITIVE <- c("cellular response to lipid", "lipid catabolic process", "phospholipid efflux", 
    "cellular response to lipoprotein particle stimulus", "lysosome organization", "positive regulation of endocytosis", 
    "phagocytosis", "microglial cell activation", "myeloid leukocyte activation")

MANUAL_NEGATIVE <- c("cholesterol biosynthetic process via lathosterol", "cholesterol biosynthetic process via desmosterol", 
    "zymosterol biosynthetic process", "cholesterol biosynthetic process", "sterol biosynthetic process", 
    "isoprenoid biosynthetic process")

MANUAL_PATHWAYS <- c(MANUAL_POSITIVE, MANUAL_NEGATIVE)

install_if_missing <- require_package

for (p in c("ggplot2", "dplyr", "ggridges", "svglite", "scales")) {
    install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggridges)
  library(svglite)
  library(scales)
})

write_csv <- write_result_csv

wrap_label <- function(x, width = 42L) {
    vapply(x, function(z) {
        paste(strwrap(z, width = width), collapse = "\n")
    }, character(1))
}

if (!file.exists(GSEA_FILE)) {
    stop("Step48 GO BP GSEA CSV not found:\n", GSEA_FILE)
}

if (!file.exists(STEP48_RDS)) {
    stop("Step48 v04 analysis object not found:\n", STEP48_RDS)
}

gsea_tab <- read.csv(GSEA_FILE, check.names = FALSE, stringsAsFactors = FALSE)

obj48 <- readRDS(STEP48_RDS)

required_gsea_cols <- c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust", "direction")

missing_gsea_cols <- setdiff(required_gsea_cols, names(gsea_tab))

if (length(missing_gsea_cols) > 0L) {
    stop("GSEA CSV lacks required column(s): ", paste(missing_gsea_cols, collapse = ", "))
}

required_obj_fields <- c("input_primary", "mapping_detail", "gene_sets_input", "GSEA")

missing_obj_fields <- setdiff(required_obj_fields, names(obj48))

if (length(missing_obj_fields) > 0L) {
    stop("Step48 v04 object lacks required field(s): ", paste(missing_obj_fields, collapse = ", "))
}

if (!"GO_BP" %in% names(obj48$GSEA)) {
    stop("Step48 v04 object lacks GSEA$GO_BP.")
}

idx <- match(MANUAL_PATHWAYS, gsea_tab$Description)

if (any(is.na(idx))) {
    stop("The following manually selected pathways are absent from GSEA_GO_BP:\n", paste(MANUAL_PATHWAYS[is.na(idx)], 
        collapse = "\n"))
}

selected <- gsea_tab[idx, , drop = FALSE]

if (anyDuplicated(selected$Description)) {
    stop("Duplicated manually selected pathway Description found in GSEA table.")
}

if (any(is.na(selected$p.adjust) | selected$p.adjust >= FDR_CUTOFF)) {
    bad <- selected[is.na(selected$p.adjust) | selected$p.adjust >= FDR_CUTOFF, c("ID", "Description", 
        "NES", "p.adjust"), drop = FALSE]
    print(bad)
    stop("At least one manually selected pathway is not significant at BH FDR < 0.05.")
}

positive_ok <- selected$NES[selected$Description %in% MANUAL_POSITIVE] > 0

negative_ok <- selected$NES[selected$Description %in% MANUAL_NEGATIVE] < 0

if (!all(positive_ok)) {
    stop("At least one manually selected positive pathway does not have NES > 0.")
}

if (!all(negative_ok)) {
    stop("At least one manually selected negative pathway does not have NES < 0.")
}

selected$minus_log10_FDR <- -log10(pmax(selected$p.adjust, .Machine$double.xmin))

gsea_go_obj <- obj48$GSEA$GO_BP

slot_names <- methods::slotNames(gsea_go_obj)

if (!"geneSets" %in% slot_names) {
    stop("GO BP gseaResult does not expose a geneSets slot.\n", "This script will not substitute leading-edge genes for complete pathway membership.")
}

gene_sets <- methods::slot(gsea_go_obj, "geneSets")

if (!is.list(gene_sets) || length(gene_sets) == 0L) {
    stop("GO BP gseaResult geneSets slot is empty.")
}

missing_ids <- setdiff(selected$ID, names(gene_sets))

if (length(missing_ids) > 0L) {
    stop("Selected GO IDs absent from complete geneSets membership:\n", paste(missing_ids, collapse = "\n"))
}

primary <- obj48$input_primary
mapping_detail <- obj48$mapping_detail

required_primary_cols <- c("gene_id_clean", "t", "logFC")

missing_primary_cols <- setdiff(required_primary_cols, names(primary))

if (length(missing_primary_cols) > 0L) {
    stop("Step48 input_primary lacks required column(s): ", paste(missing_primary_cols, collapse = ", "))
}

required_map_cols <- c("gene_id_clean", "ENTREZID")

missing_map_cols <- setdiff(required_map_cols, names(mapping_detail))

if (length(missing_map_cols) > 0L) {
    stop("Step48 mapping_detail lacks required column(s): ", paste(missing_map_cols, collapse = ", "))
}

primary_small <- unique(primary[, intersect(c("gene_id_input", "gene_symbol", "gene_symbol_restored", 
    "gene_id_clean", "t", "logFC"), names(primary)), drop = FALSE])

map_small <- unique(mapping_detail[, c("gene_id_clean", "ENTREZID"), drop = FALSE])

entrez_stats <- merge(primary_small, map_small, by = "gene_id_clean", all = FALSE, sort = FALSE)

entrez_stats <- entrez_stats[!is.na(entrez_stats$ENTREZID) & nzchar(entrez_stats$ENTREZID) & is.finite(entrez_stats$t) & 
    is.finite(entrez_stats$logFC), , drop = FALSE]

if (nrow(entrez_stats) == 0L) {
    stop("No mapped Entrez-level Primary statistics available.")
}

tie_id <- if ("gene_id_input" %in% names(entrez_stats)) {
    as.character(entrez_stats$gene_id_input)
} else {
    as.character(entrez_stats$gene_id_clean)
}

entrez_stats$.tie_id <- tie_id

entrez_stats <- entrez_stats[order(entrez_stats$ENTREZID, -abs(entrez_stats$t), entrez_stats$.tie_id), 
    , drop = FALSE]

entrez_stats <- entrez_stats[!duplicated(entrez_stats$ENTREZID), , drop = FALSE]

if (anyDuplicated(entrez_stats$ENTREZID)) {
    stop("Failed to obtain unique Entrez-level Primary statistics.")
}

frozen_gene_list <- obj48$gene_sets_input$GSEA_geneList_moderated_t

if (is.null(frozen_gene_list) || is.null(names(frozen_gene_list))) {
    stop("Step48 moderated-t GSEA gene list unavailable.")
}

if (length(frozen_gene_list) != nrow(entrez_stats)) {
    warning("Entrez-stat row count (", nrow(entrez_stats), ") differs from GSEA gene-list length (", 
        length(frozen_gene_list), "). The script will require exact Entrez-ID agreement below.")
}

if (!setequal(names(frozen_gene_list), entrez_stats$ENTREZID)) {
    stop("Reconstructed Entrez universe does not match Step48 GSEA ranked universe.")
}

ridge_list <- vector("list", nrow(selected))

membership_audit <- vector("list", nrow(selected))

for (i in seq_len(nrow(selected))) {
    pathway_id <- selected$ID[i]
    pathway_name <- selected$Description[i]
    members <- unique(as.character(gene_sets[[pathway_id]]))
    members <- members[!is.na(members) & nzchar(members)]
    z <- entrez_stats[entrez_stats$ENTREZID %in% members, , drop = FALSE]
    if (nrow(z) < 5L) {
        stop("Too few mapped genes for ridgeline density in pathway: ", pathway_name, " (n=", nrow(z), 
            ")")
    }
    z$ID <- pathway_id
    z$Description <- pathway_name
    z$NES <- selected$NES[i]
    z$pvalue_GSEA <- selected$pvalue[i]
    z$FDR_GSEA <- selected$p.adjust[i]
    z$minus_log10_FDR <- selected$minus_log10_FDR[i]
    z$GSEA_direction <- ifelse(selected$NES[i] > 0, "Positive NES", "Negative NES")
    ridge_list[[i]] <- z
    membership_audit[[i]] <- data.frame(ID = pathway_id, Description = pathway_name, frozen_setSize = selected$setSize[i], 
        complete_geneSet_members_n = length(members), mapped_logFC_members_n = nrow(z), NES = selected$NES[i], 
        pvalue = selected$pvalue[i], p.adjust = selected$p.adjust[i], stringsAsFactors = FALSE)
}

ridge_df <- dplyr::bind_rows(ridge_list)

membership_qc <- dplyr::bind_rows(membership_audit)

selected_order <- selected %>% dplyr::arrange(dplyr::desc(NES))

selected_order$Description_plot <- wrap_label(selected_order$Description, width = 42L)

plot_levels <- rev(selected_order$Description)

ridge_df$Description_factor <- factor(ridge_df$Description, levels = plot_levels)

dot_df <- selected_order
dot_df$Description_factor <- factor(dot_df$Description, levels = plot_levels)

label_lookup <- setNames(selected_order$Description_plot, selected_order$Description)

x_min_data <- min(ridge_df$logFC, na.rm = TRUE)

x_max_data <- max(ridge_df$logFC, na.rm = TRUE)

x_range <- x_max_data - x_min_data

if (!is.finite(x_range) || x_range <= 0) {
    stop("Invalid logFC range for ridgeline plot.")
}

dot_x <- x_min_data - 0.075 * x_range

plot_x_min <- dot_x - 0.05 * x_range
plot_x_max <- x_max_data + 0.04 * x_range

dot_df$dot_x <- dot_x

p <- ggplot(
  ridge_df,
  aes(
    x = logFC,
    y = Description_factor,
    group = Description_factor,
    fill = minus_log10_FDR
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#777777"
  ) +

  ggridges::geom_density_ridges(
    scale = 0.84,
    rel_min_height = 0,
    alpha = 0.90,
    color = "#4A4A4A",
    linewidth = 0.35,
    bandwidth = NULL,
    from = x_min_data,
    to = x_max_data
  ) +

  geom_point(
    data = dot_df,
    aes(
      x = dot_x,
      y = Description_factor,
      color = NES
    ),
    inherit.aes = FALSE,
    size = 3.7,
    alpha = 1
  ) +

  scale_fill_gradient(
    low = "#E2DCE5",
    high = "#E7641B",
    name = expression(-log[10] * "(FDR)")
  ) +

  scale_color_gradient2(
    low = "#3E7FB1",
    mid = "#F2F2F2",
    high = "#D73027",
    midpoint = 0,
    name = "NES"
  ) +

  scale_y_discrete(
    labels = label_lookup
  ) +

  scale_x_continuous(
    breaks = scales::pretty_breaks(
      n = 6
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  coord_cartesian(
    xlim = c(
      plot_x_min,
      plot_x_max
    ),
    clip = "off"
  ) +

  labs(
    x = expression(log[2] * "(Fold Change)"),
    y = NULL
  ) +

  theme_classic(
    base_family = "Arial",
    base_size = 10
  ) +

  theme(
    plot.title = element_blank(),
    plot.tag = element_blank(),

    axis.title.x = element_text(
      family = "Arial",
      size = 10,
      face = "plain",
      color = "black",
      margin = margin(
        t = 6
      )
    ),

    axis.text.x = element_text(
      family = "Arial",
      size = 8.5,
      color = "black"
    ),

    axis.text.y = element_text(
      family = "Arial",
      size = 8.6,
      color = "black",
      lineheight = 0.92,
      hjust = 1
    ),

    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),

    axis.line.x = element_line(
      linewidth = 0.55,
      color = "black"
    ),

    axis.ticks.x = element_line(
      linewidth = 0.45,
      color = "black"
    ),

    panel.grid = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      family = "Arial",
      size = 9,
      face = "plain"
    ),

    legend.text = element_text(
      family = "Arial",
      size = 8
    ),

    legend.key.height = grid::unit(
      4.5,
      "mm"
    ),

    plot.margin = margin(
      t = 8,
      r = 8,
      b = 6,
      l = 8,
      unit = "pt"
    )
  ) +

  guides(
    fill = guide_colorbar(
      order = 1,
      barheight = grid::unit(
        25,
        "mm"
      )
    ),
    color = guide_colorbar(
      order = 2,
      barheight = grid::unit(
        25,
        "mm"
      )
    )
  )

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))

PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

ggsave(filename = SVG_FILE, plot = p, device = svglite::svglite, width = FIG_WIDTH, height = FIG_HEIGHT, 
    units = "in", bg = "white")

ggsave(filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = FIG_WIDTH, height = FIG_HEIGHT, 
    units = "in", family = "Arial", bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", dpi = 600, 
    bg = "white")

selected_export <- selected_order[, intersect(c("direction", "ID", "Description", "setSize", "enrichmentScore", 
    "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment_SYMBOL"), names(selected_order)), 
    drop = FALSE]

write_csv(selected_export, paste0("01_", PREFIX, "_selected_GSEA_pathways.csv"))

write_csv(membership_qc, paste0("02_", PREFIX, "_geneSet_membership_QC.csv"))

ridge_export_cols <- intersect(c("ID", "Description", "ENTREZID", "gene_id_clean", "gene_id_input", "gene_symbol", 
    "gene_symbol_restored", "t", "logFC", "NES", "FDR_GSEA", "minus_log10_FDR", "GSEA_direction"), names(ridge_df))

write_csv(ridge_df[, ridge_export_cols, drop = FALSE], paste0("03_", PREFIX, "_ridge_gene_level_data.csv"))

qc <- data.frame(metric = c("GSEA_CSV_source", "Step48_RDS_source", "frozen_GSEA_ranking_metric", "ridge_x_metric", 
    "ridge_membership_definition", "selected_pathways_n", "selected_positive_NES_n", "selected_negative_NES_n", 
    "selected_all_FDR_lt_0.05", "GSEA_recomputed", "leading_edge_only", "ridge_scale", "ridge_rel_min_height", 
    "v02_analysis_changes", "panel_letter_inside_plot", "internal_plot_title"), value = c(GSEA_FILE, 
    STEP48_RDS, "limma moderated t statistic", "Step46 Primary log2FC", "complete GO gene-set membership stored in Step48 gseaResult", 
    nrow(selected), sum(selected$NES > 0), sum(selected$NES < 0), all(selected$p.adjust < FDR_CUTOFF), 
    "NO", "NO", "0.84", "0", "NONE; visual-only refinement", "NO", "NO"), stringsAsFactors = FALSE)

write_csv(qc, paste0("04_", PREFIX, "_plot_QC.csv"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, paste0("05_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

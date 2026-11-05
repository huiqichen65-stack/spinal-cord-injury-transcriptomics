source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

GO_FILE <- file.path(ROOT, "results", "50_Bulk_WGCNA_module_trait_integration", "11_50_Bulk_WGCNA_module_trait_integration_v01_key_module_ORA_GO_BP.csv")

GENE_EVIDENCE_FILE <- file.path(ROOT, "results", "50_Bulk_WGCNA_module_trait_integration", "10_50_Bulk_WGCNA_module_trait_integration_v01_key_module_gene_evidence.csv")

SUMMARY_FILE <- file.path(ROOT, "results", "50_Bulk_WGCNA_module_trait_integration", "14_50_Bulk_WGCNA_module_trait_integration_v01_key_module_enrichment_summary.csv")

OUT_DIR <- file.path(ROOT, "results", "58_Figure3H_Brown_module_GO_barplot")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "58_Figure3H_Brown_module_GO_barplot_v03"

EXPECTED_MODULE <- "brown"
EXPECTED_SELECTED_TERMS_N <- 15L

FIG_WIDTH <- 9.2
FIG_HEIGHT <- 6.4

selected_terms_def <- data.frame(category = c("Lipid handling", "Lipid handling", "Lipid handling", "Lipid handling", 
    "Lipid handling", "Lysosome / uptake", "Lysosome / uptake", "Lysosome / uptake", "Myeloid state", 
    "Myeloid state", "Myeloid state", "Myeloid state", "Repair / remodeling", "Repair / remodeling", 
    "Repair / remodeling"), Description = c("cellular response to lipid", "lipid catabolic process", 
    "response to lipoprotein particle", "lipid homeostasis", "cholesterol efflux", "phagocytosis", "lysosome organization", 
    "positive regulation of endocytosis", "myeloid leukocyte activation", "macrophage activation", "microglial cell activation", 
    "myeloid leukocyte migration", "wound healing", "tissue remodeling", "extracellular matrix organization"), 
    stringsAsFactors = FALSE)

core_genes_focus <- c("Clec7a", "Trem2", "Apoe", "Lpl", "Apoc4", "Abca1")

install_if_missing <- require_package

for (p in c("ggplot2", "dplyr", "svglite", "scales", "stringr")) {
    install_if_missing(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(scales)
  library(stringr)
})

write_csv <- write_result_csv

require_columns <- function(df, cols, df_name) {
    miss <- setdiff(cols, names(df))
    if (length(miss) > 0L) {
        stop(df_name, " lacks required column(s): ", paste(miss, collapse = ", "))
    }
}

split_slash <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(character(0))
  unique(unlist(strsplit(x, "/", fixed = TRUE), use.names = FALSE))
}

for (f in c(GO_FILE, GENE_EVIDENCE_FILE, SUMMARY_FILE)) {
  if (!file.exists(f)) {
  stop("Required Step50 input not found:\n", f)
  }
}

go_df <- read.csv(GO_FILE, check.names = FALSE, stringsAsFactors = FALSE)

gene_evidence <- read.csv(GENE_EVIDENCE_FILE, check.names = FALSE, stringsAsFactors = FALSE)

summary_df <- read.csv(SUMMARY_FILE, check.names = FALSE, stringsAsFactors = FALSE)

require_columns(go_df, c("Description", "p.adjust", "Count", "geneID_SYMBOL"), "GO BP table")

require_columns(gene_evidence, c("gene_symbol_restored", "module_color", "kME_own_module"), "Key-module gene evidence table")

require_columns(summary_df, c("database", "Description", "p.adjust", "Count", "significant_BH_0.05"), 
    "Key-module enrichment summary table")

if (!all(gene_evidence$module_color == EXPECTED_MODULE)) {
    stop("Key-module gene evidence table contains module_color values other than ", EXPECTED_MODULE, 
        ".")
}

go_sig_n <- sum(is.finite(go_df$p.adjust) & go_df$p.adjust < 0.05)

if (go_sig_n < 100L) {
    stop("Unexpectedly few significant GO BP terms in the table: ", go_sig_n)
}

summary_go_sig_n <- sum(summary_df$database == "GO_BP" & summary_df$significant_BH_0.05 %in% TRUE)

if (summary_go_sig_n != go_sig_n) {
    warning("GO BP significant-term count differs between GO table (", go_sig_n, ") and enrichment summary (", 
        summary_go_sig_n, "). Proceeding with the direct GO BP table.")
}

selected_terms <- selected_terms_def %>% left_join(go_df, by = "Description")

if (nrow(selected_terms) != EXPECTED_SELECTED_TERMS_N) {
  stop("Unexpected row count after term join.")
}

missing_terms <- selected_terms$Description[is.na(selected_terms$p.adjust)]

if (length(missing_terms) > 0L) {
    stop("Some manually selected Figure 3H GO terms were not found in the GO table: ", paste(missing_terms, 
        collapse = "; "))
}

dup_terms <- selected_terms$Description[duplicated(selected_terms$Description)]
if (length(dup_terms) > 0L) {
    stop("Duplicated GO terms detected in the manual Figure 3H selection: ", paste(unique(dup_terms), 
        collapse = "; "))
}

selected_terms <- selected_terms %>% mutate(neg_log10_padj = -log10(p.adjust), Count = as.numeric(Count), 
    core_genes_present = vapply(geneID_SYMBOL, function(x) {
        hits <- intersect(split_slash(x), core_genes_focus)
        paste(hits, collapse = "/")
    }, character(1)), n_core_genes_present = vapply(geneID_SYMBOL, function(x) {
        length(intersect(split_slash(x), core_genes_focus))
    }, integer(1)))

selected_terms <- selected_terms %>% arrange(desc(Count), desc(neg_log10_padj), Description)

selected_terms$Description_factor <- factor(selected_terms$Description, levels = rev(selected_terms$Description))

selected_terms$term_rank <- seq_len(nrow(selected_terms))

selected_terms$label_x <- pmax(1.4, selected_terms$Count * 0.022)

brown_gene_evidence_focus <- gene_evidence %>% filter(gene_symbol_restored %in% core_genes_focus) %>% 
    arrange(desc(kME_own_module), gene_symbol_restored)

term_core_support_long <- lapply(seq_len(nrow(selected_terms)), function(i) {
    genes <- intersect(split_slash(selected_terms$geneID_SYMBOL[i]), core_genes_focus)
    if (length(genes) == 0L) {
        genes <- NA_character_
    }
    data.frame(category = selected_terms$category[i], Description = selected_terms$Description[i], core_gene = genes, 
        stringsAsFactors = FALSE)
})

term_core_support_long <- do.call(rbind, term_core_support_long)

max_count <- max(selected_terms$Count, na.rm = TRUE)

p <- ggplot(
  selected_terms,
  aes(
    x = Count,
    y = Description_factor,
    fill = neg_log10_padj
  )
) +
  geom_col(
    width = 0.78,
    color = NA
  ) +

  geom_text(
    aes(
      x = label_x,
      label = Description
    ),
    hjust = 0,
    family = "Arial",
    size = 3.15,
    color = "black"
  ) +

  scale_fill_gradientn(
    colors = c(
      "#FFD84D",
      "#F6A13A",
      "#E9683F",
      "#C92E36"
    ),
    name = expression(-log[10] * "(p.adjust)")
  ) +

  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, max_count * 1.08),
    breaks = pretty(c(0, max_count), n = 6)
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  labs(
    title = "Brown module GO biological process enrichment",
    x = "Number of Gene",
    y = "Pathway"
  ) +

  theme_minimal(
    base_family = "Arial",
    base_size = 11
  ) +

  theme(
    plot.title = element_text(
      family = "Arial",
      size = 11.5,
      face = "bold",
      hjust = 0.5,
      color = "black",
      margin = margin(b = 8)
    ),

    axis.title.x = element_text(
      family = "Arial",
      size = 11,
      color = "black"
    ),

    axis.title.y = element_text(
      family = "Arial",
      size = 11,
      color = "black",
      margin = margin(r = 8)
    ),

    axis.text.y = element_blank(),

    axis.text.x = element_text(
      family = "Arial",
      size = 10,
      color = "black"
    ),

    axis.ticks.y = element_blank(),

    axis.ticks.x = element_line(
      color = "black",
      linewidth = 0.35
    ),

    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),

    panel.grid.major.x = element_line(
      color = "#E8E8E8",
      linewidth = 0.35
    ),

    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),

    axis.line = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      family = "Arial",
      size = 10.5
    ),

    legend.text = element_text(
      family = "Arial",
      size = 9.5
    ),

    plot.margin = margin(
      t = 8,
      r = 22,
      b = 10,
      l = 8,
      unit = "pt"
    )
  ) +

  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(4.4, "mm"),
      barheight = grid::unit(38, "mm"),
      frame.colour = NA
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

selected_export <- selected_terms %>% transmute(term_rank, category, Description, Count, p.adjust, neg_log10_padj, 
    core_genes_present, n_core_genes_present, geneID_SYMBOL)

write_csv(selected_export, paste0("01_", PREFIX, "_selected_GO_terms.csv"))

write_csv(brown_gene_evidence_focus, paste0("02_", PREFIX, "_core_gene_evidence.csv"))

write_csv(term_core_support_long, paste0("03_", PREFIX, "_selected_GO_term_core_gene_support_long.csv"))

plot_qc <- data.frame(metric = c("GO_result_source", "gene_evidence_source", "enrichment_summary_source", 
    "selected_module", "frozen_significant_GO_BP_terms_BH_0.05", "selected_terms_n", "selection_strategy", 
    "x_axis", "bar_fill", "term_label_location", "plot_internal_title", "panel_border", "plot_internal_panel_letter", 
    "ORA_rerun", "WGCNA_rerun"), value = c(GO_FILE, GENE_EVIDENCE_FILE, SUMMARY_FILE, EXPECTED_MODULE, 
    go_sig_n, nrow(selected_terms), "Manual curation from brown-module GO BP results", "Number of Gene (Count)", 
    "-log10(p.adjust)", "Labels start inside bars and may extend beyond short bars", "YES", "YES", "NO", 
    "NO", "NO"), stringsAsFactors = FALSE)

write_csv(plot_qc, paste0("04_", PREFIX, "_plot_QC.csv"))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT_DIR, paste0("05_", PREFIX, "_REPORT.txt")), useBytes = TRUE)

for (i in seq_len(min(5L, nrow(selected_terms)))) {
}

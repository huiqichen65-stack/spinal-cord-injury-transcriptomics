options(stringsAsFactors = FALSE)

required_pkgs <- c("monocle3", "ggplot2", "dplyr", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(dplyr)
})

ROOT <- "D:/zz"

STEP26_DIR <- file.path(ROOT, "results", "26_Monocle3_keygene_pseudotime_dynamics")

CDS_FILE <- file.path(STEP26_DIR, "26_Monocle3_partition1_pseudotime_3scRNA_v01.rds")

OUT_DIR <- file.path(ROOT, "results", "37_Figure2J_keygene_pseudotime_dynamics")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "37_Figure2J_keygene_pseudotime_dynamics_v08")

if (!file.exists(CDS_FILE)) {
    stop("Cannot find fixed Step26 partition-1 CDS:\n", CDS_FILE)
}

cds <- readRDS(CDS_FILE)

if (!inherits(cds, "cell_data_set")) {
  stop("Input RDS is not a monocle3 cell_data_set.")
}

pt <- monocle3::pseudotime(cds)

if (sum(is.finite(pt)) == 0) {
 stop("No finite pseudotime values are present in the CDS.")
}

genes_requested <- c("Trem2", "Apoe", "Lpl")

gene_meta <- as.data.frame(SummarizedExperiment::rowData(cds), stringsAsFactors = FALSE)

if ("gene_short_name" %in% colnames(gene_meta)) {
    gene_symbols <- as.character(gene_meta$gene_short_name)
    gene_ids <- rownames(cds)
    idx <- match(genes_requested, gene_symbols)
    genes_present <- genes_requested[!is.na(idx)]
    ids_present <- gene_ids[idx[!is.na(idx)]]
} else {
    genes_present <- intersect(genes_requested, rownames(cds))
    ids_present <- genes_present
}

genes_missing <- setdiff(genes_requested, genes_present)

if (length(genes_missing) > 0) {
    warning("Missing requested gene(s): ", paste(genes_missing, collapse = ", "))
}

if (length(ids_present) == 0) {
  stop("None of the requested Figure 2J genes were found.")
}

state_levels_metadata <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

display_name_map <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb⁺ macrophages", 
    `Acute-response Macrophage` = "Acute-response Macrophage")

state_cols_metadata <- c(`Homeostatic Microglia` = "#4C78A8", `Clec7a+ lipid-associated Microglia` = "#D64F4F", 
    `Activated inflammatory Microglia` = "#F2B134", `Ccl12+ acute Microglia` = "#8E6BBE", `Tnf+ inflammatory Microglia` = "#E77AA8", 
    `Proliferating Microglia` = "#4FA6A2", `MHC-II+ antigen-presenting Macrophage` = "#59A14F", `Gpnmb+ tissue-repair Macrophage` = "#D5A832", 
    `Acute-response Macrophage` = "#9C755F")

state_cols_display <- setNames(unname(state_cols_metadata[state_levels_metadata]), unname(display_name_map[state_levels_metadata]))

md <- as.data.frame(SummarizedExperiment::colData(cds), stringsAsFactors = FALSE)

if (!"myeloid_state" %in% colnames(md)) {
    stop("myeloid_state is missing from colData(cds).")
}

state_raw <- as.character(md$myeloid_state)

state_display <- unname(display_name_map[state_raw])

state_display[is.na(state_display)] <- state_raw[is.na(state_display)]

SummarizedExperiment::colData(cds)$myeloid_state_plot <- factor(state_display, levels = unname(display_name_map[state_levels_metadata]))

cds_fig <- cds[ids_present, ]

if ("gene_short_name" %in% colnames(SummarizedExperiment::rowData(cds_fig))) {
    SummarizedExperiment::rowData(cds_fig)$gene_short_name <- genes_present
}

p <- monocle3::plot_genes_in_pseudotime(cds_fig, color_cells_by = "myeloid_state_plot", min_expr = 0.5, 
    ncol = 3)

p <- p + scale_color_manual(values = state_cols_display, drop = FALSE, name = NULL) + guides(color = guide_legend(nrow = 2, 
    byrow = TRUE, override.aes = list(size = 2.4, alpha = 1))) + labs(x = "pseudotime", y = "Expression") + 
    theme_classic(base_family = "Arial", base_size = 12) + theme(text = element_text(family = "Arial", 
    color = "black"), strip.background = element_rect(fill = "#E6E6E6", color = NA), strip.text = element_text(family = "Arial", 
    size = 12.5, face = "bold", color = "black", margin = margin(2.5, 0, 2.5, 0)), axis.title.x = element_text(family = "Arial", 
    size = 12.5, face = "bold"), axis.title.y = element_text(family = "Arial", size = 12.5, face = "bold"), 
    axis.text = element_text(family = "Arial", size = 10.7, color = "black"), legend.position = "bottom", 
    legend.direction = "horizontal", legend.box = "horizontal", legend.text = element_text(family = "Arial", 
        size = 10), legend.key.height = grid::unit(3.2, "mm"), legend.key.width = grid::unit(4, "mm"), 
    panel.spacing = grid::unit(4, "mm"), plot.margin = margin(6, 6, 6, 6))

W_IN <- 9.7
H_IN <- 4.7

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

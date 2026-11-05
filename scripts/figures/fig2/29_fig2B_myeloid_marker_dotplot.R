rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- PROJECT_ROOT

IN_RDS <- file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

MARKER_STATS_FILE <- file.path(ROOT, "results", "18_myeloid_annotation_cleanup", "18_myeloid_state_markers_v01.csv")

OUT_DIR <- file.path(ROOT, "results", "29_Figure2B_myeloid_marker_dotplot")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "29_Figure2B_myeloid_marker_dotplot_v07"

FIG_W <- 11.6
FIG_H <- 5.6
DPI <- 600

COLOR_MIN <- -1
COLOR_MAX <- 1

required_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "svglite", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(svglite)
  library(scales)
})

marker_map <- data.frame(myeloid_state = rep(c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", 
    "Activated inflammatory Microglia", "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", 
    "MHC-II+ antigen-presenting Macrophage", "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage"), 
    each = 2), display_state = rep(c("Homeostatic Microglia", "Clec7a+ lipid-handling Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb⁺ macrophages", "Acute-response Macrophage"), each = 2), gene = c("P2ry12", "Tmem119", "Clec7a", 
    "Lpl", "Il1b", "Cd83", "Ccl12", "Msr1", "Tnf", "Lilrb4a", "Mki67", "Top2a", "Cd74", "H2-Aa", "Gpnmb", 
    "Atp6v0d2", "S100a8", "S100a9"), stringsAsFactors = FALSE)

state_order <- unique(marker_map$myeloid_state)
display_order <- unique(marker_map$display_state)
gene_order <- marker_map$gene

if (!file.exists(IN_RDS)) {
 stop("Step20 RDS not found:\n", IN_RDS)
}

if (!file.exists(MARKER_STATS_FILE)) {
 stop("Step18 marker table not found:\n", MARKER_STATS_FILE)
}

marker_stats <- read.csv(MARKER_STATS_FILE, check.names = FALSE, stringsAsFactors = FALSE)

required_marker_cols <- c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")

missing_marker_cols <- setdiff(required_marker_cols, colnames(marker_stats))

if (length(missing_marker_cols) > 0L) {
    stop("Step18 marker table lacks required column(s): ", paste(missing_marker_cols, collapse = ", "))
}

audit_list <- vector("list", nrow(marker_map))

for (i in seq_len(nrow(marker_map))) {
    st <- marker_map$myeloid_state[i]
    g <- marker_map$gene[i]
    hit <- marker_stats[marker_stats$cluster == st & marker_stats$gene == g, , drop = FALSE]
    if (nrow(hit) != 1L) {
        stop("Expected exactly one Step18 marker-statistics row for: ", st, " / ", g, "; observed ", 
            nrow(hit))
    }
    audit_list[[i]] <- data.frame(myeloid_state = st, display_state = marker_map$display_state[i], gene = g, 
        avg_log2FC = hit$avg_log2FC[1], pct.1 = hit$pct.1[1], pct.2 = hit$pct.2[1], delta_pct = hit$pct.1[1] - 
            hit$pct.2[1], p_val_adj = hit$p_val_adj[1], stringsAsFactors = FALSE)
}

marker_audit <- do.call(rbind, audit_list)

if (any(marker_audit$avg_log2FC <= 0)) {
  stop("At least one selected marker is not positively enriched in its intended state.")
}

obj <- readRDS(IN_RDS)

DefaultAssay(obj) <- "RNA"

obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) {
    message("JoinLayers skipped: ", conditionMessage(e))
    obj
})

if (!"myeloid_state" %in% colnames(obj[[]])) {
  stop("Metadata column 'myeloid_state' was not found.")
}

rna_layers <- SeuratObject::Layers(obj[["RNA"]])

if (!"data" %in% rna_layers) {
    stop("RNA normalized 'data' layer was not found.\n", "This manuscript plotting script will NOT rerun normalization.")
}

missing_genes <- setdiff(gene_order, rownames(obj))

if (length(missing_genes) > 0L) {
    stop("Selected marker gene(s) missing from object: ", paste(missing_genes, collapse = ", "))
}

missing_states <- setdiff(state_order, unique(as.character(obj$myeloid_state)))

if (length(missing_states) > 0L) {
    stop("Expected myeloid state(s) missing from object: ", paste(missing_states, collapse = "; "))
}

dp0 <- Seurat::DotPlot(object = obj, features = gene_order, group.by = "myeloid_state", assay = "RNA", 
    scale = TRUE, scale.by = "radius", dot.scale = 6)

dot_df <- dp0$data

required_dot_cols <- c("features.plot", "id", "pct.exp", "avg.exp", "avg.exp.scaled")

missing_dot_cols <- setdiff(required_dot_cols, colnames(dot_df))

if (length(missing_dot_cols) > 0L) {
    stop("Unexpected Seurat DotPlot output; missing column(s): ", paste(missing_dot_cols, collapse = ", "))
}

state_lookup <- setNames(unique(marker_map[, c("myeloid_state", "display_state")])$display_state, unique(marker_map[, 
    c("myeloid_state", "display_state")])$myeloid_state)

dot_df$display_state <- unname(state_lookup[as.character(dot_df$id)])

if (anyNA(dot_df$display_state)) {
  stop("Failed to map one or more final myeloid states to display labels.")
}

dot_df$gene <- factor(dot_df$features.plot, levels = gene_order)

dot_df$display_state <- factor(dot_df$display_state, levels = rev(display_order))

dot_df$avg_exp_plot <- pmax(COLOR_MIN, pmin(COLOR_MAX, dot_df$avg.exp.scaled))

pair_separators <- seq(2.5, length(gene_order) - 0.5, by = 2)

p <- ggplot(dot_df, aes(x = gene, y = display_state)) + geom_vline(xintercept = pair_separators, color = "#EFEFEF", 
    linewidth = 0.3) + geom_point(aes(size = pct.exp, color = avg_exp_plot), alpha = 1) + scale_size_continuous(name = "Percent Expressed", 
    range = c(0.45, 6.2), limits = c(0, 100), breaks = c(25, 50, 75)) + scale_color_gradient2(name = "Average Expression", 
    low = "#3B6FB6", mid = "#FFF3B0", high = "#D73027", midpoint = 0, limits = c(COLOR_MIN, COLOR_MAX), 
    breaks = c(-1, 0, 1), labels = c("-1", "0", "1"), oob = scales::squish) + scale_x_discrete(expand = expansion(add = 0.45)) + 
    scale_y_discrete(expand = expansion(add = 0.35)) + labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 1, 
    title.position = "left", title.hjust = 0, direction = "horizontal", override.aes = list(color = "black")), 
    color = guide_colorbar(order = 2, title.position = "left", title.hjust = 0, direction = "horizontal", 
        barwidth = grid::unit(1.75, "cm"), barheight = grid::unit(0.22, "cm"), ticks = TRUE)) + theme_bw(base_family = "Arial", 
    base_size = 13) + theme(text = element_text(family = "Arial", color = "black"), axis.text.x = element_text(family = "Arial", 
    face = "italic", size = 14.3, angle = 50, hjust = 1, vjust = 1, color = "black"), axis.text.y = element_text(family = "Arial", 
    size = 13.4, color = "black"), axis.ticks = element_blank(), panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.border = element_rect(color = "#666666", fill = NA, linewidth = 0.55), 
    legend.position = "top", legend.box = "horizontal", legend.box.just = "left", legend.spacing.x = grid::unit(0.1, 
        "cm"), legend.title = element_text(family = "Arial", size = 12.8, color = "black"), legend.text = element_text(family = "Arial", 
        size = 12.3, color = "black"), plot.margin = margin(t = 4, r = 7, b = 4, l = 4, unit = "pt"))

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))

PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

DOT_STATS_FILE <- file.path(OUT_DIR, paste0(PREFIX, "_dotplot_statistics.csv"))

MARKER_AUDIT_FILE <- file.path(OUT_DIR, paste0(PREFIX, "_selected_marker_audit.csv"))

MARKER_MAP_FILE <- file.path(OUT_DIR, paste0(PREFIX, "_marker_map.csv"))

ggsave(filename = SVG_FILE, plot = p, device = svglite::svglite, width = FIG_W, height = FIG_H, units = "in", 
    bg = "white")

ggsave(filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = FIG_W, height = FIG_H, units = "in", 
    family = "Arial", bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = FIG_W, height = FIG_H, units = "in", dpi = DPI, bg = "white")

write.csv(dot_df[, c("display_state", "features.plot", "pct.exp", "avg.exp", "avg.exp.scaled", "avg_exp_plot")], 
    DOT_STATS_FILE, row.names = FALSE, fileEncoding = "UTF-8")

write.csv(marker_audit, MARKER_AUDIT_FILE, row.names = FALSE, fileEncoding = "UTF-8")

write.csv(marker_map, MARKER_MAP_FILE, row.names = FALSE, fileEncoding = "UTF-8")

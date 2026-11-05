source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)

required_pkgs <- c("ComplexHeatmap", "circlize", "ggplot2", "dplyr", "tidyr", "patchwork", "svglite", 
    "stringr", "grid")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(grid)
})

ROOT <- PROJECT_ROOT

STEP27_DIR <- file.path(ROOT, "results", "27_Monocle3_trajectory_gene_modules")

MODULE_RDS <- file.path(STEP27_DIR, "27_trajectory_gene_modules_3scRNA_v01.rds")

GO_FILE <- file.path(STEP27_DIR, "27_trajectory_module_GO_BP_v01.csv")

OUT_DIR <- file.path(ROOT, "results", "38_Figure2K_pseudotime_DEG_heatmap_GO")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "38_Figure2K_pseudotime_DEG_heatmap_GO_v09")

if (!file.exists(MODULE_RDS)) {
  stop("Cannot find Step27 module RDS:\n", MODULE_RDS)
}

if (!file.exists(GO_FILE)) {
  stop("Cannot find Step27 GO file:\n", GO_FILE)
}

obj <- readRDS(MODULE_RDS)

required_obj <- c("genes", "module_summary", "bin_pseudotime", "z_matrix")

if (!all(required_obj %in% names(obj))) {
    stop("Step27 RDS is missing required components: ", paste(setdiff(required_obj, names(obj)), collapse = ", "))
}

assignment <- as.data.frame(obj$genes, stringsAsFactors = FALSE)

module_summary <- as.data.frame(obj$module_summary, stringsAsFactors = FALSE)

bin_pt <- as.data.frame(obj$bin_pseudotime, stringsAsFactors = FALSE)

pt_matrix <- as.matrix(obj$z_matrix)

if (!all(c("gene", "module") %in% colnames(assignment))) {
  stop("obj$genes must contain columns: gene and module.")
}

if (is.list(assignment$gene) || is.list(assignment$module)) {
 stop("Unexpected list-column in Step27 obj$genes ('gene' or 'module'). Please inspect the RDS.")
}

module_levels <- c("Module 1", "Module 2", "Module 3", "Module 4")

module_short <- c(`Module 1` = "M1", `Module 2` = "M2", `Module 3` = "M3", `Module 4` = "M4")

module_cols <- c(`Module 1` = "#2C7FB8", `Module 2` = "#D73027", `Module 3` = "#E6B800", `Module 4` = "#31A354")

assignment_use <- assignment %>% dplyr::filter(gene %in% rownames(pt_matrix), module %in% module_levels) %>% 
    dplyr::mutate(module = factor(module, levels = module_levels))

module_sizes <- assignment_use %>% dplyr::count(module, name = "gene_size") %>% dplyr::mutate(module = factor(module, 
    levels = module_levels)) %>% dplyr::arrange(module)

print(module_sizes)

row_module <- setNames(as.character(assignment_use$module), assignment_use$gene)

genes_use <- intersect(rownames(pt_matrix), names(row_module))

pt_matrix <- pt_matrix[genes_use, , drop = FALSE]

row_split <- factor(row_module[rownames(pt_matrix)], levels = module_levels)

pt_matrix_plot <- pmax(pmin(pt_matrix, 2), -2)

if ("pseudotime_median" %in% colnames(bin_pt)) {
  pt_values <- bin_pt$pseudotime_median
} else if ("pseudotime_mean" %in% colnames(bin_pt)) {
  pt_values <- bin_pt$pseudotime_mean
} else {
  stop("bin_pseudotime must contain pseudotime_median or pseudotime_mean.")
}

if (length(pt_values) != ncol(pt_matrix_plot)) {
  stop("Pseudotime-bin number does not match z-matrix columns.")
}

profile_df <- lapply(
  module_levels,
  function(m) {

    g <- names(row_module)[row_module == m]
    g <- intersect(g, rownames(pt_matrix_plot))

    mean_z <- colMeans(
      pt_matrix_plot[g, , drop = FALSE],
      na.rm = TRUE
    )

    sm <- stats::smooth.spline(
      x = pt_values,
      y = mean_z,
      df = min(6, max(4, floor(length(pt_values) / 8)))
    )

    pred <- predict(
      sm,
      x = pt_values
    )

    data.frame(
      module = m,
      pseudotime = pred$x,
      mean_z = pred$y,
      stringsAsFactors = FALSE
    )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    module = factor(module, levels = module_levels)
  )

marker_candidates <- list(`Module 1` = c("P2ry12", "Tmem119", "Siglech", "Sall1", "Cx3cr1"), `Module 2` = c("Trem2", 
    "Apoe", "Abca1", "Apoc4", "Clec7a"), `Module 3` = c("Lpl", "Tnf", "Il1b", "Ccl12"), `Module 4` = c("Ccnb2", 
    "Top2a", "Mki67"))

marker_genes <- unique(unlist(marker_candidates))
marker_genes <- intersect(marker_genes, rownames(pt_matrix_plot))

marker_idx <- match(marker_genes, rownames(pt_matrix_plot))
marker_idx <- marker_idx[!is.na(marker_idx)]
marker_labels <- rownames(pt_matrix_plot)[marker_idx]

left_mark <- rowAnnotation(mark = anno_mark(at = marker_idx, labels = marker_labels, which = "row", side = "left", 
    labels_gp = gpar(fontfamily = "Arial", fontsize = 12.5, fontface = "italic", col = "black"), lines_gp = gpar(col = "#4D4D4D", 
        lwd = 0.6), link_width = unit(3.8, "mm"), link_height = unit(2.2, "mm"), padding = unit(0.5, 
        "mm")), width = max_text_width(marker_labels, gp = gpar(fontfamily = "Arial", fontsize = 12.5, 
    fontface = "italic")) + unit(6, "mm"))

pt_scaled <- (pt_values - min(pt_values, na.rm = TRUE)) /
  (max(pt_values, na.rm = TRUE) - min(pt_values, na.rm = TRUE))

pseudotime_col_fun <- circlize::colorRamp2(c(0, 1), c("#F3D91B", "#4B3B8F"))

top_ha <- HeatmapAnnotation(pt = pt_scaled, col = list(pt = pseudotime_col_fun), show_annotation_name = FALSE, 
    show_legend = FALSE, simple_anno_size = unit(2, "mm"))

heat_col_fun <- circlize::colorRamp2(c(-2, -1, 0, 1, 2), c("#2166AC", "#6C6F7D", "#F1C928", "#E7771C", 
    "#B2182B"))

ht <- Heatmap(pt_matrix_plot, name = "Z-Score", col = heat_col_fun, top_annotation = top_ha, left_annotation = left_mark, 
    row_split = row_split, row_title = unname(module_short[module_levels]), row_title_rot = 90, row_title_gp = gpar(fontfamily = "Arial", 
        fontsize = 13, fontface = "bold"), show_row_names = FALSE, show_column_names = FALSE, cluster_rows = TRUE, 
    clustering_method_rows = "ward.D2", cluster_row_slices = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, 
    use_raster = FALSE, row_dend_width = unit(4.5, "mm"), row_gap = unit(1.2, "mm"), border = FALSE)

heat_grob <- grid.grabExpr(draw(ht, padding = unit(c(1, 1, 1, 1), "mm")))

p_heat <- patchwork::wrap_elements(full = heat_grob)

profile_plot_list <- vector("list", length(module_levels))
names(profile_plot_list) <- module_levels

for (m in module_levels) {
    d <- profile_df %>% dplyr::filter(module == m)
    n_gene <- module_sizes$gene_size[as.character(module_sizes$module) == m]
    profile_plot_list[[m]] <- ggplot(d, aes(x = pseudotime, y = mean_z)) + geom_line(color = module_cols[[m]], 
        linewidth = 0.9) + annotate("label", x = min(d$pseudotime), y = Inf, label = paste0("Gene Size: ", 
        n_gene), family = "Arial", fontface = "italic", size = 3.55, color = module_cols[[m]], fill = "white", 
        label.size = 0, label.padding = grid::unit(0.5, "mm"), hjust = 0, vjust = 1.05) + theme_classic(base_family = "Arial", 
        base_size = 10) + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        axis.line = element_blank(), panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.45), 
        plot.margin = margin(1.5, 2, 1.5, 2))
}

profile_heights <- module_sizes$gene_size

p_profiles <- patchwork::wrap_plots(profile_plot_list[module_levels], ncol = 1, heights = profile_heights)

selected_go <- list(`Module 1` = c("myeloid cell differentiation", "leukocyte cell-cell adhesion", "response to lipopolysaccharide", 
    "activation of innate immune response", "myeloid cell homeostasis"), `Module 2` = c("cholesterol efflux", 
    "lysosome organization", "lipid storage", "lipid homeostasis", "cholesterol transport"), `Module 3` = c("oxidative phosphorylation", 
    "cellular respiration", "aerobic respiration", "generation of precursor metabolites and energy", 
    "ATP biosynthetic process"), `Module 4` = c("chromosome segregation", "mitotic nuclear division", 
    "mitotic sister chromatid segregation", "positive regulation of cell cycle", "spindle organization"))

go_all <- read.csv(GO_FILE, stringsAsFactors = FALSE, check.names = FALSE)

go_selected <- lapply(module_levels, function(m) {
    wanted <- selected_go[[m]]
    d <- go_all %>% dplyr::filter(module == m, Description %in% wanted)
    d$Description <- factor(d$Description, levels = wanted)
    d %>% dplyr::arrange(Description) %>% dplyr::mutate(module = m, neglog10_padj = -log10(p.adjust))
}) %>% dplyr::bind_rows()

write.csv(go_selected, file.path(OUT_DIR, "38_Figure2K_selected_GO_terms_v09.csv"), row.names = FALSE)

go_panel_list <- vector("list", length(module_levels))
names(go_panel_list) <- module_levels

for (m in module_levels) {

  d <- go_selected %>%
    dplyr::filter(module == m)

  if (nrow(d) == 0) {
    go_panel_list[[m]] <- ggplot() + theme_void()
    next
  }

  d <- d %>%
    dplyr::mutate(
      Description_chr = as.character(Description),
      term = factor(
        Description_chr,
        levels = rev(Description_chr)
      )
    )

  go_text_size <- switch(m,
    "Module 1" = 4.0,
    "Module 2" = 3.5,
    "Module 3" = 4.0,
    "Module 4" = 3.4
  )
  p_text <- ggplot(
    d,
    aes(x = 0, y = term, label = Description_chr)
  ) +
    geom_text(
      hjust = 0,
      family = "Arial",
      size = go_text_size,
      color = module_cols[[m]]
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
      panel.background = element_rect(
        fill = "#F4F4F4",
        color = NA
      ),
      plot.margin = margin(2, 4, 2, 3)
    )

  max_x <- max(d$neglog10_padj, na.rm = TRUE)

  p_bar <- ggplot(
    d,
    aes(x = neglog10_padj, y = term)
  ) +
    geom_col(
      fill = module_cols[[m]],
      width = 0.38
    ) +
    annotate(
      "label",
      x = max_x * 1.16,
      y = levels(d$term)[1],
      label = module_short[[m]],
      family = "Arial",
      size = 3.0,
      color = module_cols[[m]],
      fill = "white",
      label.size = 0.25,
      hjust = 1
    ) +
    scale_x_continuous(
      limits = c(0, max_x * 1.22),
      expand = c(0, 0)
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(
        color = "#333333",
        fill = NA,
        linewidth = 0.45
      ),
      plot.margin = margin(2, 2, 2, 0)
    )

  go_panel_list[[m]] <- (p_text | p_bar) +
    plot_layout(widths = c(3.50, 1.05))
}

p_go <- patchwork::wrap_plots(go_panel_list[module_levels], ncol = 1, heights = profile_heights)

strip_list <- vector("list", length(module_levels))
names(strip_list) <- module_levels

for (m in module_levels) {
    n_gene <- module_sizes$gene_size[as.character(module_sizes$module) == m]
    strip_list[[m]] <- ggplot() + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = module_cols[[m]]) + 
        annotate("text", x = 0.5, y = 0.5, label = paste0("n:", n_gene), family = "Arial", color = "white", 
            size = 3.3, angle = 90) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + 
        theme_void() + theme(plot.margin = margin(1.5, 0.5, 1.5, 0.5))
}

p_strip <- patchwork::wrap_plots(strip_list[module_levels], ncol = 1, heights = profile_heights)

z_legend <- Legend(title = "Z-Score", col_fun = heat_col_fun, at = c(-2, -1, 0, 1, 2), direction = "vertical", 
    title_gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold"), labels_gp = gpar(fontfamily = "Arial", 
        fontsize = 9.5), legend_height = unit(32, "mm"), legend_width = unit(3, "mm"))

legend_grob <- grid.grabExpr(draw(z_legend))

p_legend <- patchwork::wrap_elements(full = legend_grob)

figure2k <- (p_profiles | p_heat | p_strip | p_go | p_legend) + plot_layout(widths = c(1.1, 2.6, 0.15, 
    4.55, 0.43))

W_IN <- 14.4
H_IN <- 9.0

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = figure2k, width = W_IN, height = H_IN, units = "in", 
    device = cairo_pdf, bg = "white", limitsize = FALSE)

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = figure2k, width = W_IN, height = H_IN, units = "in", 
    device = svglite::svglite, bg = "white", limitsize = FALSE)

ggsave(filename = paste0(OUT_BASE, ".png"), plot = figure2k, width = W_IN, height = H_IN, units = "in", 
    dpi = 600, bg = "white", limitsize = FALSE)

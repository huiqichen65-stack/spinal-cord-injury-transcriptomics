options(stringsAsFactors = FALSE)
RASTER_DPI <- 600

required_pkgs <- c("monocle3", "ggplot2", "dplyr", "viridis", "svglite", "ggrastr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(dplyr)
})

ROOT <- "D:/zz"

STEP25_DIR <- file.path(ROOT, "results", "25_Monocle3_all_myeloid_trajectory")

CDS_FILE <- file.path(STEP25_DIR, "25_Monocle3_all_myeloid_trajectory_3scRNA_v03.rds")

OUT_DIR <- file.path(ROOT, "results", "36_Figure2I_Monocle3_pseudotime_trajectory")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "36_Figure2I_Monocle3_pseudotime_trajectory_v05")

if (!file.exists(CDS_FILE)) {
    stop("Cannot find fixed Step25 Monocle3 CDS:\n", CDS_FILE)
}

cds <- readRDS(CDS_FILE)

if (!inherits(cds, "cell_data_set")) {
    stop("Input RDS is not a Monocle3 cell_data_set.")
}

if (!"UMAP" %in% names(SingleCellExperiment::reducedDims(cds))) {
    stop("UMAP reduction is missing from the CDS.")
}

pg <- monocle3::principal_graph(cds)

if (is.null(pg) || is.null(pg[["UMAP"]])) {
    stop("The CDS does not contain a learned UMAP principal graph.")
}

pt <- monocle3::pseudotime(cds)

n_cells <- ncol(cds)
n_finite <- sum(is.finite(pt))

if (n_finite == 0) {
    stop("No finite pseudotime values are present in the CDS.")
}

pt_summary <- data.frame(metric = c("total_cells", "finite_pseudotime_cells", "finite_pseudotime_percent", 
    "minimum_pseudotime", "median_pseudotime", "maximum_pseudotime"), value = c(n_cells, n_finite, 100 * 
    n_finite/n_cells, min(pt[is.finite(pt)]), median(pt[is.finite(pt)]), max(pt[is.finite(pt)])), stringsAsFactors = FALSE)

write.csv(pt_summary, paste0(OUT_BASE, "_summary.csv"), row.names = FALSE)

p <- monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "pseudotime", show_trajectory_graph = TRUE, 
    label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, 
    cell_size = 0.32, trajectory_graph_color = "black", trajectory_graph_segment_size = 0.62, graph_label_size = 2.8)

umap_xy <- SingleCellExperiment::reducedDims(cds)[["UMAP"]]

if (is.null(rownames(umap_xy))) {
  rownames(umap_xy) <- colnames(cds)
}

pt_named <- monocle3::pseudotime(cds)

if (!is.null(names(pt_named))) {
    pt_plot <- as.numeric(pt_named[rownames(umap_xy)])
} else {
    pt_plot <- as.numeric(pt_named)
}

cell_df <- data.frame(UMAP_1 = umap_xy[, 1], UMAP_2 = umap_xy[, 2], pseudotime = pt_plot, stringsAsFactors = FALSE)

if (nrow(cell_df) != ncol(cds)) {
    stop("Cell plotting table does not match CDS cell number.")
}

layer_classes <- vapply(p$layers, function(layer_obj) {
    class(layer_obj$geom)[1]
}, character(1))

point_layers <- which(vapply(p$layers, function(layer_obj) {
    inherits(layer_obj$geom, "GeomPoint")
}, logical(1)))

if (length(point_layers) == 0) {
    stop(paste0("No GeomPoint layer was detected in plot_cells(). ", "The script will stop rather than pretending rasterization succeeded."))
}

cell_point_layer <- point_layers[1]

raster_cell_layer <- ggrastr::geom_point_rast(data = cell_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = pseudotime), 
    inherit.aes = FALSE, size = 0.32, alpha = 1, stroke = 0, raster.dpi = RASTER_DPI)

p$layers[[cell_point_layer]] <- raster_cell_layer

new_layer_class <- class(p$layers[[cell_point_layer]]$geom)[1]

p <- p + scale_color_viridis_c(option = "plasma", direction = 1, name = "Pseudotime", na.value = "#D9D9D9") + 
    labs(x = "UMAP 1", y = "UMAP 2") + theme_classic(base_family = "Arial", base_size = 12) + theme(text = element_text(family = "Arial", 
    color = "black"), plot.title = element_blank(), axis.title = element_text(family = "Arial", size = 13.4, 
    color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_line(linewidth = 0.52, 
    color = "black"), legend.position = "right", legend.title = element_text(family = "Arial", size = 12.5), 
    legend.text = element_text(family = "Arial", size = 11.5), legend.key.height = grid::unit(12, "mm"), 
    plot.margin = margin(6, 8, 6, 6))

final_layer_classes <- vapply(p$layers, function(layer_obj) {
    class(layer_obj$geom)[1]
}, character(1))

W_IN <- 7.0
H_IN <- 5.6

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

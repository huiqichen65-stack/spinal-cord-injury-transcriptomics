rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- "D:/zz"

INPUT_RDS <- file.path(ROOT, "results", "45_Bulk_batch_structure_and_analysis_object", "45_Bulk_batch_structure_and_analysis_object_v02_analysis_object.rds")

OUT_DIR <- file.path(ROOT, "results", "51_Figure3A_Bulk_PCA")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- "51_Figure3A_Bulk_PCA_v02"

required_pkgs <- c("ggplot2", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0L) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

if (!file.exists(INPUT_RDS)) {
    stop("Step45 object not found:\n", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

required_fields <- c("pca_batch_adjusted", "metadata")

missing_fields <- setdiff(required_fields, names(obj))

if (length(missing_fields) > 0L) {
    stop("Step45 object is missing required field(s): ", paste(missing_fields, collapse = ", "))
}

pca <- obj$pca_batch_adjusted
meta <- obj$metadata

if (is.null(pca$x)) {
  stop("Step45 pca_batch_adjusted does not contain PCA scores.")
}

if (!"sample_id" %in% names(meta)) {
  stop("Step45 metadata lacks sample_id.")
}

scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])

scores$sample_id <- rownames(pca$x)

meta <- meta[match(scores$sample_id, meta$sample_id), , drop = FALSE]

if (anyNA(meta$sample_id) || !identical(scores$sample_id, meta$sample_id)) {
    stop("PCA scores and metadata could not be aligned.")
}

plot_df <- cbind(scores, meta[, c("batch", "timepoint", "phase"), drop = FALSE])

TIME_LEVELS <- c("Uninjured", "15 min", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "42 dpi")

plot_df$timepoint <- factor(as.character(plot_df$timepoint), levels = TIME_LEVELS)

if (nrow(plot_df) != 39L) {
    stop("Unexpected sample count: ", nrow(plot_df), " ; expected 39.")
}

if (anyNA(plot_df$timepoint)) {
  stop("Unexpected or missing timepoint labels in PCA metadata.")
}

var_explained <- 100 * pca$sdev^2/sum(pca$sdev^2)

pc1_pct <- var_explained[1]
pc2_pct <- var_explained[2]

TIME_COLORS <- c(Uninjured = "#4E79A7", `15 min` = "#76B7B2", `1 dpi` = "#59A14F", `3 dpi` = "#EDC948", 
    `7 dpi` = "#F28E2B", `14 dpi` = "#E15759", `28 dpi` = "#B07AA1", `42 dpi` = "#9C755F")

p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = timepoint)) + geom_point(shape = 16, size = 3.4, alpha = 0.92) + 
    scale_color_manual(values = TIME_COLORS, breaks = TIME_LEVELS, drop = FALSE, name = "Time point") + 
    labs(x = sprintf("PC1 (%.1f%%)", pc1_pct), y = sprintf("PC2 (%.1f%%)", pc2_pct)) + theme_classic(base_family = "Arial", 
    base_size = 10) + theme(plot.title = element_blank(), plot.tag = element_blank(), axis.title = element_text(family = "Arial", 
    size = 10, face = "plain"), axis.text = element_text(family = "Arial", size = 9, color = "black"), 
    axis.line = element_line(linewidth = 0.6, color = "black"), axis.ticks = element_line(linewidth = 0.5, 
        color = "black"), axis.ticks.length = grid::unit(1.8, "mm"), legend.position = "right", legend.title = element_text(family = "Arial", 
        size = 9, face = "plain"), legend.text = element_text(family = "Arial", size = 8.5), legend.key.height = grid::unit(4, 
        "mm"), legend.key.width = grid::unit(4, "mm"), panel.grid = element_blank(), plot.margin = margin(t = 5, 
        r = 5, b = 5, l = 5, unit = "pt")) + guides(color = guide_legend(order = 1, override.aes = list(size = 3.2, 
    shape = 16, alpha = 1)))

SVG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".svg"))

PDF_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0(PREFIX, ".png"))

CSV_FILE <- file.path(OUT_DIR, paste0(PREFIX, "_plot_data.csv"))

ggsave(filename = SVG_FILE, plot = p, device = svglite::svglite, width = 5.6, height = 4.3, units = "in")

ggsave(filename = PDF_FILE, plot = p, device = grDevices::cairo_pdf, width = 5.6, height = 4.3, units = "in", 
    family = "Arial")

ggsave(filename = PNG_FILE, plot = p, width = 5.6, height = 4.3, units = "in", dpi = 600, bg = "white")

write.csv(plot_df, CSV_FILE, row.names = FALSE, fileEncoding = "UTF-8")

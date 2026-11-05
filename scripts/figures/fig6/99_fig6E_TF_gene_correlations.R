options(stringsAsFactors = FALSE)
options(width = 180)

expected_tf_n <- 10L
expected_gene_n <- 20L
expected_pair_n <- expected_tf_n * expected_gene_n

selected_tf <- c("Mef2a", "Nr3c1", "Fli1", "Zfp148", "Bclaf1")

candidate_first <- TRUE

png_dpi <- 600

fig_width_in <- 5.5
fig_height_in <- 6.1

selected_tf_color <- "#D73027"
other_tf_color <- "#222222"

cell_border_color <- "#F2F2F2"
star_fontsize <- 8.7
row_fontsize <- 10.5
col_fontsize <- 10.0
legend_fontsize <- 10.0

step90_dir <- "D:/zz/results/90_pySCENIC_top10_TF_network_signature_correlation"

top10_file <- file.path(step90_dir, "02_90_top10_RSS_TFs_v05.csv")

key_gene_file <- file.path(step90_dir, "07_90_Figure6E_key_functional_gene_panel_v05.csv")

corr_file <- file.path(step90_dir, "08_90_Figure6E_TF_functional_gene_correlation_target_cells_v05.csv")

outdir <- "D:/zz/results/99_Figure6E_reference_style_correlation_heatmap"

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

plot_data_file <- file.path(outdir, "01_99_Figure6E_plot_data_v01.csv")

figure_pdf <- file.path(outdir, "02_99_Figure6E_reference_style_correlation_heatmap_v01.pdf")

figure_svg <- file.path(outdir, "02_99_Figure6E_reference_style_correlation_heatmap_v01.svg")

figure_png <- file.path(outdir, "02_99_Figure6E_reference_style_correlation_heatmap_v01.png")

final_qc_file <- file.path(outdir, "03_99_FINAL_QC_v01.csv")

report_file <- file.path(outdir, "04_99_REPORT_v01.txt")

required_packages <- c("ComplexHeatmap", "circlize", "grid")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

required_files <- c(top10_file, key_gene_file, corr_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
    stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

top10 <- read.csv(top10_file, stringsAsFactors = FALSE, check.names = FALSE)

key_panel <- read.csv(key_gene_file, stringsAsFactors = FALSE, check.names = FALSE)

corr <- read.csv(corr_file, stringsAsFactors = FALSE, check.names = FALSE)

required_top10_cols <- c("TF")

required_key_cols <- c("gene")

required_corr_cols <- c("TF", "key_gene", "rho", "FDR")

if (!all(required_top10_cols %in% colnames(top10))) {
    stop("Top-10 TF table is missing required column(s): ", paste(setdiff(required_top10_cols, colnames(top10)), 
        collapse = ", "))
}

if (!all(required_key_cols %in% colnames(key_panel))) {
    stop("Key-gene table is missing required column(s): ", paste(setdiff(required_key_cols, colnames(key_panel)), 
        collapse = ", "))
}

if (!all(required_corr_cols %in% colnames(corr))) {
    stop("Correlation table is missing required column(s): ", paste(setdiff(required_corr_cols, colnames(corr)), 
        collapse = ", "))
}

rss_tf_order <- unique(top10$TF)

if (length(rss_tf_order) != expected_tf_n) {
    stop("Expected exactly ", expected_tf_n, " top-10 TFs, found ", length(rss_tf_order), ".")
}

if (!all(selected_tf %in% rss_tf_order)) {
    stop("One or more selected TFs are absent from the Step 90 top-10 list: ", paste(setdiff(selected_tf, 
        rss_tf_order), collapse = ", "))
}

if (candidate_first) {
    tf_order <- c(selected_tf, setdiff(rss_tf_order, selected_tf))
} else {
    tf_order <- rss_tf_order
}

if ("selection_rank" %in% colnames(key_panel)) {
    key_panel <- key_panel[order(key_panel$selection_rank), , drop = FALSE]
}

gene_order <- unique(key_panel$gene)

if (length(gene_order) != expected_gene_n) {
    stop("Expected exactly ", expected_gene_n, " functional genes, found ", length(gene_order), ".")
}

corr <- corr[corr$TF %in% tf_order & corr$key_gene %in% gene_order, , drop = FALSE]

pair_key <- paste(corr$TF, corr$key_gene, sep = "|||")

if (anyDuplicated(pair_key)) {
    stop("Duplicated TF-gene pairs detected in correlation table.")
}

if (nrow(corr) != expected_pair_n) {
    stop("Expected ", expected_pair_n, " TF-gene pairs, found ", nrow(corr), ".")
}

if (any(!is.finite(corr$rho)) || any(!is.finite(corr$FDR))) {
    stop("Non-finite rho or FDR detected in correlation table.")
}

corr$significance <- ifelse(corr$FDR < 0.001, "***", ifelse(corr$FDR < 0.01, "**", ifelse(corr$FDR < 
    0.05, "*", "")))

corr$TF_display_order <- match(corr$TF, tf_order)

corr$gene_display_order <- match(corr$key_gene, gene_order)

corr <- corr[order(corr$gene_display_order, corr$TF_display_order), , drop = FALSE]

write.csv(corr, plot_data_file, row.names = FALSE)

rho_mat <- matrix(NA_real_, nrow = length(gene_order), ncol = length(tf_order), dimnames = list(gene_order, 
    tf_order))

star_mat <- matrix("", nrow = length(gene_order), ncol = length(tf_order), dimnames = list(gene_order, 
    tf_order))

for (i in seq_len(nrow(corr))) {
    g <- corr$key_gene[i]
    tf <- corr$TF[i]
    rho_mat[g, tf] <- corr$rho[i]
    star_mat[g, tf] <- corr$significance[i]
}

if (anyNA(rho_mat)) {
    stop("NA values remain in rho matrix after pair mapping.")
}

rho_min <- min(rho_mat, na.rm = TRUE)

rho_max <- max(rho_mat, na.rm = TRUE)

if (rho_min >= 0) {
    color_mode <- "sequential_positive"
    upper_limit <- max(0.6, ceiling(rho_max * 10)/10)
    color_fun <- circlize::colorRamp2(c(0, upper_limit * 0.33, upper_limit * 0.67, upper_limit), c("#2C7BB6", 
        "#1A9850", "#FEE08B", "#D73027"))
    legend_at <- c(0, upper_limit * 0.33, upper_limit * 0.67, upper_limit)
} else {
    color_mode <- "diverging_signed"
    abs_limit <- max(abs(c(rho_min, rho_max)))
    abs_limit <- ceiling(abs_limit * 10)/10
    color_fun <- circlize::colorRamp2(c(-abs_limit, 0, abs_limit), c("#2C7BB6", "#F7F7F7", "#D73027"))
    legend_at <- c(-abs_limit, 0, abs_limit)
}

tf_label_colors <- ifelse(tf_order %in% selected_tf, selected_tf_color, other_tf_color)

ht <- ComplexHeatmap::Heatmap(rho_mat, name = "Correlation", col = color_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
    show_row_dend = FALSE, show_column_dend = FALSE, show_row_names = TRUE, show_column_names = TRUE, 
    row_names_side = "left", column_names_side = "top", row_names_gp = grid::gpar(fontsize = row_fontsize, 
        col = "#222222"), column_names_gp = grid::gpar(fontsize = col_fontsize, col = tf_label_colors), 
    column_names_rot = 45, rect_gp = grid::gpar(col = cell_border_color, lwd = 0.7), width = grid::unit(3, 
        "in"), height = grid::unit(4.05, "in"), heatmap_legend_param = list(title = "Correlation", at = legend_at, 
        title_gp = grid::gpar(fontsize = legend_fontsize, fontface = "bold"), labels_gp = grid::gpar(fontsize = legend_fontsize), 
        grid_width = grid::unit(4, "mm"), grid_height = grid::unit(4, "mm")), cell_fun = function(j, 
        i, x, y, width, height, fill) {
        grid::grid.text(star_mat[i, j], x = x, y = y, gp = grid::gpar(fontsize = star_fontsize, col = "#111111", 
            fontface = "bold"))
    })

draw_figure <- function() {
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", padding = grid::unit(c(4, 5, 4, 4), "mm"))
}

grDevices::cairo_pdf(filename = figure_pdf, width = fig_width_in, height = fig_height_in, family = "sans")

draw_figure()

grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = figure_svg, width = fig_width_in, height = fig_height_in)
} else {
    grDevices::svg(filename = figure_svg, width = fig_width_in, height = fig_height_in, family = "sans")
}

draw_figure()

grDevices::dev.off()

if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = figure_png, width = fig_width_in, height = fig_height_in, units = "in", 
        res = png_dpi, background = "white")
} else {
    grDevices::png(filename = figure_png, width = fig_width_in, height = fig_height_in, units = "in", 
        res = png_dpi, bg = "white")
}

draw_figure()

grDevices::dev.off()

qc <- data.frame(check = c("top-10 TF count", "functional-gene count", "TF-gene pair count", "selected five present", 
    "rho matrix complete", "FDR significance used", "Selected TF labels red", "PDF written", "SVG written", 
    "600-dpi PNG written"), status = c(if (length(tf_order) == expected_tf_n) "PASS" else "FAIL", if (length(gene_order) == 
    expected_gene_n) "PASS" else "FAIL", if (nrow(corr) == expected_pair_n) "PASS" else "FAIL", if (all(selected_tf %in% 
    tf_order)) "PASS" else "FAIL", if (!anyNA(rho_mat)) "PASS" else "FAIL", "PASS", "PASS", if (file.exists(figure_pdf)) "PASS" else "FAIL", 
    if (file.exists(figure_svg)) "PASS" else "FAIL", if (file.exists(figure_png)) "PASS" else "FAIL"), 
    detail = c(paste0(length(tf_order), " TFs"), paste0(length(gene_order), " genes"), paste0(nrow(corr), 
        " pairs"), paste(selected_tf, collapse = ", "), paste0("rho range = ", sprintf("%.4f", rho_min), 
        " to ", sprintf("%.4f", rho_max)), "*** FDR<0.001; ** FDR<0.01; * FDR<0.05", paste(selected_tf, 
        collapse = ", "), figure_pdf, figure_svg, paste0(figure_png, " | ", png_dpi, " dpi")), stringsAsFactors = FALSE)

write.csv(qc, final_qc_file, row.names = FALSE)

n_fail <- sum(qc$status == "FAIL")

overall_status <- if (n_fail == 0L) {
  "PASS"
} else {
  "FAIL"
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

if (overall_status != "PASS") {
    stop("Step 99 Figure 6E failed one or more QC checks.")
}

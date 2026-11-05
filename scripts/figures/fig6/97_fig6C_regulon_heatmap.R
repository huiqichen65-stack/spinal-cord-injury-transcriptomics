source("config.R", local = TRUE)

source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

expected_regulons <- 188L
expected_states <- 9L

top_per_state <- 5L

figure6A_top10 <- c("Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1",
  "Chd1")

legacy_target_state <- "Clec7a+ lipid-associated Microglia"

formal_target_state <- "Clec7a+ lipid-handling Microglia"

state_display_map <- c(`Clec7a+ lipid-handling Microglia` = "Clec7a+ lipid MG", `Homeostatic Microglia` = "Homeostatic MG",
  `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ AP Mac", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ repair Mac",
  `Proliferating Microglia` = "Proliferating MG", `Activated inflammatory Microglia` = "Activated inflam. MG",
  `Acute-response Macrophage` = "Acute-response Mac", `Tnf+ inflammatory Microglia` = "Tnf+ inflam. MG",
  `Ccl12+ acute Microglia` = "Ccl12+ acute MG")

target_state_display <- unname(state_display_map[formal_target_state])

cluster_rows_flag <- TRUE
cluster_cols_flag <- TRUE

fig_width_in <- 9.4
fig_height_in <- 8.6
png_dpi <- 600

draw_padding <- grid::unit(c(4, 8, 36, 62), "mm")

rss_file <- paste0(project_path("results/89_pySCENIC_RSS_target_regulon_ranking/"), "05_89_Clec7a_lipid_handling_MG_regulon_RSS_ranking_v01.csv")

state_mean_auc_file <- paste0(project_path("results/88_pySCENIC_regulon_AUCell_full_atlas/"), "14_88_pySCENIC_regulon_state_mean_AUC_v02.csv")

outdir <- project_path("results/97_Figure6C_top10_regulon_state_heatmap")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

schema_file <- file.path(outdir, "01_97_state_mean_AUC_schema_audit_v11.csv")

all_raw_file <- file.path(outdir, "02_97_all188_state_mean_AUC_raw_v11.csv")

all_z_file <- file.path(outdir, "03_97_all188_state_mean_AUC_rowZ_v11.csv")

selection_long_file <- file.path(outdir, "04_97_top5_per_state_selection_long_v11.csv")

selected_union_file <- file.path(outdir, "05_97_selected_regulon_union_v11.csv")

selected_z_file <- file.path(outdir, "06_97_Figure6C_selected_regulon_rowZ_matrix_v11.csv")

row_order_file <- file.path(outdir, "07_97_Figure6C_row_order_v11.csv")

column_order_file <- file.path(outdir, "08_97_Figure6C_column_order_v11.csv")

figure_pdf <- file.path(outdir, "09_97_Figure6C_state_specific_regulon_heatmap_v11.pdf")

figure_svg <- file.path(outdir, "09_97_Figure6C_state_specific_regulon_heatmap_v11.svg")

figure_png <- file.path(outdir, "09_97_Figure6C_state_specific_regulon_heatmap_v11.png")

final_qc_file <- file.path(outdir, "10_97_FINAL_QC_v11.csv")

report_file <- file.path(outdir, "11_97_REPORT_v11.txt")

make_qc <- qc_row

normalize_regulon_name <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("`", "", x, fixed = TRUE)
    x <- gsub("\\(\\+\\)$", "", x)
    x <- gsub("\\(-\\)$", "", x)
    x <- gsub("_extended$", "", x, ignore.case = TRUE)
    x <- gsub("\\s+", "", x)
    x
}

formalize_state <- function(x) {
    x <- as.character(x)
    x[x == legacy_target_state] <- formal_target_state
    x
}

row_zscore <- function(mat) {
    z <- t(apply(mat, 1, function(x) {
        sx <- stats::sd(x, na.rm = TRUE)
        mx <- mean(x, na.rm = TRUE)
        if (!is.finite(sx) || sx == 0) {
            return(rep(0, length(x)))
        }
        (x - mx)/sx
    }))
    dimnames(z) <- dimnames(mat)
    z
}

required_packages <- c("ComplexHeatmap", "circlize")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

required_files <- c(rss_file, state_mean_auc_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

rss <- read.csv(rss_file, stringsAsFactors = FALSE, check.names = FALSE)

required_rss_cols <- c("TF", "RSS_rank", "target_RSS")

if (!all(required_rss_cols %in% colnames(rss))) {
  stop("Step89 RSS table lacks required columns.")
}

rss <- rss[order(rss$RSS_rank), , drop = FALSE]

observed_top10 <- head(rss$TF, 10)

if (!identical(observed_top10, figure6A_top10)) {
  stop("Step89 top10 mismatch.\nObserved: ", paste(observed_top10, collapse = ", "), "\nExpected: ",
    paste(figure6A_top10, collapse = ", "))
}

auc_state <- read.csv(state_mean_auc_file, stringsAsFactors = FALSE, check.names = FALSE)

schema_audit <- data.frame(column = colnames(auc_state), class = vapply(auc_state, function(x) {
  paste(class(x), collapse = "/")
}, character(1)), n_unique = vapply(auc_state, function(x) {
  length(unique(x))
}, integer(1)), stringsAsFactors = FALSE)

write.csv(schema_audit, schema_file, row.names = FALSE)

if (!"myeloid_state" %in% colnames(auc_state)) {
  stop("Verified Step88 schema requires a 'myeloid_state' column.")
}

if (nrow(auc_state) != expected_states) {
  stop("Expected ", expected_states, " myeloid-state rows; observed ", nrow(auc_state), ".")
}

state_names <- formalize_state(auc_state$myeloid_state)

if (anyNA(state_names) || any(!nzchar(state_names)) || anyDuplicated(state_names)) {
  stop("Invalid or duplicated state labels after formal-name mapping.")
}

regulon_cols <- setdiff(colnames(auc_state), "myeloid_state")

if (length(regulon_cols) != expected_regulons) {
  stop("Expected ", expected_regulons, " regulon columns; observed ", length(regulon_cols), ".")
}

tf_names <- normalize_regulon_name(regulon_cols)

if (anyDuplicated(tf_names)) {
    dup <- unique(tf_names[duplicated(tf_names)])
    stop("Duplicated normalized regulon names: ", paste(dup, collapse = ", "))
}

state_by_reg <- as.matrix(auc_state[, regulon_cols, drop = FALSE])

storage.mode(state_by_reg) <- "numeric"

if (any(!is.finite(state_by_reg))) {
  stop("Non-finite values found in Step88 state-mean AUC table.")
}

rownames(state_by_reg) <- state_names

colnames(state_by_reg) <- tf_names

all_raw <- t(state_by_reg)

if (!all(figure6A_top10 %in% rownames(all_raw))) {
  stop("One or more Figure6A top10 regulons are absent from Step88 state-mean AUC.")
}

all_z <- row_zscore(all_raw)

write.csv(data.frame(TF = rownames(all_raw), all_raw, check.names = FALSE), all_raw_file, row.names = FALSE)

write.csv(data.frame(TF = rownames(all_z), all_z, check.names = FALSE), all_z_file, row.names = FALSE)

selection_rows <- list()

for (st in colnames(all_z)) {
    ord <- order(-all_z[, st], -all_raw[, st], rownames(all_z))
    top_idx <- head(ord, top_per_state)
    selection_rows[[st]] <- data.frame(state = st, state_rank = seq_along(top_idx), TF = rownames(all_z)[top_idx], 
        row_Z = all_z[top_idx, st], mean_AUC = all_raw[top_idx, st], stringsAsFactors = FALSE)
}

selection_long <- do.call(rbind, selection_rows)

rownames(selection_long) <- NULL

write.csv(selection_long, selection_long_file, row.names = FALSE)

state_selected_tf <- unique(selection_long$TF)

selected_tf <- unique(c(state_selected_tf, figure6A_top10))

selected_tf <- selected_tf[selected_tf %in% rownames(all_z)]

selected_info <- data.frame(TF = selected_tf, selected_by_state_top10 = selected_tf %in% state_selected_tf,
  is_Figure6A_RSS_top10 = selected_tf %in% figure6A_top10, selected_state_count = vapply(selected_tf,
    function(tf) {
      length(unique(selection_long$state[selection_long$TF == tf]))
    }, integer(1)), target_state_rowZ = all_z[selected_tf, formal_target_state], target_state_mean_AUC = all_raw[selected_tf,
    formal_target_state], stringsAsFactors = FALSE)

selected_info <- selected_info[order(!selected_info$is_Figure6A_RSS_top10, -selected_info$selected_state_count,
  -selected_info$target_state_rowZ, selected_info$TF), , drop = FALSE]

write.csv(selected_info, selected_union_file, row.names = FALSE)

selected_z <- all_z[selected_tf, , drop = FALSE]

write.csv(data.frame(TF = rownames(selected_z), selected_z, check.names = FALSE), selected_z_file, row.names = FALSE)

selected_z_plot <- selected_z

missing_display_names <- setdiff(colnames(selected_z_plot), names(state_display_map))

if (length(missing_display_names) > 0) {
  stop("Missing display abbreviation(s) for state(s): ", paste(missing_display_names, collapse = ", "))
}

colnames(selected_z_plot) <- unname(state_display_map[colnames(selected_z_plot)])

if (anyNA(colnames(selected_z_plot)) || anyDuplicated(colnames(selected_z_plot))) {
  stop("Invalid or duplicated figure-only state abbreviations.")
}

library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#1A5592", "white", "#B83D3D"))

row_label_col <- ifelse(rownames(selected_z) %in% figure6A_top10, "#D73027", "black")

ht <- ComplexHeatmap::Heatmap(selected_z_plot, name = "TF activity", col = col_fun, cluster_rows = cluster_rows_flag,
  cluster_columns = cluster_cols_flag, show_row_dend = cluster_rows_flag, show_column_dend = cluster_cols_flag,
  border = TRUE, rect_gp = grid::gpar(col = "white", lwd = 0.45), row_names_gp = grid::gpar(fontsize = 9.5,
    fontface = "italic", col = row_label_col), column_names_gp = grid::gpar(fontsize = 9), column_names_rot = 35,
  column_names_centered = FALSE, heatmap_legend_param = list(title = "TF activity", at = c(-2, 0, 2),
    labels = c("-2", "0", "2"), legend_height = grid::unit(32, "mm")), use_raster = FALSE)

pdf(NULL)

ht_tmp <- ComplexHeatmap::draw(ht, heatmap_legend_side = "right", padding = draw_padding)

row_ord <- ComplexHeatmap::row_order(ht_tmp)

col_ord <- ComplexHeatmap::column_order(ht_tmp)

dev.off()

if (is.list(row_ord)) {
  row_ord <- unlist(row_ord, use.names = FALSE)
}

if (is.list(col_ord)) {
  col_ord <- unlist(col_ord, use.names = FALSE)
}

row_order_df <- data.frame(plotted_order = seq_along(row_ord), TF = rownames(selected_z_plot)[row_ord],
  is_Figure6A_RSS_top10 = rownames(selected_z_plot)[row_ord] %in% figure6A_top10, stringsAsFactors = FALSE)

column_order_df <- data.frame(plotted_order = seq_along(col_ord), state = colnames(selected_z_plot)[col_ord],
  is_target_state = colnames(selected_z_plot)[col_ord] == target_state_display, stringsAsFactors = FALSE)

write.csv(row_order_df, row_order_file, row.names = FALSE)

write.csv(column_order_df, column_order_file, row.names = FALSE)

grDevices::cairo_pdf(figure_pdf, width = fig_width_in, height = fig_height_in)

ComplexHeatmap::draw(ht, heatmap_legend_side = "right", padding = draw_padding)

dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(figure_svg, width = fig_width_in, height = fig_height_in)
} else {
    svg(figure_svg, width = fig_width_in, height = fig_height_in)
}

ComplexHeatmap::draw(ht, heatmap_legend_side = "right", padding = draw_padding)

dev.off()

png(figure_png, width = fig_width_in, height = fig_height_in, units = "in", res = png_dpi, type = "cairo")

ComplexHeatmap::draw(ht, heatmap_legend_side = "right", padding = draw_padding)

dev.off()

n_selected <- nrow(selected_z)

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step89 Figure6A top10", if (identical(observed_top10, figure6A_top10)) "PASS" else "FAIL",
  paste(observed_top10, collapse = ", "))

qc_list[[length(qc_list) + 1]] <- make_qc("Step88 matrix", if (all(dim(all_raw) == c(expected_regulons,
  expected_states))) "PASS" else "FAIL", paste(dim(all_raw), collapse = " x "))

qc_list[[length(qc_list) + 1]] <- make_qc("Per-state selection", if (nrow(selection_long) == expected_states *
  top_per_state) "PASS" else "FAIL", paste0(top_per_state, " regulons/state x ", expected_states, " states = ",
  nrow(selection_long), " state-regulon selections"))

qc_list[[length(qc_list) + 1]] <- make_qc("Selected union size", if (n_selected >= 20L && n_selected <=
  90L) "PASS" else "WARN", paste0(n_selected, " unique regulons after union + Figure6A top10"))

qc_list[[length(qc_list) + 1]] <- make_qc("Figure6A top10 retained", if (all(figure6A_top10 %in% rownames(selected_z))) "PASS" else "FAIL",
  paste0(sum(figure6A_top10 %in% rownames(selected_z)), " / 10"))

qc_list[[length(qc_list) + 1]] <- make_qc("Target-state formal label", if (formal_target_state %in% colnames(selected_z)) "PASS" else "FAIL",
  formal_target_state)

qc_list[[length(qc_list) + 1]] <- make_qc("Row scaling", "PASS", "Each of 188 regulons was Z-scored across the 9 state means before state-top5 selection.")

qc_list[[length(qc_list) + 1]] <- make_qc("Reference clustering", "PASS", paste0("cluster_rows=", cluster_rows_flag,
  "; cluster_columns=", cluster_cols_flag))

qc_list[[length(qc_list) + 1]] <- make_qc("Figure6A top10 label color", "PASS", "Figure6A RSS top10 are red; all other selected regulons are black.")

qc_list[[length(qc_list) + 1]] <- make_qc("No AUCell/RSS recalculation", "PASS", "Uses Step88 state means and Step89 RSS top10 only.")

qc_list[[length(qc_list) + 1]] <- make_qc("All state display abbreviations", if (length(state_display_map) ==
  expected_states && all(colnames(selected_z) %in% names(state_display_map)) && !anyDuplicated(unname(state_display_map))) "PASS" else "FAIL",
  paste(unname(state_display_map), collapse = " | "))

qc_list[[length(qc_list) + 1]] <- make_qc("Publication outputs", if (all(file.exists(c(figure_pdf, figure_svg,
  figure_png)))) "PASS" else "FAIL", "PDF + SVG + 600 dpi PNG")

final_qc <- do.call(rbind, qc_list)

write.csv(final_qc, final_qc_file, row.names = FALSE)

n_fail <- sum(final_qc$status == "FAIL")

n_warn <- sum(final_qc$status == "WARN")

overall_status <- if (n_fail > 0) {
  "FAIL"
} else if (n_warn > 0) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

target_selected <- selection_long[selection_long$state == formal_target_state, , drop = FALSE]

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

print(target_selected[, c("state_rank", "TF", "row_Z", "mean_AUC"), drop = FALSE], row.names = FALSE)

if (n_fail > 0) {
  stop("Step97 v03 contains blocking QC FAIL item(s).")
}

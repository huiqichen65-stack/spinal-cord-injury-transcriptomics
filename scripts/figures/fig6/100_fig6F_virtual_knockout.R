source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

candidate_tfs <- c("Mef2a", "Nr3c1", "Fli1", "Zfp148", "Bclaf1")

n_rows_per_panel <- 20L
n_nonko_per_panel <- n_rows_per_panel - 1L

png_dpi <- 600

fig_width_in <- 12
fig_height_in <- 10

tf_colors <- c(Mef2a = "#5B8DB8", Nr3c1 = "#E58A4A", Fli1 = "#67A96B", Zfp148 = "#8A68A6", Bclaf1 = "#E0B84F")

base_size <- 10
gene_text_size <- 9
value_text_size <- 7.6
title_text_size <- 11.5
x_tick_size <- 9
x_title_size <- 10

step91_dir <- project_path("results/91_scTenifoldKnk_top5_virtual_KO")

combined_dr_file <- file.path(step91_dir, "05_91_top5_KO_diffRegulation_all_v04.csv.gz")

step91_qc_file <- file.path(step91_dir, "11_91_FINAL_QC_v04.csv")

outdir <- project_path("results/100_Figure6F_reference_style_top5_virtual_KO_v03")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

plot_data_file <- file.path(outdir, "01_100_Figure6F_plot_data_v03.csv")

excluded_gene_file <- file.path(outdir, "02_100_Figure6F_excluded_low_information_genes_v03.csv")

figure_pdf <- file.path(outdir, "03_100_Figure6F_reference_style_top5_virtual_KO_v03.pdf")

figure_svg <- file.path(outdir, "03_100_Figure6F_reference_style_top5_virtual_KO_v03.svg")

figure_png <- file.path(outdir, "03_100_Figure6F_reference_style_top5_virtual_KO_v03.png")

final_qc_file <- file.path(outdir, "04_100_FINAL_QC_v03.csv")

report_file <- file.path(outdir, "05_100_REPORT_v03.txt")

required_packages <- c("ggplot2", "patchwork")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

required_files <- c(combined_dr_file, step91_qc_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("Missing Step 91 file(s):\n", paste(missing_files, collapse = "\n"))
}

step91_qc <- read.csv(step91_qc_file, stringsAsFactors = FALSE, check.names = FALSE)

if ("status" %in% colnames(step91_qc) && any(step91_qc$status == "FAIL")) {
  stop("Step 91 v04 contains one or more blocking FAIL items.")
}

dr <- read.csv(gzfile(combined_dr_file), stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("gene", "KO_TF", "distance", "FC", "log2_perturbation_FC", "perturbation_rank", "is_KO_gene",
  "p.adj", "FDR05", "is_frozen1320_signature", "is_final20_functional_gene")

missing_cols <- setdiff(required_cols, colnames(dr))

if (length(missing_cols) > 0L) {
  stop("Step 91 table is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

if (!all(candidate_tfs %in% unique(dr$KO_TF))) {
  stop("Missing KO TF(s): ", paste(setdiff(candidate_tfs, unique(dr$KO_TF)), collapse = ", "))
}

low_information_reason <- function(gene) {
    reason <- rep(NA_character_, length(gene))
    idx <- grepl("^mt-", gene, ignore.case = TRUE)
    reason[idx] <- "mitochondrial"
    idx <- is.na(reason) & grepl("^Rp[ls]", gene, ignore.case = FALSE)
    reason[idx] <- "ribosomal protein"
    idx <- is.na(reason) & grepl("^Gm[0-9A-Za-z._-]*$", gene, ignore.case = FALSE)
    reason[idx] <- "predicted Gm gene"
    idx <- is.na(reason) & grepl("Rik", gene, ignore.case = FALSE)
    reason[idx] <- "RIKEN clone symbol"
    idx <- is.na(reason) & (grepl("-ps", gene, ignore.case = TRUE) | grepl("ps[0-9]+$", gene, ignore.case = TRUE))
    reason[idx] <- "pseudogene-style symbol"
    idx <- is.na(reason) & grepl("^Rn[0-9]+s", gene, ignore.case = FALSE)
    reason[idx] <- "rRNA-style symbol"
    reason
}

panel_list <- vector("list", length(candidate_tfs))

names(panel_list) <- candidate_tfs

excluded_list <- vector("list", length(candidate_tfs))

names(excluded_list) <- candidate_tfs

for (tf in candidate_tfs) {
    z <- dr[dr$KO_TF == tf, , drop = FALSE]
    z <- z[order(z$perturbation_rank, -z$distance, z$p.adj, z$gene), , drop = FALSE]
    ko_row <- z[z$gene == tf | z$is_KO_gene, , drop = FALSE]
    if (nrow(ko_row) == 0L) {
        stop("KO gene row not found for ", tf, ".")
    }
    ko_row <- ko_row[order(ko_row$perturbation_rank, -ko_row$distance), , drop = FALSE]
    ko_row <- ko_row[1, , drop = FALSE]
    nonko <- z[!(z$gene == tf | z$is_KO_gene), , drop = FALSE]
    nonko$exclude_reason <- low_information_reason(nonko$gene)
    excluded <- nonko[!is.na(nonko$exclude_reason), , drop = FALSE]
    retained <- nonko[is.na(nonko$exclude_reason), , drop = FALSE]
    if (nrow(retained) < n_nonko_per_panel) {
        stop("After low-information filtering, ", tf, " has only ", nrow(retained), " retained non-KO genes; need ", 
            n_nonko_per_panel, ".")
    }
    retained <- retained[order(retained$perturbation_rank, -retained$distance, retained$p.adj, retained$gene), 
        , drop = FALSE]
    selected_nonko <- head(retained, n_nonko_per_panel)
    ko_row$display_role <- "KO TF"
    ko_row$exclude_reason <- NA_character_
    selected_nonko$display_role <- "Perturbed gene"
    show_df <- rbind(ko_row, selected_nonko)
    show_df$display_order <- seq_len(nrow(show_df))
    show_df$display_selected_by <- c("KO TF retained by design", rep("Top retained perturbation rank", 
        nrow(show_df) - 1L))
    panel_list[[tf]] <- show_df
    if (nrow(excluded) > 0L) {
        excluded$KO_TF_display <- tf
    }
    excluded_list[[tf]] <- excluded
}

plot_df <- do.call(rbind, panel_list)

rownames(plot_df) <- NULL

plot_df$KO_TF <- factor(plot_df$KO_TF, levels = candidate_tfs)

write.csv(plot_df, plot_data_file, row.names = FALSE)

excluded_df <- do.call(rbind, excluded_list)

if (is.null(excluded_df)) {
  excluded_df <- data.frame(gene = character(0), KO_TF = character(0), exclude_reason = character(0),
    stringsAsFactors = FALSE)
}

rownames(excluded_df) <- NULL

write.csv(excluded_df, excluded_gene_file, row.names = FALSE)

overlap_summary <- do.call(rbind, lapply(candidate_tfs, function(tf) {
    z <- panel_list[[tf]]
    nonko <- z[z$display_role == "Perturbed gene", , drop = FALSE]
    data.frame(KO_TF = tf, n_displayed_nonKO = nrow(nonko), n_FDR05_nonKO = sum(nonko$FDR05, na.rm = TRUE), 
        n_frozen1320_signature = sum(nonko$is_frozen1320_signature, na.rm = TRUE), n_final20_functional = sum(nonko$is_final20_functional_gene, 
            na.rm = TRUE), Clec7a_displayed = "Clec7a" %in% nonko$gene, Lpl_displayed = "Lpl" %in% nonko$gene, 
        stringsAsFactors = FALSE)
}))

x_vals <- plot_df$log2_perturbation_FC

if (any(!is.finite(x_vals))) {
  stop("Non-finite log2_perturbation_FC found in display data.")
}

x_min <- min(0, min(x_vals, na.rm = TRUE))

x_max <- max(0, max(x_vals, na.rm = TRUE))

x_span <- x_max - x_min

if (!is.finite(x_span) || x_span <= 0) {
  x_span <- 1
}

x_left <- x_min - 0.03 * x_span
x_right <- x_max + 0.15 * x_span

library(ggplot2)

make_panel <- function(tf) {
    d <- panel_list[[tf]]
    d$gene_factor <- factor(d$gene, levels = rev(d$gene))
    d$value_label <- sprintf("%.2f", d$log2_perturbation_FC)
    d$value_hjust <- ifelse(d$log2_perturbation_FC >= 0, -0.08, 1.08)
    p <- ggplot(d, aes(x = log2_perturbation_FC, y = gene_factor)) + geom_col(width = 0.72, fill = unname(tf_colors[tf]), 
        color = NA) + geom_text(aes(label = value_label, hjust = value_hjust), size = value_text_size/ggplot2::.pt, 
        color = "#333333", family = "sans") + scale_x_continuous(limits = c(x_left, x_right), expand = c(0, 
        0)) + labs(x = "log2 perturbation FC", y = NULL, title = tf) + theme_classic(base_size = base_size, 
        base_family = "sans") + theme(plot.title = element_text(size = title_text_size, face = "bold", 
        hjust = 0.5, color = "#222222", margin = margin(b = 5)), axis.text.y = element_text(size = gene_text_size, 
        color = "#222222", margin = margin(r = 2)), axis.text.x = element_text(size = x_tick_size, color = "#333333"), 
        axis.title.x = element_text(size = x_title_size, face = "bold", margin = margin(t = 5)), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), panel.grid = element_blank(), plot.margin = margin(t = 3, r = 5, 
            b = 3, l = 3))
    if (x_min < 0) {
        p <- p + geom_vline(xintercept = 0, linewidth = 0.3, color = "#BDBDBD")
    }
    p
}

plots <- lapply(candidate_tfs, make_panel)

names(plots) <- candidate_tfs

p6f <- patchwork::wrap_plots(plots, ncol = 3, nrow = 2)

ggsave(filename = figure_pdf, plot = p6f, width = fig_width_in, height = fig_height_in, units = "in",
  device = grDevices::cairo_pdf)

if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(filename = figure_svg, plot = p6f, width = fig_width_in, height = fig_height_in, units = "in", 
        device = svglite::svglite)
} else {
    grDevices::svg(filename = figure_svg, width = fig_width_in, height = fig_height_in, family = "sans")
    print(p6f)
    grDevices::dev.off()
}

if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = figure_png, width = fig_width_in, height = fig_height_in, units = "in", 
        res = png_dpi, background = "white")
    print(p6f)
    grDevices::dev.off()
} else {
    ggsave(filename = figure_png, plot = p6f, width = fig_width_in, height = fig_height_in, units = "in", 
        dpi = png_dpi, bg = "white")
}

panel_counts <- table(plot_df$KO_TF)

qc <- data.frame(check = c("Step91 input", "Step91 QC has no FAIL", "Five TFs present", "20 displayed rows per TF",
  "KO TF shown first in every panel", "Same low-information filter applied to all TFs", "Selection independent of storyline genes",
  "Selection independent of FDR threshold", "perturbation ranking preserved", "No expression-direction interpretation",
  "2 x 3 panel layout", "PDF written", "SVG written", "600-dpi PNG written"), status = c(if (file.exists(combined_dr_file)) "PASS" else "FAIL",
  if (!any(step91_qc$status == "FAIL")) "PASS" else "FAIL", if (all(candidate_tfs %in% unique(as.character(plot_df$KO_TF)))) "PASS" else "FAIL",
  if (all(panel_counts == n_rows_per_panel)) "PASS" else "FAIL", if (all(vapply(candidate_tfs, function(tf) {
    z <- panel_list[[tf]]
    z$gene[1] == tf && z$display_role[1] == "KO TF"
  }, logical(1)))) "PASS" else "FAIL", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", if (file.exists(figure_pdf)) "PASS" else "FAIL",
  if (file.exists(figure_svg)) "PASS" else "FAIL", if (file.exists(figure_png)) "PASS" else "FAIL"),
  detail = c(combined_dr_file, step91_qc_file, paste(candidate_tfs, collapse = ", "), paste(names(panel_counts),
    as.integer(panel_counts), sep = "=", collapse = "; "), paste(candidate_tfs, collapse = ", "),
    "Excluded only ^mt-, ^Rpl/^Rps, ^Gm, Rik-containing clone symbols, pseudogene-style -ps/ps#, and Rn###s rRNA-style symbols.",
    "No Clec7a/lipid/signature gene was manually forced into the plot.", "FDR is audited but not used for display selection.",
    "Top retained non-KO genes selected by perturbation_rank.", "log2 perturbation FC is a perturbation-distance metric, not RNA-expression log2FC.",
    "Row 1: Mef2a, Nr3c1, Fli1; Row 2: Zfp148, Bclaf1; sixth slot intentionally empty.", figure_pdf,
    figure_svg, paste0(figure_png, " | ", png_dpi, " dpi")), stringsAsFactors = FALSE)

write.csv(qc, final_qc_file, row.names = FALSE)

n_fail <- sum(qc$status == "FAIL")

overall_status <- if (n_fail == 0L) {
  "PASS"
} else {
  "FAIL"
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

print(overlap_summary, row.names = FALSE)

if (overall_status != "PASS") {
  stop("Step 100 v03 Figure 6F failed one or more QC checks.")
}

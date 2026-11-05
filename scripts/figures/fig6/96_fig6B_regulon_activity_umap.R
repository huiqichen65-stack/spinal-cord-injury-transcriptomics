source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

expected_cells <- 58533L
expected_regulons <- 188L
reduction_name <- "umap_myeloid_harmony30"

top10_tf <- c("Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1", "Chd1")

ncol_panel <- 5L
nrow_panel <- 2L

clip_lower_quantile <- 0.00
clip_upper_quantile <- 0.99

point_size <- 0.12
point_alpha <- 1

fig_width_in <- 11.0
fig_height_in <- 5.4
png_dpi <- 600

myeloid_rds <- paste0("D:/zz/results/20_AUCell_lipid_metabolism_scoring/", "20_AUCell_lipid_metabolism_scoring_3scRNA_v03.rds")

rss_file <- paste0("D:/zz/results/89_pySCENIC_RSS_target_regulon_ranking/", "05_89_Clec7a_lipid_handling_MG_regulon_RSS_ranking_v01.csv")

aucell_file <- paste0("D:/zz/results/88_pySCENIC_regulon_AUCell_full_atlas/", "11_88_pySCENIC_AUCell_full_atlas_v02.csv.gz")

outdir <- "D:/zz/results/96_Figure6B_top10_regulon_activity_UMAP"

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

input_qc_file <- file.path(outdir, "01_96_input_QC_v03.csv")

regulon_match_file <- file.path(outdir, "02_96_top10_regulon_column_matching_v03.csv")

plot_data_file <- file.path(outdir, "03_96_Figure6B_plot_data_summary_v03.csv")

figure_pdf <- file.path(outdir, "04_96_Figure6B_top10_regulon_activity_UMAP_v03.pdf")

figure_svg <- file.path(outdir, "04_96_Figure6B_top10_regulon_activity_UMAP_v03.svg")

figure_png <- file.path(outdir, "04_96_Figure6B_top10_regulon_activity_UMAP_v03.png")

final_qc_file <- file.path(outdir, "05_96_FINAL_QC_v03.csv")

report_file <- file.path(outdir, "06_96_REPORT_v03.txt")

make_qc <- qc_row

normalize_regulon_name <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x <- gsub(
    "`",
    "",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "\\(\\+\\)$",
    "",
    x
  )

  x <- gsub(
    "\\(-\\)$",
    "",
    x
  )

  x <- gsub(
    "_extended$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- gsub(
    "\\s+",
    "",
    x
  )

  x
}

resolve_top10_columns <- function(column_names, tf_names) {

  normalized_cols <- normalize_regulon_name(
    column_names
  )

  out <- lapply(
    tf_names,
    function(tf) {

      hit <- which(
        normalized_cols == tf
      )

      if (length(hit) != 1L) {

        hit2 <- which(
          tolower(normalized_cols) == tolower(tf)
        )

        hit <- unique(
          c(
            hit,
            hit2
          )
        )
      }

      data.frame(
        TF = tf,
        n_matches = length(hit),
        matched_column = if (
          length(hit) == 1L
        ) {
          column_names[hit]
        } else {
          NA_character_
        },
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    out
  )
}

required_packages <- c("SeuratObject", "ggplot2", "patchwork")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

required_files <- c(myeloid_rds, rss_file, aucell_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

rss <- read.csv(rss_file, stringsAsFactors = FALSE, check.names = FALSE)

required_rss_cols <- c("TF", "RSS_rank", "target_RSS")

if (!all(required_rss_cols %in% colnames(rss))) {
    stop("Step89 RSS file lacks required columns.")
}

rss <- rss[order(rss$RSS_rank), , drop = FALSE]

observed_top10 <- head(rss$TF, 10)

if (!identical(observed_top10, top10_tf)) {
    stop("Step89 top10 mismatch.\nObserved: ", paste(observed_top10, collapse = ", "), "\nExpected: ", 
        paste(top10_tf, collapse = ", "))
}

obj <- readRDS(myeloid_rds)

if (!reduction_name %in% names(obj@reductions)) {
    stop("reduction not found: ", reduction_name)
}

umap <- SeuratObject::Embeddings(obj, reduction = reduction_name)

if (nrow(umap) != expected_cells) {
    stop("Unexpected UMAP cell count: ", nrow(umap), " / expected ", expected_cells)
}

if (ncol(umap) < 2L) {
    stop("UMAP reduction has fewer than two dimensions.")
}

umap_df <- data.frame(cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], stringsAsFactors = FALSE)

rm(obj, umap)

gc()

auc_raw <- read.csv(gzfile(aucell_file), stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(auc_raw) != expected_cells) {
    stop("Unexpected Step88 AUC row count: ", nrow(auc_raw), " / expected ", expected_cells)
}

cell_col_candidates <- c("cell", "Cell", "cell_id", "CellID", "barcode", "Barcode")

cell_col <- intersect(cell_col_candidates, colnames(auc_raw))

if (length(cell_col) == 1L) {

  cell_col <- cell_col[1]

} else {

  match_fraction <- vapply(
    colnames(auc_raw),
    function(nm) {

      x <- as.character(
        auc_raw[[nm]]
      )

      if (
        length(unique(x)) != nrow(auc_raw)
      ) {
        return(0)
      }

      mean(
        x %in% umap_df$cell
      )
    },
    numeric(1)
  )

  best <- which.max(
    match_fraction
  )

  if (
    length(best) != 1L ||
    match_fraction[best] < 0.99
  ) {
    stop(
      "Could not confidently identify the Step88 cell-ID column."
    )
  }

  cell_col <- colnames(auc_raw)[best]
}

auc_cells <- as.character(auc_raw[[cell_col]])

if (anyDuplicated(auc_cells)) {
    stop("Duplicated cell IDs found in Step88 AUC table.")
}

if (!setequal(auc_cells, umap_df$cell)) {
    stop("Step88 AUC cells do not exactly match UMAP cells.")
}

regulon_columns <- setdiff(colnames(auc_raw), cell_col)

match_table <- resolve_top10_columns(regulon_columns, top10_tf)

if (any(match_table$n_matches != 1L)) {
    print(match_table)
    stop("Could not uniquely match all top10 TFs to Step88 AUC columns.")
}

write.csv(match_table, regulon_match_file, row.names = FALSE)

matched_cols <- match_table$matched_column
names(matched_cols) <- match_table$TF

auc_idx <- match(umap_df$cell, auc_cells)

if (anyNA(auc_idx)) {
    stop("Failed to align Step88 AUC cells to UMAP.")
}

auc_top10 <- auc_raw[auc_idx, matched_cols, drop = FALSE]

colnames(auc_top10) <- top10_tf

for (tf in top10_tf) {
    x <- suppressWarnings(as.numeric(auc_top10[[tf]]))
    if (any(!is.finite(x))) {
        stop("Non-finite AUC values detected for ", tf)
    }
    auc_top10[[tf]] <- x
}

rm(auc_raw)

gc()

pooled_values <- unlist(auc_top10[, top10_tf, drop = FALSE], use.names = FALSE)

common_limits <- as.numeric(stats::quantile(pooled_values, probs = c(clip_lower_quantile, clip_upper_quantile), 
    na.rm = TRUE, names = FALSE, type = 7))

if (any(!is.finite(common_limits)) || diff(common_limits) <= 0) {
    stop("Invalid common AUC color limits.")
}

library(ggplot2)
library(patchwork)

use_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)

plots <- vector("list", length(top10_tf))

names(plots) <- top10_tf

summary_rows <- list()

for (i in seq_along(top10_tf)) {

  tf <- top10_tf[i]

  d <- umap_df

  d$AUC <- auc_top10[[tf]]

  d <- d[
    order(
      d$AUC,
      d$cell
    ),
    ,
    drop = FALSE
  ]

  d$AUC_plot <- pmin(
    pmax(
      d$AUC,
      common_limits[1]
    ),
    common_limits[2]
  )

  p <- ggplot(
    d,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = AUC_plot
    )
  )

  if (use_ggrastr) {

    p <- p +
      ggrastr::geom_point_rast(
        size = point_size,
        alpha = point_alpha,
        raster.dpi = 600
      )

  } else {

    p <- p +
      geom_point(
        size = point_size,
        alpha = point_alpha
      )
  }

  p <- p +
    scale_color_gradientn(
      colors = c(
        "#440154",
        "#3B528B",
        "#21918C",
        "#5DC863",
        "#FDE725"
      ),
      limits = common_limits,
      oob = scales::squish,
      name = "TF activity"
    ) +
    coord_equal() +
    labs(
      title = tf,
      x = NULL,
      y = NULL
    ) +
    theme_void(
      base_size = 10.5
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "italic",
        size = 12,
        margin = margin(
          b = 2
        )
      ),
      plot.margin = margin(
        2,
        2,
        2,
        2
      ),
      legend.title = element_text(
        size = 10
      ),
      legend.text = element_text(
        size = 9
      )
    )

  plots[[tf]] <- p

  summary_rows[[tf]] <- data.frame(
    RSS_rank = i,
    TF = tf,
    matched_regulon_column = matched_cols[[tf]],
    n_cells = nrow(d),
    AUC_min = min(d$AUC),
    AUC_median = stats::median(d$AUC),
    AUC_mean = mean(d$AUC),
    AUC_max = max(d$AUC),
    common_scale_min = common_limits[1],
    common_scale_max = common_limits[2],
    stringsAsFactors = FALSE
  )
}

plot_summary <- do.call(rbind, summary_rows)

rownames(plot_summary) <- NULL

write.csv(plot_summary, plot_data_file, row.names = FALSE)

figure6b <- patchwork::wrap_plots(plots[top10_tf], ncol = ncol_panel, nrow = nrow_panel, guides = "collect", 
    byrow = TRUE) & theme(legend.position = "right")

ggsave(figure_pdf, figure6b, width = fig_width_in, height = fig_height_in, units = "in", device = cairo_pdf, 
    limitsize = FALSE)

if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(figure_svg, figure6b, width = fig_width_in, height = fig_height_in, units = "in", device = svglite::svglite, 
        limitsize = FALSE)
} else {
    svg(figure_svg, width = fig_width_in, height = fig_height_in)
    print(figure6b)
    dev.off()
}

ggsave(figure_png, figure6b, width = fig_width_in, height = fig_height_in, units = "in", dpi = png_dpi, 
    limitsize = FALSE)

input_qc <- data.frame(item = c("myeloid RDS", "Step89 RSS", "Step88 AUCell"), path = c(myeloid_rds, 
    rss_file, aucell_file), exists = file.exists(c(myeloid_rds, rss_file, aucell_file)), stringsAsFactors = FALSE)

write.csv(input_qc, input_qc_file, row.names = FALSE)

qc_list <- list()

qc_list[[length(qc_list) + 1]] <- make_qc("Step89 top10", if (identical(observed_top10, top10_tf)) "PASS" else "FAIL", 
    paste(observed_top10, collapse = ", "))

qc_list[[length(qc_list) + 1]] <- make_qc("UMAP cells", if (nrow(umap_df) == expected_cells) "PASS" else "FAIL", 
    paste0(nrow(umap_df), " / ", expected_cells))

qc_list[[length(qc_list) + 1]] <- make_qc("Step88 AUC cells match UMAP", "PASS", paste0(length(auc_cells), 
    " exact cell IDs"))

qc_list[[length(qc_list) + 1]] <- make_qc("Top10 regulon-column matching", if (all(match_table$n_matches == 
    1L)) "PASS" else "FAIL", paste0(sum(match_table$n_matches == 1L), " / 10 unique matches"))

qc_list[[length(qc_list) + 1]] <- make_qc("Panel layout", "PASS", "2 rows x 5 columns; RSS rank 1-5 top row, 6-10 bottom row.")

qc_list[[length(qc_list) + 1]] <- make_qc("Shared activity color scale", "PASS", paste0("pooled top10 quantiles ", 
    clip_lower_quantile, "-", clip_upper_quantile, ": ", signif(common_limits[1], 4), " to ", signif(common_limits[2], 
        4)))

qc_list[[length(qc_list) + 1]] <- make_qc("No AUCell recalculation", "PASS", "Reads Step88 full-atlas AUCell scores only.")

qc_list[[length(qc_list) + 1]] <- make_qc("Publication outputs", if (all(file.exists(c(figure_pdf, figure_svg, 
    figure_png)))) "PASS" else "FAIL", "PDF + SVG + 600 dpi PNG")

final_qc <- do.call(rbind, qc_list)

write.csv(final_qc, final_qc_file, row.names = FALSE)

if (any(final_qc$status == "FAIL")) {
    stop("Step96 contains blocking QC FAIL item(s).")
}

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

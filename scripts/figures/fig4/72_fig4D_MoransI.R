source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(7201)

ROOT <- "D:/zz"

STEP66_DIR <- file.path(ROOT, "results", "66_GSE195783_Clec7a_spatial_autocorrelation_neighborhood")

OUT_DIR <- file.path(ROOT, "results", "72_Figure4D_Clec7a_MoransI_section_summary")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v04"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

TIME_LABELS <- c(`3dpi` = "3 dpi", `7dpi` = "7 dpi", `14dpi` = "14 dpi", `28dpi` = "28 dpi")

required_pkgs <- c("ggplot2", "dplyr", "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step72.")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

write_csv_utf8 <- write_csv_utf8

pick_unique_col <- function(nms, exact = character(), regex = character(), label) {
    exact_hit <- intersect(exact, nms)
    if (length(exact_hit) == 1L) {
        return(exact_hit)
    }
    regex_hit <- unique(unlist(lapply(regex, function(pat) {
        grep(pat, nms, ignore.case = TRUE, value = TRUE)
    })))
    regex_hit <- setdiff(regex_hit, exact_hit)
    all_hit <- unique(c(exact_hit, regex_hit))
    if (length(all_hit) != 1L) {
        stop("Could not identify exactly one column for ", label, ".\nCandidate columns: ", ifelse(length(all_hit) == 
            0L, "<none>", paste(all_hit, collapse = " | ")), "\nAvailable columns:\n", paste(nms, collapse = " | "))
    }
    all_hit
}

normalize_timepoint <- function(x) {
    z <- trimws(as.character(x))
    z_low <- tolower(z)
    z_low <- gsub("\\s+", "", z_low)
    out <- rep(NA_character_, length(z_low))
    out[z_low %in% c("3", "3d", "3dpi", "day3", "3day", "3days")] <- "3dpi"
    out[z_low %in% c("7", "7d", "7dpi", "day7", "7day", "7days")] <- "7dpi"
    out[z_low %in% c("14", "14d", "14dpi", "day14", "14day", "14days")] <- "14dpi"
    out[z_low %in% c("28", "28d", "28dpi", "day28", "28day", "28days")] <- "28dpi"
    out
}

make_symmetric_offsets <- function(n, max_offset = 0.17) {
    if (n <= 1L) {
        return(0)
    }
    seq(-max_offset, max_offset, length.out = n)
}

INPUT_FILE <- file.path(STEP66_DIR, "02_66_Clec7a_section_global_MoransI_v02.csv")

if (!dir.exists(STEP66_DIR)) {
    stop("Step66 directory does not exist:\n", STEP66_DIR)
}

if (!file.exists(INPUT_FILE)) {
    stop("Step66 section-level Moran input is missing:\n", INPUT_FILE)
}

moran_files <- INPUT_FILE

dat0 <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

nms <- colnames(dat0)

TIME_COL <- pick_unique_col(nms = nms, exact = c("timepoint", "time", "dpi"), regex = c("^timepoint$", 
    "^time$", "^dpi$"), label = "timepoint")

SECTION_COL <- pick_unique_col(nms = nms, exact = c("section_uid_internal", "section_uid", "section_id"), 
    regex = c("section.*uid", "^section_id$", "^section$"), label = "section ID")

MORAN_COL <- pick_unique_col(nms = nms, exact = c("moran_I", "Moran_I", "morans_I", "morans_i", "moran_i"), 
    regex = c("moran.*(^|_)i$", "^I_moran$", "^moran"), label = "Moran's I")

sample_hits <- intersect(c("sample_label", "sample", "capture_area", "capture_area_label", "GSM"), nms)

SAMPLE_COL <- if (length(sample_hits) >= 1L) {
  sample_hits[[1]]
} else {
  NA_character_
}

p_hits <- grep("^p$|pvalue|p_value|p\\.value|moran.*p", nms, ignore.case = TRUE, value = TRUE)

P_COL <- if (length(p_hits) == 1L) {
  p_hits[[1]]
} else {
  NA_character_
}

fdr_hits <- grep("fdr|padj|p_adj|adj.*p|qvalue", nms, ignore.case = TRUE, value = TRUE)

FDR_COL <- if (length(fdr_hits) == 1L) {
  fdr_hits[[1]]
} else {
  NA_character_
}

plot_df <- data.frame(timepoint_raw = dat0[[TIME_COL]], section_uid_internal = as.character(dat0[[SECTION_COL]]), 
    moran_I = suppressWarnings(as.numeric(dat0[[MORAN_COL]])), stringsAsFactors = FALSE)

plot_df$timepoint <- normalize_timepoint(plot_df$timepoint_raw)

if (!is.na(SAMPLE_COL)) {
    plot_df$sample_label <- as.character(dat0[[SAMPLE_COL]])
} else {
    plot_df$sample_label <- NA_character_
}

if (!is.na(P_COL)) {
    plot_df$p_value <- suppressWarnings(as.numeric(dat0[[P_COL]]))
} else {
    plot_df$p_value <- NA_real_
}

if (!is.na(FDR_COL)) {
    plot_df$FDR <- suppressWarnings(as.numeric(dat0[[FDR_COL]]))
} else {
    plot_df$FDR <- NA_real_
}

plot_df <- plot_df %>% filter(timepoint %in% TIME_LEVELS, nzchar(section_uid_internal), is.finite(moran_I)) %>% 
    distinct(timepoint, section_uid_internal, .keep_all = TRUE)

if (nrow(plot_df) == 0L) {
    stop("No valid 3/7/14/28 dpi section-level Moran rows remain.")
}

dup_check <- plot_df %>% count(timepoint, section_uid_internal) %>% filter(n > 1L)

if (nrow(dup_check) > 0L) {
    stop("Duplicate timepoint + section rows detected.")
}

plot_df <- plot_df %>% group_by(timepoint) %>% arrange(moran_I, section_uid_internal, .by_group = TRUE) %>% 
    mutate(x_base = match(timepoint, TIME_LEVELS), x_offset = make_symmetric_offsets(dplyr::n(), max_offset = 0.16), 
        x_plot = x_base + x_offset) %>% ungroup()

summary_df <- plot_df %>% group_by(timepoint) %>% summarise(n_sections = n(), median_moran_I = median(moran_I, 
    na.rm = TRUE), mean_moran_I = mean(moran_I, na.rm = TRUE), min_moran_I = min(moran_I, na.rm = TRUE), 
    max_moran_I = max(moran_I, na.rm = TRUE), n_positive = sum(moran_I > 0, na.rm = TRUE), .groups = "drop") %>% 
    mutate(x_base = match(timepoint, TIME_LEVELS), timepoint_label = unname(TIME_LABELS[timepoint]))

PLOT_DATA_FILE <- file.path(OUT_DIR, paste0("01_72_Figure4D_section_MoransI_plot_data_", VERSION, ".csv"))

write_csv_utf8(plot_df %>% arrange(factor(timepoint, levels = TIME_LEVELS), section_uid_internal), PLOT_DATA_FILE)

SUMMARY_FILE <- file.path(OUT_DIR, paste0("02_72_Figure4D_timepoint_MoransI_summary_", VERSION, ".csv"))

write_csv_utf8(summary_df, SUMMARY_FILE)

TIME_COLORS <- c(`3dpi` = "#9A9A9A", `7dpi` = "#E2A15C", `14dpi` = "#C94844", `28dpi` = "#9E3340")

MEDIAN_LINE_COLOR <- "#7F3048"
REFERENCE_COLOR <- "#9A9A9A"

y_min_data <- min(0, min(plot_df$moran_I, na.rm = TRUE))

y_max_data <- max(plot_df$moran_I, na.rm = TRUE)

y_span <- y_max_data - y_min_data

if (!is.finite(y_span) || y_span <= 0) {
  y_span <- 0.1
}

Y_LOWER <- y_min_data - 0.04 * y_span
Y_UPPER <- y_max_data + 0.13 * y_span

summary_df$timepoint <- factor(summary_df$timepoint, levels = TIME_LEVELS)

plot_df$timepoint <- factor(plot_df$timepoint, levels = TIME_LEVELS)

p <- ggplot() +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = REFERENCE_COLOR
  ) +

  geom_line(
    data = summary_df,
    aes(
      x = x_base,
      y = median_moran_I,
      group = 1
    ),
    linewidth = 1.15,
    lineend = "round",
    color = MEDIAN_LINE_COLOR
  ) +

  geom_point(
    data = plot_df,
    aes(
      x = x_plot,
      y = moran_I,
      color = timepoint
    ),
    shape = 16,
    size = 2.05,
    alpha = 0.62
  ) +

  geom_point(
    data = summary_df,
    aes(
      x = x_base,
      y = median_moran_I,
      fill = timepoint
    ),
    shape = 21,
    size = 4.4,
    stroke = 1.0,
    color = "white"
  ) +

  geom_point(
    data = summary_df,
    aes(
      x = x_base,
      y = median_moran_I
    ),
    shape = 21,
    size = 4.4,
    stroke = 0.85,
    fill = NA,
    color = MEDIAN_LINE_COLOR
  ) +

  annotate(
    "text",
    x = 4.36,
    y = Y_UPPER - 0.035 * y_span,
    label = "Positive spatial autocorrelation",
    family = "Arial",
    size = 2.75,
    hjust = 1,
    vjust = 1,
    color = "grey38"
  ) +

  scale_color_manual(
    values = TIME_COLORS,
    guide = "none",
    drop = FALSE
  ) +

  scale_fill_manual(
    values = TIME_COLORS,
    guide = "none",
    drop = FALSE
  ) +

  scale_x_continuous(
    breaks = seq_along(TIME_LEVELS),
    labels = unname(
      TIME_LABELS[
        TIME_LEVELS
      ]
    ),
    limits = c(
      0.55,
      4.45
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_y_continuous(
    limits = c(
      Y_LOWER,
      Y_UPPER
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = NULL,
    y = "Moran's I"
  ) +

  theme_classic(
    base_family = "Arial",
    base_size = 9
  ) +

  theme(
    axis.title.y = element_text(
      size = 9.5,
      face = "plain",
      color = "black",
      margin = margin(
        r = 6
      )
    ),

    axis.text.x = element_text(
      size = 9,
      face = "bold",
      color = "black",
      margin = margin(
        t = 5
      )
    ),

    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),

    axis.line = element_line(
      linewidth = 0.55,
      color = "black"
    ),

    panel.background = element_rect(
      fill = "white",
      color = NA
    ),

    plot.background = element_rect(
      fill = "white",
      color = NA
    ),

    plot.margin = margin(
      7,
      8,
      6,
      6
    )
  )

PDF_FILE <- file.path(OUT_DIR, paste0("72_Figure4D_Clec7a_MoransI_section_summary_", VERSION, ".pdf"))

PNG_FILE <- file.path(OUT_DIR, paste0("72_Figure4D_Clec7a_MoransI_section_summary_", VERSION, ".png"))

SVG_FILE <- file.path(OUT_DIR, paste0("72_Figure4D_Clec7a_MoransI_section_summary_", VERSION, ".svg"))

ggsave(filename = PDF_FILE, plot = p, width = 4.6, height = 3.65, units = "in", device = grDevices::cairo_pdf, 
    bg = "white")

ggsave(filename = PNG_FILE, plot = p, width = 4.6, height = 3.65, units = "in", dpi = 600, bg = "white")

ggsave(filename = SVG_FILE, plot = p, width = 4.6, height = 3.65, units = "in", device = "svg", bg = "white")

expected_timepoints_ok <- setequal(unique(plot_df$timepoint), TIME_LEVELS)

all_positive <- all(plot_df$moran_I > 0)

fdr_qc_available <- any(is.finite(plot_df$FDR))

all_fdr_sig <- if (fdr_qc_available) {
    all(plot_df$FDR[is.finite(plot_df$FDR)] < 0.05)
} else {
    NA
}

final_checks <- data.frame(check = c("Step66 v02 CSV was used", "Exact section-level Moran input table was used", 
    "All 4 Figure4 timepoints are represented", "Each timepoint has at least one valid section", "Each timepoint-section combination is unique", 
    "All plotted Moran's I values are finite", "Moran's I = 0 reference line is included", "Median trajectory connects timepoint medians only", 
    "Plot background is white with no positive-zone shading", "No between-timepoint significance testing is performed", 
    "No Moran's I values are recalculated", "PDF generated", "PNG generated", "SVG generated"), pass = c(grepl("v02", 
    basename(INPUT_FILE), ignore.case = TRUE), identical(basename(INPUT_FILE), "02_66_Clec7a_section_global_MoransI_v02.csv"), 
    expected_timepoints_ok, all(summary_df$n_sections > 0), nrow(dup_check) == 0L, all(is.finite(plot_df$moran_I)), 
    TRUE, TRUE, TRUE, TRUE, TRUE, file.exists(PDF_FILE), file.exists(PNG_FILE), file.exists(SVG_FILE)), 
    stringsAsFactors = FALSE)

QC_FILE <- file.path(OUT_DIR, paste0("03_72_Figure4D_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_checks, QC_FILE)

REPORT_FILE <- file.path(OUT_DIR, paste0("04_72_Figure4D_REPORT_", VERSION, ".txt"))

summary_lines <- apply(summary_df, 1, function(z) {
    paste0("  ", z[["timepoint_label"]], ": median Moran's I=", sprintf("%.4f", as.numeric(z[["median_moran_I"]])), 
        "; mean=", sprintf("%.4f", as.numeric(z[["mean_moran_I"]])), "; range=", sprintf("%.4f", as.numeric(z[["min_moran_I"]])), 
        " to ", sprintf("%.4f", as.numeric(z[["max_moran_I"]])), "; n_sections=", z[["n_sections"]])
})

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(summary_df)
print(final_checks)

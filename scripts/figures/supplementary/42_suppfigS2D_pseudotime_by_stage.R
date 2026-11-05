source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)

required_pkgs <- c("monocle3", "SingleCellExperiment", "ggplot2", "dplyr", "ggridges", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(monocle3)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(ggridges)
})

ROOT <- PROJECT_ROOT

STEP25_DIR <- file.path(ROOT, "results", "25_Monocle3_all_myeloid_trajectory")

CDS_FILE <- file.path(STEP25_DIR, "25_Monocle3_all_myeloid_trajectory_3scRNA_v03.rds")

OUT_DIR <- file.path(ROOT, "results", "42_SuppFigure2D_pseudotime_ridge_by_time")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "42_SuppFigure2D_pseudotime_ridge_by_time_v02")

if (!file.exists(CDS_FILE)) {
    stop("Cannot find Step25 CDS:\n", CDS_FILE)
}

TIME_LEVELS <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

TIME_COLORS <- c(Control = "#7A6AAE", `1 dpi` = "#D95F59", `3 dpi` = "#E89B4A", `7 dpi` = "#E9C95D", 
    `14 dpi` = "#78B7A4", `28 dpi` = "#4C78A8")

cds <- readRDS(CDS_FILE)

if (!inherits(cds, "cell_data_set")) {
    stop("Input RDS is not a Monocle3 cell_data_set.")
}

if (!"timepoint_main" %in% colnames(colData(cds))) {
    stop("CDS does not contain colData(cds)$timepoint_main.")
}

pt <- monocle3::pseudotime(cds)

if (is.null(names(pt))) {
  names(pt) <- colnames(cds)
}

meta <- as.data.frame(colData(cds))

meta$cell <- rownames(meta)

plot_df <- data.frame(cell = colnames(cds), pseudotime = as.numeric(pt[colnames(cds)]), stringsAsFactors = FALSE) %>% 
    left_join(meta %>% select(cell, timepoint_main), by = "cell") %>% filter(is.finite(pseudotime), timepoint_main %in% 
    TIME_LEVELS) %>% mutate(timepoint_main = factor(timepoint_main, levels = rev(TIME_LEVELS)))

if (nrow(plot_df) == 0) {
    stop("No cells with finite pseudotime and valid main time point were found.")
}

summary_df <- plot_df %>% mutate(timepoint_main = factor(as.character(timepoint_main), levels = TIME_LEVELS)) %>% 
    group_by(timepoint_main) %>% summarise(n_cells = n(), median_pseudotime = median(pseudotime, na.rm = TRUE), 
    mean_pseudotime = mean(pseudotime, na.rm = TRUE), q25_pseudotime = quantile(pseudotime, 0.25, na.rm = TRUE), 
    q75_pseudotime = quantile(pseudotime, 0.75, na.rm = TRUE), .groups = "drop")

write.csv(summary_df, file.path(OUT_DIR, "42_SuppFigure2D_pseudotime_summary_by_time_v02.csv"), row.names = FALSE)

print(summary_df)

p <- ggplot(
  plot_df,
  aes(
    x = pseudotime,
    y = timepoint_main,
    fill = timepoint_main
  )
) +

  ggridges::geom_density_ridges(
    scale = 1.22,
    rel_min_height = 0.035,
    linewidth = 0.45,
    color = "#444444",
    alpha = 0.78
  ) +

  geom_segment(
    data = summary_df %>%
      mutate(
        timepoint_main = factor(
          as.character(timepoint_main),
          levels = rev(TIME_LEVELS)
        ),
        y_num = as.numeric(
          timepoint_main
        )
      ),
    aes(
      x = median_pseudotime,
      xend = median_pseudotime,
      y = y_num - 0.04,
      yend = y_num + 0.34
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +

  scale_fill_manual(
    values = TIME_COLORS,
    guide = "none"
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +

  labs(
    x = "Monocle3 pseudotime",
    y = NULL
  ) +

  theme_classic(
    base_family = "Arial",
    base_size = 9
  ) +

  theme(
    text = element_text(
      family = "Arial",
      color = "black"
    ),

    axis.title.x = element_text(
      family = "Arial",
      size = 9.2,
      margin = margin(
        t = 5
      )
    ),

    axis.text.x = element_text(
      family = "Arial",
      size = 7.8,
      color = "black"
    ),

    axis.text.y = element_text(
      family = "Arial",
      size = 8.2,
      color = "black"
    ),

    axis.line = element_line(
      linewidth = 0.48,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.42,
      color = "black"
    ),

    panel.grid = element_blank(),

    plot.margin = margin(
      6, 8, 6, 6
    )
  )

W_IN <- 6.8
H_IN <- 4.8

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

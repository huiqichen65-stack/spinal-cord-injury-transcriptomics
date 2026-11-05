source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)

required_pkgs <- c("ggplot2", "dplyr", "gghalves", "svglite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gghalves)
})

set.seed(1234)

ROOT <- PROJECT_ROOT

AUCELL_DIR <- file.path(ROOT, "results", "20_AUCell_lipid_metabolism_scoring")

CELL_SCORE_FILE <- file.path(AUCELL_DIR, "20_AUCell_cell_scores_v03.csv")

OUT_DIR <- file.path(ROOT, "results", "34_Figure2G_Clec7a_AUCell")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(OUT_DIR, "34_Figure2G_Clec7a_AUCell_v12")

if (!file.exists(CELL_SCORE_FILE)) {
    stop("Cannot find fixed AUCell result file:\n", CELL_SCORE_FILE)
}

state_levels_metadata <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

display_name_map <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb+ macrophages", 
    `Acute-response Macrophage` = "Acute-response Macrophage")

state_cols_metadata <- c(`Homeostatic Microglia` = "#4C78A8", `Clec7a+ lipid-associated Microglia` = "#D64F4F", 
    `Activated inflammatory Microglia` = "#F2B134", `Ccl12+ acute Microglia` = "#8E6BBE", `Tnf+ inflammatory Microglia` = "#E77AA8", 
    `Proliferating Microglia` = "#4FA6A2", `MHC-II+ antigen-presenting Macrophage` = "#59A14F", `Gpnmb+ tissue-repair Macrophage` = "#D5A832", 
    `Acute-response Macrophage` = "#9C755F")

state_levels_display <- unname(display_name_map[state_levels_metadata])

display_label_wrapped <- c(`Homeostatic Microglia` = "Homeostatic\nMicroglia", `Clec7a+ lipid-handling Microglia` = "Clec7a+ lipid-handling\nMicroglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory\nMicroglia", `Ccl12+ acute Microglia` = "Ccl12+ acute\nMicroglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory\nMicroglia", `Proliferating Microglia` = "Proliferating\nMicroglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting\nMacrophage", `Gpnmb+ macrophages` = "Gpnmb+\nmacrophages", 
    `Acute-response Macrophage` = "Acute-response\nMacrophage")

state_cols_display <- setNames(unname(state_cols_metadata[state_levels_metadata]), state_levels_display)

target_metadata <- "Clec7a+ lipid-associated Microglia"
target_display  <- "Clec7a+ lipid-handling Microglia"

df_raw <- read.csv(CELL_SCORE_FILE, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("myeloid_state", "AUC_FA_uptake_transport")

missing_cols <- setdiff(required_cols, colnames(df_raw))

if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), "\nAvailable columns:\n", 
        paste(colnames(df_raw), collapse = ", "))
}

df <- df_raw %>% transmute(state_metadata = as.character(myeloid_state), score = as.numeric(AUC_FA_uptake_transport)) %>% 
    filter(!is.na(state_metadata), nzchar(state_metadata), is.finite(score))

df <- df %>% filter(state_metadata %in% state_levels_metadata)

if (nrow(df) == 0) {
  stop("No usable AUCell rows remained.")
}

df$state_display <- unname(display_name_map[df$state_metadata])

state_summary <- df %>% group_by(state_metadata, state_display) %>% summarise(n_cells = n(), mean_score = mean(score, 
    na.rm = TRUE), median_score = median(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), q1 = quantile(score, 
    0.25, na.rm = TRUE), q3 = quantile(score, 0.75, na.rm = TRUE), .groups = "drop")

plot_order <- state_summary %>% arrange(desc(median_score), desc(mean_score)) %>% pull(state_display)

state_summary <- state_summary %>% mutate(plot_order = match(state_display, plot_order)) %>% arrange(plot_order)

write.csv(state_summary, paste0(OUT_BASE, "_state_summary.csv"), row.names = FALSE)

target_scores <- df$score[df$state_metadata == target_metadata]

other_scores <- df$score[df$state_metadata != target_metadata]

wilcox_p <- wilcox.test(target_scores, other_scores)$p.value

stats_lines <- c("============================================================", "Figure 2G — FA uptake/transport AUCell raincloud", 
    "============================================================", paste0("Input: ", CELL_SCORE_FILE), 
    paste0("Metadata target: ", target_metadata), paste0("Figure display target: ", target_display), 
    "", paste0("Target cells: ", length(target_scores)), paste0("Other cells: ", length(other_scores)), 
    paste0("Target median: ", signif(median(target_scores), 5)), paste0("Other median: ", signif(median(other_scores), 
        5)), paste0("Wilcoxon P: ", formatC(wilcox_p, format = "e", digits = 3)), "", "Plot order by median score:", 
    paste(plot_order, collapse = " > "))

writeLines(stats_lines, paste0(OUT_BASE, "_stats.txt"), useBytes = TRUE)

df$state_display <- factor(df$state_display, levels = plot_order)

state_index <- setNames(seq_along(plot_order), plot_order)

df$x_center <- unname(state_index[as.character(df$state_display)])

MAX_POINTS_PER_STATE <- 450L

df_points <- df %>% group_by(state_display) %>% group_modify(~{
    n_take <- min(MAX_POINTS_PER_STATE, nrow(.x))
    dplyr::slice_sample(.x, n = n_take, replace = FALSE)
}) %>% ungroup()

df_points$x_point <- df_points$x_center + runif(nrow(df_points), min = -0.38, max = -0.1)

df$x_plot <- df$x_center

p <- ggplot() +

  gghalves::geom_half_violin(
    data = df,
    aes(
      x = x_plot,
      y = score,
      group = state_display,
      fill = state_display
    ),
    side = "r",
    width = 0.78,
    scale = "width",
    trim = TRUE,
    alpha = 0.38,
    color = NA
  ) +

  geom_point(
    data = df_points,
    aes(
      x = x_point,
      y = score,
      color = state_display
    ),
    size = 0.70,
    alpha = 0.56,
    shape = 16
  ) +

  geom_boxplot(
    data = df,
    aes(
      x = x_plot,
      y = score,
      group = state_display,
      fill = state_display
    ),
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.78,
    color = "#333333",
    linewidth = 0.42
  ) +

  stat_summary(
    data = df,
    aes(
      x = x_plot,
      y = score,
      group = state_display
    ),
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.55,
    stroke = 0.35,
    fill = "white",
    color = "#222222"
  ) +

  scale_fill_manual(
    values = state_cols_display,
    drop = FALSE
  ) +

  scale_color_manual(
    values = state_cols_display,
    drop = FALSE
  ) +

  scale_x_continuous(
    breaks = seq_along(
      plot_order
    ),
    labels = unname(
      display_label_wrapped[plot_order]
    ),
    expand = expansion(
      mult = c(
        0.025,
        0.08
      )
    )
  ) +

  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(
      mult = c(
        0.00,
        0.10
      )
    )
  ) +

  labs(
    x = NULL,
    y = "AUCell lipid metabolism score"
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  theme_classic(
    base_family = "Arial",
    base_size = 12
  ) +

  theme(
    text = element_text(
      family = "Arial",
      color = "black"
    ),

    axis.title.y = element_text(
      size = 13.2,
      margin = margin(
        r = 7
      )
    ),

    axis.text.y = element_text(
      size = 11.5,
      color = "black"
    ),

    axis.text.x = element_text(
      size = 11.1,
      color = "black",
      angle = 22,
      hjust = 1,
      vjust = 1,
      lineheight = 0.95,
      margin = margin(t = 6)
    ),

    axis.line = element_line(
      linewidth = 0.55,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),

    panel.grid.major.y = element_line(
      color = "#D9D9D9",
      linewidth = 0.35,
      linetype = "dashed"
    ),

    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    legend.position = "none",

    plot.margin = margin(
      8, 22, 78, 32
    )
  )

W_IN <- 11.2
H_IN <- 6.4

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

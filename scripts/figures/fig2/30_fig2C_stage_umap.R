options(stringsAsFactors = FALSE)

required_pkgs <- c("Seurat", "ggplot2", "dplyr", "scales", "svglite", "ggrastr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggrastr)
})

ROOT <- "D:/zz"
RES_DIR <- file.path(ROOT, "results")
STEP_DIR <- file.path(RES_DIR, "30_Figure2C_myeloid_timepoint_UMAP")
dir.create(STEP_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(STEP_DIR, "30_Figure2C_myeloid_timepoint_UMAP_v12")

POINT_SIZE_RASTER <- 0.32
POINT_ALPHA_RASTER <- 1
RASTER_DPI <- 600

INPUT_RDS <- NA_character_

EXPECTED_CELLS <- 58533L

has_required_umap <- function(x) {
  inherits(x, "Seurat") && "umap_myeloid_harmony30" %in% Reductions(x)
}

if (is.na(INPUT_RDS)) {
  rds_files <- list.files(RES_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

  priority <- c("18_myeloid_annotation_cleanup", "19_myeloid_state_robustness_temporal", "20_AUCell_lipid_metabolism_scoring",
    "21_Clec7a_DEG", "24_CytoTRACE2_myeloid", "25_Monocle3_all_myeloid_trajectory")

  score_file <- function(f) {
    nm <- gsub("\\\\", "/", f)
    hit <- which(vapply(priority, function(p) grepl(p, nm, fixed = TRUE), logical(1)))
    if (length(hit) == 0)
      999L else hit[1]
  }

  ord <- order(vapply(rds_files, score_file, integer(1)))
  candidates <- rds_files[ord]

  obj <- NULL
  chosen_rds <- NULL

  for (f in candidates) {
    tmp <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(tmp) && has_required_umap(tmp) && ncol(tmp) == EXPECTED_CELLS) {
      obj <- tmp
      chosen_rds <- f
      rm(tmp)
      break
    }
    rm(tmp)
    invisible(gc())
  }

  if (is.null(obj)) {
    stop("Could not find the 58,533-cell myeloid Seurat object with ", "'umap_myeloid_harmony30' under D:/zz/results.\n",
      "Do NOT use the 59,409-cell Step 15 object for Figure 2.\n", "Please set INPUT_RDS manually to the Step 18/19 object.")
  }
} else {
  if (!file.exists(INPUT_RDS)) {
    stop("INPUT_RDS does not exist: ", INPUT_RDS)
  }
  obj <- readRDS(INPUT_RDS)
  chosen_rds <- INPUT_RDS

  if (!has_required_umap(obj)) {
    stop("The selected RDS is not a Seurat object with reduction ", "'umap_myeloid_harmony30'.")
  }

  if (ncol(obj) != EXPECTED_CELLS) {
    stop("Wrong myeloid object: selected object has ", ncol(obj), " cells, but Figure 2 object must have ",
      EXPECTED_CELLS, " cells.")
  }
}

um <- Embeddings(obj, reduction = "umap_myeloid_harmony30")

if (ncol(um) < 2) {
  stop("UMAP reduction has fewer than 2 dimensions.")
}

plot_df <- data.frame(cell = rownames(um), UMAP_1 = um[, 1], UMAP_2 = um[, 2], stringsAsFactors = FALSE)

meta <- obj@meta.data
meta$cell <- rownames(meta)

all_time_levels <- c("Control", "4h", "1dpi", "3dpi", "7dpi", "14dpi", "28dpi", "38dpi")

display_levels <- c("Control", "1dpi", "3dpi", "7dpi", "14dpi", "28dpi")

normalize_time <- function(x) {
  z <- trimws(as.character(x))
  z2 <- tolower(z)
  z2 <- gsub("[ _.-]", "", z2)

  out <- rep(NA_character_, length(z2))

  out[grepl("^(control|ctrl|sham|naive|normal|vehiclecontrol|uninjured)$", z2)] <- "Control"

  out[grepl("^(4h|4hr|4hrs|4hour|4hours)$", z2)] <- "4h"

  day_map <- c(`1` = "1dpi", `3` = "3dpi", `7` = "7dpi", `14` = "14dpi", `28` = "28dpi", `38` = "38dpi")

  for (d in names(day_map)) {
    idx <- grepl(paste0("^(", d, "d|d", d, "|", d, "dpi|day", d, "|days", d, "|dpi", d, ")$"), z2)
    out[idx] <- day_map[[d]]
  }

  out
}

candidate_cols <- c("timepoint", "time_point", "time", "Time", "dpi", "group", "Group", "condition",
  "Condition", "orig.ident", "sample", "sample_id", "Sample")
candidate_cols <- intersect(candidate_cols, colnames(meta))

best_col <- NULL
best_n <- -1
best_time <- NULL

for (cc in candidate_cols) {
  tmp <- normalize_time(meta[[cc]])
  n_ok <- sum(!is.na(tmp))
  n_levels <- length(unique(tmp[!is.na(tmp)]))

  score <- n_ok + 1e+05 * n_levels
  if (score > best_n) {
    best_n <- score
    best_col <- cc
    best_time <- tmp
  }
}

if (is.null(best_col) || sum(!is.na(best_time)) == 0) {
    for (cc in colnames(meta)) {
        tmp <- normalize_time(meta[[cc]])
        n_ok <- sum(!is.na(tmp))
        n_levels <- length(unique(tmp[!is.na(tmp)]))
        score <- n_ok + 1e+05 * n_levels
        if (score > best_n) {
            best_n <- score
            best_col <- cc
            best_time <- tmp
        }
    }
}

if (is.null(best_col) || sum(!is.na(best_time)) == 0) {
  stop("Could not identify the time-point metadata column.\n", "Expected labels equivalent to Control, 1dpi, 3dpi, 7dpi, 14dpi, 28dpi.\n",
    "Please inspect colnames(obj@meta.data) and set the time-point column manually.")
}

meta$Figure2C_timepoint <- factor(best_time, levels = all_time_levels)

print(table(meta[[best_col]], useNA = "ifany"))
print(table(meta$Figure2C_timepoint, useNA = "ifany"))

plot_df <- plot_df %>% left_join(meta[, c("cell", "Figure2C_timepoint"), drop = FALSE], by = "cell") %>% 
    filter(!is.na(Figure2C_timepoint))

if (nrow(plot_df) == 0) {
  stop("No cells remained after time-point parsing.")
}

print(display_levels)

plot_df <- plot_df %>% filter(Figure2C_timepoint %in% display_levels) %>% mutate(Figure2C_timepoint = factor(as.character(Figure2C_timepoint), 
    levels = display_levels))

if (nrow(plot_df) == 0) {
  stop("No cells remained after restricting to display_levels.")
}

time_cols <- c(Control = "#A7A9AC", `1dpi` = "#5DA5A4", `3dpi` = "#8FC17A", `7dpi` = "#E6C84F", `14dpi` = "#F29E4C",
  `28dpi` = "#E36A4A")

x_rng <- range(plot_df$UMAP_1, na.rm = TRUE)
y_rng <- range(plot_df$UMAP_2, na.rm = TRUE)

p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) + ggrastr::geom_point_rast(aes(color = Figure2C_timepoint),
  size = POINT_SIZE_RASTER, alpha = POINT_ALPHA_RASTER, stroke = 0, raster.dpi = RASTER_DPI) + scale_color_manual(values = time_cols,
  breaks = display_levels, drop = FALSE, name = NULL) + coord_equal(expand = FALSE) + theme_void(base_family = "Arial") +
  theme(text = element_text(family = "Arial", color = "black"), legend.position = "right", legend.justification = "center",
    legend.text = element_text(size = 14, family = "Arial"), legend.key.height = grid::unit(0.42,
      "cm"), legend.key.width = grid::unit(0.42, "cm"), plot.margin = margin(4, 8, 4, 4)) + guides(color = guide_legend(override.aes = list(size = 3,
  alpha = 1), ncol = 1))

dx <- diff(x_rng)
dy <- diff(y_rng)

x0 <- x_rng[1] + 0.035 * dx
y0 <- y_rng[1] + 0.045 * dy
axis_len_x <- 0.115 * dx
axis_len_y <- 0.115 * dy

p <- p + annotate("segment", x = x0, xend = x0 + axis_len_x, y = y0, yend = y0, linewidth = 0.55, color = "black",
  lineend = "square", arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed")) + annotate("segment",
  x = x0, xend = x0, y = y0, yend = y0 + axis_len_y, linewidth = 0.55, color = "black", lineend = "square",
  arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed")) + annotate("text", x = x0 +
  axis_len_x/2, y = y0 - 0.028 * dy, label = "UMAP 1", family = "Arial", size = 3.65, color = "black") +
  annotate("text", x = x0 - 0.022 * dx, y = y0 + axis_len_y/2, label = "UMAP 2", family = "Arial",
    size = 3.65, angle = 90, color = "black")

W_IN <- 7.7
H_IN <- 4.9

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite,
  bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf,
  bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600,
  bg = "white")

ggsave(filename = paste0(OUT_BASE, ".tif"), plot = p, device = "tiff", width = W_IN, height = H_IN, units = "in",
  dpi = 600, compression = "lzw", bg = "white")

audit <- plot_df %>% count(Figure2C_timepoint, name = "n_cells") %>% mutate(percent = 100 * n_cells/sum(n_cells))

write.csv(audit, paste0(OUT_BASE, "_cell_counts.csv"), row.names = FALSE)

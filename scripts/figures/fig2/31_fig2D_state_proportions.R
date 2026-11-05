options(stringsAsFactors = FALSE)

required_pkgs <- c("Seurat", "ggplot2", "dplyr", "tidyr", "scales", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

ROOT <- "D:/zz"
RES_DIR <- file.path(ROOT, "results")
STEP_DIR <- file.path(RES_DIR, "31_Figure2D_myeloid_state_temporal_proportions")
dir.create(STEP_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(STEP_DIR, "31_Figure2D_myeloid_state_temporal_proportions_v04")

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
            999L
        else hit[1]
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
        stop("Could not find the 58,533-cell myeloid Seurat object.\n", "Do NOT use the 59,409-cell Step 15 object.\n", 
            "Please set INPUT_RDS manually to the Step 18/19 object.")
    }
} else {
    if (!file.exists(INPUT_RDS)) {
        stop("INPUT_RDS does not exist: ", INPUT_RDS)
    }
    obj <- readRDS(INPUT_RDS)
    chosen_rds <- INPUT_RDS
    if (!has_required_umap(obj)) {
        stop("Selected RDS does not contain reduction 'umap_myeloid_harmony30'.")
    }
    if (ncol(obj) != EXPECTED_CELLS) {
        stop("Wrong object: ", ncol(obj), " cells found; Figure 2 object must contain ", EXPECTED_CELLS, 
            " cells.")
    }
}

meta <- obj@meta.data
meta$cell <- rownames(meta)

all_time_levels <- c("Control", "4h", "1dpi", "3dpi", "7dpi", "14dpi", "28dpi", "38dpi")

display_levels <- c("Control", "1dpi", "3dpi", "7dpi", "14dpi", "28dpi")

normalize_time <- function(x) {
  z <- trimws(as.character(x))
  z2 <- tolower(z)
  z2 <- gsub("[ _.-]", "", z2)

  out <- rep(NA_character_, length(z2))

  out[grepl(
    "^(control|ctrl|sham|naive|normal|vehiclecontrol|uninjured)$",
    z2
  )] <- "Control"

  out[grepl(
    "^(4h|4hr|4hrs|4hour|4hours)$",
    z2
  )] <- "4h"

  day_map <- c(
    "1" = "1dpi",
    "3" = "3dpi",
    "7" = "7dpi",
    "14" = "14dpi",
    "28" = "28dpi",
    "38" = "38dpi"
  )

  for (d in names(day_map)) {
    idx <- grepl(
      paste0(
        "^(",
        d, "d|d", d, "|", d, "dpi|day", d,
        "|days", d, "|dpi", d,
        ")$"
      ),
      z2
    )
    out[idx] <- day_map[[d]]
  }

  out
}

candidate_cols <- c("timepoint", "time_point", "time", "Time", "dpi", "group", "Group", "condition", 
    "Condition", "orig.ident", "sample", "sample_id", "Sample")
candidate_cols <- intersect(candidate_cols, colnames(meta))

best_col <- NULL
best_score <- -Inf
best_time <- NULL

for (cc in candidate_cols) {
    tmp <- normalize_time(meta[[cc]])
    n_ok <- sum(!is.na(tmp))
    n_levels <- length(unique(tmp[!is.na(tmp)]))
    score <- n_ok + 1e+05 * n_levels
    if (score > best_score) {
        best_score <- score
        best_col <- cc
        best_time <- tmp
    }
}

if (is.null(best_col) || sum(!is.na(best_time)) == 0) {
  stop("Could not identify the time-point metadata column.")
}

meta$Figure2D_timepoint <- factor(best_time, levels = all_time_levels)

print(table(meta[[best_col]], useNA = "ifany"))
print(table(meta$Figure2D_timepoint, useNA = "ifany"))

expected_states_metadata <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", 
    "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage", 
    "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")

candidate_state_cols <- c("myeloid_state", "myeloid_states", "myeloid_annotation", "myeloid_celltype", 
    "celltype", "cell_type", "annotation")
candidate_state_cols <- intersect(candidate_state_cols, colnames(meta))

best_state_col <- NULL
best_match <- -1

for (cc in candidate_state_cols) {
    vals <- unique(as.character(meta[[cc]]))
    n_match <- sum(expected_states_metadata %in% vals)
    if (n_match > best_match) {
        best_match <- n_match
        best_state_col <- cc
    }
}

if (is.null(best_state_col) || best_match < 7) {
    for (cc in colnames(meta)) {
        vals <- unique(as.character(meta[[cc]]))
        n_match <- sum(expected_states_metadata %in% vals)
        if (n_match > best_match) {
            best_match <- n_match
            best_state_col <- cc
        }
    }
}

if (is.null(best_state_col) || best_match < 7) {
    stop("Could not confidently identify the 9-state annotation column.\n", "Please inspect colnames(obj@meta.data) and set it manually.")
}

meta$Figure2D_state_metadata <- as.character(meta[[best_state_col]])

meta2 <- meta %>% filter(Figure2D_state_metadata %in% expected_states_metadata, Figure2D_timepoint %in% 
    display_levels)

if (nrow(meta2) == 0) {
 stop("No cells remained after filtering to states and display time points.")
}

display_name_map <- c(`Homeostatic Microglia` = "Homeostatic Microglia", `Clec7a+ lipid-associated Microglia` = "Clec7a+ lipid-handling Microglia", 
    `Activated inflammatory Microglia` = "Activated inflammatory Microglia", `Ccl12+ acute Microglia` = "Ccl12+ acute Microglia", 
    `Tnf+ inflammatory Microglia` = "Tnf+ inflammatory Microglia", `Proliferating Microglia` = "Proliferating Microglia", 
    `MHC-II+ antigen-presenting Macrophage` = "MHC-II+ antigen-presenting Macrophage", `Gpnmb+ tissue-repair Macrophage` = "Gpnmb⁺ macrophages", 
    `Acute-response Macrophage` = "Acute-response Macrophage")

state_display_levels <- unname(display_name_map[expected_states_metadata])

meta2$Figure2D_state <- unname(display_name_map[meta2$Figure2D_state_metadata])

meta2$Figure2D_state <- factor(meta2$Figure2D_state, levels = state_display_levels)

meta2$Figure2D_timepoint <- factor(as.character(meta2$Figure2D_timepoint), levels = display_levels)

prop_df <- meta2 %>% count(Figure2D_timepoint, Figure2D_state, name = "n_cells") %>% tidyr::complete(Figure2D_timepoint = factor(display_levels, 
    levels = display_levels), Figure2D_state = factor(state_display_levels, levels = state_display_levels), 
    fill = list(n_cells = 0)) %>% group_by(Figure2D_timepoint) %>% mutate(total_cells = sum(n_cells), 
    proportion = n_cells/total_cells, percent = 100 * proportion) %>% ungroup()

print(prop_df %>% filter(Figure2D_state == "Clec7a+ lipid-handling Microglia") %>% select(Figure2D_timepoint, 
    n_cells, total_cells, percent))

state_cols <- c(`Homeostatic Microglia` = "#4C78A8", `Clec7a+ lipid-handling Microglia` = "#D64F4F", 
    `Activated inflammatory Microglia` = "#F2B134", `Ccl12+ acute Microglia` = "#8E6BBE", `Tnf+ inflammatory Microglia` = "#E77AA8", 
    `Proliferating Microglia` = "#4FA6A2", `MHC-II+ antigen-presenting Macrophage` = "#59A14F", `Gpnmb⁺ macrophages` = "#D5A832", 
    `Acute-response Macrophage` = "#9C755F")

p <- ggplot(prop_df, aes(x = Figure2D_timepoint, y = proportion, fill = Figure2D_state)) + geom_col(width = 0.7, 
    linewidth = 0.18, color = "white") + scale_fill_manual(values = state_cols, breaks = state_display_levels, 
    drop = FALSE, name = NULL) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) + 
    labs(x = NULL, y = "Proportion of myeloid cells") + theme_classic(base_family = "Arial", base_size = 13) + 
    theme(text = element_text(family = "Arial", color = "black"), axis.title.y = element_text(size = 14, 
        margin = margin(r = 7)), axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 12.4, 
        color = "black"), axis.line = element_line(linewidth = 0.55, color = "black"), axis.ticks = element_line(linewidth = 0.45, 
        color = "black"), legend.position = "right", legend.text = element_text(size = 13, family = "Arial"), 
        legend.key.height = grid::unit(0.42, "cm"), legend.key.width = grid::unit(0.42, "cm"), plot.margin = margin(5, 
            7, 5, 5)) + guides(fill = guide_legend(ncol = 1, byrow = TRUE))

W_IN <- 8.0
H_IN <- 5.2

ggsave(filename = paste0(OUT_BASE, ".svg"), plot = p, width = W_IN, height = H_IN, units = "in", device = svglite::svglite, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".pdf"), plot = p, width = W_IN, height = H_IN, units = "in", device = cairo_pdf, 
    bg = "white")

ggsave(filename = paste0(OUT_BASE, ".png"), plot = p, width = W_IN, height = H_IN, units = "in", dpi = 600, 
    bg = "white")

write.csv(prop_df, paste0(OUT_BASE, "_proportions.csv"), row.names = FALSE)

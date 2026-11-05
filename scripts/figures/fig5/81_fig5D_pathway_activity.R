rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(grid)
})

ROOT <- "D:/zz"

STEP77_DIR <- file.path(ROOT, "results", "77_Global_CellChat_Clec7a_centered")

IN_CELLCHAT <- file.path(STEP77_DIR, "04_77_Global_CellChat_object_v05.rds")

IN_OUT_PATHWAY <- file.path(STEP77_DIR, "09_77_Clec7a_outgoing_pathways_v05.csv")

IN_IN_PATHWAY <- file.path(STEP77_DIR, "10_77_Clec7a_incoming_pathways_v05.csv")

OUT_DIR <- file.path(ROOT, "results", "81_Figure5D_ReferenceFig5C_signaling_role_heatmap")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PDF <- file.path(OUT_DIR, "81_Figure5D_ReferenceFig5C_signaling_role_heatmap_v02.pdf")

OUT_PNG <- file.path(OUT_DIR, "81_Figure5D_ReferenceFig5C_signaling_role_heatmap_v02.png")

OUT_SVG <- file.path(OUT_DIR, "81_Figure5D_ReferenceFig5C_signaling_role_heatmap_v02.svg")

OUT_SELECTED <- file.path(OUT_DIR, "81_Figure5D_selected_pathways_v02.csv")

OUT_QC <- file.path(OUT_DIR, "81_Figure5D_FINAL_QC_v02.csv")

TARGET_GROUP <- "Clec7a+ lipid-handling Microglia"
DISPLAY_TARGET <- "Clec7a+ lipid-handling MG"

TOP_N_OUT <- 8
TOP_N_IN <- 8

required_files <- c(IN_CELLCHAT, IN_OUT_PATHWAY, IN_IN_PATHWAY)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

cellchat <- readRDS(IN_CELLCHAT)

if (!inherits(cellchat, "CellChat")) {
  stop("Input object is not a CellChat object.")
}

if (is.null(cellchat@netP$centr) || length(cellchat@netP$centr) == 0) {
  stop("cellchat@netP$centr is absent. ", "This plotting script will NOT recompute centrality.")
}

if (!(TARGET_GROUP %in% levels(cellchat@idents))) {
  stop("Target group not found in CellChat object: ", TARGET_GROUP)
}

out_tbl <- readr::read_csv(IN_OUT_PATHWAY, show_col_types = FALSE)

in_tbl <- readr::read_csv(IN_IN_PATHWAY, show_col_types = FALSE)

needed_cols <- c("source", "target", "pathway_name", "prob")

if (!all(needed_cols %in% colnames(out_tbl))) {
  stop("Outgoing pathway table does not contain all required columns: ", paste(needed_cols, collapse = ", "))
}

if (!all(needed_cols %in% colnames(in_tbl))) {
  stop("Incoming pathway table does not contain all required columns: ", paste(needed_cols, collapse = ", "))
}

out_rank <- out_tbl %>% filter(source == TARGET_GROUP, !is.na(prob), prob > 0) %>% group_by(pathway_name) %>% 
    summarise(summed_probability = sum(prob), n_partner_edges = n(), .groups = "drop") %>% arrange(desc(summed_probability), 
    pathway_name) %>% slice_head(n = TOP_N_OUT) %>% mutate(selected_from = "Outgoing")

in_rank <- in_tbl %>% filter(target == TARGET_GROUP, !is.na(prob), prob > 0) %>% group_by(pathway_name) %>% 
    summarise(summed_probability = sum(prob), n_partner_edges = n(), .groups = "drop") %>% arrange(desc(summed_probability), 
    pathway_name) %>% slice_head(n = TOP_N_IN) %>% mutate(selected_from = "Incoming")

pathways_selected <- unique(c(out_rank$pathway_name, in_rank$pathway_name))

available_pathways <- names(cellchat@netP$centr)

pathways_show <- pathways_selected[pathways_selected %in% available_pathways]

missing_selected <- setdiff(pathways_selected, available_pathways)

if (length(missing_selected) > 0) {
  warning("Selected pathway(s) absent from cellchat@netP$centr and omitted: ", paste(missing_selected,
    collapse = ", "))
}

if (length(pathways_show) == 0) {
  stop("No selected pathways are available in cellchat@netP$centr.")
}

selected_export <- bind_rows(out_rank, in_rank) %>% mutate(included_in_heatmap = pathway_name %in% pathways_show)

write_csv(selected_export, OUT_SELECTED)

cellchat_plot <- cellchat

old_levels <- levels(cellchat_plot@idents)

new_levels <- old_levels
new_levels[new_levels == TARGET_GROUP] <- DISPLAY_TARGET

levels(cellchat_plot@idents) <- new_levels

ht_out <- netAnalysis_signalingRole_heatmap(cellchat_plot, signaling = pathways_show, pattern = "outgoing",
  slot.name = "netP", title = "Outgoing interaction strength", font.size = 11, font.size.title = 13,
  cluster.rows = FALSE, cluster.cols = FALSE)

ht_in <- netAnalysis_signalingRole_heatmap(cellchat_plot, signaling = pathways_show, pattern = "incoming",
  slot.name = "netP", title = "Incoming interaction strength", font.size = 11, font.size.title = 13,
  cluster.rows = FALSE, cluster.cols = FALSE)

draw_figure5d <- function() {
    ComplexHeatmap::draw(ht_out + ht_in, ht_gap = grid::unit(5, "mm"), merge_legends = TRUE)
}

grDevices::cairo_pdf(filename = OUT_PDF, width = 14.5, height = 7.2)

draw_figure5d()

dev.off()

png(filename = OUT_PNG, width = 14.5, height = 7.2, units = "in", res = 600, bg = "white")

draw_figure5d()

dev.off()

grDevices::svg(filename = OUT_SVG, width = 14.5, height = 7.2, bg = "white")

draw_figure5d()

dev.off()

qc <- data.frame(item = c("Step77 CellChat object exists", "Outgoing target pathway table exists", "Incoming target pathway table exists",
  "Input is CellChat object", "netP centrality available", "Target group present", "At least one pathway selected",
  "All displayed pathways available in netP centrality", "Display target label changed to MG", "PDF generated",
  "PNG generated", "SVG generated"), pass = c(file.exists(IN_CELLCHAT), file.exists(IN_OUT_PATHWAY),
  file.exists(IN_IN_PATHWAY), inherits(cellchat, "CellChat"), !is.null(cellchat@netP$centr) && length(cellchat@netP$centr) >
    0, TARGET_GROUP %in% old_levels, length(pathways_show) > 0, all(pathways_show %in% available_pathways),
  DISPLAY_TARGET %in% levels(cellchat_plot@idents), file.exists(OUT_PDF), file.exists(OUT_PNG), file.exists(OUT_SVG)),
  stringsAsFactors = FALSE)

write.csv(qc, OUT_QC, row.names = FALSE, fileEncoding = "UTF-8")

print(qc)

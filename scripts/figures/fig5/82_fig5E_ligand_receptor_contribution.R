rm(list = ls())
source("config.R", local = TRUE)


options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
})

ROOT <- PROJECT_ROOT

STEP77_DIR <- file.path(ROOT, "results", "77_Global_CellChat_Clec7a_centered")

IN_CELLCHAT <- file.path(STEP77_DIR, "04_77_Global_CellChat_object_v05.rds")

IN_OUT_PATHWAY <- file.path(STEP77_DIR, "09_77_Clec7a_outgoing_pathways_v05.csv")

IN_IN_PATHWAY <- file.path(STEP77_DIR, "10_77_Clec7a_incoming_pathways_v05.csv")

OUT_DIR <- file.path(ROOT, "results", "82_Figure5E_ReferenceFig5D_LR_contribution")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PDF <- file.path(OUT_DIR, "82_Figure5E_ReferenceFig5D_LR_contribution_v02.pdf")

OUT_PNG <- file.path(OUT_DIR, "82_Figure5E_ReferenceFig5D_LR_contribution_v02.png")

OUT_SVG <- file.path(OUT_DIR, "82_Figure5E_ReferenceFig5D_LR_contribution_v02.svg")

OUT_DATA <- file.path(OUT_DIR, "82_Figure5E_LR_contribution_data_v02.csv")

OUT_PATHWAYS <- file.path(OUT_DIR, "82_Figure5E_selected_pathways_v02.csv")

OUT_QC <- file.path(OUT_DIR, "82_Figure5E_FINAL_QC_v02.csv")

TARGET_GROUP <- "Clec7a+ lipid-handling Microglia"

TOP_N_OUT <- 8
TOP_N_IN  <- 8

LR_P_THRESHOLD <- 0.05

required_files <- c(IN_CELLCHAT, IN_OUT_PATHWAY, IN_IN_PATHWAY)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

cellchat <- readRDS(IN_CELLCHAT)

if (!inherits(cellchat, "CellChat")) {
  stop("Input object is not a CellChat object.")
}

if (!(TARGET_GROUP %in% levels(cellchat@idents))) {
    stop("Target group not found in CellChat object: ", TARGET_GROUP)
}

out_tbl <- readr::read_csv(IN_OUT_PATHWAY, show_col_types = FALSE)

in_tbl <- readr::read_csv(IN_IN_PATHWAY, show_col_types = FALSE)

needed_cols <- c("source", "target", "pathway_name", "prob")

if (!all(needed_cols %in% colnames(out_tbl))) {
    stop("Outgoing pathway table lacks required columns: ", paste(needed_cols, collapse = ", "))
}

if (!all(needed_cols %in% colnames(in_tbl))) {
    stop("Incoming pathway table lacks required columns: ", paste(needed_cols, collapse = ", "))
}

out_rank <- out_tbl %>% filter(source == TARGET_GROUP, !is.na(prob), prob > 0) %>% group_by(pathway_name) %>% 
    summarise(summed_probability = sum(prob), n_partner_edges = n(), .groups = "drop") %>% arrange(desc(summed_probability), 
    pathway_name) %>% slice_head(n = TOP_N_OUT) %>% mutate(selected_from = "Outgoing")

in_rank <- in_tbl %>% filter(target == TARGET_GROUP, !is.na(prob), prob > 0) %>% group_by(pathway_name) %>% 
    summarise(summed_probability = sum(prob), n_partner_edges = n(), .groups = "drop") %>% arrange(desc(summed_probability), 
    pathway_name) %>% slice_head(n = TOP_N_IN) %>% mutate(selected_from = "Incoming")

pathways_show <- unique(c(out_rank$pathway_name, in_rank$pathway_name))

if (length(pathways_show) == 0) {
  stop("No pathways were selected.")
}

write_csv(bind_rows(out_rank, in_rank), OUT_PATHWAYS)

extract_lr_contribution <- function(pathway, direction) {
    cat("  ", direction, " | ", pathway, "\n", sep = "")
    res <- tryCatch({
        if (direction == "Source") {
            netAnalysis_contribution(cellchat, signaling = pathway, sources.use = TARGET_GROUP, thresh = LR_P_THRESHOLD, 
                return.data = TRUE, title = NULL)
        }
        else {
            netAnalysis_contribution(cellchat, signaling = pathway, targets.use = TARGET_GROUP, thresh = LR_P_THRESHOLD, 
                return.data = TRUE, title = NULL)
        }
    }, error = function(e) {
        message("    SKIP: ", direction, " / ", pathway, " -> ", conditionMessage(e))
        return(NULL)
    })
    if (is.null(res) || is.null(res$LR.contribution) || nrow(res$LR.contribution) == 0) {
        return(NULL)
    }
    dat <- as.data.frame(res$LR.contribution)
    if (!all(c("name", "contribution") %in% colnames(dat))) {
        stop("Unexpected netAnalysis_contribution() output for pathway ", pathway, ". Columns found: ", 
            paste(colnames(dat), collapse = ", "))
    }
    dat <- dat %>% transmute(direction = direction, pathway = pathway, lr_pair = as.character(name), 
        contribution = as.numeric(contribution)) %>% filter(is.finite(contribution), contribution > 0)
    if (nrow(dat) == 0) {
        return(NULL)
    }
    dat
}

source_df <- bind_rows(lapply(pathways_show, extract_lr_contribution, direction = "Source"))

target_df <- bind_rows(lapply(pathways_show, extract_lr_contribution, direction = "Target"))

if (nrow(source_df) == 0) {
    stop("No Clec7a source LR-contribution rows were obtained.")
}

if (nrow(target_df) == 0) {
    stop("No Clec7a target LR-contribution rows were obtained.")
}

lr_df <- bind_rows(source_df, target_df)

write_csv(lr_df, OUT_DATA)

source_df <- source_df %>% mutate(pathway = factor(pathway, levels = pathways_show))

target_df <- target_df %>% mutate(pathway = factor(pathway, levels = pathways_show))

prepare_plot_df <- function(dat) {
    dat %>% group_by(pathway) %>% arrange(desc(contribution), lr_pair, .by_group = TRUE) %>% mutate(lr_index = row_number(), 
        x_id = paste0(as.character(pathway), "__", sprintf("%03d", lr_index))) %>% ungroup()
}

source_plot_df <- prepare_plot_df(source_df)

target_plot_df <- prepare_plot_df(target_df)

label_map_source <- setNames(source_plot_df$lr_pair, source_plot_df$x_id)

label_map_target <- setNames(target_plot_df$lr_pair, target_plot_df$x_id)

make_row_plot <- function(dat, label_map, y_title) {
    ggplot(dat, aes(x = x_id, y = contribution, fill = pathway)) + geom_col(width = 0.78) + facet_grid(cols = vars(pathway), 
        scales = "free_x", space = "free_x", drop = TRUE) + scale_x_discrete(labels = label_map) + labs(x = NULL, 
        y = y_title) + theme_classic(base_size = 14) + theme(legend.position = "none", strip.background = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold", colour = "black", margin = margin(b = 2)), 
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 13.5, colour = "black"), 
        axis.text.y = element_text(size = 13, colour = "black"), axis.title.y = element_text(size = 15, 
            face = "bold"), panel.spacing.x = grid::unit(0.8, "mm"), plot.margin = margin(9, 12, 32, 
            14))
}

p_source <- make_row_plot(source_plot_df, label_map_source, "Source Contribution")

p_target <- make_row_plot(target_plot_df, label_map_target, "Target Contribution")

invisible(ggplot_build(p_source))

invisible(ggplot_build(p_target))

p_final <- cowplot::plot_grid(p_source, p_target, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 
    1))

invisible(ggplot_build(p_final))

grDevices::cairo_pdf(filename = OUT_PDF, width = 18, height = 11)

print(p_final)

dev.off()

png(filename = OUT_PNG, width = 18, height = 11, units = "in", res = 600, bg = "white")

print(p_final)

dev.off()

grDevices::svg(filename = OUT_SVG, width = 18, height = 11, bg = "white")

print(p_final)

dev.off()

qc <- data.frame(item = c("Step77 CellChat object exists", "Outgoing pathway table exists", "Incoming pathway table exists", 
    "Input is CellChat object", "Target group present", "Figure5D pathway rule reproduced", "Source LR contribution extracted", 
    "Target LR contribution extracted", "Source ggplot validated", "Target ggplot validated", "Contribution CSV generated", 
    "PDF generated", "PNG generated", "SVG generated"), pass = c(file.exists(IN_CELLCHAT), file.exists(IN_OUT_PATHWAY), 
    file.exists(IN_IN_PATHWAY), inherits(cellchat, "CellChat"), TARGET_GROUP %in% levels(cellchat@idents), 
    length(pathways_show) > 0, nrow(source_df) > 0, nrow(target_df) > 0, TRUE, TRUE, file.exists(OUT_DATA), 
    file.exists(OUT_PDF), file.exists(OUT_PNG), file.exists(OUT_SVG)), stringsAsFactors = FALSE)

write.csv(qc, OUT_QC, row.names = FALSE, fileEncoding = "UTF-8")

print(qc)

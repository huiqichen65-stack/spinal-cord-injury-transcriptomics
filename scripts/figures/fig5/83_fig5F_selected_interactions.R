rm(list = ls())
source("config.R", local = TRUE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(readr)
})

ROOT <- PROJECT_ROOT

STEP77_DIR <- file.path(ROOT, "results", "77_Global_CellChat_Clec7a_centered")

IN_CELLCHAT <- file.path(STEP77_DIR, "04_77_Global_CellChat_object_v05.rds")

OUT_DIR <- file.path(ROOT, "results", "83_Figure5F_Clec7a_selected6_LR_bubble")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PDF <- file.path(OUT_DIR, "83_Figure5F_Clec7a_selected6_LR_bubble_v03.pdf")

OUT_SVG <- file.path(OUT_DIR, "83_Figure5F_Clec7a_selected6_LR_bubble_v03.svg")

OUT_PNG <- file.path(OUT_DIR, "83_Figure5F_Clec7a_selected6_LR_bubble_v03.png")

OUT_LR <- file.path(OUT_DIR, "83_Figure5F_selected6_LR_v03.csv")

OUT_QC <- file.path(OUT_DIR, "83_Figure5F_FINAL_QC_v03.csv")

TARGET_GROUP <- "Clec7a+ lipid-handling Microglia"
DISPLAY_TARGET <- "Clec7a+ lipid-handling MG"

P_THRESHOLD <- 0.05

SELECTED_LR <- c("APOE_TREM2_TYROBP", "PSAP_GPR37L1", "SPP1_CD44", "FN1_CD44", "APP_TREM2_TYROBP", "IL34_CSF1R")

LR_DISPLAY <- c(APOE_TREM2_TYROBP = "Apoe - (Trem2+Tyrobp)", PSAP_GPR37L1 = "Psap - Gpr37l1", SPP1_CD44 = "Spp1 - Cd44", 
    FN1_CD44 = "Fn1 - Cd44", APP_TREM2_TYROBP = "App - (Trem2+Tyrobp)", IL34_CSF1R = "Il34 - Csf1r")

PAIR_LR_USE <- data.frame(interaction_name = SELECTED_LR, stringsAsFactors = FALSE)

write.csv(data.frame(interaction_name = SELECTED_LR, display_name = unname(LR_DISPLAY[SELECTED_LR]), 
    stringsAsFactors = FALSE), OUT_LR, row.names = FALSE, fileEncoding = "UTF-8")

if (!file.exists(IN_CELLCHAT)) {
    stop("CellChat object not found:\n", IN_CELLCHAT)
}

cellchat <- readRDS(IN_CELLCHAT)

if (!inherits(cellchat, "CellChat")) {
  stop("Input object is not a CellChat object.")
}

cell_levels <- levels(cellchat@idents)

if (!(TARGET_GROUP %in% cell_levels)) {
    stop("Target group not found: ", TARGET_GROUP)
}

target_idx <- which(cell_levels == TARGET_GROUP)

available_lr <- dimnames(cellchat@net$prob)[[3]]

missing_lr <- setdiff(SELECTED_LR, available_lr)

if (length(missing_lr) > 0) {
    stop("Selected LR pair(s) missing from cellchat@net$prob:\n", paste(missing_lr, collapse = "\n"))
}

cellchat_plot <- cellchat

prob <- cellchat_plot@net$prob
pval <- cellchat_plot@net$pval

n_group <- dim(prob)[1]

keep_mask <- matrix(FALSE, nrow = n_group, ncol = n_group)

keep_mask[target_idx, ] <- TRUE

keep_mask[, target_idx] <- TRUE

for (lr in SELECTED_LR) {

  m_prob <- prob[
    ,
    ,
    lr
  ]

  m_pval <- pval[
    ,
    ,
    lr
  ]

  m_prob[
    is.na(m_pval) |
      m_pval > P_THRESHOLD
  ] <- 0

  m_prob[
    !keep_mask
  ] <- 0

  m_pval[
    !keep_mask
  ] <- 1

  prob[
    ,
    ,
    lr
  ] <- m_prob

  pval[
    ,
    ,
    lr
  ] <- m_pval
}

cellchat_plot@net$prob <- prob
cellchat_plot@net$pval <- pval

display_levels <- cell_levels

display_levels[display_levels == TARGET_GROUP] <- "Clec7a MG"

display_levels[display_levels == "Other myeloid"] <- "Other myel."

display_levels[display_levels == "Neutrophil"] <- "Neutro."

display_levels[display_levels == "Astrocyte"] <- "Astro."

display_levels[display_levels == "Oligodendrocyte"] <- "Oligo."

display_levels[display_levels == "OPC"] <- "OPC"

display_levels[display_levels == "Neuron"] <- "Neuron"

display_levels[display_levels == "Ependymal"] <- "Epend."

display_levels[display_levels == "Endothelial"] <- "Endo."

display_levels[display_levels == "Pericyte/VSMC"] <- "Peri./VSMC"

display_levels[display_levels == "Fibroblast/Meningeal"] <- "Fibro."

levels(cellchat_plot@idents) <- display_levels

for (slot_name in c("prob", "pval")) {
    arr <- cellchat_plot@net[[slot_name]]
    dn <- dimnames(arr)
    dn[[1]] <- display_levels
    dn[[2]] <- display_levels
    dimnames(arr) <- dn
    cellchat_plot@net[[slot_name]] <- arr
}

p <- CellChat::netVisual_bubble(cellchat_plot, sources.use = seq_along(display_levels), targets.use = seq_along(display_levels), 
    pairLR.use = PAIR_LR_USE, remove.isolate = TRUE, thresh = P_THRESHOLD, angle.x = 55, vjust.x = 1, 
    hjust.x = 1, font.size = 12, font.size.title = 13, dot.size.min = 2.5, dot.size.max = 8.5, title.name = NULL)

p <- p + theme(plot.margin = margin(12, 22, 20, 14), axis.text.x = element_text(size = 11, angle = 55, 
    hjust = 1, vjust = 1, colour = "black"), axis.text.y = element_text(size = 12, face = "bold", colour = "black"), 
    axis.title = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11), 
    legend.position = "right")

FIG_W <- 13.8
FIG_H <- 6.8

ggsave(filename = OUT_PDF, plot = p, width = FIG_W, height = FIG_H, units = "in", device = cairo_pdf, 
    bg = "white", limitsize = FALSE)

ggsave(filename = OUT_SVG, plot = p, width = FIG_W, height = FIG_H, units = "in", device = grDevices::svg, 
    bg = "white", limitsize = FALSE)

ggsave(filename = OUT_PNG, plot = p, width = FIG_W, height = FIG_H, units = "in", dpi = 600, bg = "white", 
    limitsize = FALSE)

selected_positive <- vapply(SELECTED_LR, function(lr) {
    sum(cellchat_plot@net$prob[, , lr], na.rm = TRUE) > 0
}, logical(1))

qc <- data.frame(item = c("CellChat object exists", "Input is CellChat object", "Target group present", 
    "All 6 selected LR pairs present", paste0("Positive Clec7a-centered communication: ", SELECTED_LR), 
    "Bubble plot generated", "PDF generated", "SVG generated", "PNG generated"), pass = c(file.exists(IN_CELLCHAT), 
    inherits(cellchat, "CellChat"), TARGET_GROUP %in% cell_levels, length(missing_lr) == 0, unname(selected_positive), 
    inherits(p, "ggplot"), file.exists(OUT_PDF), file.exists(OUT_SVG), file.exists(OUT_PNG)), stringsAsFactors = FALSE)

write.csv(qc, OUT_QC, row.names = FALSE, fileEncoding = "UTF-8")

print(qc)

rm(list = ls())
source("config.R", local = TRUE)


options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
})

ROOT <- PROJECT_ROOT

IN_FILE <- file.path(ROOT, "results", "77_Global_CellChat_Clec7a_centered", "04_77_Global_CellChat_object_v05.rds")

OUT_DIR <- file.path(ROOT, "results", "80_Figure5C_Clec7a_incoming_outgoing_strength")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PDF <- file.path(OUT_DIR, "80_Figure5C_Clec7a_incoming_outgoing_strength_v02.pdf")

OUT_PNG <- file.path(OUT_DIR, "80_Figure5C_Clec7a_incoming_outgoing_strength_v02.png")

OUT_SVG <- file.path(OUT_DIR, "80_Figure5C_Clec7a_incoming_outgoing_strength_v02.svg")

OUT_QC <- file.path(OUT_DIR, "80_Figure5C_FINAL_QC_v02.csv")

if (!file.exists(IN_FILE)) {
  stop("CellChat object not found:\n", IN_FILE)
}

cellchat <- readRDS(IN_FILE)

if (!inherits(cellchat, "CellChat")) {
  stop("Input object is not a CellChat object.")
}

if (is.null(cellchat@net$weight)) {
  stop("cellchat@net$weight is missing.")
}

TARGET_GROUP <- "Clec7a+ lipid-handling Microglia"
DISPLAY_TARGET <- "Clec7a+ lipid-handling MG"

weight_mat <- cellchat@net$weight

if (!(TARGET_GROUP %in% rownames(weight_mat))) {
  stop("Target group not found in cellchat@net$weight: ", TARGET_GROUP)
}

group_size <- as.numeric(table(cellchat@idents))

names(group_size) <- names(table(cellchat@idents))

out_weight_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))

out_weight_mat[TARGET_GROUP, ] <- weight_mat[TARGET_GROUP, ]

in_weight_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))

in_weight_mat[, TARGET_GROUP] <- weight_mat[, TARGET_GROUP]

replace_plot_label <- function(mat) {
    rn <- rownames(mat)
    cn <- colnames(mat)
    if (!is.null(rn)) {
        rn[rn == TARGET_GROUP] <- DISPLAY_TARGET
    }
    if (!is.null(cn)) {
        cn[cn == TARGET_GROUP] <- DISPLAY_TARGET
    }
    rownames(mat) <- rn
    colnames(mat) <- cn
    mat
}

out_weight_plot <- replace_plot_label(out_weight_mat)

in_weight_plot <- replace_plot_label(in_weight_mat)

draw_figure5c <- function() {
    par(mfrow = c(1, 2), xpd = TRUE, mar = c(1, 1, 3, 1))
    netVisual_circle(out_weight_plot, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, 
        edge.weight.max = max(weight_mat, na.rm = TRUE), title.name = "Clec7a outgoing | interaction strength", 
        vertex.label.cex = 0.95, margin = 0.3)
    netVisual_circle(in_weight_plot, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, 
        edge.weight.max = max(weight_mat, na.rm = TRUE), title.name = "Clec7a incoming | interaction strength", 
        vertex.label.cex = 0.95, margin = 0.3)
}

grDevices::cairo_pdf(filename = OUT_PDF, width = 14, height = 6.5)

draw_figure5c()

dev.off()

png(filename = OUT_PNG, width = 14, height = 6.5, units = "in", res = 600, bg = "white")

draw_figure5c()

dev.off()

grDevices::svg(filename = OUT_SVG, width = 14, height = 6.5, bg = "white")

draw_figure5c()

dev.off()

qc <- data.frame(item = c("Step77 CellChat object exists", "Input is CellChat object", "cellchat@net$weight available",
  "Target group found", "Outgoing target matrix generated", "Incoming target matrix generated", "Displayed target label shortened",
  "Figure5C PDF generated", "Figure5C PNG generated", "Figure5C SVG generated"), pass = c(file.exists(IN_FILE),
  inherits(cellchat, "CellChat"), !is.null(cellchat@net$weight), TARGET_GROUP %in% rownames(weight_mat),
  sum(out_weight_mat, na.rm = TRUE) > 0, sum(in_weight_mat, na.rm = TRUE) > 0, DISPLAY_TARGET %in%
    rownames(out_weight_plot), file.exists(OUT_PDF), file.exists(OUT_PNG), file.exists(OUT_SVG)),
  stringsAsFactors = FALSE)

write.csv(qc, OUT_QC, row.names = FALSE, fileEncoding = "UTF-8")

print(qc)

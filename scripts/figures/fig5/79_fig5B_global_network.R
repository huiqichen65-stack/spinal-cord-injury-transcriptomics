rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
})

ROOT <- PROJECT_ROOT
IN_FILE <- file.path( ROOT, "results", "77_Global_CellChat_Clec7a_centered", "04_77_Global_CellChat_object_v05.rds" )

OUT_DIR <- file.path( ROOT, "results", "79_Figure5B_Global_CellChat_network" )

dir.create( OUT_DIR, recursive = TRUE, showWarnings = FALSE )

OUT_PDF <- file.path( OUT_DIR, "79_Figure5B_Global_interaction_number_strength_v07.pdf" )

OUT_PNG <- file.path( OUT_DIR, "79_Figure5B_Global_interaction_number_strength_v07.png" )

OUT_SVG <- file.path( OUT_DIR, "79_Figure5B_Global_interaction_number_strength_v07.svg" )

OUT_QC <- file.path( OUT_DIR, "79_Figure5B_FINAL_QC_v07.csv" )

if (!file.exists(IN_FILE)) {
  stop( "CellChat object not found:\n", IN_FILE )
}

cellchat <- readRDS( IN_FILE )

if (!inherits(
  cellchat,
  "CellChat"
)) {
  stop( "Input object is not a CellChat object." )
}

if (
  is.null(cellchat@net$count) ||
  is.null(cellchat@net$weight)
) {
  stop( "cellchat@net$count or cellchat@net$weight is missing. ", "Do not rerun Figure 5B until the Step77 final object is complete." )
}

groupSize <- as.numeric( table( cellchat@idents ) )

count_plot <- cellchat@net$count
weight_plot <- cellchat@net$weight

OLD_LABEL <- "Clec7a+ lipid-handling Microglia"
NEW_LABEL <- "Clec7a+ state"

LEGACY_DISPLAY_LABEL <- "Clec7a+ lipid-handling MG"

replace_plot_label <- function(mat, old_label, new_label) {
  rn <- rownames(mat)
  cn <- colnames(mat)

  if (!is.null(rn)) rn[rn == old_label] <- new_label
  if (!is.null(cn)) cn[cn == old_label] <- new_label

  rownames(mat) <- rn
  colnames(mat) <- cn
  mat
}

count_plot <- replace_plot_label(count_plot, OLD_LABEL, NEW_LABEL)
weight_plot <- replace_plot_label(weight_plot, OLD_LABEL, NEW_LABEL)

cat("[1/3] Drawing PDF...\n")

grDevices::cairo_pdf( filename = OUT_PDF, width = 13, height = 6.5 )

par( mfrow = c(1, 2), xpd = TRUE, mar = c(2, 2.5, 3.2, 2.5) )

netVisual_circle( count_plot, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions", vertex.label.cex = 1.02 )

netVisual_circle( weight_plot, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength", vertex.label.cex = 1.02 )

dev.off()
cat("PDF finished.\n")

cat("[2/3] Drawing PNG...\n")

png( filename = OUT_PNG, width = 13, height = 6.5, units = "in", res = 600, bg = "white" )

par( mfrow = c(1, 2), xpd = TRUE, mar = c(2, 2.5, 3.2, 2.5) )

netVisual_circle( count_plot, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions", vertex.label.cex = 1.02 )

netVisual_circle( weight_plot, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength", vertex.label.cex = 1.02 )

dev.off()
cat("PNG finished.\n")

cat("[3/3] Drawing SVG...\n")

grDevices::svg( filename = OUT_SVG, width = 13, height = 6.5, bg = "white" )

par( mfrow = c(1, 2), xpd = TRUE, mar = c(2, 2.5, 3.2, 2.5) )

netVisual_circle( count_plot, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions", vertex.label.cex = 1.02 )

netVisual_circle( weight_plot, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength", vertex.label.cex = 1.02 )

dev.off()
cat("SVG finished.\n")

qc <- data.frame( item = c( "Step77 CellChat object exists", "Input is CellChat object", "cellchat@net$count available", "cellchat@net$weight available",
  "Old label present in original matrix", "New Clec7a+ state label present in plotting matrix", "Legacy MG display label absent from plotting matrix", "Figure5B PDF generated",
  "Figure5B PNG generated", "Figure5B SVG generated" ), pass = c( file.exists(IN_FILE), inherits(cellchat, "CellChat"), !is.null(cellchat@net$count), !is.null(cellchat@net$weight),
  OLD_LABEL %in% rownames(cellchat@net$count), NEW_LABEL %in% rownames(count_plot), !(LEGACY_DISPLAY_LABEL %in% rownames(count_plot)) && !(LEGACY_DISPLAY_LABEL %in% colnames(count_plot)),
  file.exists(OUT_PDF), file.exists(OUT_PNG), file.exists(OUT_SVG) ), stringsAsFactors = FALSE )

write.csv( qc, OUT_QC, row.names = FALSE, fileEncoding = "UTF-8" )

print(qc)

cat( "\n============================================================\n", "Figure 5B v07 Reviewer 2 terminology revision completed.\n",
  "Only label changed: ", OLD_LABEL, " -> ", NEW_LABEL, "\n", "PDF: ", OUT_PDF, "\n", "PNG: ", OUT_PNG, "\n", "SVG: ", OUT_SVG, "\n",
  "============================================================\n", sep = "" )

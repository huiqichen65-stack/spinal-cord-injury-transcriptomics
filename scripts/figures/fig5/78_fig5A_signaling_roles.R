rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

ROOT <- "D:/zz"

IN_FILE <- file.path(ROOT, "results", "77_Global_CellChat_Clec7a_centered", "04_77_Global_CellChat_object_v05.rds")

OUT_DIR <- file.path(ROOT, "results", "78_Figure5A_CellChat_signaling_role_scatter")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PDF <- file.path(OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_v02.pdf")

OUT_PNG <- file.path(OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_v02.png")

OUT_SVG <- file.path(OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_v02.svg")

OUT_QC <- file.path(OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_QC_v02.csv")

if (!file.exists(IN_FILE)) {
  stop("CellChat object not found:\n", IN_FILE)
}

cellchat <- readRDS(IN_FILE)

if (!inherits(cellchat, "CellChat")) {
  stop("Input object is not a CellChat object.")
}

if (is.null(cellchat@netP$centr) || length(cellchat@netP$centr) == 0) {
  stop("CellChat netP centrality is absent. ", "Please check the Step77 final object; this plotting script ",
    "will not rerun centrality automatically.")
}

p <- netAnalysis_signalingRole_scatter(cellchat, signaling = NULL, slot.name = "netP", x.measure = "outdeg",
  y.measure = "indeg", xlabel = "Outgoing interaction strength", ylabel = "Incoming interaction strength",
  title = "Global Intercellular Signaling", dot.size = c(2, 6), label.size = 4, dot.alpha = 0.8, font.size = 12.5,
  font.size.title = 14, do.label = TRUE, show.legend = TRUE, show.axes = TRUE)

p <- p + ggplot2::theme(axis.title = ggplot2::element_text(size = 13), axis.text = ggplot2::element_text(size = 11.5),
  legend.title = ggplot2::element_text(size = 11.5), legend.text = ggplot2::element_text(size = 10.5),
  plot.margin = ggplot2::margin(12, 24, 12, 14))

ggsave(filename = OUT_PDF, plot = p, width = 6.2, height = 5.3, units = "in", device = cairo_pdf, bg = "white")

ggsave(filename = OUT_PNG, plot = p, width = 6.2, height = 5.3, units = "in", dpi = 600, bg = "white")

ggsave(filename = OUT_SVG, plot = p, width = 6.2, height = 5.3, units = "in", device = grDevices::svg,
  bg = "white")

qc <- data.frame(item = c("Step77 CellChat object exists", "Input is CellChat object", "netP centrality available",
  "Figure5A PDF generated", "Figure5A PNG generated", "Figure5A SVG generated"), pass = c(file.exists(IN_FILE),
  inherits(cellchat, "CellChat"), !is.null(cellchat@netP$centr) && length(cellchat@netP$centr) > 0,
  file.exists(OUT_PDF), file.exists(OUT_PNG), file.exists(OUT_SVG)), stringsAsFactors = FALSE)

write.csv(qc, OUT_QC, row.names = FALSE, fileEncoding = "UTF-8")

print(qc)

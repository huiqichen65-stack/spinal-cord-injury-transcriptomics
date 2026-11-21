rm(list = ls())
source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

ROOT <- PROJECT_ROOT
IN_FILE <- file.path( ROOT, "results", "77_Global_CellChat_Clec7a_centered", "04_77_Global_CellChat_object_v05.rds" )

OUT_DIR <- file.path( ROOT, "results", "78_Figure5A_CellChat_signaling_role_scatter" )

dir.create( OUT_DIR, recursive = TRUE, showWarnings = FALSE )

OUT_PDF <- file.path( OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_v05.pdf" )

OUT_PNG <- file.path( OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_v05.png" )

OUT_SVG <- file.path( OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_v05.svg" )

OUT_QC <- file.path( OUT_DIR, "78_Figure5A_CellChat_signaling_role_scatter_QC_v05.csv" )

if (!file.exists(IN_FILE)) {
  stop( "CellChat object not found:\n", IN_FILE )
}

cellchat <- readRDS( IN_FILE )

if (!inherits(cellchat, "CellChat")) {
  stop( "Input object is not a CellChat object." )
}

if (
  is.null(cellchat@netP$centr) ||
  length(cellchat@netP$centr) == 0
) {
  stop( "CellChat netP centrality is absent. ", "Please check the Step77 final object; this plotting script ", "will not rerun centrality automatically." )
}

p <- netAnalysis_signalingRole_scatter( cellchat, signaling = NULL, slot.name = "netP", x.measure = "outdeg", y.measure = "indeg", xlabel = "Outgoing interaction strength",
  ylabel = "Incoming interaction strength", title = "Global Intercellular Signaling", dot.size = c(2, 6), label.size = 4.0, dot.alpha = 0.8, font.size = 12.5, font.size.title = 14,
  do.label = TRUE, show.legend = TRUE, show.axes = TRUE )

p <- p + ggplot2::theme( axis.title = ggplot2::element_text(size = 13), axis.text = ggplot2::element_text(size = 11.5), legend.title = ggplot2::element_text(size = 11.5),
  legend.text = ggplot2::element_text(size = 10.5), plot.margin = ggplot2::margin(12, 24, 12, 14) )

old_focal_labels <- c( "Clec7a+ lipid-handling Microglia", "Clec7a+ lipid-associated Microglia" )

new_focal_label <- "Clec7a+ lipid-handling myeloid state"

replace_label_value <- function(x) {

  if (is.factor(x)) {
    x_chr <- as.character(x)
    for (old_value in old_focal_labels) {
      x_chr <- gsub( old_value, new_focal_label, x_chr, fixed = TRUE )
    }
    return(factor(x_chr, levels = unique(x_chr)))
  }

  if (is.character(x)) {
    for (old_value in old_focal_labels) {
      x <- gsub( old_value, new_focal_label, x, fixed = TRUE )
    }
  }

  x
}

replace_labels_in_df <- function(d) {
  if (is.null(d) || !is.data.frame(d)) return(d)

  for (nm in names(d)) {
    if (is.character(d[[nm]]) || is.factor(d[[nm]])) {
      d[[nm]] <- replace_label_value(d[[nm]])
    }
  }

  d
}

replace_labels_in_list <- function(x) {
  if (is.null(x)) return(NULL)

  if (is.character(x) || is.factor(x)) {
    return(replace_label_value(x))
  }

  if (is.list(x)) {

    for (i in seq_along(x)) {
      value_i <- replace_labels_in_list(x[[i]])
      x[i] <- list(value_i)
    }
  }

  x
}

p$data <- replace_labels_in_df(p$data)

for (i in seq_along(p$layers)) {

  if (is.data.frame(p$layers[[i]]$data)) {
    p$layers[[i]]$data <- replace_labels_in_df( p$layers[[i]]$data )
  }

  if (!is.null(p$layers[[i]]$aes_params)) {
    p$layers[[i]]$aes_params <- replace_labels_in_list( p$layers[[i]]$aes_params )
  }

  if (!is.null(p$layers[[i]]$geom_params)) {
    p$layers[[i]]$geom_params <- replace_labels_in_list( p$layers[[i]]$geom_params )
  }

  if (!is.null(p$layers[[i]]$stat_params)) {
    p$layers[[i]]$stat_params <- replace_labels_in_list( p$layers[[i]]$stat_params )
  }
}

p$labels <- replace_labels_in_list(p$labels)

collect_strings <- function(x) {
  out <- character(0)

  if (is.factor(x)) {
    return(as.character(x))
  }

  if (is.character(x)) {
    return(x)
  }

  if (is.data.frame(x)) {
    for (nm in names(x)) {
      out <- c(out, collect_strings(x[[nm]]))
    }
    return(out)
  }

  if (is.list(x)) {
    for (i in seq_along(x)) {
      out <- c(out, collect_strings(x[[i]]))
    }
  }

  out
}

plot_strings <- c(
  collect_strings(p$data),
  unlist(lapply(p$layers, function(z) {
    c( collect_strings(z$data), collect_strings(z$aes_params), collect_strings(z$geom_params), collect_strings(z$stat_params) )
  }), use.names = FALSE),
  collect_strings(p$labels)
)

old_focal_label_remaining <- any( vapply( old_focal_labels, function(z) any(grepl(z, plot_strings, fixed = TRUE)), logical(1) ) )

new_focal_label_present <- any( grepl( new_focal_label, plot_strings, fixed = TRUE ) )

if (old_focal_label_remaining) {
  stop( "Reviewer 2 label replacement failed: an old Clec7a+ Microglia label ", "is still present in the final ggplot object." )
}

if (!new_focal_label_present) {
  warning( "The new Clec7a+ myeloid-state label was not detected in the final ", "ggplot object. Please visually inspect the exported figure." )
}

ggsave( filename = OUT_PDF, plot = p, width = 6.2, height = 5.3, units = "in", device = cairo_pdf, bg = "white" )

ggsave( filename = OUT_PNG, plot = p, width = 6.2, height = 5.3, units = "in", dpi = 600, bg = "white" )

ggsave( filename = OUT_SVG, plot = p, width = 6.2, height = 5.3, units = "in", device = grDevices::svg, bg = "white" )

qc <- data.frame( item = c( "Step77 CellChat object exists", "Input is CellChat object", "netP centrality available", "Reviewer2 focal-state display label present",
  "Old focal-state Microglia display label absent", "Figure5A PDF generated", "Figure5A PNG generated", "Figure5A SVG generated" ), pass = c( file.exists(IN_FILE),
  inherits(cellchat, "CellChat"), !is.null(cellchat@netP$centr) && length(cellchat@netP$centr) > 0, new_focal_label_present, !old_focal_label_remaining, file.exists(OUT_PDF),
  file.exists(OUT_PNG), file.exists(OUT_SVG) ), stringsAsFactors = FALSE )

write.csv( qc, OUT_QC, row.names = FALSE, fileEncoding = "UTF-8" )

print(qc)

cat( "\n============================================================\n", "Figure 5A v05 Reviewer 2 label-fix completed.\n", "PDF: ", OUT_PDF, "\n", "PNG: ", OUT_PNG, "\n",
  "SVG: ", OUT_SVG, "\n", "============================================================\n", sep = "" )

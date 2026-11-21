source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE)
options(width = 180)

topN <- 10L

expected_regulons <- 188L

expected_top10 <- c( "Mef2a", "Nr3c1", "Pura", "Ikzf1", "Bach1", "Fli1", "Zfp148", "Ep300", "Bclaf1", "Chd1" )

target_state <- "Clec7a+ lipid-handling myeloid state"

point_size <- 2.4
label_size <- 4.25

fig_width_in <- 5.0
fig_height_in <- 4.0
png_dpi <- 600

rss_file <- file.path(PROJECT_ROOT, "results", "89_pySCENIC_RSS_target_regulon_ranking", "05_89_Clec7a_lipid_handling_MG_regulon_RSS_ranking_v01.csv")

outdir <- file.path(PROJECT_ROOT, "results", "95_Figure6A_RSS_specificity_ranking")

if (!dir.exists(outdir)) {
  dir.create( outdir, recursive = TRUE )
}

plot_data_file <- file.path( outdir, "01_95_Figure6A_RSS_plot_data_v02.csv" )

top10_file <- file.path( outdir, "02_95_Figure6A_top10_RSS_regulons_v02.csv" )

figure_pdf <- file.path( outdir, "03_95_Figure6A_RSS_specificity_ranking_v02.pdf" )

figure_svg <- file.path( outdir, "03_95_Figure6A_RSS_specificity_ranking_v02.svg" )

figure_png <- file.path( outdir, "03_95_Figure6A_RSS_specificity_ranking_v02.png" )

qc_file <- file.path( outdir, "04_95_FINAL_QC_v02.csv" )

report_file <- file.path( outdir, "05_95_REPORT_v02.txt" )

required_packages <- c( "ggplot2", "ggrepel" )

missing_packages <- required_packages[ !vapply( required_packages, requireNamespace, logical(1), quietly = TRUE ) ]

if (length(missing_packages) > 0) {
  stop( "Missing package(s): ", paste( missing_packages, collapse = ", " ) )
}

if (!file.exists(rss_file)) {
  stop( "Frozen Step89 RSS file not found:\n", rss_file )
}

rss <- read.csv( rss_file, stringsAsFactors = FALSE, check.names = FALSE )

required_cols <- c( "regulon", "TF", "target_RSS", "RSS_rank" )

if (!all(required_cols %in% colnames(rss))) {
  stop( "Step89 RSS table lacks required columns: ", paste( setdiff( required_cols, colnames(rss) ), collapse = ", " ) )
}

rss <- rss[ is.finite(rss$target_RSS) & is.finite(rss$RSS_rank), , drop = FALSE ]

rss <- rss[ order( rss$RSS_rank, -rss$target_RSS, rss$TF ), , drop = FALSE ]

if (nrow(rss) != expected_regulons) {
  stop( "Unexpected regulon count: ", nrow(rss), " / expected ", expected_regulons )
}

if (anyDuplicated(rss$TF)) {
  stop( "Duplicated TF names detected in frozen RSS table." )
}

if (
  !identical( rss$RSS_rank, seq_len( nrow(rss) ) )
) {
  stop( "RSS_rank is not a complete 1..N sequence." )
}

observed_top10 <- rss$TF[ seq_len( topN ) ]

if (!identical(observed_top10, expected_top10)) {
  stop( "Frozen RSS top10 mismatch.\nObserved: ", paste( observed_top10, collapse = ", " ), "\nExpected: ", paste( expected_top10, collapse = ", " ) )
}

df <- data.frame( Regulons = rss$RSS_rank, SpecificityScore = rss$target_RSS, TF = rss$TF, regulon = rss$regulon, stringsAsFactors = FALSE )

df$point_class <- ifelse( df$Regulons <= topN, "Top10", "Other" )

df$label <- df$TF

df$label[ df$Regulons > topN ] <- NA_character_

write.csv( df, plot_data_file, row.names = FALSE )

write.csv( rss[ seq_len(topN), , drop = FALSE ], top10_file, row.names = FALSE )

library(ggplot2)
library(ggrepel)

x_max <- max( df$Regulons )

p <- ggplot( df, aes( x = Regulons, y = SpecificityScore ) ) +
  geom_point( aes( color = point_class ), size = point_size ) +
  ggrepel::geom_text_repel( aes( label = label ), size = label_size, fontface = "italic", box.padding = 0.55, point.padding = 0.30, min.segment.length = 0, max.overlaps = Inf, seed = 9501,
    segment.color = "grey55", segment.linewidth = 0.35 ) +
  scale_color_manual( values = c( "Top10" = "red", "Other" = "gray" ), guide = "none" ) +
  scale_x_continuous( limits = c( 0, x_max + 5 ), breaks = c( 1, 50, 100, 150, x_max ), expand = expansion( mult = c( 0, 0 ) ) ) +
  labs( title = target_state, x = "Regulons Rank", y = "Specificity Score" ) +
  theme_bw( base_size = 13 ) +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text( hjust = 0.5, face = "bold", size = 13 ), axis.title = element_text( face = "bold"
    ), axis.text = element_text( color = "black" ), panel.border = element_rect( color = "black", linewidth = 0.6 ) )

ggsave( figure_pdf, p, width = fig_width_in, height = fig_height_in, units = "in", device = cairo_pdf )

if (requireNamespace("svglite", quietly = TRUE)) {

  ggsave( figure_svg, p, width = fig_width_in, height = fig_height_in, units = "in", device = svglite::svglite )

} else {

  svg( figure_svg, width = fig_width_in, height = fig_height_in )

  print(p)

  dev.off()
}

ggsave( figure_png, p, width = fig_width_in, height = fig_height_in, units = "in", dpi = png_dpi )

qc <- data.frame( check = c( "Frozen Step89 file exists", "Regulon count", "RSS rank sequence", "Frozen top10 identity", "Top10 labels", "Reviewer2 focal-state display terminology",
  "No RSS recalculation", "Publication outputs" ), status = c( "PASS", if ( nrow(rss) == expected_regulons ) "PASS" else "FAIL", if ( identical( rss$RSS_rank, seq_len(nrow(rss)) )
  ) "PASS" else "FAIL", if ( identical( observed_top10, expected_top10 ) ) "PASS" else "FAIL", if ( sum(!is.na(df$label)) == topN ) "PASS" else "FAIL", if ( identical( target_state,
  "Clec7a+ lipid-handling myeloid state" ) ) "PASS" else "FAIL", "PASS", if ( all( file.exists( c( figure_pdf, figure_svg, figure_png ) ) ) ) "PASS" else "FAIL" ), detail = c( rss_file,
  paste0( nrow(rss), " / ", expected_regulons ), "RSS_rank = 1..188", paste( observed_top10, collapse = ", " ), paste0( sum(!is.na(df$label)), " labels" ), target_state,
  "Figure95 reads frozen target_RSS/RSS_rank only.", "PDF + SVG + 600 dpi PNG" ), stringsAsFactors = FALSE )

write.csv( qc, qc_file, row.names = FALSE )

if (any(qc$status == "FAIL")) {
  stop( "Step95 contains blocking QC FAIL item(s)." )
}

report_lines <- c( "============================================================", "Step 95: Figure 6A RSS specificity ranking",
  "============================================================", "", "STATUS: PASS", "", "REFERENCE CODE LOGIC", "All regulons are ranked by specificity score.",
  "Top regulons are colored red and labeled; remaining regulons are gray.", "X-axis = Regulons Rank; Y-axis = Specificity Score.", "", "PROJECT ADAPTATION", paste0( "Target state: ",
  target_state ), paste0( "Frozen regulons: ", nrow(rss) ), paste0( "Highlighted top N: ", topN ), paste0( "Top10: ", paste( observed_top10, collapse = ", " ) ), "", "IMPORTANT",
  "RSS was NOT recalculated.", "The figure reads the frozen Step89 RSS ranking directly.", "", "============================================================" )

writeLines( report_lines, report_file )

cat("\n")
cat("============================================================\n")
cat("STEP 95 COMPLETE\n")
cat("============================================================\n")
cat("Top10 RSS TFs:\n")
print( rss[ seq_len(topN), c( "RSS_rank", "TF", "regulon", "target_RSS" ), drop = FALSE ], row.names = FALSE )
cat("\nOutputs:\n")
cat(figure_png, "\n")
cat(figure_pdf, "\n")
cat(figure_svg, "\n")
cat(qc_file, "\n")
cat(report_file, "\n")
cat("============================================================\n")

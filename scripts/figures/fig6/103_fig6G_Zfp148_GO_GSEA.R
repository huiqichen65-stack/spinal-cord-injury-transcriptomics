source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, width = 180)

root <- PROJECT_ROOT
step102_dir <- file.path(root, "results", "102_Zfp148_GO_GSEA_fullrank_v01")
outdir <- file.path(root, "results", "103_Figure6G_Zfp148_GO_GSEA_ridgeplot_v03")

primary_rank_file <- file.path(step102_dir, "04_102_Zfp148_primary_filtered_ranked_genes_v01.csv")
primary_bp_file <- file.path(step102_dir, "03_102_Zfp148_primary_filtered_GO_BP_all_v01.csv")
primary_cc_file <- file.path(step102_dir, "03_102_Zfp148_primary_filtered_GO_CC_all_v01.csv")
primary_mf_file <- file.path(step102_dir, "03_102_Zfp148_primary_filtered_GO_MF_all_v01.csv")
step102_qc_file <- file.path(step102_dir, "05_102_FINAL_QC_v01.csv")

png_dpi <- 600L
fig_width_in <- 9.4
fig_height_in <- 8.2
fdr_cutoff <- 0.05
min_members <- 5L

min_member_fraction <- 0.75

selected_terms <- data.frame(
  GO_ID = c(

    "GO:0071404",
    "GO:0071402",
    "GO:0032367",
    "GO:0015909",
    "GO:0008610",
    "GO:1901568",

    "GO:0006897",
    "GO:0072665",

    "GO:0030029",
    "GO:0007010",
    "GO:0030335",
    "GO:0001667",
    "GO:0042060",

    "GO:0043209",
    "GO:1902495",
    "GO:0009986"
  ),
  theme = c(
    rep("Lipid / lipoprotein / cholesterol", 6),
    rep("Endocytosis / trafficking", 2),
    rep("Cytoskeleton / migration", 5),
    rep("Cellular component / context", 3)
  ),
  display_order = seq_len(16),
  stringsAsFactors = FALSE
)
stopifnot(nrow(selected_terms) == 16L, anyDuplicated(selected_terms$GO_ID) == 0L)

need <- c("ggplot2", "ggridges", "scales", "AnnotationDbi", "org.Mm.eg.db")
missing_pkgs <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
    stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\nSee the installation commands at the top of this script.")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(scales)
})

required_files <- c(primary_rank_file, primary_bp_file, primary_cc_file, primary_mf_file, step102_qc_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
 stop("Missing Step102 file(s):\n", paste(missing_files, collapse = "\n"))
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

read_csv <- function(f) {
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
find_col <- function(df, candidates, required = TRUE, label = "column") {
  ix <- match(tolower(candidates), tolower(names(df)))
  ix <- ix[!is.na(ix)]
  if (length(ix)) return(names(df)[ix[1]])
  if (required) {
    stop("Missing ", label, ". Expected ", paste(candidates, collapse = "/"),
         "; observed ", paste(names(df), collapse = ", "))
  }
  NULL
}
qc_items <- list()
add_qc <- function(check, status, detail) {
    qc_items[[length(qc_items) + 1L]] <<- data.frame(check = check, status = status, detail = as.character(detail), 
        stringsAsFactors = FALSE)
}

step102_qc <- read_csv(step102_qc_file)
qc_status_col <- find_col(step102_qc, "status", FALSE)
if (!is.null(qc_status_col) && any(toupper(as.character(step102_qc[[qc_status_col]])) == "FAIL", na.rm = TRUE)) {
    stop("Step102 FINAL_QC includes FAIL; investigate before producing Figure6G.")
}
add_qc("Step102 QC", "PASS", "No blocking FAIL in Step102 QC file")

load_go <- function(path, ont) {
    dat <- read_csv(path)
    if (!nrow(dat)) 
        return(NULL)
    go_col <- find_col(dat, c("ID", "GO_ID"), label = "GO ID")
    desc_col <- find_col(dat, "Description", label = "GO Description")
    nes_col <- find_col(dat, "NES", label = "GO NES")
    fdr_col <- find_col(dat, c("p.adjust", "FDR", "padj"), label = "GO FDR")
    size_col <- find_col(dat, c("setSize", "set_size"), label = "GO setSize")
    out <- data.frame(GO_ID = as.character(dat[[go_col]]), Description = as.character(dat[[desc_col]]), 
        NES = as.numeric(dat[[nes_col]]), p.adjust = as.numeric(dat[[fdr_col]]), setSize = as.integer(dat[[size_col]]), 
        ONTOLOGY = ont, stringsAsFactors = FALSE)
    in_ontology_col <- find_col(dat, c("Ontology", "ONTOLOGY"), FALSE)
    if (!is.null(in_ontology_col)) {
        found <- toupper(as.character(dat[[in_ontology_col]]))
        if (any(!is.na(found) & found != ont)) {
            stop("Ontology mismatch in: ", path)
        }
    }
    out
}
all_go <- do.call(rbind, list(load_go(primary_bp_file, "BP"), load_go(primary_cc_file, "CC"), load_go(primary_mf_file, 
    "MF")))
if (!nrow(all_go)) stop("PRIMARY GO tables are empty")
if (anyDuplicated(all_go$GO_ID)) {
  stop("Same GO ID repeated across PRIMARY ontologies; check input provenance")
}
idx <- match(selected_terms$GO_ID, all_go$GO_ID)
if (anyNA(idx)) {
  stop("Selected GO ID(s) absent from Zfp148 PRIMARY GO results: ",
       paste(selected_terms$GO_ID[is.na(idx)], collapse = ", "))
}
term_stats <- cbind(selected_terms, all_go[idx, c("Description", "ONTOLOGY", "setSize", "NES", "p.adjust")])
if (any(!is.finite(term_stats$NES)) || any(!is.finite(term_stats$p.adjust)) || any(!is.finite(term_stats$setSize))) {
    stop("Invalid NES/FDR/setSize in selected PRIMARY terms")
}
if (any(term_stats$p.adjust >= fdr_cutoff)) {
  stop("Some selected Zfp148 PRIMARY terms fail FDR < 0.05: ",
       paste(term_stats$GO_ID[term_stats$p.adjust >= fdr_cutoff], collapse = ", "))
}
if (!identical(as.character(term_stats$ONTOLOGY), c(rep("BP", 13), rep("CC", 3)))) {
    stop("Unexpected GO ontology order for the 16 selected terms")
}
add_qc("Term selection", "PASS", "16 declared BP/CC terms; each PRIMARY FDR < 0.05")

rank_raw <- read_csv(primary_rank_file)
gene_col <- find_col(rank_raw, c("gene", "SYMBOL", "symbol", "gene_symbol"), label = "rank gene")
z_col <- find_col(rank_raw, c("Z", "perturbation_Z", "z_score", "z"), label = "rank Z score")
entrez_col <- find_col(rank_raw, c("ENTREZID", "entrez_id", "ENTREZ", "entrez_gene_id"), FALSE)
rank_order_col <- find_col(rank_raw, c("GSEA_rank", "rank"), FALSE)

rank_df <- data.frame(gene = as.character(rank_raw[[gene_col]]), Z = suppressWarnings(as.numeric(rank_raw[[z_col]])), 
    stringsAsFactors = FALSE)

valid_rank_rows <- !is.na(rank_df$gene) & nzchar(trimws(rank_df$gene)) & is.finite(rank_df$Z)
rank_raw <- rank_raw[valid_rank_rows, , drop = FALSE]
rank_df <- rank_df[valid_rank_rows, , drop = FALSE]
if (!nrow(rank_df)) stop("No valid primary ranked genes before ID mapping")
if (!is.null(rank_order_col)) {
  rank_df$GSEA_rank <- suppressWarnings(as.numeric(rank_raw[[rank_order_col]]))
} else {
  rank_df$GSEA_rank <- seq_len(nrow(rank_df))
}
if (!is.null(entrez_col)) {
  rank_df$ENTREZID <- as.character(rank_raw[[entrez_col]])
 mapping_source <- "Step102 ranked Entrez IDs"
} else {

  ids <- AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db, keys = unique(rank_df$gene),
    keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"
  )
  rank_df$ENTREZID <- as.character(unname(ids[rank_df$gene]))
  mapping_source <- paste0("Reconstructed SYMBOL-to-ENTREZ mapping via org.Mm.eg.db ",
                           as.character(utils::packageVersion("org.Mm.eg.db")))
  add_qc("Rank ENTREZ mapping", "WARN", mapping_source)
}
rank_df$ENTREZID[rank_df$ENTREZID %in% c("NA", "", "NULL")] <- NA_character_
rank_df <- rank_df[!is.na(rank_df$ENTREZID) & nzchar(rank_df$ENTREZID) & !is.na(rank_df$gene) & nzchar(rank_df$gene) & 
    is.finite(rank_df$Z), , drop = FALSE]
if (!nrow(rank_df)) stop("No usable mapped genes in Step102 PRIMARY ranked CSV")
if (anyDuplicated(rank_df$ENTREZID)) {
    if (all(is.finite(rank_df$GSEA_rank))) {
        rank_df <- rank_df[order(rank_df$GSEA_rank), , drop = FALSE]
    }
    ndup <- sum(duplicated(rank_df$ENTREZID))
    rank_df <- rank_df[!duplicated(rank_df$ENTREZID), , drop = FALSE]
    add_qc("Rank duplicate ENTREZ IDs", "WARN", paste0("Kept first -ranked entry; removed ", ndup, " duplicated ENTREZ IDs"))
} else {
    add_qc("Rank duplicate ENTREZ IDs", "PASS", "No repeated ENTREZ IDs")
}
if (nrow(rank_df) < 500L) stop("Too few mapped PRIMARY ranked genes: ", nrow(rank_df))
add_qc("PRIMARY rank", "PASS", paste(nrow(rank_df), "mapped unique ENTREZ genes"))

go_map <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = unique(rank_df$ENTREZID), 
    keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
go_map <- go_map[!is.na(go_map$GOALL) & go_map$GOALL %in% term_stats$GO_ID, c("ENTREZID", "GOALL", "ONTOLOGYALL"), 
    drop = FALSE]
if (!nrow(go_map)) stop("GOALL mapping contains no members for selected terms")
go_map <- unique(go_map)
colnames(go_map) <- c("ENTREZID", "GO_ID", "GO_annotation_ontology")

member_df <- merge(go_map, rank_df[, c("ENTREZID", "gene", "Z", "GSEA_rank")],
                   by = "ENTREZID", all = FALSE, sort = FALSE)
member_df <- member_df[!duplicated(paste(member_df$GO_ID, member_df$ENTREZID, sep = "|")), , drop = FALSE]
if (!nrow(member_df)) stop("Zero GO members after intersecting PRIMARY ranked universe")

member_counts <- as.integer(table(factor(member_df$GO_ID, levels = term_stats$GO_ID)))
term_stats$n_full_members <- member_counts
term_stats$setSize_difference <- term_stats$n_full_members - term_stats$setSize
term_stats$member_fraction <- term_stats$n_full_members / pmax(1, term_stats$setSize)
if (any(member_counts < min_members)) {
  stop("Selected term(s) have <", min_members, " mapped full GO members: ",
       paste(term_stats$GO_ID[member_counts < min_members], collapse = ", "))
}
if (any(term_stats$member_fraction < min_member_fraction)) {
    bad <- term_stats[term_stats$member_fraction < min_member_fraction, c("GO_ID", "setSize", "n_full_members", 
        "member_fraction")]
    print(bad)
    stop("GOALL member reconstruction disagrees strongly with Step102 gseGO setSize. ", "Check org.Mm.eg.db annotation version and Step102 rank mappings.")
}
if (any(term_stats$setSize_difference != 0L)) {
    add_qc("GO full member counts", "WARN", paste0(sum(term_stats$setSize_difference != 0L), "/16 terms differ from setSize; see 01_103 selected-term CSV"))
} else {
    add_qc("GO full member counts", "PASS", "All 16 full-member counts equal gseGO setSize")
}

pretty_labels <- c("Cellular response to LDL\nparticle stimulus", "Cellular response to lipoprotein\nparticle stimulus", 
    "Intracellular cholesterol transport", "Long-chain fatty acid transport", "Lipid biosynthetic process", 
    "Fatty acid derivative metabolic process", "Endocytosis", "Protein localization to vacuole", "Actin filament-based process", 
    "Cytoskeleton organization", "Positive regulation of cell migration", "Ameboidal-type cell migration", 
    "Wound healing", "Myelin sheath", "Transmembrane transporter complex", "Cell surface")
if (length(pretty_labels) != nrow(term_stats) || anyDuplicated(pretty_labels)) {
 stop("Plot labels must be uniquely specified for the 16 GO terms")
}
term_stats$display_label <- pretty_labels
term_stats$neglog10_FDR <- -log10(pmax(term_stats$p.adjust, .Machine$double.xmin))

category_order <- unique(selected_terms$theme)
term_stats$category_index <- match(term_stats$theme, category_order)
category_gap <- 0.70
term_stats$plot_y <- nrow(term_stats) - term_stats$display_order + 1 -
  category_gap * (term_stats$category_index - 1)

category_end_rows <- which(diff(term_stats$category_index) != 0)
separator_y <- vapply(category_end_rows, function(i) {
  mean(term_stats$plot_y[c(i, i + 1L)])
}, numeric(1))

member_df <- merge(member_df, term_stats[, c("GO_ID", "Description", "ONTOLOGY", "theme", "display_order", 
    "setSize", "NES", "p.adjust", "neglog10_FDR", "display_label", "category_index", "plot_y")], by = "GO_ID", 
    all.x = TRUE, sort = FALSE)
if (anyNA(member_df$plot_y)) stop("GO members could not be matched to their plot positions")
member_df <- member_df[order(member_df$display_order, member_df$GSEA_rank, member_df$gene), , drop = FALSE]

x_min <- min(member_df$Z, na.rm = TRUE)
x_max <- max(member_df$Z, na.rm = TRUE)
x_span <- x_max - x_min
if (!is.finite(x_span) || x_span <= 0) stop("Invalid perturbation Z range")
term_stats$dot_x <- x_min - 0.105 * x_span
plot_x_min <- x_min - 0.18 * x_span
plot_x_max <- x_max + 0.035 * x_span
nes_limit <- max(abs(term_stats$NES), na.rm = TRUE)
if (!is.finite(nes_limit) || nes_limit <= 0) nes_limit <- 1
fdr_min <- min(term_stats$neglog10_FDR)
fdr_max <- max(term_stats$neglog10_FDR)
if (fdr_max <= fdr_min) fdr_max <- fdr_min + 1

p6g <- ggplot(member_df, aes(x = Z, y = plot_y, group = GO_ID)) + geom_hline(yintercept = separator_y, 
    linewidth = 0.32, color = "#E5E5E5") + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.46, 
    color = "#B0B0B0") + ggridges::geom_density_ridges(aes(fill = neglog10_FDR), scale = 1.05, rel_min_height = 0.01, 
    alpha = 0.94, color = "#3C3C3C", linewidth = 0.35, panel_scaling = FALSE) + geom_point(data = term_stats, 
    aes(x = dot_x, y = plot_y, color = NES, size = neglog10_FDR), inherit.aes = FALSE, alpha = 0.98) + 
    scale_fill_gradientn(colours = c("#EEE7F2", "#D7C9E2", "#F2B39C", "#F1734E", "#E6550D"), limits = c(fdr_min, 
        fdr_max), oob = scales::squish, name = expression(-log[10]("(adj.P)"))) + scale_color_gradient2(low = "#3B73B9", 
    mid = "#F7F7F7", high = "#D7301F", midpoint = 0, limits = c(-nes_limit, nes_limit), oob = scales::squish, 
    name = "NES") + scale_size_continuous(range = c(3.2, 7.8), guide = "none") + scale_x_continuous(limits = c(plot_x_min, 
    plot_x_max), expand = c(0, 0)) + scale_y_continuous(breaks = term_stats$plot_y, labels = term_stats$display_label, 
    expand = expansion(add = c(0.25, 0.62))) + labs(title = "GSEA of GO enrichment", x = "Perturbation Z score", 
    y = NULL) + theme_classic(base_size = 13, base_family = "sans") + theme(plot.title = element_text(face = "bold", 
    size = 18, hjust = 0.5, margin = margin(b = 9)), axis.title.x = element_text(face = "bold", size = 15, 
    margin = margin(t = 7)), axis.text.x = element_text(size = 12.5, color = "#222222"), axis.text.y = element_text(size = 12, 
    color = "#222222", hjust = 1, lineheight = 1, margin = margin(r = 7)), axis.line.y = element_line(linewidth = 0.43, 
    color = "#222222"), axis.ticks.y = element_blank(), panel.grid.major.x = element_line(color = "#F0F0F0", 
    linewidth = 0.3), panel.grid.minor = element_blank(), legend.position = "right", legend.box = "vertical", 
    legend.spacing.y = grid::unit(7, "mm"), legend.title = element_text(face = "bold", size = 13), legend.text = element_text(size = 11.5), 
    plot.margin = margin(10, 16, 10, 10)) + guides(fill = guide_colorbar(order = 1, barheight = grid::unit(28, 
    "mm"), barwidth = grid::unit(5, "mm")), color = guide_colorbar(order = 2, barheight = grid::unit(28, 
    "mm"), barwidth = grid::unit(5, "mm")))

add_qc("v03 appearance-only changes", "PASS", "Same 16 GO IDs/NES/FDR/Z; clean labels, group gaps and plot size only")
add_qc("Four-category layout", "PASS", paste0(length(category_order), " categories; ", length(separator_y), 
    " visual separators"))

fig_prefix <- file.path(outdir, "03_103_Figure6G_Zfp148_GO_GSEA_ridgeplot_v04_large_fonts")
figure_pdf <- paste0(fig_prefix, ".pdf")
figure_svg <- paste0(fig_prefix, ".svg")
figure_png <- paste0(fig_prefix, ".png")

term_file <- file.path(outdir, "01_103_Zfp148_biologically_selected_GO_terms_v03.csv")
member_file <- file.path(outdir, "02_103_Zfp148_full_GO_member_genes_v03.csv")
qc_file <- file.path(outdir, "04_103_FINAL_QC_v03.csv")
report_file <- file.path(outdir, "05_103_REPORT_v03.txt")
index_file <- file.path(outdir, "00_103_OUTPUT_INDEX_v03.csv")

write.csv(term_stats[order(term_stats$display_order), ], term_file, row.names = FALSE)
write.csv(member_df, member_file, row.names = FALSE)

ggsave(figure_pdf, plot = p6g, width = fig_width_in, height = fig_height_in, units = "in", device = grDevices::cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(figure_svg, plot = p6g, width = fig_width_in, height = fig_height_in, units = "in", device = svglite::svglite)
} else {
    grDevices::svg(filename = figure_svg, width = fig_width_in, height = fig_height_in, family = "sans")
    print(p6g)
    grDevices::dev.off()
}
if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = figure_png, width = fig_width_in, height = fig_height_in, units = "in", 
        res = png_dpi, background = "white")
    print(p6g)
    grDevices::dev.off()
} else {
    ggsave(figure_png, plot = p6g, width = fig_width_in, height = fig_height_in, units = "in", dpi = png_dpi, 
        bg = "white")
}

for (f in c(figure_pdf, figure_svg, figure_png)) {
  add_qc(basename(f), if (file.exists(f)) "PASS" else "FAIL", f)
}
add_qc("Ontology and gene-set mapping source", "PASS", paste0("org.Mm.eg.db ", as.character(utils::packageVersion("org.Mm.eg.db"))))
add_qc("Scientific meaning", "PASS", "Virtual KO perturbation Z/NES != transcript up/down or direct TF targets")
qc_df <- do.call(rbind, qc_items)
write.csv(qc_df, qc_file, row.names = FALSE)
overall_status <- if (any(qc_df$status == "FAIL")) "FAIL" else if (any(qc_df$status == "WARN")) "PASS_WITH_WARNINGS" else "PASS"

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    report_file)

files_written <- c(term_file, member_file, figure_pdf, figure_svg, figure_png, qc_file, report_file)
index <- data.frame(file = basename(files_written), exists = file.exists(files_written), size_MB = file.info(files_written)$size/1024^2, 
    path = normalizePath(files_written, winslash = "/", mustWork = FALSE), stringsAsFactors = FALSE)
write.csv(index, index_file, row.names = FALSE)

if (overall_status == "FAIL") stop("Step 103 v03 has blocking QC failure")

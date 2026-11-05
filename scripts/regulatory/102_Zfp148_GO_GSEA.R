source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, width = 160)

project_root <- PROJECT_ROOT
primary_tf <- "Zfp148"
compare_with_mef2a <- TRUE
fdr_cutoff <- 0.05
min_gs_size <- 10L
max_gs_size <- 500L
ridge_top_n <- 16L
png_dpi <- 450L
seed <- 10201L

technical_patterns <- c("^mt-", "^Rpl", "^Rps", "^Gm[0-9]+", "Rik$", "-ps[0-9]*$", "^Rn45s$", "^Mir[0-9]+",
  "^Hist", "^Xist$", "^AY[0-9]+", "^AI[0-9]+", "^AU[0-9]+")

input_file <- file.path(project_root, "results", "91_scTenifoldKnk_top5_virtual_KO", "05_91_top5_KO_diffRegulation_all_v04.csv.gz")
outdir <- file.path(project_root, "results", "102_Zfp148_GO_GSEA_fullrank_v01")
if (!file.exists(input_file)) stop("Cannot find Step 91 input: ", input_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_packages <- c("clusterProfiler", "org.Mm.eg.db", "AnnotationDbi")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
needed <- c(bioc_packages, "ggplot2")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages after attempted install: ", paste(missing, collapse = ", "))

db <- org.Mm.eg.db::org.Mm.eg.db

write_table <- function(x, basename) {
  f <- file.path(outdir, basename)
  utils::write.csv(x, f, row.names = FALSE, na = "")
  message("Saved: ", f)
  invisible(f)
}
write_lines <- function(x, basename) {
  f <- file.path(outdir, basename)
  writeLines(x, f, useBytes = TRUE)
  message("Saved: ", f)
  invisible(f)
}
is_technical <- function(g) {
  g[is.na(g)] <- ""
  Reduce(`|`, lapply(technical_patterns, grepl, x = g, ignore.case = TRUE, perl = TRUE))
}

raw <- utils::read.csv(gzfile(input_file), check.names = FALSE)
required <- c("gene", "KO_TF", "Z", "distance", "p.adj", "is_KO_gene")
if (!all(required %in% names(raw))) {
  stop("Input lacks required columns: ", paste(setdiff(required, names(raw)), collapse = ", "))
}
tfs <- if (compare_with_mef2a) c("Zfp148", "Mef2a") else "Zfp148"
if (!all(tfs %in% unique(raw$KO_TF))) stop("Step 91 input missing a requested KO TF")
if (anyDuplicated(raw[, c("KO_TF", "gene")])) {
  stop("Unexpected duplicated (KO_TF, gene) in Step 91; audit upstream data")
}

raw <- raw[raw$KO_TF %in% tfs, , drop = FALSE]
raw$gene <- trimws(as.character(raw$gene))
raw <- raw[!is.na(raw$gene) & nzchar(raw$gene) & is.finite(raw$Z), , drop = FALSE]
raw <- raw[!(raw$is_KO_gene %in% TRUE) & raw$gene != raw$KO_TF, , drop = FALSE]
if (any(table(raw$KO_TF) < 1000)) stop("Too few genes for at least one selected TF")
raw$technical_gene <- is_technical(raw$gene)

all_syms <- sort(unique(raw$gene))
entrez_map <- suppressMessages(AnnotationDbi::mapIds(db, keys = all_syms, keytype = "SYMBOL", column = "ENTREZID",
  multiVals = "first"))
map_df <- data.frame(gene = names(entrez_map), ENTREZID = as.character(unname(entrez_map)), stringsAsFactors = FALSE)
map_df$mapped <- !is.na(map_df$ENTREZID) & nzchar(map_df$ENTREZID)
map_df$technical_gene <- is_technical(map_df$gene)
write_table(map_df, "01_102_gene_ID_mapping_audit_v01.csv")

raw$ENTREZID <- map_df$ENTREZID[match(raw$gene, map_df$gene)]
write_table(aggregate(cbind(total = rep(1L, nrow(raw)), mapped = as.integer(!is.na(raw$ENTREZID)), technical = as.integer(raw$technical_gene)) ~
  KO_TF, data = raw, FUN = sum), "02_102_input_gene_counts_v01.csv")

empty_go <- function() {
  data.frame(ID = character(), Description = character(), setSize = integer(), enrichmentScore = numeric(),
    NES = numeric(), pvalue = numeric(), p.adjust = numeric(), qvalues = numeric(), rank = integer(),
    leading_edge = character(), core_enrichment = character(), Ontology = character(), KO_TF = character(),
    Analysis = character(), stringsAsFactors = FALSE)
}

make_rank <- function(df, analysis) {
  x <- df[!is.na(df$ENTREZID) & nzchar(df$ENTREZID), , drop = FALSE]
  if (identical(analysis, "primary_filtered")) {
    x <- x[!x$technical_gene, , drop = FALSE]
  }
  x <- x[order(-abs(x$Z), x$gene), , drop = FALSE]
  n_collisions <- sum(duplicated(x$ENTREZID))
  x <- x[!duplicated(x$ENTREZID), , drop = FALSE]
  x <- x[order(-x$Z, x$ENTREZID), , drop = FALSE]
  x$rank <- seq_len(nrow(x))
  if (nrow(x) < 1000L)
    stop("Too few Entrez-mapped genes: ", analysis)
  scores <- setNames(as.numeric(x$Z), x$ENTREZID)
  if (anyNA(scores) || anyDuplicated(names(scores)))
    stop("Invalid GO-GSEA ranking")
  attr(x, "n_collisions") <- n_collisions
  list(df = x, ranks = scores)
}

run_go <- function(ranks, tf, analysis) {
  all_out <- list()
  for (ont in c("BP", "CC", "MF")) {
    message("Running gseGO: TF=", tf, " analysis=", analysis, " ontology=", ont)
    set.seed(seed)
    gs <- tryCatch(clusterProfiler::gseGO(geneList = ranks, OrgDb = db, keyType = "ENTREZID", ont = ont,
      minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = 1, pAdjustMethod = "BH",
      by = "fgsea", eps = 1e-10, seed = TRUE, verbose = FALSE), error = function(e) stop("gseGO failed (",
      tf, "/", analysis, "/", ont, "): ", conditionMessage(e)))
    tab <- if (is.null(gs))
      empty_go() else as.data.frame(gs)
    if (nrow(tab) == 0L)
      tab <- empty_go()
    if (nrow(tab) > 0L) {
      tab$Ontology <- ont
      tab$KO_TF <- tf
      tab$Analysis <- analysis
      tab <- tab[order(tab$p.adjust, tab$pvalue), , drop = FALSE]
      tab$BH_FDR_within_ontology <- tab$p.adjust
    } else {
      tab$BH_FDR_within_ontology <- numeric(0)
    }
    filebase <- paste0("03_102_", tf, "_", analysis, "_GO_", ont)
    write_table(tab, paste0(filebase, "_all_v01.csv"))
    sig <- tab[!is.na(tab$p.adjust) & tab$p.adjust < fdr_cutoff, , drop = FALSE]
    write_table(sig, paste0(filebase, "_FDR05_v01.csv"))
    all_out[[ont]] <- tab
  }
  do.call(rbind, all_out)
}

rank_data <- list()
go_results <- list()
qclist <- list()
for (tf in tfs) {
  d <- raw[raw$KO_TF == tf, , drop = FALSE]
  for (analysis in c("primary_filtered", "fullnetwork_sensitivity")) {
    rr <- make_rank(d, analysis)
    key <- paste(tf, analysis, sep = "__")
    rank_data[[key]] <- rr$df
    write_table(rr$df, paste0("04_102_", tf, "_", analysis, "_ranked_genes_v01.csv"))
    go_results[[key]] <- run_go(rr$ranks, tf, analysis)
    a <- go_results[[key]]
    qclist[[key]] <- data.frame(KO_TF = tf, Analysis = analysis, tested_nonKO_genes = nrow(d), mapped_unique_Entrez = nrow(rr$df),
      Entrez_duplicate_collisions = attr(rr$df, "n_collisions"), tested_BP = sum(a$Ontology ==
        "BP"), tested_CC = sum(a$Ontology == "CC"), tested_MF = sum(a$Ontology == "MF"), significant_BP = sum(a$Ontology ==
        "BP" & !is.na(a$p.adjust) & a$p.adjust < fdr_cutoff), significant_CC = sum(a$Ontology ==
        "CC" & !is.na(a$p.adjust) & a$p.adjust < fdr_cutoff), significant_MF = sum(a$Ontology ==
        "MF" & !is.na(a$p.adjust) & a$p.adjust < fdr_cutoff), stringsAsFactors = FALSE)
  }
}
qc <- do.call(rbind, qclist)
rownames(qc) <- NULL
write_table(qc, "05_102_FINAL_QC_v01.csv")

lipid_pattern <- paste(c("lipid", "lipoprotein", "cholesterol", "sterol", "lysosom", "sphingo", "fatty acid",
  "phospholipid", "triglyceride", "eicosanoid", "arachidon", "endocyt", "lipid droplet"), collapse = "|")
all_go <- do.call(rbind, go_results)
rownames(all_go) <- NULL
if (nrow(all_go) > 0L) {
  all_go$is_lipid_or_lysosome_related <- grepl(lipid_pattern, all_go$Description, ignore.case = TRUE)
  write_table(all_go[all_go$is_lipid_or_lysosome_related, , drop = FALSE][order(all_go[all_go$is_lipid_or_lysosome_related,
    "p.adjust"], na.last = TRUE), , drop = FALSE], "06_102_lipid_related_GO_all_v01.csv")
} else {
  write_table(all_go, "06_102_lipid_related_GO_all_v01.csv")
}

expand_leading_edge <- function(tab, ranked, significant_only = TRUE) {
  if (significant_only)
    tab <- tab[!is.na(tab$p.adjust) & tab$p.adjust < fdr_cutoff, , drop = FALSE]
  if (!nrow(tab))
    return(data.frame())
  z <- setNames(ranked$Z, ranked$ENTREZID)
  sym <- setNames(ranked$gene, ranked$ENTREZID)
  x <- lapply(seq_len(nrow(tab)), function(i) {
    ids <- unlist(strsplit(as.character(tab$core_enrichment[i]), "/", fixed = TRUE))
    ids <- ids[nzchar(ids) & ids %in% names(z)]
    if (!length(ids))
      return(NULL)
    data.frame(KO_TF = tab$KO_TF[i], Analysis = tab$Analysis[i], Ontology = tab$Ontology[i], GO_ID = tab$ID[i],
      GO_term = tab$Description[i], NES = tab$NES[i], GO_FDR = tab$p.adjust[i], gene = unname(sym[ids]),
      ENTREZID = ids, perturbation_Z = unname(z[ids]), stringsAsFactors = FALSE)
  })
  x <- Filter(Negate(is.null), x)
  if (length(x))
    do.call(rbind, x) else data.frame()
}
for (tf in tfs) {
  key <- paste(tf, "primary_filtered", sep = "__")
  led <- expand_leading_edge(go_results[[key]], rank_data[[key]])
  write_table(led, paste0("07_102_", tf, "_primary_GO_leading_edge_FDR05_v01.csv"))
}

if (compare_with_mef2a) {
  z <- go_results[["Zfp148__primary_filtered"]]
  m <- go_results[["Mef2a__primary_filtered"]]
  keep <- c("ID", "Description", "Ontology", "NES", "pvalue", "p.adjust", "setSize", "core_enrichment")
  z <- z[, keep, drop = FALSE]
  m <- m[, keep, drop = FALSE]
  names(z)[-c(1, 2, 3)] <- paste0("Zfp148_", names(z)[-c(1, 2, 3)])
  names(m)[-c(1, 2, 3)] <- paste0("Mef2a_", names(m)[-c(1, 2, 3)])
  compare <- merge(z, m, by = c("ID", "Description", "Ontology"), all = TRUE, sort = FALSE)
  compare$Zfp148_FDR05 <- !is.na(compare$Zfp148_p.adjust) & compare$Zfp148_p.adjust < fdr_cutoff
  compare$Mef2a_FDR05 <- !is.na(compare$Mef2a_p.adjust) & compare$Mef2a_p.adjust < fdr_cutoff
  compare$is_lipid_or_lysosome_related <- grepl(lipid_pattern, compare$Description, ignore.case = TRUE)
  compare <- compare[order(compare$Zfp148_p.adjust, compare$Mef2a_p.adjust, na.last = TRUE), , drop = FALSE]
  write_table(compare, "08_102_Zfp148_vs_Mef2a_all_GO_primary_v01.csv")
  write_table(compare[compare$is_lipid_or_lysosome_related, , drop = FALSE], "09_102_Zfp148_vs_Mef2a_lipid_related_GO_v01.csv")
}

get_term_members <- function(go_ids, term_tab, ranked) {
  sym <- setNames(ranked$gene, ranked$ENTREZID)
  zs <- setNames(ranked$Z, ranked$ENTREZID)
  members <- tryCatch(suppressMessages(AnnotationDbi::select(db, keys = unique(go_ids), keytype = "GOALL",
    columns = c("ENTREZID", "ONTOLOGYALL"))), error = function(e) NULL)
  out <- vector("list", length(go_ids))
  for (i in seq_along(go_ids)) {
    tid <- go_ids[i]
    tt <- term_tab[term_tab$ID == tid, , drop = FALSE][1, , drop = FALSE]
    if (!is.null(members) && "GOALL" %in% names(members) && "ENTREZID" %in% names(members)) {
      ids <- unique(as.character(members$ENTREZID[!is.na(members$GOALL) & members$GOALL == tid]))
      source <- "GOALL_membership"
    } else {
      ids <- character()
      source <- "core_enrichment_fallback"
    }
    ids <- intersect(ids, names(zs))
    if (length(ids) < 3L) {
      ids <- unlist(strsplit(as.character(tt$core_enrichment), "/", fixed = TRUE))
      ids <- intersect(ids, names(zs))
      source <- "core_enrichment_fallback"
    }
    if (!length(ids))
      next
    core <- unlist(strsplit(as.character(tt$core_enrichment), "/", fixed = TRUE))
    out[[i]] <- data.frame(GO_ID = tid, Ontology = tt$Ontology, Description = tt$Description, NES = tt$NES,
      p.adjust = tt$p.adjust, ENTREZID = ids, gene = unname(sym[ids]), Z = unname(zs[ids]), in_leading_edge = ids %in%
        core, membership_source = source, stringsAsFactors = FALSE)
  }
  out <- Filter(Negate(is.null), out)
  if (!length(out))
    data.frame() else do.call(rbind, out)
}

select_plot_terms <- function(tab) {
  if (!nrow(tab))
    return(list(data = tab, status = "NO_RESULTS"))
  sig <- tab[!is.na(tab$p.adjust) & tab$p.adjust < fdr_cutoff, , drop = FALSE]
  if (nrow(sig)) {
    x <- head(sig[order(sig$p.adjust, -abs(sig$NES)), , drop = FALSE], ridge_top_n)
    return(list(data = x, status = "FDR_SIGNIFICANT"))
  }
  list(data = head(tab[order(tab$p.adjust, tab$pvalue), , drop = FALSE], ridge_top_n), status = "EXPLORATORY_NO_FDR_SIGNIFICANT_TERMS")
}

plot_ridges <- function(tf, tab, ranked, basename) {
  selected <- select_plot_terms(tab)
  terms <- selected$data
  terms$plot_status <- rep(selected$status, nrow(terms))
  write_table(terms, paste0(basename, "_selected_terms_v01.csv"))
  if (nrow(terms) == 0L) {
    message("No GO terms available for ridge plot: ", tf)
    return(invisible(NULL))
  }
  members <- get_term_members(terms$ID, terms, ranked)
  write_table(members, paste0(basename, "_member_gene_Z_v01.csv"))
  if (!nrow(members)) {
    message("No term gene memberships for plot: ", tf)
    return(invisible(NULL))
  }
  terms <- terms[terms$ID %in% members$GO_ID, , drop = FALSE]
  terms <- terms[!duplicated(terms$ID), , drop = FALSE]
  terms$y_index <- rev(seq_len(nrow(terms)))
  terms$logFDR <- -log10(pmax(terms$p.adjust, 1e-15))
  lab <- paste0("GO", terms$Ontology, "_", toupper(gsub("[^[:alnum:]]+", "_", terms$Description)))

  lab <- ifelse(nchar(lab) > 57L, paste0(substr(lab, 1L, 54L), "..."), lab)
  if (anyDuplicated(lab))
    lab <- make.unique(lab, sep = "_")
  terms$label <- lab
  members$y_index <- terms$y_index[match(members$GO_ID, terms$ID)]
  members$logFDR <- terms$logFDR[match(members$GO_ID, terms$ID)]
  members$label <- terms$label[match(members$GO_ID, terms$ID)]
  x0 <- min(ranked$Z, na.rm = TRUE)
  x1 <- max(ranked$Z, na.rm = TRUE)
  dx <- max(x1 - x0, 1)
  dotx <- x0 - 0.07 * dx
  density_rows <- list()
  for (i in seq_len(nrow(terms))) {
    d0 <- members[members$GO_ID == terms$ID[i], , drop = FALSE]
    zs <- d0$Z[is.finite(d0$Z)]
    if (length(zs) >= 3L && length(unique(zs)) >= 2L) {
      den <- stats::density(zs, n = 256L, from = x0, to = x1, na.rm = TRUE)
      h <- den$y/max(den$y)
      density_rows[[length(density_rows) + 1L]] <- data.frame(x = den$x, ymin = terms$y_index[i],
        ymax = terms$y_index[i] + h * 0.8, ID = terms$ID[i], logFDR = terms$logFDR[i], stringsAsFactors = FALSE)
    }
  }
  if (!length(density_rows)) {
    message("Selected terms have too few varying member Z scores for density: ", tf)
    return(invisible(NULL))
  }
  dr <- do.call(rbind, density_rows)
  title <- paste0(tf, " virtual KO | GO-GSEA of perturbation Z")
  subtitle <- if (selected$status == "FDR_SIGNIFICANT") {
    "Top GO terms passing BH-FDR < 0.05 (BP / CC / MF tested separately)"
  } else "EXPLORATORY: NO GO TERMS PASS BH-FDR < 0.05"
  p <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey65",
    linewidth = 0.35) + ggplot2::geom_ribbon(data = dr, ggplot2::aes(x = x, ymin = ymin, ymax = ymax,
    group = ID, fill = logFDR), alpha = 0.82, color = "grey38", linewidth = 0.25) + ggplot2::geom_point(data = terms,
    ggplot2::aes(x = dotx, y = y_index, color = NES, size = logFDR), inherit.aes = FALSE) + ggplot2::scale_color_gradient2(low = "#416FA9",
    mid = "#E6E8EE", high = "#D95536", midpoint = 0, name = "NES (Z ranking)") + ggplot2::scale_fill_gradient(low = "#D5D8EB",
    high = "#E76C4A", name = "-log10(BH-FDR)") + ggplot2::scale_size_continuous(range = c(2.5, 5.5),
    name = "-log10(BH-FDR)") + ggplot2::scale_y_continuous(breaks = terms$y_index, labels = terms$label,
    expand = ggplot2::expansion(mult = c(0.02, 0.07))) + ggplot2::scale_x_continuous(limits = c(dotx -
    0.025 * dx, x1 + 0.04 * dx)) + ggplot2::labs(title = title, subtitle = subtitle, x = "scTenifoldKnk perturbation Z score (not expression fold-change)",
    y = NULL, caption = "Density: GO member gene Z scores. NES reflects ranked perturbation, not activation/repression.") +
    ggplot2::theme_classic(base_size = 11, base_family = "sans") + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = "grey92",
    linewidth = 0.3), plot.title = ggplot2::element_text(face = "bold", size = 14), plot.subtitle = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 8.5, color = "grey20"), axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(8, 8, 8, 8), legend.position = "right", plot.caption = ggplot2::element_text(size = 8,
      hjust = 0))
  pdf_file <- file.path(outdir, paste0(basename, ".pdf"))
  png_file <- file.path(outdir, paste0(basename, ".png"))
  h <- max(5.2, min(11, 3.1 + 0.4 * nrow(terms)))
  ggplot2::ggsave(pdf_file, p, width = 13, height = h, device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(png_file, p, width = 13, height = h, dpi = png_dpi, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(file.path(outdir, paste0(basename, ".svg")), p, width = 13, height = h, device = svglite::svglite)
  }
  message("Ridge plot saved: ", basename)
  invisible(p)
}

for (tf in tfs) {
  key <- paste(tf, "primary_filtered", sep = "__")
  plot_ridges(tf, go_results[[key]], rank_data[[key]], paste0("10_102_", tf, "_Figure6G_GO_GSEA_ridgeplot_v01"))
}

write_lines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    "11_102_REPORT_v01.txt")
outputs <- list.files(outdir, full.names = FALSE)
write_table(data.frame(file = outputs, stringsAsFactors = FALSE), "00_102_OUTPUT_INDEX_v01.csv")

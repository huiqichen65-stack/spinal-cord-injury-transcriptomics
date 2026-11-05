source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)
set.seed(1234L)

ROOT <- PROJECT_ROOT
SCRIPT_NAME <- "129_microglia_pseudotime_sensitivity_v01.R"
OUT <- file.path(ROOT, "results", "129_microglia_pseudotime_sensitivity")
INPUT_RDS <- file.path(ROOT, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
INPUT_IDENTITY <- file.path(ROOT, "results", "106_identity_programs", "02_106_per_cell_identity_scores_v01.csv.gz")
INPUT_ORIGINAL_PT <- file.path(ROOT, "results", "25_Monocle3_all_myeloid_trajectory", "25_Monocle3_cell_pseudotime_v03.csv")
V <- "v01"
SEED <- 1234L
N_PCS <- 30L
N_CORES <- 1L
MIN_CELLS_MAIN <- 500L
MIN_CELLS_REFIT <- 300L
MIN_ANCHOR_CELLS_PER_PARTITION <- 10L
MIN_SHARED_FOR_COMPARISON <- 100L
MIN_OVERLAP_JACCARD <- 0.25

MG_MIN_SCORE <- 0
MG_DOMINANCE_MARGIN <- 0.25

CYCLE_TOP_QUANTILE <- 0.9
MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE <- 10L
SAVE_MAIN_CDS <- TRUE
RUN_LEAVE_ONE_SOURCE_OUT <- TRUE
RUN_WITHIN_SOURCE <- TRUE
RUN_CYCLE_EXCLUSION <- TRUE
SOURCE_ORDER <- c("17702045", "GSE162610", "GSE196928")
STAGE_ORDER <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")
STAGE_DAY <- c(Control = 0, `4 h` = 4/24, `1 dpi` = 1, `3 dpi` = 3, `7 dpi` = 7, `14 dpi` = 14, `28 dpi` = 28,
  `38 dpi` = 38)
MG_STATES <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia",
  "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia")
HOME <- "Homeostatic Microglia"
FOCUS <- "Clec7a+ lipid-associated Microglia"

pkgs <- c("Seurat", "SeuratObject", "Matrix", "monocle3", "SingleCellExperiment", "igraph", "ggplot2")
missing_pkg <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) stop("Missing R packages: ", paste(missing_pkg, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(monocle3)
  library(SingleCellExperiment)
  library(igraph)
  library(ggplot2)
})
for (f in c(INPUT_RDS, INPUT_IDENTITY)) if (!file.exists(f)) stop("Required input missing: ", f)
if (!dir.exists(file.path(ROOT, "script"))) stop("the project script directory directory missing: check ROOT")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "129_microglia_pseudotime_sensitivity_v01.log")
if (file.exists(LOG)) unlink(LOG)
logit <- function(...) {
  x <- paste0(..., collapse = "")
  cat(x, "\n")
  cat(x, "\n", file = LOG, append = TRUE)
}
path_csv <- function(stem, compressed = FALSE) {
  file.path(OUT, paste0(stem, "_", V, if (compressed)
    ".csv.gz" else ".csv"))
}
write_csv <- function(x, stem, compressed = FALSE) {
  f <- path_csv(stem, compressed)
  if (compressed) {
    con <- gzfile(f, "wt")
    on.exit(close(con), add = TRUE)
    utils::write.csv(x, con, row.names = FALSE, na = "")
  } else {
    utils::write.csv(x, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
  logit("Wrote ", basename(f), ": ", nrow(x), " rows")
  invisible(f)
}
rbind_list <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0, x)
  if (!length(x))
    return(data.frame())
  ans <- do.call(rbind, x)
  rownames(ans) <- NULL
  ans
}
stat_pt <- function(z) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  c(n_finite = length(z), median = if (length(z)) stats::median(z) else NA_real_, q25 = if (length(z)) unname(stats::quantile(z,
    0.25)) else NA_real_, q75 = if (length(z)) unname(stats::quantile(z, 0.75)) else NA_real_)
}
scale_rank <- function(v) {
  v <- as.numeric(v)
  good <- which(is.finite(v))
  out <- rep(NA_real_, length(v))
  if (length(good) >= 2L)
    out[good] <- (rank(v[good], ties.method = "average") - 1)/(length(good) - 1)
  out
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < MIN_SHARED_FOR_COMPARISON || length(unique(x[ok])) < 3L || length(unique(y[ok])) <
    3L)
    return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
logit("Step 129 v01 started ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logit("INPUT_RDS: ", INPUT_RDS)
logit("INPUT_IDENTITY: ", INPUT_IDENTITY)
logit("Output dir: ", OUT)
logit("Monocle3 version: ", as.character(utils::packageVersion("monocle3")))

obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat") || !"RNA" %in% Seurat::Assays(obj)) stop("Step18 input not RNA Seurat")
md <- obj[[]]
need_md <- c("source_dataset", "timepoint_main", "myeloid_state")
if (length(setdiff(need_md, names(md)))) stop("Missing Step18 metadata: ", paste(setdiff(need_md, names(md)),
  collapse = ", "))
md$cell_id <- rownames(md)
if (anyDuplicated(md$cell_id)) stop("Duplicate Step18 cell IDs")
for (z in need_md) md[[z]] <- as.character(md[[z]])
if (anyNA(md$source_dataset) || anyNA(md$timepoint_main) || anyNA(md$myeloid_state)) stop("Missing required source/stage/state")
if (any(!md$source_dataset %in% SOURCE_ORDER)) stop("Unexpected source: ", paste(setdiff(unique(md$source_dataset),
  SOURCE_ORDER), collapse = ", "))
if (any(!md$timepoint_main %in% STAGE_ORDER)) stop("Unexpected stage metadata")
if (!"sample_id" %in% names(md)) md$sample_id <- NA_character_
if (!"treatment" %in% names(md)) md$treatment <- NA_character_
md$sample_id <- as.character(md$sample_id)
md$treatment <- as.character(md$treatment)
if (any(grepl("PLX5622", md$treatment, ignore.case = TRUE), na.rm = TRUE)) {
  stop("PLX5622 samples present in retained Step18 cohort: resolve treatment confounding first")
}
identity <- utils::read.csv(gzfile(INPUT_IDENTITY), check.names = FALSE, stringsAsFactors = FALSE)
need_scores <- c("cell_id", "source_dataset", "timepoint_main", "myeloid_state", "eligible_identity_scoring",
  "MG_score", "MDM_score", "BAM_score")
if (length(setdiff(need_scores, names(identity)))) stop("Step106 missing columns: ", paste(setdiff(need_scores,
  names(identity)), collapse = ", "))
if (anyDuplicated(identity$cell_id)) stop("Duplicate Step106 cell IDs")
if (!setequal(identity$cell_id, md$cell_id)) stop("Step106 and Step18 cell identifiers differ: abort rather than silently subset")
identity <- identity[match(md$cell_id, identity$cell_id), , drop = FALSE]
for (z in c("source_dataset", "timepoint_main", "myeloid_state")) {
  if (!identical(as.character(identity[[z]]), as.character(md[[z]])))
    stop("Step106 vs Step18 metadata mismatch: ", z)
}
for (z in c("MG_score", "MDM_score", "BAM_score")) identity[[z]] <- suppressWarnings(as.numeric(identity[[z]]))
md$MG_score <- identity$MG_score
md$MDM_score <- identity$MDM_score
md$BAM_score <- identity$BAM_score
md$is_MG_label <- md$myeloid_state %in% MG_STATES
if (sum(md$is_MG_label) < MIN_CELLS_MAIN) stop("Fewer annotated MG than expected")
md$score_available <- is.finite(md$MG_score) & is.finite(md$MDM_score) & is.finite(md$BAM_score) & !is.na(identity$eligible_identity_scoring) &
  identity$eligible_identity_scoring
md$high_confidence <- md$is_MG_label & md$score_available & md$MG_score >= MG_MIN_SCORE & md$MG_score -
  md$MDM_score >= MG_DOMINANCE_MARGIN & md$MG_score - md$BAM_score >= MG_DOMINANCE_MARGIN
md$selection_reason <- ifelse(!md$is_MG_label, "not_annotated_MG", ifelse(!md$score_available, "identity_score_missing_or_RNA_zero",
  ifelse(md$MG_score < MG_MIN_SCORE, "MG_score_below_zero", ifelse(md$MG_score - md$MDM_score < MG_DOMINANCE_MARGIN,
    "MG_vs_MDM_margin_fail", ifelse(md$MG_score - md$BAM_score < MG_DOMINANCE_MARGIN, "MG_vs_BAM_margin_fail",
      "HC_MG_consistent")))))

Seurat::DefaultAssay(obj) <- "RNA"
obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) stop("RNA JoinLayers required (no silent fallback): ",
  conditionMessage(e)))
cts <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) stop("RNA counts not accessible: ",
  conditionMessage(e)))
if (!inherits(cts, "sparseMatrix")) stop("RNA counts not sparse: stop to avoid memory explosion")
if (!setequal(colnames(cts), md$cell_id)) stop("RNA count cells differ from Step18 metadata")
cts <- cts[, md$cell_id, drop = FALSE]
libsize <- as.numeric(Matrix::colSums(cts))
md$positive_RNA <- is.finite(libsize) & libsize > 0
md$high_confidence <- md$high_confidence & md$positive_RNA
md$selection_reason[!md$positive_RNA & md$is_MG_label] <- "RNA_zero_or_nonfinite"
if (sum(md$high_confidence) < MIN_CELLS_MAIN) stop("Too few HC microglia; do not silently relax thresholds")
if (!any(md$high_confidence & md$myeloid_state == FOCUS)) stop("No retained Clec7a+ cells; stop")
logit("Step18 cells: ", nrow(md), "; MG annotated: ", sum(md$is_MG_label), "; HC candidate: ", sum(md$high_confidence),
  "; HC Clec7a+: ", sum(md$high_confidence & md$myeloid_state == FOCUS))

audit <- data.frame(audit_item = c("Step18_total", "MG_label_total", "identity_score_available_MG", "zero_RNA_total",
  "high_confidence_MG", "HC_Clec7a", "HC_MG_score_min", "HC_MG_vs_MDM_margin_min", "HC_MG_vs_BAM_margin_min"),
  value = as.character(c(nrow(md), sum(md$is_MG_label), sum(md$is_MG_label & md$score_available), sum(!md$positive_RNA),
    sum(md$high_confidence), sum(md$high_confidence & md$myeloid_state == FOCUS), MG_MIN_SCORE, MG_DOMINANCE_MARGIN,
    MG_DOMINANCE_MARGIN)))
write_csv(audit, "01_129_input_and_cohort_audit")
selection_export <- md[, c("cell_id", "source_dataset", "timepoint_main", "sample_id", "treatment", "myeloid_state",
  "is_MG_label", "score_available", "MG_score", "MDM_score", "BAM_score", "positive_RNA", "high_confidence",
  "selection_reason"), drop = FALSE]
write_csv(selection_export, "02_129_cell_inclusion_audit", compressed = TRUE)

coverage <- as.data.frame(table(source_dataset = md$source_dataset, timepoint_main = factor(md$timepoint_main,
  levels = STAGE_ORDER), myeloid_state = md$myeloid_state, high_confidence = md$high_confidence), stringsAsFactors = FALSE)
names(coverage)[4] <- "high_confidence"
coverage <- coverage[coverage$Freq > 0, , drop = FALSE]
write_csv(coverage, "03_129_source_stage_state_coverage")

cc_env <- new.env(parent = emptyenv())
utils::data("cc.genes.updated.2019", package = "Seurat", envir = cc_env)
if (!exists("cc.genes.updated.2019", envir = cc_env, inherits = FALSE)) {
  stop("Seurat::cc.genes.updated.2019 is unavailable; cannot audit cell cycle")
}
cc <- get("cc.genes.updated.2019", envir = cc_env)
if (!all(c("s.genes", "g2m.genes") %in% names(cc))) stop("Unexpected cell-cycle gene set format")
map_gene <- function(z) {
  hit <- match(toupper(z), toupper(rownames(cts)))
  rownames(cts)[hit]
}
s_genes <- unique(stats::na.omit(map_gene(cc$s.genes)))
g_genes <- unique(stats::na.omit(map_gene(cc$g2m.genes)))
if (length(s_genes) < MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE || length(g_genes) < MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE) {
  stop("Too few cell-cycle genes matched to mouse RNA features")
}
cycle_genes <- unique(c(s_genes, g_genes))
manifest <- data.frame(program = c(rep("S", length(cc$s.genes)), rep("G2M", length(cc$g2m.genes))), reference_gene = c(cc$s.genes,
  cc$g2m.genes), matched_RNA_gene = c(map_gene(cc$s.genes), map_gene(cc$g2m.genes)))
write_csv(manifest, "04_129_cell_cycle_gene_manifest")
hc_md <- md[md$high_confidence, , drop = FALSE]
hc_cells <- hc_md$cell_id
hc_lib <- libsize[match(hc_cells, md$cell_id)]
cycle_counts <- cts[cycle_genes, hc_cells, drop = FALSE]

cycle_log <- as.matrix(cycle_counts)
cycle_log <- log1p(sweep(cycle_log, 2, hc_lib/10000, "/"))
rm(cycle_counts)
cycle_S <- rep(NA_real_, length(hc_cells))
cycle_G2M <- rep(NA_real_, length(hc_cells))
cycle_nS <- integer(length(hc_cells))
cycle_nG <- integer(length(hc_cells))
for (src in SOURCE_ORDER) {
  ii <- which(hc_md$source_dataset == src)
  if (!length(ii))
    next
  nz <- rowSums(cycle_log[, ii, drop = FALSE] > 0)
  eligible <- names(nz)[nz >= 10L]
  usableS <- intersect(s_genes, eligible)
  usableG <- intersect(g_genes, eligible)
  cycle_nS[ii] <- length(usableS)
  cycle_nG[ii] <- length(usableG)
  if (length(usableS) < MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE || length(usableG) < MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE) {
    logit("Cell-cycle scores unavailable for ", src, ": S=", length(usableS), ", G2M=", length(usableG))
    next
  }
  geneZ <- function(genes) {
    x <- cycle_log[genes, ii, drop = FALSE]
    mu <- rowMeans(x)
    sdv <- apply(x, 1, stats::sd)
    ok <- is.finite(sdv) & sdv > 1e-08
    if (sum(ok) < MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE)
      return(rep(NA_real_, length(ii)))
    x <- sweep(x[ok, , drop = FALSE], 1, mu[ok], "-")
    x <- sweep(x, 1, sdv[ok], "/")
    colMeans(pmax(pmin(x, 3), -3))
  }
  cycle_S[ii] <- geneZ(usableS)
  cycle_G2M[ii] <- geneZ(usableG)
}
rm(cycle_log)
invisible(gc(verbose = FALSE))
hc_md$S_program_z <- cycle_S
hc_md$G2M_program_z <- cycle_G2M
hc_md$cycle_program_score <- pmax(cycle_S, cycle_G2M)
hc_md$cycle_score_available <- is.finite(hc_md$cycle_program_score)
hc_md$cycle_source_q90 <- NA_real_
hc_md$high_cycle_program <- FALSE
for (src in SOURCE_ORDER) {
  ii <- which(hc_md$source_dataset == src & hc_md$cycle_score_available)
  if (length(ii) < 50L)
    next
  cut <- unname(stats::quantile(hc_md$cycle_program_score[ii], CYCLE_TOP_QUANTILE))
  hc_md$cycle_source_q90[ii] <- cut
  hc_md$high_cycle_program[ii] <- hc_md$cycle_program_score[ii] >= cut
}
hc_md$annotated_proliferating <- hc_md$myeloid_state == "Proliferating Microglia"
hc_md$cycle_exclusion <- hc_md$annotated_proliferating | hc_md$high_cycle_program
cycle_ok_all_sources <- all(vapply(SOURCE_ORDER, function(src) {
  x <- hc_md[hc_md$source_dataset == src, , drop = FALSE]
  nrow(x) > 0 && mean(x$cycle_score_available) > 0.95 && all(cycle_nS[hc_md$source_dataset == src] >=
    MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE) && all(cycle_nG[hc_md$source_dataset == src] >= MIN_CYCLE_GENES_PER_PROGRAM_PER_SOURCE)
}, logical(1)))
cycle_export <- hc_md[, c("cell_id", "source_dataset", "timepoint_main", "myeloid_state", "S_program_z",
  "G2M_program_z", "cycle_program_score", "cycle_score_available", "cycle_source_q90", "high_cycle_program",
  "annotated_proliferating", "cycle_exclusion"), drop = FALSE]
write_csv(cycle_export, "05_129_cell_cycle_scores", compressed = TRUE)
logit("Cell-cycle scoring complete; fraction available: ", round(mean(hc_md$cycle_score_available), 3),
  "; removed in cycle branch: ", sum(hc_md$cycle_exclusion), "; all sources sufficiently scored: ",
  cycle_ok_all_sources)

rm(obj, identity)
invisible(gc(verbose = FALSE))

choose_roots <- function(node_df, graph_obj, scheme) {
  parts <- sort(unique(node_df$partition))
  rows <- list()
  roots <- character()
  for (part in parts) {
    x <- node_df[node_df$partition == part, , drop = FALSE]
    ax <- x[x$myeloid_state == HOME & x$timepoint_main == "Control", , drop = FALSE]
    chosen <- NA_character_
    rule <- "no_applicable_reference"
    n_candidate <- 0L
    if (scheme == "control_homeostatic") {
      if (nrow(ax) >= MIN_ANCHOR_CELLS_PER_PARTITION) {
        tab <- sort(table(ax$principal_node), decreasing = TRUE)
        chosen <- names(tab)[1]
        n_candidate <- nrow(ax)
        rule <- "Control_plus_Homeostatic_at_least_10"
      }
    } else if (scheme == "all_homeostatic") {
      y <- x[x$myeloid_state == HOME, , drop = FALSE]
      if (nrow(y) >= MIN_ANCHOR_CELLS_PER_PARTITION) {
        tab <- sort(table(y$principal_node), decreasing = TRUE)
        chosen <- names(tab)[1]
        n_candidate <- nrow(y)
        rule <- "All_stage_Homeostatic_at_least_10"
      }
    } else if (scheme == "earliest_in_partition") {
      ranks <- STAGE_DAY[x$timepoint_main]
      if (any(is.finite(ranks))) {
        y <- x[is.finite(ranks) & ranks == min(ranks, na.rm = TRUE), , drop = FALSE]
        tab <- sort(table(y$principal_node), decreasing = TRUE)
        chosen <- names(tab)[1]
        n_candidate <- nrow(y)
        rule <- paste0("Earliest_available_stage_", paste(unique(y$timepoint_main), collapse = ";"))
      }
    } else if (scheme == "farthest_from_control_anchor") {
      if (nrow(ax) >= MIN_ANCHOR_CELLS_PER_PARTITION) {
        tab <- sort(table(ax$principal_node), decreasing = TRUE)
        anchor <- names(tab)[1]
        nodes <- intersect(unique(x$principal_node), igraph::V(graph_obj)$name)
        dist <- as.numeric(igraph::distances(graph_obj, v = anchor, to = nodes, weights = NA)[1,
          ])
        if (any(is.finite(dist) & dist > 0)) {
          chosen <- nodes[which.max(replace(dist, !is.finite(dist), -Inf))]
          n_candidate <- nrow(ax)
          rule <- "Topology_stress_test_not_a_biological_root"
        }
      }
    } else stop("Unknown root scheme")
    if (!is.na(chosen))
      roots <- c(roots, chosen)
    rows[[length(rows) + 1L]] <- data.frame(partition = part, scheme = scheme, root_node = chosen,
      selection_rule = rule, n_partition_cells = nrow(x), n_control_homeostatic = nrow(ax), n_candidate_cells = n_candidate,
      n_root_node_cells = if (is.na(chosen))
        0L else sum(x$principal_node == chosen), status = if (is.na(chosen))
        "UNROOTED_NO_DEFENSIBLE_ANCHOR" else "ROOT_SELECTED")
  }
  list(roots = unique(roots), audit = rbind_list(rows))
}

extract_nodes <- function(cds) {
  graph <- monocle3::principal_graph(cds)[["UMAP"]]
  if (is.null(graph) || igraph::vcount(graph) < 1)
    stop("Principal graph missing")
  part <- cds@clusters[["UMAP"]]$partitions
  if (is.null(names(part)))
    names(part) <- colnames(cds)
  part <- as.character(part[colnames(cds)])
  proj <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(proj))
    stop("Cell-to-graph node projection missing")
  node_raw <- as.character(as.matrix(proj[colnames(cds), , drop = FALSE])[, 1])
  v <- as.character(igraph::V(graph)$name)
  if (all(node_raw %in% v)) {
    node <- node_raw
  } else {
    k <- suppressWarnings(as.integer(node_raw))
    if (anyNA(k) || any(k < 1L | k > length(v)))
      stop("Graph node IDs cannot be resolved")
    node <- v[k]
  }
  if (anyNA(part) || anyNA(node))
    stop("Partition/node assignment incomplete")
  list(node = node, part = part, graph = graph)
}

fit_trajectory <- function(branch_id, cells, matrix_counts, cellmd) {
  logit("Fitting ", branch_id, ", n=", length(cells))
  if (length(cells) < MIN_CELLS_REFIT)
    stop("Too few cells in branch")
  m <- matrix_counts[, cells, drop = FALSE]
  totals <- as.numeric(Matrix::colSums(m))
  if (any(!is.finite(totals) | totals <= 0))
    stop("Zero-RNA cells unexpectedly retained")
  meta <- cellmd[match(cells, cellmd$cell_id), , drop = FALSE]
  if (anyNA(meta$cell_id))
    stop("Cell metadata mismatch for branch ", branch_id)
  rownames(meta) <- meta$cell_id
  gmeta <- data.frame(gene_short_name = rownames(m), row.names = rownames(m))
  cds <- monocle3::new_cell_data_set(m, cell_metadata = meta, gene_metadata = gmeta)
  cds <- monocle3::estimate_size_factors(cds)
  if (any(!is.finite(monocle3::size_factors(cds)) | monocle3::size_factors(cds) <= 0))
    stop("Invalid Monocle3 size factors")
  set.seed(SEED)
  cds <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = N_PCS, norm_method = "log")
  set.seed(SEED)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", umap.fast_sgd = FALSE,
    cores = N_CORES)
  set.seed(SEED)
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
  set.seed(SEED)
  cds <- monocle3::learn_graph(cds, use_partition = TRUE, learn_graph_control = list(euclidean_distance_ratio = 0.8))
  nd <- extract_nodes(cds)
  node_df <- data.frame(cell_id = colnames(cds), partition = nd$part, principal_node = nd$node, myeloid_state = as.character(SummarizedExperiment::colData(cds)$myeloid_state),
    timepoint_main = as.character(SummarizedExperiment::colData(cds)$timepoint_main), source_dataset = as.character(SummarizedExperiment::colData(cds)$source_dataset))
  schemes <- c("control_homeostatic", "all_homeostatic", "earliest_in_partition", "farthest_from_control_anchor")
  root_rows <- list()
  pt_rows <- list()
  primary_cds <- NULL
  for (scheme in schemes) {
    root_info <- choose_roots(node_df, nd$graph, scheme)
    root_info$audit$branch <- branch_id
    root_rows[[scheme]] <- root_info$audit
    if (length(root_info$roots) == 0) {
      pt <- rep(NA_real_, nrow(node_df))
    } else {

      ordered <- monocle3::order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_info$roots)
      pt_raw <- monocle3::pseudotime(ordered)
      pt <- as.numeric(pt_raw[node_df$cell_id])
      pt[!is.finite(pt)] <- NA_real_
      rooted <- root_info$audit$partition[!is.na(root_info$audit$root_node)]
      pt[!node_df$partition %in% rooted] <- NA_real_
      if (scheme == "control_homeostatic")
        primary_cds <- ordered
    }
    t <- node_df[, c("cell_id", "partition", "principal_node", "myeloid_state", "timepoint_main",
      "source_dataset"), drop = FALSE]
    t$branch <- branch_id
    t$root_scheme <- scheme
    t$pseudotime <- pt
    t$rank_within_partition <- NA_real_
    for (part in unique(t$partition)) {
      ii <- which(t$partition == part)
      t$rank_within_partition[ii] <- scale_rank(t$pseudotime[ii])
    }
    pt_rows[[scheme]] <- t
  }
  list(pt = rbind_list(pt_rows), roots = rbind_list(root_rows), main_cds = primary_cds, cds = cds)
}

branch_pt <- list()
branch_roots <- list()
branch_status <- list()
main_fit <- NULL
fit_one <- function(id, cells, essential = FALSE) {
  t0 <- Sys.time()
  x <- tryCatch(fit_trajectory(id, cells, cts, hc_md), error = function(e) e)
  if (inherits(x, "error")) {
    logit("BRANCH ERROR: ", id, ": ", conditionMessage(x))
    branch_status[[id]] <<- data.frame(branch = id, n_requested = length(cells), n_rooted_primary = NA_integer_,
      n_partition_primary = NA_integer_, status = "FAILED", notes = conditionMessage(x), elapsed_sec = as.numeric(difftime(Sys.time(),
        t0, units = "secs")))
    if (essential)
      stop("Primary HC_MG trajectory failed: ", conditionMessage(x))
    return(invisible(NULL))
  }
  branch_pt[[id]] <<- x$pt
  branch_roots[[id]] <<- x$roots
  x_main <- x$pt[x$pt$root_scheme == "control_homeostatic", , drop = FALSE]
  branch_status[[id]] <<- data.frame(branch = id, n_requested = length(cells), n_rooted_primary = sum(is.finite(x_main$pseudotime)),
    n_partition_primary = length(unique(x_main$partition)), status = "COMPLETE", notes = if (all(!is.finite(x_main$pseudotime)))
      "No Control+Homeostatic anchor: no primary clock" else "", elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (identical(id, "HC_main")) {
    main_fit <<- list(cds = x$main_cds, graph = x$cds)
    if (SAVE_MAIN_CDS && !is.null(x$main_cds)) {
      saveRDS(x$main_cds, file.path(OUT, "129_HC_microglia_main_v01.rds"), compress = TRUE)
      logit("Saved primary rooted Monocle3 CDS (main branch only)")
    }
  }
  invisible(TRUE)
}
fit_one("HC_main", hc_cells, essential = TRUE)
if (RUN_CYCLE_EXCLUSION) {
  if (!cycle_ok_all_sources) {
    branch_status[["HC_cycle_excluded"]] <- data.frame(branch = "HC_cycle_excluded", n_requested = NA_integer_,
      n_rooted_primary = NA_integer_, n_partition_primary = NA_integer_, status = "SKIPPED", notes = "Cell-cycle gene coverage insufficient in >=1 source; no unverifiable sensitivity run",
      elapsed_sec = NA_real_)
  } else {
    keep <- hc_md$cell_id[!hc_md$cycle_exclusion]
    fit_one("HC_cycle_excluded", keep)
  }
}
if (RUN_LEAVE_ONE_SOURCE_OUT) for (src in SOURCE_ORDER) {
  keep <- hc_md$cell_id[hc_md$source_dataset != src]
  fit_one(paste0("HC_leave_out_", src), keep)
}
if (RUN_WITHIN_SOURCE) for (src in SOURCE_ORDER) {
  keep <- hc_md$cell_id[hc_md$source_dataset == src]
  fit_one(paste0("HC_source_only_", src), keep)
}
write_csv(rbind_list(branch_status), "06_129_branch_run_status")
roots_all <- rbind_list(branch_roots)
write_csv(roots_all, "07_129_partition_and_root_audit")
pt_all <- rbind_list(branch_pt)
write_csv(pt_all, "08_129_cell_pseudotime_all_roots", compressed = TRUE)

summary_rows <- list()
unique_combos <- unique(pt_all[, c("branch", "root_scheme", "partition"), drop = FALSE])
for (combo_i in seq_len(nrow(unique_combos))) {
  combo <- unique_combos[combo_i, , drop = FALSE]
  x <- pt_all[pt_all$branch == combo$branch & pt_all$root_scheme == combo$root_scheme & pt_all$partition ==
    combo$partition, , drop = FALSE]
  grouping <- list(state = "myeloid_state", stage = "timepoint_main", source = "source_dataset")
  for (group in names(grouping)) {
    spl <- split(seq_len(nrow(x)), x[[grouping[[group]]]])
    for (nm in names(spl)) {
      z <- x$pseudotime[spl[[nm]]]
      st <- stat_pt(z)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(branch = combo$branch, root_scheme = combo$root_scheme,
        partition = combo$partition, grouping = group, group_label = nm, n_cells = length(z),
        n_finite = unname(st["n_finite"]), median_pseudotime = unname(st["median"]), q25_pseudotime = unname(st["q25"]),
        q75_pseudotime = unname(st["q75"]))
    }
  }
}
write_csv(rbind_list(summary_rows), "09_129_partition_specific_pseudotime_summary")

root_cmp <- list()
for (branch in names(branch_pt)) {
  tb <- branch_pt[[branch]]
  base <- tb[tb$root_scheme == "control_homeostatic", , drop = FALSE]
  alts <- setdiff(unique(tb$root_scheme), "control_homeostatic")
  for (alt_scheme in alts) {
    alt <- tb[tb$root_scheme == alt_scheme, , drop = FALSE]
    stopifnot(identical(base$cell_id, alt$cell_id), identical(base$partition, alt$partition))
    for (part in unique(base$partition)) {
      ii <- which(base$partition == part)
      good <- is.finite(base$pseudotime[ii]) & is.finite(alt$pseudotime[ii])
      root_cmp[[length(root_cmp) + 1L]] <- data.frame(branch = branch, partition = part, baseline_root = "control_homeostatic",
        alternative_root = alt_scheme, n_partition = length(ii), n_matched_finite = sum(good),
        spearman_within_same_graph = safe_cor(base$pseudotime[ii], alt$pseudotime[ii]), note = if (sum(good) <
          MIN_SHARED_FOR_COMPARISON)
          "UNESTIMABLE_or_too_few_rooted_cells" else "Descriptive root sensitivity; not a time-direction test")
    }
  }
}
write_csv(rbind_list(root_cmp), "10_129_root_choice_sensitivity")

base <- branch_pt[["HC_main"]]
base <- base[base$root_scheme == "control_homeostatic", , drop = FALSE]
overlap_rows <- list()
refit_rows <- list()
for (branch in setdiff(names(branch_pt), "HC_main")) {
  alt <- branch_pt[[branch]]
  alt <- alt[alt$root_scheme == "control_homeostatic", , drop = FALSE]
  for (pa in unique(base$partition)) for (pb in unique(alt$partition)) {
    ca <- base$cell_id[base$partition == pa]
    cb <- alt$cell_id[alt$partition == pb]
    nshare <- sum(ca %in% cb)
    ja <- nshare/(length(ca) + length(cb) - nshare)
    overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(branch = branch, main_partition = pa,
      refit_partition = pb, n_main_partition = length(ca), n_refit_partition = length(cb), n_shared = nshare,
      jaccard = ja)
  }
  for (pa in unique(base$partition)) {
    pairs <- Filter(function(z) z$branch == branch && z$main_partition == pa, overlap_rows)
    z <- rbind_list(pairs)
    if (!nrow(z))
      next
    z <- z[order(-z$jaccard, -z$n_shared), , drop = FALSE]
    best <- z[1, , drop = FALSE]
    a <- base[base$partition == pa, , drop = FALSE]
    b <- alt[alt$partition == best$refit_partition, , drop = FALSE]
    ix <- match(a$cell_id, b$cell_id)
    shared <- which(!is.na(ix))
    x <- a$rank_within_partition[shared]
    y <- b$rank_within_partition[ix[shared]]
    finite <- is.finite(x) & is.finite(y)
    eligible <- best$jaccard >= MIN_OVERLAP_JACCARD && sum(finite) >= MIN_SHARED_FOR_COMPARISON
    refit_rows[[length(refit_rows) + 1L]] <- data.frame(branch = branch, main_partition = pa, matched_refit_partition = best$refit_partition,
      jaccard = best$jaccard, n_shared = length(shared), n_both_rooted = sum(finite), rank_spearman = if (eligible)
        safe_cor(x, y) else NA_real_, status = if (eligible)
        "DESCRIPTIVE_COMPARABLE_PARTITIONS" else "NONCOMPARABLE_OR_INSUFFICIENT_ANCHOR", note = "Ranks are recomputed within their own partitions; not global pseudotime")
  }
}
write_csv(rbind_list(overlap_rows), "11_129_branch_partition_overlap")
write_csv(rbind_list(refit_rows), "12_129_refit_sensitivity")

orig_audit <- data.frame(status = "OPTIONAL_REFERENCE_NOT_PRESENT", n_original = NA_integer_, n_shared = NA_integer_,
  n_rooted_shared = NA_integer_, spearman_raw = NA_real_, notes = "No comparison attempted")
if (file.exists(INPUT_ORIGINAL_PT)) {
  orig <- utils::read.csv(INPUT_ORIGINAL_PT, stringsAsFactors = FALSE)
  if (all(c("cell", "partition", "pseudotime") %in% names(orig)) && !anyDuplicated(orig$cell)) {
    j <- match(base$cell_id, orig$cell)
    good <- which(!is.na(j) & is.finite(base$pseudotime))

    original_partition <- as.character(orig$partition[j])
    best_part <- if (length(good))
      names(sort(table(original_partition[good]), decreasing = TRUE))[1] else NA_character_
    idx <- good[original_partition[good] == best_part]
    co <- if (length(idx))
      safe_cor(base$pseudotime[idx], orig$pseudotime[j[idx]]) else NA_real_
    orig_audit <- data.frame(status = "DESCRIPTIVE_ONLY_DIFFERENT_GRAPHS", n_original = nrow(orig),
      n_shared = sum(!is.na(j)), n_rooted_shared = length(idx), spearman_raw = co, notes = paste0("Original partition ",
        best_part, " selected for within-original-partition descriptive check; do not infer shared clock"))
  } else orig_audit$status <- "OPTIONAL_REFERENCE_INVALID_SKIPPED"
}
write_csv(orig_audit, "13_129_original_vs_restricted_audit")

assoc_rows <- list()
for (branch in names(branch_pt)) {
  bx <- branch_pt[[branch]]
  bx <- bx[bx$root_scheme == "control_homeostatic", , drop = FALSE]
  cyc_ix <- match(bx$cell_id, hc_md$cell_id)
  bx$cycle_score <- hc_md$cycle_program_score[cyc_ix]
  bx$stage_day <- unname(STAGE_DAY[bx$timepoint_main])
  for (part in unique(bx$partition)) for (src in c("ALL_SOURCES_DESCRIPTIVE", SOURCE_ORDER)) {
    z <- bx[bx$partition == part & (src == "ALL_SOURCES_DESCRIPTIVE" | bx$source_dataset == src),
      , drop = FALSE]
    if (!nrow(z))
      next
    n_joint <- sum(is.finite(z$pseudotime) & is.finite(z$cycle_score))
    assoc_rows[[length(assoc_rows) + 1L]] <- data.frame(branch = branch, partition = part, source_dataset = src,
      n_cells = nrow(z), n_pseudotime_finite = sum(is.finite(z$pseudotime)), n_cycle_and_pseudotime_finite = n_joint,
      n_distinct_stages = length(unique(z$timepoint_main[is.finite(z$pseudotime)])), spearman_cycle_vs_pseudotime = safe_cor(z$pseudotime,
        z$cycle_score), spearman_stage_day_vs_pseudotime = safe_cor(z$pseudotime, z$stage_day),
      interpretation = "Cell-level descriptive Spearman ONLY; pooled stages may be source-confounded")
  }
}
write_csv(rbind_list(assoc_rows), "17_129_cycle_stage_source_diagnostics")

pdf_path <- file.path(OUT, "14_129_trajectory_diagnostics_v01.pdf")
grDevices::pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
tryCatch({
  if (!is.null(main_fit$graph)) {
    coords <- as.data.frame(SingleCellExperiment::reducedDim(main_fit$graph, "UMAP"))
    names(coords)[1:2] <- c("UMAP_1", "UMAP_2")
    coords$cell_id <- rownames(coords)
    d <- merge(coords, base, by = "cell_id", sort = FALSE)
    d$timepoint_main <- factor(d$timepoint_main, levels = STAGE_ORDER)
    print(ggplot2::ggplot(d, ggplot2::aes(UMAP_1, UMAP_2, color = myeloid_state)) + ggplot2::geom_point(size = 0.16,
      alpha = 0.65) + ggplot2::facet_wrap(~partition, scales = "free") + ggplot2::theme_bw(base_size = 10) +
      ggplot2::labs(title = "129 HC microglia-only: de novo Monocle3 graph partitions", subtitle = "Partitions are separate; no common cross-partition clock",
        color = "State"))
    print(ggplot2::ggplot(d, ggplot2::aes(UMAP_1, UMAP_2, color = source_dataset)) + ggplot2::geom_point(size = 0.16,
      alpha = 0.65) + ggplot2::facet_wrap(~partition, scales = "free") + ggplot2::theme_bw(base_size = 10) +
      ggplot2::labs(title = "129 HC microglia-only: source contribution", color = "Source"))
    print(ggplot2::ggplot(d, ggplot2::aes(UMAP_1, UMAP_2, color = pseudotime)) + ggplot2::geom_point(size = 0.16,
      alpha = 0.65) + ggplot2::facet_wrap(~partition, scales = "free") + ggplot2::scale_color_viridis_c(na.value = "grey85") +
      ggplot2::theme_bw(base_size = 10) + ggplot2::labs(title = "129 HC microglia-only: anchored pseudotime",
      subtitle = "Within-partition distances only; unrooted partitions shown grey"))
    dd <- d[is.finite(d$pseudotime), , drop = FALSE]
    if (nrow(dd)) {
      print(ggplot2::ggplot(dd, ggplot2::aes(x = timepoint_main, y = pseudotime, fill = source_dataset)) +
        ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::facet_wrap(~partition, scales = "free_y") +
        ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(title = "Stage distributions within rooted partitions (descriptive only)",
        x = "Injury stage", y = "Partition-specific pseudotime"))
      print(ggplot2::ggplot(dd, ggplot2::aes(x = myeloid_state, y = pseudotime, fill = source_dataset)) +
        ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::facet_wrap(~partition, scales = "free_y") +
        ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(title = "State distributions within rooted partitions",
        x = "State", y = "Partition-specific pseudotime"))
    }
    cycle_plot <- hc_md[, c("cell_id", "cycle_program_score", "cycle_exclusion"), drop = FALSE]
    cp <- merge(d, cycle_plot, by = "cell_id", sort = FALSE)
    cp <- cp[is.finite(cp$pseudotime) & is.finite(cp$cycle_program_score), , drop = FALSE]
    if (nrow(cp))
      print(ggplot2::ggplot(cp, ggplot2::aes(pseudotime, cycle_program_score)) + ggplot2::geom_bin2d(bins = 35) +
        ggplot2::facet_wrap(~partition, scales = "free") + ggplot2::theme_bw(base_size = 10) +
        ggplot2::labs(title = "Cell-cycle program score along HC pseudotime", subtitle = "Cell-level descriptive association; not biological replicates"))
  }
}, error = function(e) logit("Plot warning: ", conditionMessage(e)))
grDevices::dev.off()
logit("Wrote ", basename(pdf_path))

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, "15_129_methods_and_limitations_v01.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUT, "16_129_sessionInfo_v01.txt"))
all_out <- list.files(OUT, full.names = TRUE, recursive = FALSE)
info <- file.info(all_out)
index <- data.frame(filename = basename(all_out), file_bytes = as.numeric(info$size), purpose = ifelse(grepl("\\.rds$",
  all_out, ignore.case = TRUE), "Optional main-CDS object; not required for interpreting CSV evidence",
  "Step129 versioned audit/result"))
write_csv(index, "00_129_output_index")
logit("Step129 COMPLETE. Outputs in: ", OUT)

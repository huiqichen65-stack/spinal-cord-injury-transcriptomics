source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)
set.seed(1234L)
ROOT <- PROJECT_ROOT
V <- "v04"
OUT <- file.path(ROOT, "results", "133_Harmony_matched_trajectory")
IN18 <- file.path(ROOT, "results", "18_myeloid_annotation_cleanup", "18_myeloid_annotation_cleanup_3scRNA_v01.rds")
IN129 <- file.path(ROOT, "results", "129_microglia_pseudotime_sensitivity")
IN_COHORT <- file.path(IN129, "02_129_cell_inclusion_audit_v01.csv.gz")
IN_CYCLE <- file.path(IN129, "05_129_cell_cycle_scores_v01.csv.gz")
IN_OLDPT <- file.path(IN129, "08_129_cell_pseudotime_all_roots_v01.csv.gz")
SEED <- 1234L
N_HVG <- 3000L
N_PCA_COMPUTE <- 50L
N_DIM <- 30L
HARMONY_ITER <- 20L
UMAP_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.3
MIN_CELLS <- 300L
MIN_ROOT_ANCHOR <- 10L
MIN_MATCHED <- 100L
JACCARD_CUTOFF <- 0.25
RUN_CYCLE_EXCLUSION <- TRUE
RUN_LEAVE_ONE_SOURCE_OUT <- TRUE
RUN_INHERITED_HARMONY_UMAP <- TRUE
SAVE_HARMONY_CDS <- FALSE
RUN_PDF_PLOTS <- TRUE
SOURCE_ORDER <- c("17702045", "GSE162610", "GSE196928")
STAGE_ORDER <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")
STAGE_DAY <- c(Control = 0, `4 h` = 4/24, `1 dpi` = 1, `3 dpi` = 3, `7 dpi` = 7, `14 dpi` = 14, `28 dpi` = 28,
  `38 dpi` = 38)
HOME <- "Homeostatic Microglia"
FOCUS <- "Clec7a+ lipid-associated Microglia"
ROOT_SCHEMES <- c("control_homeostatic", "all_homeostatic", "earliest_in_partition", "farthest_from_anchor")

ORIGINAL_LIBRARY <- file.path(ROOT, "Rlibs", "Step15_original_versions")
if (dir.exists(ORIGINAL_LIBRARY)) {
  .libPaths(c(ORIGINAL_LIBRARY, .libPaths()))
}
STRICT_VERSION_REQUIREMENTS <- c(harmony = "1.2.4", Seurat = "5.4.0", uwot = "0.2.4")
for (pkg in names(STRICT_VERSION_REQUIREMENTS)) {
  expect <- unname(STRICT_VERSION_REQUIREMENTS[[pkg]])
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Step133 v04 requires ", pkg, " ", expect, "; package is missing. See installation instructions at top of script.",
      call. = FALSE)
  }
  actual <- as.character(utils::packageVersion(pkg))
  if (!identical(actual, expect)) {
    stop(paste0("Step133 v04 STOPPED before analysis: ", pkg, " version ", actual, " is active; original Step15 requires ",
      expect, ".\n", "Do NOT swap max.iter.harmony to max_iter. Install original version into ",
      ORIGINAL_LIBRARY, " and restart R (Ctrl+Shift+F10 in RStudio).\n", "Setup instructions are at the top of the v04 script."),
      call. = FALSE)
  }
}
needed <- c("Seurat", "SeuratObject", "harmony", "Matrix", "monocle3", "SingleCellExperiment", "SummarizedExperiment",
  "igraph", "ggplot2")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
for (p in c(IN18, IN_COHORT, IN_CYCLE)) if (!file.exists(p)) stop("Required DIRECT project file missing: ",
  p)
if (!dir.exists(file.path(ROOT, "script"))) stop("the project script directory/ missing; check ROOT")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "133_Harmony_matched_trajectory_v04.log")
if (file.exists(LOG)) unlink(LOG)
logit <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n")
  cat(s, "\n", file = LOG, append = TRUE)
}
write_out <- function(df, stem, gz = FALSE) {
  if (!is.data.frame(df))
    stop("Output not a data.frame: ", stem)
  f <- file.path(OUT, paste0(stem, "_", V, if (gz)
    ".csv.gz" else ".csv"))
  if (gz) {
    con <- gzfile(f, "wt")
    on.exit(close(con), add = TRUE)
    utils::write.csv(df, con, row.names = FALSE, na = "")
  } else utils::write.csv(df, f, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  logit("WROTE ", basename(f), " | ", nrow(df), " rows")
  invisible(f)
}
rbind_frames <- function(xs) {
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (!length(xs))
    return(data.frame())
  ans <- do.call(rbind, xs)
  rownames(ans) <- NULL
  ans
}
stat_vec <- function(x) {
  x <- x[is.finite(x)]
  c(n = length(x), median = if (length(x)) stats::median(x) else NA_real_, q25 = if (length(x)) as.numeric(stats::quantile(x,
    0.25)) else NA_real_, q75 = if (length(x)) as.numeric(stats::quantile(x, 0.75)) else NA_real_)
}
safe_rho <- function(x, y) {
  good <- is.finite(x) & is.finite(y)
  if (sum(good) < MIN_MATCHED || length(unique(x[good])) < 3L || length(unique(y[good])) < 3L)
    return(NA_real_)
  suppressWarnings(stats::cor(x[good], y[good], method = "spearman"))
}
rank_in_part <- function(x) {
  good <- which(is.finite(x))
  r <- rep(NA_real_, length(x))
  if (length(good) >= 2)
    r[good] <- (rank(x[good], ties.method = "average") - 1)/(length(good) - 1)
  r
}
get_partition_node <- function(cds) {
  graph <- monocle3::principal_graph(cds)[["UMAP"]]
  if (is.null(graph) || igraph::vcount(graph) == 0)
    stop("No Monocle3 principal graph")
  part <- cds@clusters[["UMAP"]]$partitions
  if (is.null(names(part)))
    names(part) <- colnames(cds)
  part <- as.character(part[colnames(cds)])
  raw <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(raw))
    stop("Missing projected principal nodes")
  if (is.null(rownames(raw)))
    stop("Principal-node projection has no cell IDs")
  raw <- as.character(as.matrix(raw[colnames(cds), , drop = FALSE])[, 1])
  vertices <- as.character(igraph::V(graph)$name)
  if (all(raw %in% vertices))
    node <- raw else {
    v <- suppressWarnings(as.integer(raw))
    if (anyNA(v) || any(v < 1L | v > length(vertices)))
      stop("Cannot map projected node IDs")
    node <- vertices[v]
  }
  if (anyNA(part) || anyNA(node))
    stop("Missing partition or principal node")
  list(part = part, node = node, graph = graph)
}
select_roots <- function(nodes, graph, scheme) {
  rootlist <- character()
  rows <- list()
  for (part in sort(unique(nodes$partition))) {
    x <- nodes[nodes$partition == part, , drop = FALSE]
    anchor <- x[x$myeloid_state == HOME & x$timepoint_main == "Control", , drop = FALSE]
    chosen <- NA_character_
    candidates <- 0L
    rule <- "NO_DEFENSIBLE_ANCHOR"
    if (scheme == "control_homeostatic") {

      if (nrow(anchor) >= MIN_ROOT_ANCHOR) {
        chosen <- names(sort(table(anchor$principal_node), decreasing = TRUE))[1]
        candidates <- nrow(anchor)
        rule <- "Control_and_Homeostatic_at_least_10"
      } else {
        days <- STAGE_DAY[x$timepoint_main]
        if (any(is.finite(days))) {
          y <- x[is.finite(days) & days == min(days, na.rm = TRUE), , drop = FALSE]
          chosen <- names(sort(table(y$principal_node), decreasing = TRUE))[1]
          candidates <- nrow(y)
          rule <- paste0("Fallback_earliest_available_", paste(unique(y$timepoint_main), collapse = ";"),
          "_original_Step25_rule")
        }
      }
    } else if (scheme == "all_homeostatic") {
      y <- x[x$myeloid_state == HOME, , drop = FALSE]
      if (nrow(y) >= MIN_ROOT_ANCHOR) {
        chosen <- names(sort(table(y$principal_node), decreasing = TRUE))[1]
        candidates <- nrow(y)
        rule <- "All_stage_Homeostatic_at_least_10"
      }
    } else if (scheme == "earliest_in_partition") {
      days <- STAGE_DAY[x$timepoint_main]
      if (any(is.finite(days))) {
        y <- x[is.finite(days) & days == min(days, na.rm = TRUE), , drop = FALSE]
        chosen <- names(sort(table(y$principal_node), decreasing = TRUE))[1]
        candidates <- nrow(y)
        rule <- paste0("Earliest_available_", paste(unique(y$timepoint_main), collapse = ";"),
          "_NOT_independent_biological_clock")
      }
    } else if (scheme == "farthest_from_anchor") {
      if (nrow(anchor) >= MIN_ROOT_ANCHOR) {
        original <- names(sort(table(anchor$principal_node), decreasing = TRUE))[1]
        dest <- intersect(unique(x$principal_node), igraph::V(graph)$name)
        d <- as.numeric(igraph::distances(graph, v = original, to = dest, weights = NA)[1, ])
        if (any(is.finite(d) & d > 0)) {
          chosen <- dest[which.max(replace(d, !is.finite(d), -Inf))]
          candidates <- nrow(anchor)
          rule <- "Topology_stress_test_NOT_biological"
        }
      }
    } else stop("Unknown root scheme: ", scheme)
    if (!is.na(chosen))
      rootlist <- c(rootlist, chosen)
    rows[[length(rows) + 1L]] <- data.frame(partition = part, root_scheme = scheme, root_node = chosen,
      selection_rule = rule, n_partition = nrow(x), n_control_homeostatic = nrow(anchor), n_candidates = candidates,
      status = if (is.na(chosen))
        "UNROOTED" else "ROOTED", stringsAsFactors = FALSE)
  }
  list(roots = unique(rootlist), audit = rbind_frames(rows))
}

logit("Step133 started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logit("Inputs: ", IN18, " | ", IN_COHORT, " | ", IN_CYCLE)
logit("Output: ", OUT)
for (pkg in c("Seurat", "harmony", "monocle3", "uwot")) {
  if (requireNamespace(pkg, quietly = TRUE))
    logit(pkg, " version: ", as.character(utils::packageVersion(pkg)))
}
logit("R version: ", R.version.string)
logit("Local isolated library: ", if (dir.exists(ORIGINAL_LIBRARY)) ORIGINAL_LIBRARY else "NOT_FOUND (global library used)")
logit("Harmony 1.2.4 verified; strict package version preflight PASSED")

obj <- readRDS(IN18)
if (!inherits(obj, "Seurat") || !"RNA" %in% Seurat::Assays(obj)) stop("Step18 must be a Seurat object with RNA assay")
md <- obj[[]]
req_md <- c("source_dataset", "timepoint_main", "myeloid_state")
if (length(setdiff(req_md, names(md)))) stop("Step18 metadata missing: ", paste(setdiff(req_md, names(md)),
  collapse = ", "))
md$cell_id <- rownames(md)
for (x in req_md) md[[x]] <- as.character(md[[x]])
if (anyDuplicated(md$cell_id)) stop("Duplicate Step18 cell IDs")
inc <- utils::read.csv(gzfile(IN_COHORT), check.names = FALSE, stringsAsFactors = FALSE)
req_inc <- c("cell_id", "source_dataset", "timepoint_main", "myeloid_state", "high_confidence")
if (length(setdiff(req_inc, names(inc)))) stop("Step129 inclusion file missing columns")
if (anyDuplicated(inc$cell_id) || !setequal(inc$cell_id, md$cell_id)) stop("Step129 inclusion IDs do not match current Step18 object; refuse cohort drift")
inc <- inc[match(md$cell_id, inc$cell_id), , drop = FALSE]
for (x in req_md) if (!identical(as.character(inc[[x]]), md[[x]])) stop("Step18 metadata mismatch with Step129 inclusion table: ",
  x)
inc$high_confidence <- as.logical(inc$high_confidence)
if (anyNA(inc$high_confidence)) stop("Missing high_confidence flags in Step129 input")
hc_cells <- md$cell_id[inc$high_confidence]
if (length(hc_cells) < MIN_CELLS || !any(md$cell_id %in% hc_cells & md$myeloid_state == FOCUS)) stop("High-confidence cohort too small or missing Clec7a+ cells")
cc <- utils::read.csv(gzfile(IN_CYCLE), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("cell_id", "source_dataset", "timepoint_main", "myeloid_state", "cycle_exclusion", "cycle_score_available") %in%
  names(cc))) stop("Step129 cell-cycle input missing required columns")
if (anyDuplicated(cc$cell_id) || !setequal(cc$cell_id, hc_cells)) stop("Step129 cycle score file not identical to fixed high-confidence cohort")
cc <- cc[match(hc_cells, cc$cell_id), , drop = FALSE]
hc_md <- md[match(hc_cells, md$cell_id), , drop = FALSE]
if (anyNA(hc_md$cell_id)) stop("HC metadata ID mapping failure")
for (x in req_md) if (!identical(as.character(cc[[x]]), as.character(hc_md[[x]]))) stop("Step129 cycle score vs Step18 metadata mismatch: ",
  x)
cc$cycle_exclusion <- as.logical(cc$cycle_exclusion)
cc$cycle_score_available <- as.logical(cc$cycle_score_available)
if (anyNA(cc$cycle_exclusion) || anyNA(cc$cycle_score_available)) stop("Missing flags in Step129 cell-cycle table")
if (RUN_CYCLE_EXCLUSION && !all(cc$cycle_score_available)) stop("Cycle sensitivity requested but some HC cells lack cycle scores; do not silently use partial thresholds")
hc_md$cycle_exclusion <- cc$cycle_exclusion
hc_md$cycle_score_available <- cc$cycle_score_available
hc_md$cycle_program_score <- if ("cycle_program_score" %in% names(cc)) cc$cycle_program_score else NA_real_
if (any(!hc_md$source_dataset %in% SOURCE_ORDER)) stop("Unexpected source_dataset")
if (any(!hc_md$timepoint_main %in% STAGE_ORDER)) stop("Unexpected timepoint")
if ("treatment" %in% names(hc_md) && any(grepl("PLX5622", hc_md$treatment, ignore.case = TRUE), na.rm = TRUE)) stop("PLX5622 detected in Step18 metadata; review treatment inclusion")
if (!"umap_myeloid_harmony30" %in% names(obj@reductions) && RUN_INHERITED_HARMONY_UMAP) stop("Step18 lacks original umap_myeloid_harmony30 reduction")
old_umap <- NULL
if (RUN_INHERITED_HARMONY_UMAP) {
  old_umap <- Seurat::Embeddings(obj, reduction = "umap_myeloid_harmony30")
  if (!all(hc_cells %in% rownames(old_umap)) || ncol(old_umap) < 2L)
    stop("Original Harmony UMAP does not contain every HC cell")
  old_umap <- old_umap[hc_cells, 1:2, drop = FALSE]
  if (any(!is.finite(old_umap)))
    stop("Non-finite inherited UMAP coordinates")
}
Seurat::DefaultAssay(obj) <- "RNA"
obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = "RNA"), error = function(e) stop("Step18 RNA JoinLayers failed: ",
  conditionMessage(e)))
obj <- subset(obj, cells = hc_cells)
cts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
if (!inherits(cts, "sparseMatrix") || !setequal(colnames(cts), hc_cells)) stop("RNA counts must be sparse and match HC cohort")
cts <- cts[, hc_cells, drop = FALSE]
if (any(Matrix::colSums(cts) <= 0)) stop("Zero-count HC cells")
if (!all(vapply(SOURCE_ORDER, function(s) sum(hc_md$source_dataset == s) > 0, logical(1)))) stop("One of the three study sources absent from fixed cohort")
cohort_audit <- data.frame(item = c("source_object", "cohort_rule", "n_HC", "n_HC_Clec7a", "n_sources",
  "n_cycle_excluded", "n_cycle_retained", "n_RNA_genes", "initial_umap_reduction", "hvg", "pca_computed",
  "pca_used", "seed", "Harmony_batch_covariate", "Harmony_iterations", "UMAP_neighbors", "UMAP_min_dist"),
  value = as.character(c(IN18, IN_COHORT, length(hc_cells), sum(hc_md$myeloid_state == FOCUS), length(unique(hc_md$source_dataset)),
    sum(hc_md$cycle_exclusion), sum(!hc_md$cycle_exclusion), nrow(cts), "umap_myeloid_harmony30",
    N_HVG, N_PCA_COMPUTE, N_DIM, SEED, "source_dataset", HARMONY_ITER, UMAP_NEIGHBORS, UMAP_MIN_DIST)))
write_out(cohort_audit, "01_133_cohort_input_audit")
coverage <- as.data.frame(table(source_dataset = hc_md$source_dataset, timepoint_main = factor(hc_md$timepoint_main,
  levels = STAGE_ORDER), myeloid_state = hc_md$myeloid_state, cycle_exclusion = hc_md$cycle_exclusion),
  stringsAsFactors = FALSE)
coverage <- coverage[coverage$Freq > 0L, , drop = FALSE]
write_out(coverage, "09_133_source_stage_cycle_coverage")
logit("Fixed HC cells=", length(hc_cells), "; HC Clec7a+=", sum(hc_md$myeloid_state == FOCUS), "; cycle-excluded=",
  sum(hc_md$cycle_exclusion))

make_paired_umaps <- function(cells, cohort_label) {
  logit("Shared preprocessing for ", cohort_label, "; n=", length(cells))
  if (length(cells) < MIN_CELLS)
    stop("Insufficient cells for paired embedding")
  so <- subset(obj, cells = cells)
  Seurat::DefaultAssay(so) <- "RNA"
  set.seed(SEED)
  so <- Seurat::NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
  hvg <- Seurat::VariableFeatures(so)
  if (length(hvg) < N_DIM + 2L)
    stop("Too few variable genes for Step133 PCA")
  so <- Seurat::ScaleData(so, features = hvg, verbose = FALSE)
  set.seed(SEED)
  so <- Seurat::RunPCA(so, features = hvg, npcs = N_PCA_COMPUTE, reduction.name = "pca_hc133", reduction.key = "HCPC_",
    seed.use = SEED, verbose = FALSE)
  n_dim_available <- ncol(Seurat::Embeddings(so, "pca_hc133"))
  if (n_dim_available < N_DIM)
    stop("Fewer than 30 PCA components available")
  set.seed(SEED)
  so <- Seurat::RunUMAP(so, reduction = "pca_hc133", dims = seq_len(N_DIM), seed.use = SEED, reduction.name = "umap_hc133_unadjusted",
    reduction.key = "HCU_", verbose = FALSE)
  if (length(unique(so$source_dataset)) < 2L)
    stop("Harmony cannot test source integration with fewer than 2 sources")
  set.seed(SEED)
  so <- harmony::RunHarmony(object = so, group.by.vars = "source_dataset", reduction.use = "pca_hc133",
    dims.use = seq_len(N_DIM), reduction.save = "harmony_hc133", max.iter.harmony = HARMONY_ITER,
    plot_convergence = FALSE, verbose = FALSE)
  if (!"harmony_hc133" %in% names(so@reductions))
    stop("Harmony reduction missing")
  if (ncol(Seurat::Embeddings(so, "harmony_hc133")) < N_DIM)
    stop("Harmony embedding has insufficient dimensions")
  set.seed(SEED)
  so <- Seurat::RunUMAP(so, reduction = "harmony_hc133", dims = seq_len(N_DIM), seed.use = SEED, reduction.name = "umap_hc133_harmony",
    reduction.key = "HCH_", verbose = FALSE)
  extract <- function(red) {
    xy <- Seurat::Embeddings(so, red)
    if (!all(cells %in% rownames(xy)) || ncol(xy) < 2L)
      stop("UMAP missing cells in ", red)
    xy <- xy[cells, 1:2, drop = FALSE]
    if (any(!is.finite(xy)))
      stop("Non-finite UMAP coordinates in ", red)
    xy
  }
  result <- list(raw = extract("umap_hc133_unadjusted"), harmony = extract("umap_hc133_harmony"), hvg_count = length(hvg))
  rm(so)
  invisible(gc(verbose = FALSE))
  result
}

make_base_cds <- function(cells) {
  m <- cts[, cells, drop = FALSE]
  z <- hc_md[match(cells, hc_md$cell_id), , drop = FALSE]
  if (anyNA(z$cell_id))
    stop("Unknown cells in CDS")
  rownames(z) <- z$cell_id
  gene <- data.frame(gene_short_name = rownames(m), row.names = rownames(m))
  cds <- monocle3::new_cell_data_set(m, cell_metadata = z, gene_metadata = gene)
  cds <- monocle3::estimate_size_factors(cds)
  set.seed(SEED)
  cds <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = N_DIM, norm_method = "log")
  set.seed(SEED)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", umap.fast_sgd = FALSE,
    cores = 1L)
  cds
}

all_pt <- list()
all_roots <- list()
status_rows <- list()
partition_diag <- list()
summaries <- list()
root_sensitivity <- list()
focus_position <- list()

plot_branch <- function(cds, branch, main_only = FALSE) {
  if (!RUN_PDF_PLOTS || !main_only)
    return(invisible(NULL))
  for (col in c("source_dataset", "timepoint_main", "myeloid_state")) {
    p <- tryCatch(monocle3::plot_cells(cds, color_cells_by = col, show_trajectory_graph = TRUE, label_cell_groups = FALSE,
      label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.12, graph_label_size = 2) +
      ggplot2::ggtitle(paste(branch, col, sep = " | ")), error = function(e) e)
    if (inherits(p, "error")) {
      logit("Monocle plot fallback for ", branch, " ", col, ": ", conditionMessage(p))
      xy <- SingleCellExperiment::reducedDim(cds, "UMAP")
      d <- as.data.frame(xy[, 1:2, drop = FALSE])
      d$label <- as.character(SummarizedExperiment::colData(cds)[[col]])
      names(d)[1:2] <- c("UMAP1", "UMAP2")
      if (nrow(d) > 9000L)
        d <- d[sample.int(nrow(d), 9000L), , drop = FALSE]
      p <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, color = label)) + ggplot2::geom_point(size = 0.13,
        alpha = 0.7) + ggplot2::theme_classic() + ggplot2::ggtitle(paste(branch, col, "(graph plot unavailable)"))
    }
    print(p)
  }
}

fit_umap_graph <- function(branch, base_cds, xy, main_only = FALSE) {
  logit("Monocle3 fit: ", branch, " n=", nrow(xy))
  cell <- colnames(base_cds)
  if (!setequal(cell, rownames(xy)))
    stop("UMAP vs CDS cell IDs mismatch")
  xy <- xy[cell, 1:2, drop = FALSE]
  cds <- base_cds
  SingleCellExperiment::reducedDim(cds, "UMAP") <- as.matrix(xy)
  set.seed(SEED)
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
  set.seed(SEED)
  cds <- monocle3::learn_graph(cds, use_partition = TRUE, learn_graph_control = list(euclidean_distance_ratio = 0.8))
  nd <- get_partition_node(cds)
  m <- as.data.frame(SummarizedExperiment::colData(cds))
  nodes <- data.frame(cell_id = cell, partition = nd$part, principal_node = nd$node, myeloid_state = as.character(m[cell,
    "myeloid_state"]), timepoint_main = as.character(m[cell, "timepoint_main"]), source_dataset = as.character(m[cell,
    "source_dataset"]), stringsAsFactors = FALSE)

  for (part in unique(nodes$partition)) {
    x <- nodes[nodes$partition == part, , drop = FALSE]
    counts <- table(factor(x$source_dataset, levels = SOURCE_ORDER))
    p <- counts[counts > 0L]/sum(counts)
    entropy <- if (length(p))
      -sum(p * log(p))/log(length(SOURCE_ORDER)) else NA_real_
    partition_diag[[length(partition_diag) + 1L]] <<- data.frame(branch = branch, partition = part,
      n_cells = nrow(x), n_sources = sum(counts > 0L), dominant_source = names(counts)[which.max(counts)],
      dominant_source_fraction = max(counts)/nrow(x), normalized_source_entropy = entropy, n_17702045 = unname(counts[1]),
      n_GSE162610 = unname(counts[2]), n_GSE196928 = unname(counts[3]))
  }
  branch_pts <- list()
  root_tables <- list()
  for (scheme in ROOT_SCHEMES) {
    xroot <- select_roots(nodes, nd$graph, scheme)
    xroot$audit$branch <- branch
    root_tables[[scheme]] <- xroot$audit
    if (!length(xroot$roots))
      pt <- rep(NA_real_, length(cell)) else {
      ordered <- monocle3::order_cells(cds, reduction_method = "UMAP", root_pr_nodes = xroot$roots)
      raw_pt <- monocle3::pseudotime(ordered)
      pt <- as.numeric(raw_pt[cell])
      pt[!is.finite(pt)] <- NA_real_
      rooted <- xroot$audit$partition[!is.na(xroot$audit$root_node)]
      pt[!nodes$partition %in% rooted] <- NA_real_
      if (scheme == "control_homeostatic" && main_only && branch == "HC_Harmony_refit" && SAVE_HARMONY_CDS) {
        saveRDS(ordered, file.path(OUT, "133_HC_Harmony_main_v04.rds"), compress = TRUE)
        logit("Saved Step133 HC Harmony CDS")
      }
    }
    t <- nodes
    t$branch <- branch
    t$root_scheme <- scheme
    t$pseudotime <- pt
    t$rank_within_partition <- NA_real_
    for (part in unique(t$partition)) {
      ix <- which(t$partition == part)
      t$rank_within_partition[ix] <- rank_in_part(t$pseudotime[ix])
    }
    branch_pts[[scheme]] <- t
  }
  root_out <- rbind_frames(root_tables)
  pt_out <- rbind_frames(branch_pts)

  reference <- branch_pts[["control_homeostatic"]]
  for (alt in setdiff(ROOT_SCHEMES, "control_homeostatic")) {
    other <- branch_pts[[alt]]
    if (!identical(reference$cell_id, other$cell_id) || !identical(reference$partition, other$partition))
      stop("Root test mismatched cells")
    for (part in unique(reference$partition)) {
      ix <- which(reference$partition == part)
      valid <- is.finite(reference$pseudotime[ix]) & is.finite(other$pseudotime[ix])
      root_sensitivity[[length(root_sensitivity) + 1L]] <<- data.frame(branch = branch, partition = part,
        baseline = "control_homeostatic", alternative = alt, n_partition = length(ix), n_both_rooted = sum(valid),
        spearman_same_graph = safe_rho(reference$pseudotime[ix], other$pseudotime[ix]), interpretation = "Same graph; within-partition root diagnostic only")
    }
  }
  primary <- branch_pts[["control_homeostatic"]]
  for (part in unique(primary$partition)) {
    x <- primary[primary$partition == part, , drop = FALSE]
    foc <- x[x$myeloid_state == FOCUS, , drop = FALSE]
    focus_stats <- stat_vec(foc$rank_within_partition)
    other <- x[x$myeloid_state != FOCUS, , drop = FALSE]
    other_stats <- stat_vec(other$rank_within_partition)
    focus_position[[length(focus_position) + 1L]] <<- data.frame(branch = branch, partition = part,
      n_partition = nrow(x), n_Clec7a = nrow(foc), n_Clec7a_rooted = unname(focus_stats["n"]),
      n_other_rooted = unname(other_stats["n"]), Clec7a_median_within_partition_rank = unname(focus_stats["median"]),
      other_state_median_within_partition_rank = unname(other_stats["median"]), Clec7a_q25_rank = unname(focus_stats["q25"]),
      Clec7a_q75_rank = unname(focus_stats["q75"]), note = if (unname(focus_stats["n"]) == 0)
        "UNROOTED_OR_NO_CLEC7A" else "Descriptive within-partition position, NOT developmental stage")
    for (variable in c("myeloid_state", "timepoint_main", "source_dataset")) {
      inds <- split(seq_len(nrow(x)), x[[variable]])
      for (value in names(inds)) {
        pts <- x$pseudotime[inds[[value]]]
        stats <- stat_vec(pts)
        summaries[[length(summaries) + 1L]] <<- data.frame(branch = branch, partition = part,
          group_type = variable, group_value = value, n_cells = length(pts), n_finite = unname(stats["n"]),
          median_pseudotime = unname(stats["median"]), q25_pseudotime = unname(stats["q25"]),
          q75_pseudotime = unname(stats["q75"]))
      }
    }
  }
  tryCatch(plot_branch(cds, branch, main_only), error = function(e) logit("Diagnostic plot error (trajectory remains valid) for ",
    branch, ": ", conditionMessage(e)))
  all_pt[[branch]] <<- pt_out
  all_roots[[branch]] <<- root_out
  rooted <- sum(is.finite(primary$pseudotime))
  nfocus <- sum(primary$myeloid_state == FOCUS)
  rooted_focus <- sum(is.finite(primary$pseudotime) & primary$myeloid_state == FOCUS)
  list(n_cells = length(cell), n_partition = length(unique(nodes$partition)), n_rooted_primary = rooted,
    n_focus = nfocus, n_focus_rooted = rooted_focus)
}

run_pair <- function(cells, cohort_name, fatal = FALSE, inherited = FALSE) {
  t0 <- Sys.time()
  attempted <- c(paste0(cohort_name, "_unadjusted"), paste0(cohort_name, "_Harmony_refit"))
  if (inherited)
    attempted <- c(attempted, paste0(cohort_name, "_Harmony_inherited"))
  current <- paste0(cohort_name, "_paired_preprocessing")
  completed <- character()
  tryCatch({
    em <- make_paired_umaps(cells, cohort_name)
    cds <- make_base_cds(cells)
    for (b in attempted) {
      current <- b
      xy <- if (grepl("_Harmony_inherited$", b))
        old_umap[cells, , drop = FALSE] else if (grepl("_unadjusted$", b))
        em$raw else em$harmony
      x <- fit_umap_graph(b, cds, xy, main_only = (cohort_name == "HC"))
      completed <- c(completed, b)
      status_rows[[b]] <<- data.frame(branch = b, n_requested = length(cells), n_rooted_primary = x$n_rooted_primary,
        n_partition = x$n_partition, n_focus = x$n_focus, n_focus_rooted = x$n_focus_rooted,
        status = "COMPLETE", notes = "", elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      logit("DONE ", b, " partitions=", x$n_partition, " rooted=", x$n_rooted_primary, " focus_rooted=",
        x$n_focus_rooted, "/", x$n_focus)
    }
    rm(cds, em)
    invisible(gc(verbose = FALSE))
  }, error = function(e) {
    msg <- paste0("Step133 branch ", current, " failed: ", conditionMessage(e))
    logit(msg)
    status_rows[[current]] <<- data.frame(branch = current, n_requested = length(cells), n_rooted_primary = NA_integer_,
      n_partition = NA_integer_, n_focus = NA_integer_, n_focus_rooted = NA_integer_, status = "FAILED",
      notes = msg, elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    if (!current %in% attempted) {
      for (b in attempted) status_rows[[b]] <<- data.frame(branch = b, n_requested = length(cells),
        n_rooted_primary = NA_integer_, n_partition = NA_integer_, n_focus = NA_integer_, n_focus_rooted = NA_integer_,
        status = "NOT_RUN_PREPROCESSING_FAILED", notes = msg, elapsed_sec = NA_real_)
    }
    for (b in setdiff(attempted, c(completed, current))) status_rows[[b]] <<- data.frame(branch = b,
      n_requested = length(cells), n_rooted_primary = NA_integer_, n_partition = NA_integer_, n_focus = NA_integer_,
      n_focus_rooted = NA_integer_, status = "NOT_RUN_DEPENDENCY_FAILED", notes = msg, elapsed_sec = NA_real_)
    if (fatal) {
      if (RUN_PDF_PLOTS && grDevices::dev.cur() > 1L)
        grDevices::dev.off()
      stop(msg, call. = FALSE)
    }
  })

  write_out(rbind_frames(status_rows), "02_133_branch_status")
  invisible(NULL)
}

if (RUN_PDF_PLOTS) {
  grDevices::pdf(file.path(OUT, "12_133_trajectory_diagnostics_v04.pdf"), width = 11, height = 8.5,
    onefile = TRUE)
}
run_pair(hc_cells, "HC", fatal = TRUE, inherited = RUN_INHERITED_HARMONY_UMAP)
if (RUN_CYCLE_EXCLUSION) {
  retained <- hc_cells[!hc_md$cycle_exclusion]
  run_pair(retained, "HC_cycle_excluded")
}
if (RUN_LEAVE_ONE_SOURCE_OUT) for (src in SOURCE_ORDER) {
  retained <- hc_cells[hc_md$source_dataset != src]
  run_pair(retained, paste0("HC_LOO_", src))
}
if (RUN_PDF_PLOTS && grDevices::dev.cur() > 1L) grDevices::dev.off()

write_out(rbind_frames(partition_diag), "03_133_partition_source_diagnostics")
write_out(rbind_frames(all_roots), "04_133_partition_root_audit")
write_out(rbind_frames(summaries), "05_133_partition_state_stage_source_summary")
write_out(rbind_frames(root_sensitivity), "06_133_root_choice_sensitivity")
write_out(rbind_frames(focus_position), "15_133_Clec7a_partition_position")
all_cell <- rbind_frames(all_pt)
write_out(all_cell, "11_133_cell_pseudotime_all_branches", gz = TRUE)

main_primary <- function(b) {
  if (!b %in% names(all_pt))
    return(NULL)
  tb <- all_pt[[b]]
  tb[tb$root_scheme == "control_homeostatic", , drop = FALSE]
}
pairs <- list(c("HC_unadjusted", "HC_Harmony_refit", "Harmony_effect_fixed_HC_cohort"), c("HC_Harmony_refit",
  "HC_Harmony_inherited", "Inherited_all_myeloid_vs_HC_Harmony"), c("HC_unadjusted", "HC_Harmony_inherited",
  "Inherited_vs_unadjusted"), c("HC_unadjusted", "HC_cycle_excluded_unadjusted", "Cycle_exclusion_raw"),
  c("HC_Harmony_refit", "HC_cycle_excluded_Harmony_refit", "Cycle_exclusion_Harmony"), c("HC_cycle_excluded_unadjusted",
    "HC_cycle_excluded_Harmony_refit", "Harmony_effect_cycle_excluded"))
for (src in SOURCE_ORDER) {
  p <- paste0("HC_LOO_", src)
  pairs <- c(pairs, list(c(paste0(p, "_unadjusted"), paste0(p, "_Harmony_refit"), paste0("Harmony_effect_LOO_",
    src)), c("HC_unadjusted", paste0(p, "_unadjusted"), paste0("Leave_one_source_out_raw_", src)),
    c("HC_Harmony_refit", paste0(p, "_Harmony_refit"), paste0("Leave_one_source_out_Harmony_", src))))
}
part_overlap <- list()
rank_cmp <- list()
for (pr in pairs) {
  aa <- main_primary(pr[1])
  bb <- main_primary(pr[2])
  if (is.null(aa) || is.null(bb))
    next
  combo <- pr[3]
  pa <- unique(aa$partition)
  pb <- unique(bb$partition)
  overlap_now <- list()
  for (a in pa) for (b in pb) {
    ca <- aa$cell_id[aa$partition == a]
    cb <- bb$cell_id[bb$partition == b]
    shared <- sum(ca %in% cb)
    j <- shared/(length(ca) + length(cb) - shared)
    overlap_now[[length(overlap_now) + 1L]] <- data.frame(comparison = combo, branch_A = pr[1], branch_B = pr[2],
      partition_A = a, partition_B = b, n_A = length(ca), n_B = length(cb), n_shared = shared,
      jaccard = j)
  }
  ov <- rbind_frames(overlap_now)
  part_overlap[[length(part_overlap) + 1L]] <- ov
  for (a in pa) {
    z <- ov[ov$partition_A == a, , drop = FALSE]
    z <- z[order(-z$jaccard, -z$n_shared), , drop = FALSE]
    best <- z[1, , drop = FALSE]
    xa <- aa[aa$partition == a, , drop = FALSE]
    xb <- bb[bb$partition == best$partition_B, , drop = FALSE]
    match_b <- match(xa$cell_id, xb$cell_id)
    ind <- which(!is.na(match_b))
    v1 <- xa$rank_within_partition[ind]
    v2 <- xb$rank_within_partition[match_b[ind]]
    good <- is.finite(v1) & is.finite(v2)
    eligible <- is.finite(best$jaccard) && best$jaccard >= JACCARD_CUTOFF && sum(good) >= MIN_MATCHED
    rank_cmp[[length(rank_cmp) + 1L]] <- data.frame(comparison = combo, branch_A = pr[1], branch_B = pr[2],
      partition_A = a, best_partition_B = best$partition_B, jaccard = best$jaccard, n_shared = length(ind),
      n_both_rooted = sum(good), spearman_within_matched_partitions = if (eligible)
        safe_rho(v1, v2) else NA_real_, interpretable = eligible, note = if (eligible)
        "Conditional shared-cell rank comparison; not a common biological clock" else "Unmatched/unrooted/too few cells: no valid rank comparison")
  }
}
write_out(rbind_frames(part_overlap), "07_133_pair_partition_overlap")
write_out(rbind_frames(rank_cmp), "08_133_pair_rank_sensitivity")

old_comp <- data.frame(status = "OPTIONAL_129_PSEUDOTIME_NOT_FOUND", n_shared = NA_integer_, n_same_partition_label = NA_integer_,
  note = "No Step129 original per-cell reference read")
if (file.exists(IN_OLDPT)) {
  old <- utils::read.csv(gzfile(IN_OLDPT), stringsAsFactors = FALSE)
  req <- c("cell_id", "branch", "root_scheme", "partition", "pseudotime")
  if (all(req %in% names(old))) {
    old <- old[old$branch == "HC_main" & old$root_scheme == "control_homeostatic", , drop = FALSE]
    old$partition <- as.character(old$partition)
    now <- main_primary("HC_unadjusted")
    matched <- match(now$cell_id, old$cell_id)
    old_comp <- data.frame(status = "DESCRIPTIVE_ONLY_DIFFERENT_PCA_PIPELINES", n_shared = sum(!is.na(matched)),
      n_same_partition_label = NA_integer_, note = paste("Step129 used Monocle3 PCA to build raw UMAP; Step133 uses",
        "the same Seurat PCA for the paired no-Harmony/Harmony contrast.", "Partition numeric labels are not comparable; use overlap tables."))
  } else old_comp <- data.frame(status = "INCOMPATIBLE_129_REFERENCE_COLUMNS", n_shared = NA_integer_,
    n_same_partition_label = NA_integer_, note = "Step129 primary reference is present but lacks expected columns")
}
write_out(old_comp, "10_133_original129_comparison")

limitations <- c("Step133 v04 | METHODS AND LIMITATIONS", paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "All inputs are direct project files, not archives; no old numbered R script is sourced.",
  "Software: fixed original Step15 harmony 1.2.4, uwot 0.2.4 and Seurat 5.4.0; strict preflight before analysis.",
  "RunHarmony: source_dataset on the first 30 PC, max.iter.harmony=20, as original Step15.", "RunUMAP: Seurat 5.4.0 default uwot backend and original Step15 default neighbors/min.dist, seed=1234.",
  "v04 re-runs full analysis; Step18 and Step129 inclusion/cycle tables are inputs; no prior Step133 result reuse.",
  "v04 output index includes only v04 files; older v01/v02/v03 files in the folder remain untouched.",
  "Cohort: high_confidence=TRUE cell IDs from Step129 inclusion CSV.gz; RNA counts from Step18 Seurat RDS.",
  "Cell-cycle sensitivity: exactly the Step129 cycle_exclusion flags; no threshold re-fitting or tuning to pseudotime.",
  "Main paired embeddings: SAME MG-only cells and SAME Seurat NormalizeData, 3000 HVGs, 50-PC PCA.",
  "No-Harmony: first 30 PCA dimensions -> uwot UMAP. Harmony: same PCA, source_dataset correction (20 iterations) -> uwot UMAP.",
  "Both UMAPs: Seurat 5.4.0 default n.neighbors=30/min.dist=0.3 and seed=1234. Graphs relearned identically.",
  "Inherited Harmony: original all-myeloid umap_myeloid_harmony30 from Step18, restricted to SAME HC cells, then NEW Monocle graph.",
  "The inherited embedding is contextual, not an MG-only Harmony re-fit and not evidence for correct batch correction.",
  "Additional matched paired re-fits: cell-cycle-excluded cohort and each leave-one-source-out cohort.",
  "Primary root rule restored to original Step25 logic: use the node enriched for Control + Homeostatic Microglia when >=10 such cells are present in a partition; otherwise use the node containing the largest number of cells from the earliest available time point in that partition.",
  "Root alternatives all-stage Homeostatic, earliest available stage, and farthest topological root are same-graph sensitivity checks only; the earliest-stage fallback is part of the primary rule only when the preferred Control + Homeostatic anchor is insufficient.",
  "DO NOT interpret pseudotime from differently rooted partitions as the same clock; all summary tables retain branch+partition.",
  "Rank comparison across re-fitted graphs is only attempted on matching cells within best-Jaccard matched partitions.",
  "Partition matches are not guaranteed to be one-to-one; rank correlations conditional on both partitions having defensible roots.",
  "Harmony may remove legitimate biological variation because study source, injury time and cell state are partly confounded.",
  "An apparent increase in source mixing after Harmony is NOT evidence for biological preservation or progression.",
  "Cell-cycle exclusion alters both selected cell composition and learned graph, so effect is not isolatable as pure cell cycle.",
  "No stage-time monotonicity, lineage conversion, mature-cell potency or directional maturation is assumed or established.",
  "All correlations are DESCRIPTIVE at the cell level, not specimen-level independent replicate statistics.",
  "Comparing Step129 and Step133 unadjusted results is contextual only: PCA preprocessing differs between pipelines.",
  "If de novo Harmony branches fail, inspect 02_133_branch_status and log; do not claim negative/positive robustness.")
writeLines(limitations, file.path(OUT, "13_133_methods_and_limitations_v04.txt"), useBytes = TRUE)
sink(file.path(OUT, "14_133_sessionInfo_v04.txt"))
print(utils::sessionInfo())
sink()
write_out(rbind_frames(status_rows), "02_133_branch_status")
files <- list.files(OUT, full.names = TRUE)
files <- files[file.info(files)$isdir == FALSE]

files <- files[grepl(paste0("_", V, "."), basename(files), fixed = TRUE)]
index <- data.frame(filename = basename(files), file_bytes = file.info(files)$size, purpose = "Step133 audit / trajectory / sensitivity result",
  stringsAsFactors = FALSE)
index <- rbind(index, data.frame(filename = "00_133_OUTPUT_INDEX_v04.csv", file_bytes = NA_real_, purpose = "Self-index; size varies after writing"))
write_out(index, "00_133_OUTPUT_INDEX")
logit("STEP133 COMPLETE; ", length(all_pt), " completed trajectory branches")
logit("Completed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
ROOT <- PROJECT_ROOT
STEM <- "139_stage_refined_CellChat_spatial_validation"
VERSION <- "v03"
OUT <- file.path(ROOT, "results", STEM)
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
RDSDIR <- file.path(OUT, "RDS")
CKPT <- file.path(OUT, "checkpoints_v03")
for (d in c(OUT, TAB, FIG, RDSDIR, CKPT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, paste0(STEM, "_", VERSION, ".log"))
cat("Step139 v03 START ", format(Sys.time()), "\n", file = LOG)
log_it <- function(...) {
  m <- paste0(..., collapse = "")
  message(m)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", m, "\n", file = LOG, append = TRUE, sep = "")
}
assert <- function(cond, what) if (!isTRUE(cond)) stop(what, call. = FALSE)
write_tab <- function(x, name) {
  assert(is.data.frame(x), paste("Not data.frame:", name))
  p <- file.path(TAB, paste0(name, "_", VERSION, ".csv"))
  utils::write.csv(x, p, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  log_it("SAVED ", basename(p), " [", nrow(x), " rows]")
  invisible(p)
}
write_gz <- function(x, name) {
  p <- file.path(TAB, paste0(name, "_", VERSION, ".csv.gz"))
  con <- gzfile(p, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "")
  log_it("SAVED ", basename(p), " [", nrow(x), " rows]")
  invisible(p)
}
read_table <- function(p) {
  assert(file.exists(p), paste("Missing input:", p))
  utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
}
need_cols <- function(x, cols, nm) {
  missing <- setdiff(cols, names(x))
  assert(!length(missing), paste(nm, "missing:", paste(missing, collapse = ", ")))
}
relative <- function(x) substring(gsub("\\\\", "/", x), nchar(gsub("\\\\", "/", OUT)) + 2L)
stage_norm <- function(x) {
  y <- gsub("[[:space:]_.-]", "", tolower(as.character(x)))
  y[y %in% c("14d", "14dpi", "14day", "14days")] <- "14dpi"
  y[y %in% c("28d", "28dpi", "28day", "28days")] <- "28dpi"
  y[!y %in% c("14dpi", "28dpi")] <- NA_character_
  y
}
STAGES <- c("14dpi", "28dpi")
TARGET_FROZEN <- "Clec7a+ lipid-associated Microglia"
TARGET <- "Clec7a+ lipid-handling Microglia"
PARENTS <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid")
STATES <- c("Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia",
  "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia", "Proliferating Microglia", "MHC-II+ antigen-presenting Macrophage",
  "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage")
NONMY <- c("Neutrophil", "T/NK cell", "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", "Ependymal", "Endothelial",
  "Pericyte/VSMC", "Fibroblast/Meningeal")
AXES <- data.frame(axis = c("APOE", "PSAP", "SPP1", "APP"), ligand = c("Apoe", "Psap", "Spp1", "App"),
  receptor_genes = c("Trem2;Tyrobp", "Gpr37l1", "Cd44", "Trem2;Tyrobp"), canonical_interaction = c("APOE_TREM2_TYROBP",
    "PSAP_GPR37L1", "SPP1_CD44", "APP_TREM2_TYROBP"))
FOCUS_GENES <- unique(c(AXES$ligand, unlist(strsplit(AXES$receptor_genes, ";", fixed = TRUE))))
NBOOT <- 100L
MIN_CELLS <- 10L
SEED <- 13902L
N_WORKERS <- 1L
FORCE_RESTART <- FALSE
SPATIAL_WEIGHT_THRESHOLDS <- c(0.05, 0.1)
SPATIAL_NEIGHBOR_FACTOR <- 1.55
RUN_GSE281205_VALIDATION <- TRUE
set.seed(SEED)
IN <- list(atlas = file.path(ROOT, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds"),
  myeloid = file.path(ROOT, "results", "19_myeloid_state_robustness_temporal_dynamics", "19_myeloid_state_robustness_temporal_dynamics_3scRNA_v02.rds"),
  step77 = file.path(ROOT, "results", "77_Global_CellChat_Clec7a_centered", "02_77_Global_CellChat_input_provenance_v05.csv"),
  spatial = file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation", "12_64_GSE195783_spatial_region_annotated_v03.rds"),
  rctd = file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping", "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv.gz"),
  gse281205_object = file.path(ROOT, "results", "130_GSE281205_myeloid_reclustering", "RDS", "130_GSE281205_myeloid_LogNorm_PCA_Harmony_clusters_v02.rds"),
  gse281205_manifest = file.path(ROOT, "results", "131_GSE281205_myeloid_subtype_annotation_and_microglia_groups",
    "04_131_per_cell_subtype_and_analysis_manifest_v01.csv"), gse281205_features = file.path(ROOT,
    "results", "124_GSE281205_lineage_review_myeloid_extraction", "12_124_gene_feature_mapping_v01.csv"))
write_tab(data.frame(role = names(IN), path = unlist(IN), exists = file.exists(unlist(IN)), size_bytes = as.numeric(file.info(unlist(IN))$size)),
  "00_139_INPUT_AUDIT")
for (nm in c("atlas", "myeloid", "spatial", "rctd")) assert(file.exists(IN[[nm]]), paste("Required local input missing:",
  IN[[nm]]))

init_step <- function(label, expression) {
  log_it("INIT_START: ", label)
  tryCatch({
    ans <- force(expression)
    log_it("INIT_OK: ", label)
    ans
  }, error = function(e) {
    msg <- paste0("INIT_FAILED [", label, "]: ", conditionMessage(e))
    log_it(msg)
    stack <- vapply(sys.calls(), function(call) paste(deparse(call, width.cutoff = 160), collapse = " "),
      "")
    writeLines(c(format(Sys.time()), msg, "CALL_STACK:", paste(seq_along(stack), stack, sep = " | ")),
      con = file.path(OUT, paste0("139_STARTUP_ERROR_TRACE_", VERSION, ".txt")), useBytes = TRUE)
    stop(msg, call. = FALSE)
  })
}
log_it("R VERSION: ", R.version.string)
req <- c("Matrix", "SeuratObject", "Seurat", "future", "CellChat")
for (pkg in req) {

  if (identical(pkg, "CellChat")) {
    init_step("future sequential backend (before CellChat)", future::plan(future::sequential))
  }
  ok <- init_step(paste0("requireNamespace(", pkg, ")"), requireNamespace(pkg, quietly = TRUE))
  assert(isTRUE(ok), paste("Failed to load package", pkg, "; inspect startup log"))
  log_it("PACKAGE_VERSION ", pkg, " = ", as.character(utils::packageVersion(pkg)))
}

options(future.globals.maxSize = 8 * 1024^3)
DB <- init_step("CellChat mouse database", CellChat::subsetDB(CellChat::CellChatDB.mouse))
DBGENES <- init_step("CellChat database gene extraction", CellChat::extractGene(DB))
assert(length(DBGENES) > 100L, "CellChat database returned too few genes")
log_it("INIT_READY: all dependencies/database loaded; beginning biological input RDS reads")
get_layer <- function(obj, assay, layer) {
  out <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = layer), error = function(e) NULL)
  if (is.null(out))
    out <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, slot = layer), error = function(e) NULL)
  assert(!is.null(out), paste("Missing", assay, layer, "layer"))
  out
}
check_data <- function(obj, assay) {
  if (inherits(obj[[assay]], "Assay5")) {
    layers <- SeuratObject::Layers(obj[[assay]])
    if (sum(grepl("^counts(\\.|$)", layers)) > 1L || sum(grepl("^data(\\.|$)", layers)) > 1L) {
      log_it("Joining split ", assay, " layers in TEMPORARY subset")
      obj <- SeuratObject::JoinLayers(obj, assay = assay)
    }
    if (!"data" %in% SeuratObject::Layers(obj[[assay]])) {
      log_it("RNA normalized layer absent; LogNormalize temporary subset only")
      obj <- Seurat::NormalizeData(obj, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000,
        verbose = FALSE)
    }
  } else {
    datum <- tryCatch(get_layer(obj, assay, "data"), error = function(e) NULL)
    if (is.null(datum) || !nrow(datum) || !ncol(datum))
      obj <- Seurat::NormalizeData(obj, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000,
        verbose = FALSE)
  }
  obj
}
log_it("STEP139 v03 | local input RDS/CSV (no ZIP) | stage-specific refined CellChat")

full <- readRDS(IN$atlas)
my <- readRDS(IN$myeloid)
assert(inherits(full, "Seurat") && inherits(my, "Seurat"), "Step09/19 must be Seurat objects")
assert("RNA" %in% SeuratObject::Assays(full), "Step09 RNA assay missing")
fm <- full[[]]
mm <- my[[]]
need_cols(fm, c("major_celltype", "time", "source_dataset", "sample_id"), "Step09 meta")
need_cols(mm, "myeloid_state", "Step19 meta")
assert(identical(rownames(fm), colnames(full)), "Full metadata/cell order differs")
assert(identical(rownames(mm), colnames(my)), "Myeloid metadata/cell order differs")
assert(!anyDuplicated(rownames(fm)) && !anyDuplicated(rownames(mm)), "Duplicate cell IDs")
assert(all(rownames(mm) %in% rownames(fm)), "Step19 IDs not exactly present in Step09; STOP")
assert(all(unique(as.character(mm$myeloid_state)) %in% STATES), paste("Unexpected Step19 labels:", paste(setdiff(unique(as.character(mm$myeloid_state)),
  STATES), collapse = ";")))
assert(TARGET_FROZEN %in% as.character(mm$myeloid_state), "Clec7a+ state missing")
map <- setNames(as.character(mm$myeloid_state), rownames(mm))
fm$cell_id_139 <- rownames(fm)
fm$stage_139 <- stage_norm(fm$time)
fm$major_139 <- as.character(fm$major_celltype)
fm$source_139 <- as.character(fm$source_dataset)
fm$sample_139 <- as.character(fm$sample_id)
fm$is_refined_139 <- fm$cell_id_139 %in% names(map)
fm$removed_139 <- fm$major_139 %in% PARENTS & !fm$is_refined_139
fm$group_139 <- fm$major_139
idx <- which(fm$is_refined_139)
fm$group_139[idx] <- map[fm$cell_id_139[idx]]
fm$group_139[fm$group_139 == TARGET_FROZEN] <- TARGET
selected <- !is.na(fm$stage_139) & !fm$removed_139 & !is.na(fm$group_139) & nzchar(trimws(fm$group_139))
md <- fm[selected, , drop = FALSE]
assert(all(STAGES %in% unique(md$stage_139)), "14dpi or 28dpi missing")
assert(all(md$group_139 %in% c(setdiff(STATES, TARGET_FROZEN), TARGET, NONMY)), paste("Unrecognized stage communication group:",
  paste(setdiff(unique(md$group_139), c(setdiff(STATES, TARGET_FROZEN), TARGET, NONMY)), collapse = ";")))
assert(all(!is.na(md$source_139) & nzchar(md$source_139)), "Missing source_dataset")
assert(all(!is.na(md$sample_139) & nzchar(md$sample_139)), "Missing sample_id")
LEVELS <- c(TARGET, setdiff(STATES, TARGET_FROZEN), NONMY)
write_tab(data.frame(source = rep(c("Step09", "Step19"), c(ncol(fm), ncol(mm))), metadata_field = c(colnames(fm),
  colnames(mm))), "01_139_METADATA_FIELDS")
manifest <- as.data.frame(table(stage = factor(md$stage_139, levels = STAGES), source_dataset = md$source_139,
  sample_id = md$sample_139, refined_group = md$group_139), stringsAsFactors = FALSE)
manifest <- manifest[manifest$Freq > 0, , drop = FALSE]
names(manifest)[5] <- "n_cells"
write_tab(manifest, "02_139_STAGE_SOURCE_LIBRARY_REFINED_GROUP_COUNTS")
write_gz(data.frame(cell = rownames(md), stage = md$stage_139, source_dataset = md$source_139, sample_id = md$sample_139,
  parent_lineage = md$major_139, refined_group = md$group_139), "03_139_STAGE_CELL_MANIFEST")
refcheck <- data.frame(stage = STAGES, n_step139 = vapply(STAGES, function(s) sum(md$stage_139 == s),
  0L), n_step77 = NA_integer_, stage_matches_reference = NA, stringsAsFactors = FALSE)
if (file.exists(IN$step77)) {
  ref <- read_table(IN$step77)
  need_cols(ref, c("time", "n_cells"), "Step77 provenance")
  ref$stage <- stage_norm(ref$time)
  refcheck$n_step77 <- vapply(STAGES, function(s) as.integer(sum(ref$n_cells[ref$stage == s], na.rm = TRUE)),
    0L)
  refcheck$stage_matches_reference <- refcheck$n_step139 == refcheck$n_step77
  write_tab(refcheck, "04_139_COHORT_STEP77_RECONCILIATION")
  assert(all(refcheck$stage_matches_reference), "Cohort does not match Step77! Inspect 04 reconciliation before running CellChat")
} else {
  log_it("WARNING: Step77 provenance absent, count reconciliation unavailable")
  write_tab(refcheck, "04_139_COHORT_STEP77_RECONCILIATION")
}
unit_names <- intersect(c("animal_id", "mouse_id", "specimen_id", "donor_id", "biological_sample_id",
  "animal", "tissue_pool_id", "pool_id"), names(fm))
feas <- do.call(rbind, lapply(STAGES, function(s) {
  z <- md[md$stage_139 == s, , drop = FALSE]
  data.frame(stage = s, n_cells = nrow(z), source_datasets = paste(sort(unique(z$source_139)), collapse = ";"),
    n_sources = length(unique(z$source_139)), sample_or_library_ids = paste(sort(unique(z$sample_139)),
      collapse = ";"), n_library_ids = length(unique(z$sample_139)), candidate_biological_unit_fields = paste(unit_names,
      collapse = ";"), n_verified_independent_animals = NA_integer_, same_stage_cross_dataset_replication = if (length(unique(z$source_139)) >
      1)
      "SOURCE_STRATIFICATION_POSSIBLE_NEEDS_VALIDATION" else "NOT_IDENTIFIABLE_ONE_SOURCE", specimen_level_CellChat_replication = "NOT_TESTED_WITHOUT_VERIFIED_ANIMAL_MAP",
    stringsAsFactors = FALSE)
}))
write_tab(feas, "05_139_SOURCE_SAMPLE_REPLICATION_FEASIBILITY")
rm(my, mm, fm, map, idx, manifest)
gc(verbose = FALSE)

checkpoint_file <- function(s, step) file.path(CKPT, paste0("139_", s, "_", step, "_", VERSION, ".rds"))
final_file <- function(s) file.path(RDSDIR, paste0("139_", s, "_refined_CellChat_", VERSION, ".rds"))
if (FORCE_RESTART) {
  f <- c(list.files(CKPT, full.names = TRUE, pattern = "_v03\\.rds$"), vapply(STAGES, final_file, ""))
  unlink(f[file.exists(f)], force = TRUE)
  log_it("Deleted ONLY Step139 v03 checkpoints and final CellChat objects")
}
safe_save <- function(o, path) {
  tmp <- paste0(path, ".writing")
  saveRDS(o, tmp, compress = FALSE)
  if (file.exists(path))
    unlink(path)
  assert(file.rename(tmp, path), paste("Could not finalize:", path))
}
load_cached <- function(s, expected_cells) {
  paths <- c(final_file(s), checkpoint_file(s, "04_network"), checkpoint_file(s, "03_prob"), checkpoint_file(s,
    "02_over"), checkpoint_file(s, "01_subset"))
  levels <- c(4L, 4L, 3L, 2L, 1L)
  for (i in seq_along(paths)) {
    if (!file.exists(paths[i]))
      next
    o <- tryCatch(readRDS(paths[i]), error = function(e) NULL)
    if (is.null(o) || !inherits(o, "CellChat"))
      next
    same <- tryCatch({
      saved_cells <- rownames(o@meta)
      saved_ids <- as.character(o@idents)
      assert(length(saved_cells) == length(saved_ids), "Checkpoint group length mismatch")
      saved_groups_by_cell <- setNames(saved_ids, saved_cells)
      setequal(saved_cells, expected_cells) && identical(unname(saved_groups_by_cell[expected_cells]),
        unname(as.character(md[expected_cells, "group_139"])))
    }, error = function(e) FALSE)
    if (!isTRUE(same))
      next
    log_it("RESUME ", s, " from ", basename(paths[i]), " (level ", levels[i], ")")
    return(list(object = o, level = levels[i]))
  }
  list(object = NULL, level = 0L)
}
empty_lr <- function() data.frame(source = character(), target = character(), ligand = character(), receptor = character(),
  interaction_name = character(), pathway_name = character(), prob = numeric(), pval = numeric())
empty_pw <- function() data.frame(source = character(), target = character(), pathway_name = character(),
  prob = numeric(), pval = numeric())
fetch_chat_table <- function(ch, kind) {
  z <- tryCatch(as.data.frame(CellChat::subsetCommunication(ch, slot.name = kind)), error = function(e) {
    log_it("WARNING: subsetCommunication ", kind, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(z) || !nrow(z))
    z <- if (kind == "net")
      empty_lr() else empty_pw()
  need_cols(z, if (kind == "net")
    c("source", "target", "ligand", "receptor", "prob", "pval") else c("source", "target", "pathway_name", "prob", "pval"), "CellChat table")
  z$source <- as.character(z$source)
  z$target <- as.character(z$target)
  z
}

expr_summary <- function(mat, zmd, stage) {
  cell_ids <- colnames(mat)
  assert(identical(cell_ids, rownames(zmd)), "RNA matrix/meta cell mismatch")
  keys <- paste(zmd$source_139, zmd$sample_139, zmd$group_139, sep = "@@")
  ind <- split(seq_along(keys), keys)
  rows <- lapply(ind, function(k) {
    x <- zmd[k[1], , drop = FALSE]
    do.call(rbind, lapply(FOCUS_GENES, function(g) {
      ok <- g %in% rownames(mat)
      v <- if (ok)
        as.numeric(mat[g, k, drop = TRUE]) else rep(NA_real_, length(k))
      data.frame(stage = stage, source_dataset = x$source_139, sample_id = x$sample_139, group = x$group_139,
        gene = g, feature_present = ok, n_cells = length(k), n_detected = if (ok)
          sum(v > 0) else NA_integer_, fraction_detected = if (ok)
          mean(v > 0) else NA_real_, mean_normalized_logexpr = if (ok)
          mean(v) else NA_real_)
    }))
  })
  do.call(rbind, rows)
}
complex_summary <- function(mat, zmd, stage) {
  ind <- split(seq_len(nrow(zmd)), paste(zmd$sample_139, zmd$group_139, sep = "@@"))
  do.call(rbind, lapply(ind, function(k) {
    do.call(rbind, lapply(seq_len(nrow(AXES)), function(j) {
      genes <- strsplit(AXES$receptor_genes[j], ";", fixed = TRUE)[[1]]
      ok <- all(genes %in% rownames(mat))
      if (ok) {
        det <- Matrix::colSums(mat[genes, k, drop = FALSE] > 0) == length(genes)
        fraction <- mean(det)
      } else fraction <- NA_real_
      data.frame(stage = stage, source_dataset = zmd$source_139[k[1]], sample_id = zmd$sample_139[k[1]],
        group = zmd$group_139[k[1]], axis = AXES$axis[j], receptor_genes = AXES$receptor_genes[j],
        all_genes_available = ok, n_cells = length(k), receptor_complex_codetection_fraction = fraction)
    }))
  }))
}
axis_of <- function(ligand) {
  match(tolower(as.character(ligand)), tolower(AXES$ligand))
}
all_results <- list()
all_expr <- list()
all_complex <- list()
for (s in STAGES) {
  zmd <- md[md$stage_139 == s, , drop = FALSE]
  counts <- table(factor(zmd$group_139, levels = LEVELS))
  groups <- names(counts)[counts >= MIN_CELLS]
  assert(TARGET %in% groups, paste(s, "Clec7a target <", MIN_CELLS, " cells"))
  write_tab(data.frame(stage = s, group = names(counts), n_cells = as.integer(counts), included_CellChat = as.integer(counts) >=
    MIN_CELLS), paste0("06_139_", s, "_GROUP_SIZE_GATE"))
  run_md <- zmd[zmd$group_139 %in% groups, , drop = FALSE]
  cp <- load_cached(s, rownames(run_md))

  log_it(s, " loading temporary stage RNA to compute expression support")
  stage_obj <- base::subset(full, cells = rownames(run_md))
  Seurat::DefaultAssay(stage_obj) <- "RNA"
  stage_obj <- check_data(stage_obj, "RNA")
  dat_all <- get_layer(stage_obj, "RNA", "data")
  genes <- intersect(rownames(dat_all), union(DBGENES, FOCUS_GENES))
  assert(length(intersect(genes, DBGENES)) > 100, paste(s, "fewer than 100 CellChatDB genes present"))
  dat <- dat_all[genes, , drop = FALSE]
  order_idx <- match(colnames(dat), rownames(run_md))
  assert(!anyNA(order_idx), paste(s, "unmatched Seurat columns"))
  run_md <- run_md[order_idx, , drop = FALSE]
  assert(identical(colnames(dat), rownames(run_md)), paste(s, "column ordering error"))
  all_expr[[s]] <- expr_summary(dat, run_md, s)
  all_complex[[s]] <- complex_summary(dat, run_md, s)
  rm(stage_obj, dat_all)
  gc(verbose = FALSE)
  if (cp$level == 0L) {
    meta <- data.frame(communication_group = factor(run_md$group_139, levels = groups), row.names = rownames(run_md))
    assert(identical(colnames(dat), rownames(meta)), "CellChat input order mismatch")
    log_it(s, ": create refined CellChat (", nrow(meta), " cells; ", length(groups), " groups)")
    chat <- CellChat::createCellChat(object = dat[intersect(rownames(dat), DBGENES), , drop = FALSE],
      meta = meta, group.by = "communication_group")
    chat@DB <- DB
  } else chat <- cp$object
  rm(dat, run_md)
  gc(verbose = FALSE)
  if (cp$level < 1L) {
    chat <- CellChat::subsetData(chat)
    safe_save(chat, checkpoint_file(s, "01_subset"))
  }
  if (cp$level < 2L) {
    log_it(s, ": overexpressed genes/interactions")
    chat <- CellChat::identifyOverExpressedGenes(chat)
    chat <- CellChat::identifyOverExpressedInteractions(chat)
    safe_save(chat, checkpoint_file(s, "02_over"))
  }
  if (cp$level < 3L) {
    log_it(s, ": compute communication probability (100 bootstrap)")
    this_seed <- SEED + match(s, STAGES)
    set.seed(this_seed)
    log_it("RUN_COMPUTE_PROB_START: ", s)
    chat <- CellChat::computeCommunProb(chat, type = "triMean", raw.use = TRUE, population.size = FALSE,
      nboot = NBOOT, seed.use = this_seed)
    log_it("RUN_COMPUTE_PROB_OK: ", s)
    chat <- CellChat::filterCommunication(chat, min.cells = MIN_CELLS)
    safe_save(chat, checkpoint_file(s, "03_prob"))
  }
  if (cp$level < 4L) {
    chat <- CellChat::computeCommunProbPathway(chat)
    chat <- CellChat::aggregateNet(chat)
    chat <- CellChat::netAnalysis_computeCentrality(chat, slot.name = "netP")
    safe_save(chat, checkpoint_file(s, "04_network"))
    safe_save(chat, final_file(s))
  } else if (!file.exists(final_file(s)))
    safe_save(chat, final_file(s))
  lr <- fetch_chat_table(chat, "net")
  pw <- fetch_chat_table(chat, "netP")
  lr$stage <- s
  pw$stage <- s
  write_tab(lr, paste0("07_139_", s, "_ALL_SIGNIFICANT_LR"))
  write_tab(pw, paste0("08_139_", s, "_ALL_SIGNIFICANT_PATHWAYS"))
  focal <- lr[lr$source == TARGET | lr$target == TARGET, , drop = FALSE]
  write_tab(focal, paste0("09_139_", s, "_CLEC7A_CENTERED_REFINED_LR"))
  x <- focal[!is.na(axis_of(focal$ligand)), , drop = FALSE]
  x$axis <- AXES$axis[axis_of(x$ligand)]
  x$is_predeclared_canonical_pair <- toupper(as.character(x$interaction_name)) %in% AXES$canonical_interaction
  write_tab(x, paste0("10_139_", s, "_APOE_PSAP_SPP1_APP_LINKS"))

  if (!is.null(chat@net$weight) && length(chat@net$weight)) {
    W <- chat@net$weight
    strengths <- data.frame(stage = s, group = rownames(W), inferred_outgoing_sum = as.numeric(rowSums(W)),
      inferred_incoming_sum = as.numeric(colSums(W)), stringsAsFactors = FALSE)
    write_tab(strengths, paste0("11_139_", s, "_NETWORK_STRENGTH"))
  }

  pdfp <- file.path(FIG, paste0("139_", s, "_refined_communications_", VERSION, ".pdf"))
  grDevices::pdf(pdfp, width = 12, height = 9)
  par(mar = c(4, 18, 4, 2))
  o <- order(as.integer(counts))
  barplot(as.integer(counts)[o], names.arg = names(counts)[o], horiz = TRUE, las = 1, cex.names = 0.68,
    col = "grey55", main = paste(s, "Refined groups"), xlab = "Number of retained cells")
  if (nrow(x)) {
    x$label <- paste0(x$axis, " | ", x$source, " -> ", x$target)
    xx <- x[order(-x$prob), , drop = FALSE]
    xx <- head(xx, 25L)
    par(mar = c(5, 20, 4, 2))
    barplot(rev(xx$prob), names.arg = rev(xx$label), horiz = TRUE, las = 1, cex.names = 0.55, col = "grey35",
      main = paste(s, "Significant Clec7a-centered four-axis links"), xlab = "CellChat communication probability")
  } else {
    plot.new()
    text(0.5, 0.5, paste(s, "no significant four-axis links"))
  }
  grDevices::dev.off()
  all_results[[s]] <- list(lr = lr, pw = pw, focus = x, groups = groups, source = unique(zmd$source_139),
    samples = unique(zmd$sample_139))
  rm(chat, cp, lr, pw, focal, x)
  gc(verbose = FALSE)
  log_it("STAGE FINISHED: ", s)
}
rm(full)
gc(verbose = FALSE)
write_tab(do.call(rbind, all_expr), "12_139_PER_SOURCE_LIBRARY_GROUP_GENE_EXPRESSION")
write_tab(do.call(rbind, all_complex), "13_139_PER_SOURCE_LIBRARY_GROUP_RECEPTOR_CODETECTION")

axis_summary <- do.call(rbind, lapply(STAGES, function(s) {
  z <- all_results[[s]]$focus
  do.call(rbind, lapply(seq_len(nrow(AXES)), function(j) {
    x <- z[z$axis == AXES$axis[j], , drop = FALSE]
    cn <- if (nrow(x))
      x[x$is_predeclared_canonical_pair, , drop = FALSE] else x
    data.frame(stage = s, source_dataset = paste(all_results[[s]]$source, collapse = ";"), axis = AXES$axis[j],
      ligand = AXES$ligand[j], receptor_genes = AXES$receptor_genes[j], n_Clec7a_centered_significant_links_all_receptors = nrow(x),
      n_predeclared_canonical_significant_links = nrow(cn), n_canonical_outgoing = if (nrow(cn))
        sum(cn$source == TARGET) else 0L, n_canonical_incoming = if (nrow(cn))
        sum(cn$target == TARGET) else 0L, status = if (nrow(cn))
        "SIGNIFICANT_CANONICAL_CELLCHAT_LINK" else "CANONICAL_LINK_NOT_REPORTED_SIGNIFICANT", stringsAsFactors = FALSE)
  }))
}))
write_tab(axis_summary, "14_139_FOUR_AXIS_STAGE_STATUS")
lrkey <- function(z) if (!nrow(z)) character() else paste(z$source, z$target, z$interaction_name, sep = "||")
a <- all_results[["14dpi"]]$focus
b <- all_results[["28dpi"]]$focus
unionkey <- union(lrkey(a), lrkey(b))
if (length(unionkey)) {
  j14 <- match(unionkey, lrkey(a))
  j28 <- match(unionkey, lrkey(b))
  grab <- function(z, idx, col) {
    if (!col %in% names(z) || !length(idx))
      return(rep(NA, length(idx)))
    out <- z[[col]][idx]
    out
  }
  ab <- strsplit(unionkey, "||", fixed = TRUE)
  cross <- data.frame(source = vapply(ab, function(z) z[1L], ""), target = vapply(ab, function(z) z[2L],
    ""), interaction_name = vapply(ab, function(z) z[3L], ""), pathway_14dpi = grab(a, j14, "pathway_name"),
    pathway_28dpi = grab(b, j28, "pathway_name"), prob_14dpi = grab(a, j14, "prob"), prob_28dpi = grab(b,
      j28, "prob"), pval_14dpi = grab(a, j14, "pval"), pval_28dpi = grab(b, j28, "pval"), significant_14dpi = !is.na(j14),
    significant_28dpi = !is.na(j28), stringsAsFactors = FALSE)
  cross$group_available_14dpi <- cross$source %in% all_results[["14dpi"]]$groups & cross$target %in%
    all_results[["14dpi"]]$groups
  cross$group_available_28dpi <- cross$source %in% all_results[["28dpi"]]$groups & cross$target %in%
    all_results[["28dpi"]]$groups
  cross$interpretation <- "Descriptive, confounded between stages and sources; not same-stage biological replication. Missing significant row does not imply absent expression."
} else cross <- data.frame(source = character(), target = character(), interaction_name = character(),
  significant_14dpi = logical(), significant_28dpi = logical(), interpretation = character())
write_tab(cross, "15_139_FOUR_AXIS_CROSS_STAGE_LINK_DESCRIPTIVE_OVERLAP")

ex <- do.call(rbind, all_expr)
co <- do.call(rbind, all_complex)
canonical_links <- do.call(rbind, lapply(STAGES, function(s) {
  z <- all_results[[s]]$focus
  z[z$is_predeclared_canonical_pair, , drop = FALSE]
}))
if (nrow(canonical_links)) {
  evidence <- do.call(rbind, lapply(seq_len(nrow(canonical_links)), function(i) {
    q <- canonical_links[i, , drop = FALSE]
    j <- match(q$axis, AXES$axis)
    libs <- unique(ex[ex$stage == q$stage, c("source_dataset", "sample_id"), drop = FALSE])
    do.call(rbind, lapply(seq_len(nrow(libs)), function(k) {
      src <- ex[ex$stage == q$stage & ex$source_dataset == libs$source_dataset[k] & ex$sample_id ==
        libs$sample_id[k] & ex$group == q$source & ex$gene == AXES$ligand[j], , drop = FALSE]
      tar <- co[co$stage == q$stage & co$source_dataset == libs$source_dataset[k] & co$sample_id ==
        libs$sample_id[k] & co$group == q$target & co$axis == q$axis, , drop = FALSE]
      data.frame(stage = q$stage, source_dataset = libs$source_dataset[k], sample_or_library = libs$sample_id[k],
        axis = q$axis, sender = q$source, receiver = q$target, interaction_name = q$interaction_name,
        n_sender_cells = if (nrow(src))
          src$n_cells[1] else NA_integer_, ligand_fraction_detected_sender = if (nrow(src))
          src$fraction_detected[1] else NA_real_, ligand_mean_logexpr_sender = if (nrow(src))
          src$mean_normalized_logexpr[1] else NA_real_, n_receiver_cells = if (nrow(tar))
          tar$n_cells[1] else NA_integer_, receptor_ALL_subunits_codetected_receiver = if (nrow(tar))
          tar$receptor_complex_codetection_fraction[1] else NA_real_, same_stage_expression_compatible = if (nrow(src) && nrow(tar))
          isTRUE(src$fraction_detected[1] > 0 && tar$receptor_complex_codetection_fraction[1] >
          0) else NA, note = "Cell-level same-stage RNA supports expression only; CellChat is computational; not independently replicated animals.")
    }))
  }))
} else evidence <- data.frame(stage = character(), axis = character(), sender = character(), receiver = character(),
  same_stage_expression_compatible = logical(), note = character())
write_tab(evidence, "16_139_CANONICAL_LR_SAME_STAGE_SENDER_RECEIVER_EXPRESSION")

log_it("SPATIAL: reading Step65 v06 section-aware RCTD table")
spot <- read_table(IN$rctd)
need_cols(spot, c("GSM", "sample_label", "timepoint", "array_x", "array_y", "section_uid_internal", "section_eligible",
  TARGET), "Step65 RCTD")
spot$stage_139 <- stage_norm(spot$timepoint)
spot <- spot[spot$stage_139 %in% STAGES, , drop = FALSE]
assert(nrow(spot) > 0, "No eligible 14/28 spatial spots in Step65")

valid_section <- !is.na(spot$section_uid_internal) & nzchar(trimws(spot$section_uid_internal))
valid_eligible <- !is.na(spot$section_eligible) & tolower(as.character(spot$section_eligible)) %in% c("true",
  "t", "1")
spot <- spot[valid_section & valid_eligible, , drop = FALSE]
assert(all(STAGES %in% unique(spot$stage_139)), "Spatial 14/28 dpi missing after section QC")
spot$key_139 <- paste(spot$GSM, spot$array_x, spot$array_y, sep = "||")
assert(!anyDuplicated(spot$key_139), "Step65 spatial duplicate coordinate keys")
log_it("SPATIAL: reading Step64 Seurat Spatial assay (only focal gene rows will be used)")
sp_obj <- readRDS(IN$spatial)
assert(inherits(sp_obj, "Seurat") && "Spatial" %in% SeuratObject::Assays(sp_obj), "Step64 Spatial Seurat assay missing")
sp_meta <- sp_obj[[]]
need_cols(sp_meta, c("GSM", "array_x", "array_y", "sample_label"), "Step64 spatial metadata")
sp_data <- get_layer(sp_obj, "Spatial", "data")
sp_pos <- match(colnames(sp_data), rownames(sp_meta))
assert(!anyNA(sp_pos), "Step64 Spatial matrix columns cannot be matched to metadata")
sp_meta <- sp_meta[sp_pos, , drop = FALSE]
sp_meta$key_139 <- paste(sp_meta$GSM, sp_meta$array_x, sp_meta$array_y, sep = "||")
assert(!anyDuplicated(sp_meta$key_139), "Step64 spatial duplicate coordinate keys")
jsp <- match(spot$key_139, sp_meta$key_139)
assert(!anyNA(jsp), paste("Step65 coords not found in Step64 spatial assay:", sum(is.na(jsp))))
assert(all(as.character(spot$sample_label) == as.character(sp_meta$sample_label[jsp])), "Step64/65 capture labels disagree after coordinate join")

sg <- intersect(FOCUS_GENES, rownames(sp_data))
sp_expr <- as.matrix(sp_data[sg, jsp, drop = FALSE])
colnames(sp_expr) <- spot$key_139
rm(sp_obj, sp_meta, sp_data)
gc(verbose = FALSE)
missing_spatial <- setdiff(FOCUS_GENES, sg)
write_tab(data.frame(gene = FOCUS_GENES, present_Step64_Spatial = FOCUS_GENES %in% sg), "17_139_SPATIAL_FOCAL_GENE_FEATURE_AVAILABILITY")
if (length(missing_spatial)) log_it("WARNING spatial missing genes: ", paste(missing_spatial, collapse = ";"))

all_weight_cols <- setdiff(names(spot), c("barcode", "GSM", "sample_label", "timepoint", "spatial_region_v03",
  "array_x", "array_y", "is_major_component", "section_position_internal", "section_uid_internal",
  "section_status", "section_eligible", "stage_139", "key_139"))
map_qc <- data.frame(stage_group = LEVELS, exists_in_Step65_RCTD = LEVELS %in% all_weight_cols, stringsAsFactors = FALSE)
write_tab(map_qc, "18_139_CELLCHAT_GROUP_TO_SPATIAL_RCTD_MAPPING")

links_spatial <- canonical_links[, intersect(c("stage", "axis", "source", "target", "ligand", "receptor",
  "interaction_name", "prob", "pval"), names(canonical_links)), drop = FALSE]
write_tab(links_spatial, "19_139_CANONICAL_LR_INPUT_FOR_SPATIAL_AUDIT")

spatial_gene_detection <- do.call(rbind, lapply(split(seq_len(nrow(spot)), spot$sample_label), function(k) {
  do.call(rbind, lapply(FOCUS_GENES, function(g) {
    v <- if (g %in% rownames(sp_expr))
      sp_expr[g, k] else rep(NA_real_, length(k))
    data.frame(stage = spot$stage_139[k[1]], sample_label = spot$sample_label[k[1]], gene = g, n_spots = length(k),
      present = g %in% rownames(sp_expr), fraction_spots_detected = if (g %in% rownames(sp_expr))
        mean(v > 0) else NA_real_, mean_Spatial_logexpr = if (g %in% rownames(sp_expr))
        mean(v) else NA_real_)
  }))
}))
write_tab(spatial_gene_detection, "20_139_SAME_STAGE_CAPTURE_AREA_SPATIAL_GENE_DETECTION")

build_neighborhood <- function(pos) {
  n <- nrow(pos)
  neighbors <- matrix(FALSE, nrow = n, ncol = n)
  sec <- split(seq_len(n), as.character(pos$section_uid_internal))
  nr <- list()
  for (u in names(sec)) {
    k <- sec[[u]]
    if (length(k) < 2L)
      next
    xy <- cbind(as.numeric(pos$array_x[k]), as.numeric(pos$array_y[k]))
    d <- as.matrix(stats::dist(xy))
    diag(d) <- Inf
    nearest <- apply(d, 1, min, na.rm = TRUE)
    radius <- SPATIAL_NEIGHBOR_FACTOR * stats::median(nearest[is.finite(nearest)])
    if (!is.finite(radius) || radius <= 0)
      next
    adj <- (d <= radius)
    diag(adj) <- FALSE
    neighbors[k, k] <- adj
    nr[[length(nr) + 1L]] <- data.frame(section_uid = u, sample_label = pos$sample_label[k[1]], stage = pos$stage_139[k[1]],
      n_spots = length(k), median_nn_distance = stats::median(nearest), radius = radius, n_undirected_edges = sum(adj)/2)
  }
  list(matrix = neighbors, qc = if (length(nr)) do.call(rbind, nr) else data.frame(section_uid = character(),
    sample_label = character(), stage = character(), n_spots = integer(), median_nn_distance = numeric(),
    radius = numeric(), n_undirected_edges = numeric()))
}
areas <- split(seq_len(nrow(spot)), as.character(spot$sample_label))
area_nets <- lapply(areas, function(k) build_neighborhood(spot[k, , drop = FALSE]))
write_tab(do.call(rbind, lapply(area_nets, function(z) z[["qc"]])), "21_139_SPATIAL_WITHIN_SECTION_NEIGHBOR_GRAPH_QC")
compat <- list()
if (nrow(links_spatial)) {
  for (u in seq_len(nrow(links_spatial))) {
    q <- links_spatial[u, , drop = FALSE]
    ax <- AXES[match(q$axis, AXES$axis), , drop = FALSE]
    receptor <- strsplit(ax$receptor_genes, ";", fixed = TRUE)[[1]]
    for (nm in names(areas)) {
      k <- areas[[nm]]
      if (spot$stage_139[k[1]] != q$stage)
        next
      area <- spot[k, , drop = FALSE]
      avail <- all(c(q$source, q$target) %in% names(area)) && all(c(ax$ligand, receptor) %in% rownames(sp_expr))
      for (cut in SPATIAL_WEIGHT_THRESHOLDS) {
        if (!avail) {
          row <- data.frame(stage = q$stage, sample_label = nm, axis = q$axis, sender = q$source,
          receiver = q$target, interaction_name = q$interaction_name, source_RCTD_weight_cutoff = cut,
          n_spots = length(k), n_source_weight_and_ligand_positive = NA_integer_, n_target_weight_and_all_receptors_positive = NA_integer_,
          n_same_spot_source_target_gene_weight_positive = NA_integer_, n_source_positive_spots_with_neighbor_receiver = NA_integer_,
          fraction_source_positive_with_neighbor_receiver = NA_real_, support_flag = "NOT_EVALUABLE_GROUP_OR_GENE_MISSING",
          limitation = "Visium mixed spots + RCTD inferred weights; spatial compatibility is descriptive, not direct communication")
        } else {
          src <- as.numeric(area[[q$source]]) >= cut & (sp_expr[ax$ligand, k] > 0)
          tar <- as.numeric(area[[q$target]]) >= cut & (colSums(sp_expr[receptor, k, drop = FALSE] >
          0) == length(receptor))
          src[is.na(src)] <- FALSE
          tar[is.na(tar)] <- FALSE
          both <- src & tar
          adj <- area_nets[[nm]]$matrix

          near_receiver <- as.vector(adj %*% as.numeric(tar)) > 0
          near_n <- sum(src & near_receiver)
          row <- data.frame(stage = q$stage, sample_label = nm, axis = q$axis, sender = q$source,
          receiver = q$target, interaction_name = q$interaction_name, source_RCTD_weight_cutoff = cut,
          n_spots = length(k), n_source_weight_and_ligand_positive = sum(src), n_target_weight_and_all_receptors_positive = sum(tar),
          n_same_spot_source_target_gene_weight_positive = sum(both), n_source_positive_spots_with_neighbor_receiver = near_n,
          fraction_source_positive_with_neighbor_receiver = if (sum(src))
            near_n/sum(src) else NA_real_, support_flag = if (sum(src) == 0 || sum(tar) == 0)
            "NO_POSITIVE_SOURCE_OR_RECEIVER_SPOTS_AT_THRESHOLD" else if (sum(both) || near_n)
            "DESCRIPTIVE_SAME_SPOT_OR_NEIGHBOR_COMPATIBILITY" else "NO_SAME_SPOT_OR_NEIGHBOR_COMPATIBILITY_OBSERVED", limitation = "Mixed Visium spots; RCTD and gene expression from same assay; NO experimental proof, NO spot-level inference")
        }
        compat[[length(compat) + 1L]] <- row
      }
    }
  }
}
spat_compat <- if (length(compat)) do.call(rbind, compat) else data.frame(stage = character(), sample_label = character(),
  axis = character(), sender = character(), receiver = character(), support_flag = character(), limitation = character())
write_tab(spat_compat, "22_139_LR_SAME_STAGE_SPATIAL_COMPATIBILITY_BY_CAPTURE_AREA")

if (nrow(spat_compat) && "n_source_positive_spots_with_neighbor_receiver" %in% names(spat_compat)) {
  ss <- split(seq_len(nrow(spat_compat)), paste(spat_compat$stage, spat_compat$axis, spat_compat$sender,
    spat_compat$receiver, spat_compat$interaction_name, spat_compat$source_RCTD_weight_cutoff, sep = "@@"))
  summary_sp <- do.call(rbind, lapply(ss, function(k) {
    z <- spat_compat[k, , drop = FALSE]
    data.frame(stage = z$stage[1], axis = z$axis[1], sender = z$sender[1], receiver = z$receiver[1],
      interaction_name = z$interaction_name[1], cutoff = z$source_RCTD_weight_cutoff[1], n_capture_areas = nrow(z),
      n_capture_areas_with_any_neighbor_support = sum(z$n_source_positive_spots_with_neighbor_receiver >
        0, na.rm = TRUE), n_capture_areas_with_any_same_spot_support = sum(z$n_same_spot_source_target_gene_weight_positive >
        0, na.rm = TRUE), n_capture_areas_not_evaluable = sum(is.na(z$n_source_positive_spots_with_neighbor_receiver)),
      support_areas = paste(z$sample_label[!is.na(z$n_source_positive_spots_with_neighbor_receiver) &
        z$n_source_positive_spots_with_neighbor_receiver > 0], collapse = ";"), note = "Capture-area-level DESCRIPTIVE consistency only; may share animals; no spatial-contact or independent-replicate inference")
  }))
} else summary_sp <- data.frame(stage = character(), axis = character(), sender = character(), receiver = character(),
  note = character())
write_tab(summary_sp, "23_139_LR_SPATIAL_CONSISTENCY_ACROSS_CAPTURE_AREAS")

axis_spatial_status <- merge(axis_summary, if (nrow(summary_sp)) aggregate(n_capture_areas_with_any_neighbor_support ~
  stage + axis, summary_sp, sum) else data.frame(stage = character(), axis = character(), n_capture_areas_with_any_neighbor_support = integer()),
  by = c("stage", "axis"), all.x = TRUE, sort = FALSE)
axis_spatial_status$spatial_status <- ifelse(axis_spatial_status$n_predeclared_canonical_significant_links ==
  0L, "NO_SIGNIFICANT_CANONICAL_CELLCHAT_LINK_DO_NOT_ASSERT_SPATIAL_SIGNAL", ifelse(is.na(axis_spatial_status$n_capture_areas_with_any_neighbor_support),
  "NO_TESTABLE_SPATIAL_LINK", "SPATIAL_AUDIT_REPORTED_SEPARATELY"))
write_tab(axis_spatial_status, "24_139_FOUR_AXIS_COMBINED_EVIDENCE_STATUS")

if (nrow(summary_sp)) {
  p <- file.path(FIG, paste0("139_four_axis_spatial_compatibility_", VERSION, ".pdf"))
  grDevices::pdf(p, width = 11, height = 8)
  for (s in STAGES) {
    z <- summary_sp[summary_sp$stage == s & summary_sp$cutoff == SPATIAL_WEIGHT_THRESHOLDS[1], ,
      drop = FALSE]
    if (nrow(z)) {
      z <- z[order(-z$n_capture_areas_with_any_neighbor_support), , drop = FALSE]
      z <- head(z, 25L)
      labs <- paste(z$axis, z$sender, "to", z$receiver, sep = " | ")
      par(mar = c(4, 20, 4, 2))
      barplot(rev(z$n_capture_areas_with_any_neighbor_support), names.arg = rev(labs), horiz = TRUE,
        las = 1, cex.names = 0.5, xlim = c(0, 2.1), col = "grey45", xlab = "Capture areas with neighboring-spot support (0-2)",
        main = paste(s, "Spatial compatibility (not proof of signaling)"))
    } else {
      plot.new()
      text(0.5, 0.5, paste(s, "no canonical links to evaluate"))
    }
  }
  grDevices::dev.off()
}
rm(spot, sp_expr, area_nets, areas, compat)
gc(verbose = FALSE)

external_status <- "NOT_RUN_USER_OPTION"
if (RUN_GSE281205_VALIDATION) {
  files_281 <- unlist(IN[c("gse281205_object", "gse281205_manifest", "gse281205_features")])
  if (!all(file.exists(files_281))) {
    external_status <- "NOT_RUN_MISSING_STEP124_130_131_INPUT"
    log_it("GSE281205 external per-GSM expression not run: missing file(s)")
  } else {
    log_it("GSE281205: optional per-GSM four-axis gene expression validation")

    external_status <- tryCatch({
      d <- read_table(IN$gse281205_manifest)
      feat <- read_table(IN$gse281205_features)
      need_cols(d, c("cell_id", "gsm", "dpi", "primary_group_131", "included_primary"), "Step131 per-cell manifest")
      need_cols(feat, c("gene_symbol", "gene_id"), "Step124 mapping")
      assert(!anyDuplicated(d$cell_id), "Step131 duplicated cell IDs")
      valid_groups <- !is.na(d$included_primary) & tolower(as.character(d$included_primary)) %in%
        c("true", "t", "1")
      d <- d[valid_groups & d$dpi %in% c(14, 28), , drop = FALSE]
      ids <- setNames(seq_len(nrow(d)), d$cell_id)
      mapgene <- lapply(FOCUS_GENES, function(g) unique(as.character(feat$gene_id[feat$gene_symbol ==
        g])))
      names(mapgene) <- FOCUS_GENES
      gmap <- data.frame(gene = FOCUS_GENES, n_id_matches = lengths(mapgene), gene_id = vapply(mapgene,
        function(z) if (length(z) == 1L)
          z else NA_character_, ""))
      write_tab(gmap, "25_139_GSE281205_GENE_MAPPING_AUDIT")

      obj <- readRDS(IN$gse281205_object)
      assert(inherits(obj, "Seurat"), "Step130 GSE281205 is not Seurat")
      assert("DecontX" %in% SeuratObject::Assays(obj), "Step130 DecontX assay missing")
      counts <- get_layer(obj, "DecontX", "counts")
      ix <- match(d$cell_id, colnames(counts))
      assert(!anyNA(ix), "GSE281205 Step131 cell IDs mismatch Step130 RDS")
      mgenes <- gmap$gene_id
      valid <- !is.na(mgenes) & mgenes %in% rownames(counts)
      gmap$feature_found <- valid
      write_tab(gmap, "26_139_GSE281205_FEATURE_AND_CELL_ID_AUDIT")
      ind <- split(seq_len(nrow(d)), paste(d$dpi, d$gsm, d$primary_group_131, sep = "@@"))
      e <- do.call(rbind, lapply(ind, function(k) {
        do.call(rbind, lapply(seq_along(FOCUS_GENES), function(i) {
          v <- if (valid[i])
          as.numeric(counts[mgenes[i], ix[k], drop = TRUE]) else rep(NA_real_, length(k))
          data.frame(stage = paste0(d$dpi[k[1]], "dpi"), gsm = d$gsm[k[1]], group_131 = d$primary_group_131[k[1]],
          gene = FOCUS_GENES[i], feature_identified = valid[i], n_cells = length(k), fraction_DecontX_corrected_counts_ge1 = if (valid[i])
            mean(v >= 1) else NA_real_, mean_DecontX_corrected_counts = if (valid[i])
            mean(v) else NA_real_, note = "Separate GSE281205 cohort: Step131 working identities; DecontX corrected counts; NOT replicated CellChat")
        }))
      }))
      write_tab(e, "27_139_GSE281205_PER_GSM_GROUP_FOUR_AXIS_EXPRESSION")

      both_genes <- c("Trem2", "Tyrobp")
      if (all(both_genes %in% gmap$gene[valid])) {
        rr <- gmap$gene_id[match(both_genes, gmap$gene)]
        ce <- do.call(rbind, lapply(ind, function(k) {
          co <- Matrix::colSums(counts[rr, ix[k], drop = FALSE] >= 1) == 2L
          data.frame(stage = paste0(d$dpi[k[1]], "dpi"), gsm = d$gsm[k[1]], group_131 = d$primary_group_131[k[1]],
          n_cells = length(k), Trem2_Tyrobp_same_cell_codetection_fraction = mean(co), note = "Both subunits detected in the same single cell; descriptive")
        }))
        write_tab(ce, "28_139_GSE281205_TREM2_TYROBP_CODETECTION")
      }
      rm(d, feat, obj, counts, e, ix)
      gc(verbose = FALSE)
      "COMPLETED_DESCRIPTIVE_GSM_EXPRESSION_NOT_CELLCHAT_REPLICATION"
    }, error = function(e) {
      log_it("WARNING GSE281205 external validation failed: ", conditionMessage(e))
      paste0("FAILED_", conditionMessage(e))
    })
  }
}
write_tab(data.frame(component = "GSE281205 per-GSM external expression", status = external_status, interpretation = "Different cohort/working states. Per-GSM gene detection is corroboration only, not independently replicated CellChat inference."),
  "29_139_GSE281205_EXTERNAL_VALIDATION_STATUS")

write_tab(data.frame(issue = c("14dpi within-dataset CellChat", "28dpi within-dataset CellChat", "Same-stage cross-source CellChat replicates",
  "Independent animal CellChat replicates", "Nine refined myeloid states", "APOE/PSAP/SPP1/APP CellChat",
  "Within-stage ligand and same-cell receptor complex expression", "Same-stage spatial RCTD/Visium compatibility",
  "External GSE281205 GSM-level myeloid expression", "Autocrine or direct binding validated"), status = c("PERFORMED",
  "PERFORMED", "NOT_IDENTIFIABLE_IN_STEP77_LATE_STAGE_COHORT", "NOT_VERIFIED_FROM_SOURCE_METADATA",
  "PERFORMED_WITH_MIN_10_CELLS_PER_STAGE", "PERFORMED", "PERFORMED", "PERFORMED_DESCRIPTIVE_CAPTURE_AREA",
  external_status, "NO_EXPERIMENTAL_VALIDATION"), scientific_note = c("14dpi is Figshare 17702045",
  "28dpi is GSE196928", "Different source at each stage: stage and source confounded", "GSM, library, section, portion and individual cells are not verified animals",
  "No Other myeloid pooling; subtypes <10 cells removed at their own stage", "Bootstrap tests cell-group model, not specimen-level reproducibility",
  "L/R co-detection supports transcript presence, not receptor function", "Spot mixture, RCTD and RNA share data; no pooled-spot p-values",
  "Independent sample evidence only for the genes and groups actually observed", "Within-population APOE model output does not establish autocrine signaling")),
  "30_139_VALIDATION_STATUS")
writeLines("Numerical results are provided in the accompanying output tables.", 
    con = file.path(OUT, "139_METHODS_AND_LIMITATIONS_v03.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), con = file.path(OUT, "139_sessionInfo_v03.txt"), useBytes = TRUE)
qcnames <- c("INPUTS_ATLAS_MYELOID_SPATIAL_RCTD", "STEP77_COHORT_CHECK_IF_AVAILABLE", "REFINED_CELLCHAT_TWO_STAGE_RDS",
  "GROUPS_REFINED_NOT_COLLAPSED", "FOUR_AXIS_STATUS_ROWS_8", "SAME_STAGE_EXPRESSION_TABLE", "SPATIAL_GENE_FEATURE_AUDIT",
  "SPATIAL_CAPTURE_AREA_GENE_DETECTION", "SPATIAL_COMPATIBILITY_OUTPUT", "SAMPLE_IDENTIFIABILITY_EXPLICIT",
  "NO_ZIP_DEPENDENCY")
qc <- data.frame(check = qcnames, pass = c(all(file.exists(unlist(IN[c("atlas", "myeloid", "spatial",
  "rctd")]))), if (file.exists(IN$step77)) all(refcheck$stage_matches_reference) else NA, all(file.exists(vapply(STAGES,
  final_file, ""))), all(!vapply(all_results, function(z) "Other myeloid" %in% z$groups, FALSE)), nrow(axis_summary) ==
  8L, nrow(ex) > 0 && nrow(co) > 0, file.exists(file.path(TAB, "17_139_SPATIAL_FOCAL_GENE_FEATURE_AVAILABILITY_v03.csv")),
  nrow(spatial_gene_detection) > 0, file.exists(file.path(TAB, "22_139_LR_SAME_STAGE_SPATIAL_COMPATIBILITY_BY_CAPTURE_AREA_v03.csv")),
  all(is.na(feas$n_verified_independent_animals)), TRUE))
write_tab(qc, "31_139_FINAL_QC")
allfiles <- list.files(OUT, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
fileidx <- data.frame(relative_path = relative(allfiles), size_bytes = as.numeric(file.info(allfiles)$size),
  stringsAsFactors = FALSE)
fileidx <- fileidx[order(fileidx$relative_path), , drop = FALSE]
write_tab(fileidx, "00_139_OUTPUT_INDEX")
log_it("SUCCESS STEP139 v03 | See output: ", OUT)
log_it("Outputs include CSV/TXT/LOG/PDF; RDS/checkpoints are not required.")
future::plan(future::sequential)

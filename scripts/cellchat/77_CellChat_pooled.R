rm(list = ls())

options(stringsAsFactors = FALSE)
options(warn = 1)

set.seed(7701)

options(future.globals.maxSize = 8 * 1024^3)

ROOT <- "D:/zz"

FULL_ATLAS_FILE <- file.path(ROOT, "results", "09_Figure1A_major_celltype_UMAP", "09_Figure1A_major_celltype_UMAP_object_v04.rds")

MYELOID_FILE <- file.path(ROOT, "results", "19_myeloid_state_robustness_temporal_dynamics", "19_myeloid_state_robustness_temporal_dynamics_3scRNA_v02.rds")

OUT_DIR <- file.path(ROOT, "results", "77_Global_CellChat_Clec7a_centered")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v05"

TARGET_GROUP <- "Clec7a+ lipid-handling Microglia"
OTHER_MYELOID_GROUP <- "Other myeloid"

TARGET_GROUP_FROZEN <- "Clec7a+ lipid-associated Microglia"

MYELOID_PARENT_TYPES <- c("Microglia", "Macrophage/Monocyte", "APC/DC-like myeloid")

MIN_CELLS <- 10L
NBOOT <- 100L
SEED_USE <- 7701L

N_WORKERS <- 3L

FORCE_RESTART <- FALSE

CHECKPOINT_DIR <- file.path(OUT_DIR, "_checkpoints_v05")

dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

CHECKPOINT_FILES <- c(subsetData = file.path(CHECKPOINT_DIR, "01_after_subsetData_v05.rds"), overexpression = file.path(CHECKPOINT_DIR, 
    "02_after_overexpression_v05.rds"), communication_probability = file.path(CHECKPOINT_DIR, "03_after_computeCommunProb_filtered_v05.rds"), 
    final_network = file.path(CHECKPOINT_DIR, "04_after_pathway_aggregate_centrality_v05.rds"))

if (isTRUE(FORCE_RESTART)) {
    old_ckpt <- CHECKPOINT_FILES[file.exists(CHECKPOINT_FILES)]
    if (length(old_ckpt) > 0L) {
        unlink(old_ckpt, force = TRUE)
    }
}

required_pkgs <- c("Seurat", "SeuratObject", "CellChat", "Matrix", "dplyr", "tidyr", "ggplot2", "patchwork", 
    "future")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them before rerunning Step77.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(future)
})

if (N_WORKERS <= 1L) {
    future::plan("sequential")
} else {
    future::plan("multisession", workers = N_WORKERS)
}

message("future.globals.maxSize = ", round(getOption("future.globals.maxSize")/1024^3, 1), " GiB")

write_csv_utf8 <- function(x, path) {
    write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

write_csv_gz <- function(x, path) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(x, file = con, row.names = FALSE, na = "")
}

save_cellchat_checkpoint <- function(cellchat, path, stage_label) {
    message("Saving checkpoint: ", stage_label)
    saveRDS(cellchat, file = path, compress = FALSE)
    if (!file.exists(path)) {
        stop("Checkpoint was not created: ", path)
    }
    message("Checkpoint saved: ", path)
    invisible(path)
}

load_latest_valid_checkpoint <- function(
  checkpoint_files
) {

  idx <- rev(
    seq_along(
      checkpoint_files
    )
  )

  for (i in idx) {

    path <- checkpoint_files[
      i
    ]

    if (!file.exists(path)) {
      next
    }

    message(
      "Found checkpoint candidate: ",
      path
    )

    obj <- tryCatch(
      readRDS(path),
      error = function(e) {
        warning(
          "Could not read checkpoint; trying an earlier one.\n",
          path,
          "\nReason: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(obj)) {

      message(
        "Resuming from checkpoint stage ",
        i,
        ": ",
        names(checkpoint_files)[i]
      )

      return(
        list(
          cellchat = obj,
          stage = i,
          stage_name = names(checkpoint_files)[i],
          path = path
        )
      )
    }
  }

  list(
    cellchat = NULL,
    stage = 0L,
    stage_name = "none",
    path = NA_character_
  )
}

sort_prob_desc <- function(x) {
    if (is.data.frame(x) && nrow(x) > 0L && "prob" %in% colnames(x)) {
        x <- x[order(x$prob, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    }
    x
}

summarize_lr_by_partner <- function(df, direction = c("outgoing", "incoming")) {
    direction <- match.arg(direction)
    if (is.null(df) || nrow(df) == 0L) {
        return(data.frame(direction = character(0), partner = character(0), n_LR_pairs = integer(0), 
            total_probability = numeric(0), max_probability = numeric(0), stringsAsFactors = FALSE))
    }
    partner_col <- if (direction == "outgoing") {
        "target"
    }
    else {
        "source"
    }
    out <- df %>% mutate(partner = .data[[partner_col]]) %>% group_by(partner) %>% summarise(n_LR_pairs = n(), 
        total_probability = sum(prob, na.rm = TRUE), max_probability = max(prob, na.rm = TRUE), .groups = "drop") %>% 
        mutate(direction = direction, .before = 1) %>% arrange(desc(total_probability))
    out
}

summarize_pathway_by_partner <- function(df, direction = c("outgoing", "incoming")) {
    direction <- match.arg(direction)
    if (is.null(df) || nrow(df) == 0L) {
        return(data.frame(direction = character(0), partner = character(0), pathway_name = character(0), 
            probability = numeric(0), stringsAsFactors = FALSE))
    }
    partner_col <- if (direction == "outgoing") {
        "target"
    }
    else {
        "source"
    }
    pathway_col <- if ("pathway_name" %in% colnames(df)) {
        "pathway_name"
    }
    else if ("pathway" %in% colnames(df)) {
        "pathway"
    }
    else {
        stop("Could not identify pathway-name column in CellChat pathway table.")
    }
    df %>% transmute(direction = direction, partner = .data[[partner_col]], pathway_name = .data[[pathway_col]], 
        probability = prob) %>% arrange(desc(probability))
}

for (f in c(FULL_ATLAS_FILE, MYELOID_FILE)) {
    if (!file.exists(f)) {
        stop("Required input file not found:\n", f)
    }
}

full_obj <- readRDS(FULL_ATLAS_FILE)

myeloid_obj <- readRDS(MYELOID_FILE)

if (!inherits(full_obj, "Seurat")) {
    stop("Full atlas input is not a Seurat object.")
}

if (!inherits(myeloid_obj, "Seurat")) {
    stop("Myeloid input is not a Seurat object.")
}

full_md <- full_obj[[]]
my_md <- myeloid_obj[[]]

required_full_meta <- c("major_celltype")

missing_full_meta <- setdiff(required_full_meta, colnames(full_md))

if (length(missing_full_meta) > 0L) {
    stop("Full atlas missing metadata field(s): ", paste(missing_full_meta, collapse = ", "))
}

if (!"myeloid_state" %in% colnames(my_md)) {
    stop("myeloid object is missing myeloid_state.")
}

full_cells <- Cells(full_obj)

my_cells <- Cells(myeloid_obj)

missing_myeloid_cells <- setdiff(my_cells, full_cells)

if (length(missing_myeloid_cells) > 0L) {
    stop(length(missing_myeloid_cells), " myeloid cells are absent from the full atlas. ", "Exact cell-name transfer cannot proceed.")
}

my_state <- as.character(my_md$myeloid_state)

names(my_state) <- rownames(my_md)

if (!TARGET_GROUP_FROZEN %in% my_state) {
    stop("target state not found: ", TARGET_GROUP_FROZEN)
}

full_md$cell <- rownames(full_md)

full_md$major_celltype_chr <- as.character(full_md$major_celltype)

full_md$in_refined_myeloid <- (full_md$cell %in% my_cells)

full_md$is_original_myeloid_parent <- (full_md$major_celltype_chr %in% MYELOID_PARENT_TYPES)

full_md$removed_during_myeloid_cleanup <- (full_md$is_original_myeloid_parent & !full_md$in_refined_myeloid)

full_md$communication_group <- full_md$major_celltype_chr

idx_refined <- which(full_md$in_refined_myeloid)

refined_state_for_full <- my_state[full_md$cell[idx_refined]]

full_md$communication_group[idx_refined] <- ifelse(refined_state_for_full == TARGET_GROUP_FROZEN, TARGET_GROUP, 
    OTHER_MYELOID_GROUP)

formal_keep <- (!full_md$removed_during_myeloid_cleanup & !is.na(full_md$communication_group) & nzchar(trimws(full_md$communication_group)))

comm_md <- full_md[formal_keep, , drop = FALSE]

preferred_nonmyeloid_order <- c("Neutrophil", "T/NK cell", "Astrocyte", "Oligodendrocyte", "OPC", "Neuron", 
    "Ependymal", "Endothelial", "Pericyte/VSMC", "Fibroblast/Meningeal")

observed_groups <- unique(as.character(comm_md$communication_group))

other_groups <- setdiff(observed_groups, c(TARGET_GROUP, OTHER_MYELOID_GROUP))

ordered_other_groups <- c(preferred_nonmyeloid_order[preferred_nonmyeloid_order %in% other_groups], sort(setdiff(other_groups, 
    preferred_nonmyeloid_order)))

GROUP_LEVELS <- c(TARGET_GROUP, OTHER_MYELOID_GROUP, ordered_other_groups)

comm_md$communication_group <- factor(as.character(comm_md$communication_group), levels = GROUP_LEVELS)

group_counts <- comm_md %>% count(communication_group, name = "n_cells", .drop = FALSE) %>% mutate(eligible_ge10 = n_cells >= 
    MIN_CELLS)

if (any(group_counts$n_cells < MIN_CELLS)) {
    bad_groups <- as.character(group_counts$communication_group[group_counts$n_cells < MIN_CELLS])
    warning("The following pooled group(s) have < ", MIN_CELLS, " cells and will be removed before CellChat: ", 
        paste(bad_groups, collapse = " | "))
    GROUP_LEVELS <- setdiff(GROUP_LEVELS, bad_groups)
    comm_md <- comm_md[as.character(comm_md$communication_group) %in% GROUP_LEVELS, , drop = FALSE]
    comm_md$communication_group <- factor(as.character(comm_md$communication_group), levels = GROUP_LEVELS)
    group_counts <- comm_md %>% count(communication_group, name = "n_cells", .drop = FALSE) %>% mutate(eligible_ge10 = n_cells >= 
        MIN_CELLS)
}

if (!TARGET_GROUP %in% as.character(group_counts$communication_group)) {
    stop("Target group was lost during pooled-group QC.")
}

if (group_counts$n_cells[group_counts$communication_group == TARGET_GROUP] < MIN_CELLS) {
    stop("Target group contains fewer than ", MIN_CELLS, " cells.")
}

write_csv_utf8(group_counts, file.path(OUT_DIR, "01_77_Global_CellChat_group_counts_v05.csv"))

provenance_cols <- intersect(c("source_dataset", "timepoint_main", "time", "sample_id"), colnames(comm_md))

if (length(provenance_cols) > 0L) {
    prov_df <- comm_md %>% mutate(communication_group = as.character(communication_group)) %>% group_by(across(all_of(c(provenance_cols, 
        "communication_group")))) %>% summarise(n_cells = n(), .groups = "drop")
    write_csv_utf8(prov_df, file.path(OUT_DIR, "02_77_Global_CellChat_input_provenance_v05.csv"))
}

metadata_out <- data.frame(cell = rownames(comm_md), communication_group = as.character(comm_md$communication_group), 
    major_celltype = as.character(comm_md$major_celltype_chr), refined_myeloid = comm_md$in_refined_myeloid, 
    stringsAsFactors = FALSE)

write_csv_gz(metadata_out, file.path(OUT_DIR, "03_77_Global_CellChat_cell_metadata_v05.csv.gz"))

keep_cells <- rownames(comm_md)

obj_comm <- subset(full_obj, cells = keep_cells)

rm(full_obj)

gc(verbose = FALSE)

DefaultAssay(obj_comm) <- "RNA"

if (inherits(obj_comm[["RNA"]], "Assay5")) {
    layer_names <- Layers(obj_comm[["RNA"]])
    if (length(layer_names) > 1L) {
        obj_comm <- JoinLayers(obj_comm, assay = "RNA")
    }
}

obj_comm$communication_group <- comm_md[Cells(obj_comm), "communication_group"]

if (anyNA(obj_comm$communication_group)) {
    stop("NA communication_group detected after Seurat subset.")
}

obj_comm$communication_group <- factor(as.character(obj_comm$communication_group), levels = GROUP_LEVELS)

rna_layers_now <- if (inherits(obj_comm[["RNA"]], "Assay5")) {
    Layers(obj_comm[["RNA"]])
} else {
    c("counts", "data", "scale.data")
}

has_data_layer <- any(rna_layers_now == "data")

if (!has_data_layer) {
    obj_comm <- NormalizeData(obj_comm, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, 
        verbose = FALSE)
}

CellChatDB.use <- subsetDB(CellChatDB.mouse)

db_gene_set <- extractGene(CellChatDB.use)

features_use <- intersect(rownames(obj_comm), db_gene_set)

if (length(features_use) < 100L) {
    stop("Too few CellChatDB genes matched the RNA assay: ", length(features_use))
}

data_input <- tryCatch(GetAssayData(obj_comm, assay = "RNA", layer = "data")[features_use, , drop = FALSE], 
    error = function(e) {
        GetAssayData(obj_comm, assay = "RNA", slot = "data")[features_use, , drop = FALSE]
    })

meta_input <- data.frame(communication_group = obj_comm$communication_group, row.names = Cells(obj_comm), 
    stringsAsFactors = FALSE)

if (!identical(colnames(data_input), rownames(meta_input))) {
    stop("Cell order mismatch between expression matrix and metadata.")
}

rm(obj_comm, myeloid_obj, full_md, my_md)

gc(verbose = FALSE)

resume_info <- load_latest_valid_checkpoint(CHECKPOINT_FILES)

resume_stage <- resume_info$stage

if (resume_stage > 0L) {
    cellchat <- resume_info$cellchat
    message("The script will skip those completed stages.")
} else {
    cellchat <- createCellChat(object = data_input, meta = meta_input, group.by = "communication_group")
    cellchat@DB <- CellChatDB.use
}

if (resume_stage < 1L) {
    cellchat <- subsetData(cellchat)
    save_cellchat_checkpoint(cellchat, CHECKPOINT_FILES["subsetData"], "after subsetData")
    resume_stage <- 1L
} else {
    message("SKIP: subsetData already completed.")
}

if (resume_stage < 2L) {
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    save_cellchat_checkpoint(cellchat, CHECKPOINT_FILES["overexpression"], "after identifyOverExpressedGenes + identifyOverExpressedInteractions")
    resume_stage <- 2L
} else {
    message("SKIP: over-expression analysis already completed.")
}

if (resume_stage < 3L) {
    set.seed(SEED_USE)
    cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = FALSE, 
        nboot = NBOOT, seed.use = SEED_USE)
    cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
    save_cellchat_checkpoint(cellchat, CHECKPOINT_FILES["communication_probability"], "after computeCommunProb + filterCommunication")
    resume_stage <- 3L
} else {
    message("SKIP: communication probabilities already completed.")
}

if (resume_stage < 4L) {
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    save_cellchat_checkpoint(cellchat, CHECKPOINT_FILES["final_network"], "after pathway + aggregateNet + centrality")
    resume_stage <- 4L
} else {
    message("SKIP: pathway/network/centrality analysis already completed.")
}

CELLCHAT_RDS <- file.path(OUT_DIR, "04_77_Global_CellChat_object_v05.rds")

saveRDS(cellchat, CELLCHAT_RDS, compress = FALSE)

all_lr <- subsetCommunication(cellchat, slot.name = "net")

all_lr <- sort_prob_desc(all_lr)

write_csv_gz(all_lr, file.path(OUT_DIR, "05_77_Global_CellChat_all_LR_interactions_v05.csv.gz"))

all_pathways <- subsetCommunication(cellchat, slot.name = "netP")

all_pathways <- sort_prob_desc(all_pathways)

write_csv_gz(all_pathways, file.path(OUT_DIR, "06_77_Global_CellChat_all_pathway_interactions_v05.csv.gz"))

target_idx <- which(levels(cellchat@idents) == TARGET_GROUP)

if (length(target_idx) != 1L) {
    stop("Could not uniquely identify target CellChat identity.")
}

clec_out_lr <- subsetCommunication(cellchat, slot.name = "net", sources.use = target_idx)

clec_out_lr <- sort_prob_desc(clec_out_lr)

clec_in_lr <- subsetCommunication(cellchat, slot.name = "net", targets.use = target_idx)

clec_in_lr <- sort_prob_desc(clec_in_lr)

clec_out_path <- subsetCommunication(cellchat, slot.name = "netP", sources.use = target_idx)

clec_out_path <- sort_prob_desc(clec_out_path)

clec_in_path <- subsetCommunication(cellchat, slot.name = "netP", targets.use = target_idx)

clec_in_path <- sort_prob_desc(clec_in_path)

write_csv_utf8(clec_out_lr, file.path(OUT_DIR, "07_77_Clec7a_outgoing_LR_v05.csv"))

write_csv_utf8(clec_in_lr, file.path(OUT_DIR, "08_77_Clec7a_incoming_LR_v05.csv"))

write_csv_utf8(clec_out_path, file.path(OUT_DIR, "09_77_Clec7a_outgoing_pathways_v05.csv"))

write_csv_utf8(clec_in_path, file.path(OUT_DIR, "10_77_Clec7a_incoming_pathways_v05.csv"))

out_lr_partner <- summarize_lr_by_partner(clec_out_lr, direction = "outgoing")

in_lr_partner <- summarize_lr_by_partner(clec_in_lr, direction = "incoming")

partner_lr_summary <- bind_rows(out_lr_partner, in_lr_partner)

write_csv_utf8(partner_lr_summary, file.path(OUT_DIR, "11_77_Clec7a_partner_LR_summary_v05.csv"))

out_path_partner <- summarize_pathway_by_partner(clec_out_path, direction = "outgoing")

in_path_partner <- summarize_pathway_by_partner(clec_in_path, direction = "incoming")

partner_pathway_summary <- bind_rows(out_path_partner, in_path_partner)

write_csv_utf8(partner_pathway_summary, file.path(OUT_DIR, "12_77_Clec7a_partner_pathway_summary_v05.csv"))

count_mat <- cellchat@net$count
weight_mat <- cellchat@net$weight

target_partner_matrix <- data.frame(partner = rownames(weight_mat), outgoing_count = as.numeric(count_mat[TARGET_GROUP, 
    , drop = TRUE]), incoming_count = as.numeric(count_mat[, TARGET_GROUP, drop = TRUE]), outgoing_weight = as.numeric(weight_mat[TARGET_GROUP, 
    , drop = TRUE]), incoming_weight = as.numeric(weight_mat[, TARGET_GROUP, drop = TRUE]), stringsAsFactors = FALSE)

write_csv_utf8(target_partner_matrix, file.path(OUT_DIR, "13_77_Clec7a_partner_network_matrix_v05.csv"))

PDF_FILE <- file.path(OUT_DIR, "14_77_Global_CellChat_diagnostics_v05.pdf")

grDevices::cairo_pdf(filename = PDF_FILE, width = 11, height = 8.5, onefile = TRUE)

p_counts <- ggplot(group_counts, aes(x = reorder(as.character(communication_group), n_cells), y = n_cells)) + 
    geom_col(width = 0.72) + coord_flip() + labs(title = "Global pooled CellChat input", subtitle = "All datasets and all available time points pooled", 
    x = NULL, y = "Cell count") + theme_classic(base_size = 10)

print(p_counts)

group_size <- as.numeric(table(cellchat@idents))

names(group_size) <- names(table(cellchat@idents))

par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, 
    title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, 
    title.name = "Interaction strength")

out_count_mat <- matrix(0, nrow = nrow(count_mat), ncol = ncol(count_mat), dimnames = dimnames(count_mat))

out_count_mat[TARGET_GROUP, ] <- count_mat[TARGET_GROUP, ]

in_count_mat <- matrix(0, nrow = nrow(count_mat), ncol = ncol(count_mat), dimnames = dimnames(count_mat))

in_count_mat[, TARGET_GROUP] <- count_mat[, TARGET_GROUP]

par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(out_count_mat, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, 
    edge.weight.max = max(count_mat, na.rm = TRUE), title.name = "Clec7a outgoing | interaction number")

netVisual_circle(in_count_mat, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, edge.weight.max = max(count_mat, 
    na.rm = TRUE), title.name = "Clec7a incoming | interaction number")

out_weight_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))

out_weight_mat[TARGET_GROUP, ] <- weight_mat[TARGET_GROUP, ]

in_weight_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))

in_weight_mat[, TARGET_GROUP] <- weight_mat[, TARGET_GROUP]

par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(out_weight_mat, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, 
    edge.weight.max = max(weight_mat, na.rm = TRUE), title.name = "Clec7a outgoing | interaction strength")

netVisual_circle(in_weight_mat, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, 
    edge.weight.max = max(weight_mat, na.rm = TRUE), title.name = "Clec7a incoming | interaction strength")

partner_plot_df <- target_partner_matrix %>% filter(partner != TARGET_GROUP) %>% dplyr::select(partner, 
    outgoing_weight, incoming_weight) %>% pivot_longer(cols = c(outgoing_weight, incoming_weight), names_to = "direction", 
    values_to = "weight") %>% mutate(direction = recode(direction, outgoing_weight = "Clec7a -> partner", 
    incoming_weight = "partner -> Clec7a"))

p_partner <- ggplot(partner_plot_df, aes(x = reorder(partner, weight), y = weight, fill = direction)) + 
    geom_col(position = "dodge", width = 0.72) + coord_flip() + labs(title = "Clec7a-centered global communication strength", 
    x = NULL, y = "Aggregated CellChat interaction strength", fill = NULL) + theme_classic(base_size = 10)

print(p_partner)

p_role <- netAnalysis_signalingRole_scatter(cellchat)

print(p_role)

dev.off()

qc <- data.frame(check = c("Full atlas input exists", "myeloid input exists", "Every refined myeloid cell exists in full atlas", 
    "Target state is present", "Other refined myeloid states collapsed to Other myeloid", "No timepoint stratification used", 
    "No dataset stratification used", "No downsampling used", "population.size = FALSE", "CellChatDB non-protein signaling excluded", 
    "Target outgoing LR table generated", "Target incoming LR table generated", "Target outgoing pathway table generated", 
    "Target incoming pathway table generated", "CellChat object saved", "Diagnostic PDF generated"), 
    pass = c(file.exists(FULL_ATLAS_FILE), file.exists(MYELOID_FILE), length(missing_myeloid_cells) == 
        0L, TARGET_GROUP %in% levels(cellchat@idents), OTHER_MYELOID_GROUP %in% levels(cellchat@idents), 
        TRUE, TRUE, TRUE, TRUE, TRUE, file.exists(file.path(OUT_DIR, "07_77_Clec7a_outgoing_LR_v05.csv")), 
        file.exists(file.path(OUT_DIR, "08_77_Clec7a_incoming_LR_v05.csv")), file.exists(file.path(OUT_DIR, 
            "09_77_Clec7a_outgoing_pathways_v05.csv")), file.exists(file.path(OUT_DIR, "10_77_Clec7a_incoming_pathways_v05.csv")), 
        file.exists(CELLCHAT_RDS), file.exists(PDF_FILE)), stringsAsFactors = FALSE)

write_csv_utf8(qc, file.path(OUT_DIR, "15_77_Global_CellChat_FINAL_QC_v05.csv"))

REPORT_FILE <- file.path(OUT_DIR, "16_77_Global_CellChat_REPORT_v05.txt")

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(OUT_DIR, "17_77_sessionInfo_v05.txt"))

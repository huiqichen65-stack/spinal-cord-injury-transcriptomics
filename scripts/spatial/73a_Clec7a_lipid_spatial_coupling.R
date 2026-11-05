source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 1200)

set.seed(7301)

ROOT <- "D:/zz"

STEP64_DIR <- file.path(ROOT, "results", "64_GSE195783_spatial_cluster_marker_region_annotation")

STEP65_DIR <- file.path(ROOT, "results", "65_GSE195783_RCTD_myeloid_state_mapping")

SPATIAL_OBJECT <- file.path(STEP64_DIR, "12_64_GSE195783_spatial_region_annotated_v03.rds")

RCTD_FILE <- file.path(STEP65_DIR, "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv.gz")

OUT_DIR <- file.path(ROOT, "results", "73_GSE195783_Clec7a_lipid_spatial_coupling")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v03"

TARGET_STATE <- "Clec7a+ lipid-handling Microglia"

TIME_LEVELS <- c("3dpi", "7dpi", "14dpi", "28dpi")

SAMPLE_ORDER <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1", 
    "28dpi_rep2")

EXPECTED_CAPTURE_AREAS <- 8L
EXPECTED_MAJOR_SECTIONS <- 32L
MIN_SECTION_SPOTS <- 100L

GRAPH_RADIUS_MULTIPLIER <- 1.55

FORMAL_ANALYSES <- c("direct_same_spot", "local_neighbor_lipid")

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "ggplot2", "patchwork", "spdep", 
    "scales")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, then rerun Step73.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(spdep)
  library(scales)
})

write_csv_utf8 <- write_csv_utf8

write_csv_gz_base <- write_gzip_csv_utf8

read_csv_gz_base <- function(filename) {
    read.csv(gzfile(filename), stringsAsFactors = FALSE, check.names = FALSE)
}

pick_first_existing_col <- function(df, candidates, label) {
    hit <- candidates[candidates %in% colnames(df)]
    if (length(hit) == 0L) {
        stop("Could not find column for ", label, ".\nCandidates: ", paste(candidates, collapse = " | "), 
            "\nAvailable columns:\n", paste(colnames(df), collapse = " | "))
    }
    hit[[1]]
}

safe_spearman <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) < 10L) {
        return(NA_real_)
    }
    if (!is.finite(sd(x)) || !is.finite(sd(y)) || sd(x) == 0 || sd(y) == 0) {
        return(NA_real_)
    }
    suppressWarnings(cor(x, y, method = "spearman"))
}

fisher_z <- function(r) {
    if (!is.finite(r)) {
        return(NA_real_)
    }
    r <- max(min(r, 0.999999), -0.999999)
    atanh(r)
}

inv_fisher_z <- function(z) {
    ifelse(is.finite(z), tanh(z), NA_real_)
}

wilcox_greater_zero <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3L) {
        return(NA_real_)
    }
    if (all(abs(x) < 1e-12)) {
        return(1)
    }
    suppressWarnings(wilcox.test(x, mu = 0, alternative = "greater", exact = FALSE, conf.int = FALSE)$p.value)
}

build_visium_graph <- function(coords, radius_multiplier = GRAPH_RADIUS_MULTIPLIER) {
    if (nrow(coords) < 3L) {
        stop("Too few coordinates to build a spatial graph.")
    }
    if (anyDuplicated(data.frame(x = coords[, 1], y = coords[, 2])) > 0L) {
        stop("Duplicated array coordinates detected within a tissue section.")
    }
    kn1 <- spdep::knearneigh(coords, k = 1, longlat = FALSE)
    nn1 <- spdep::knn2nb(kn1, row.names = rownames(coords))
    nn_distance <- vapply(seq_along(nn1), function(i) {
        ids <- nn1[[i]]
        if (length(ids) == 0L) {
            return(NA_real_)
        }
        sqrt(sum((coords[i, ] - coords[ids[1], ])^2))
    }, FUN.VALUE = numeric(1))
    nn_distance <- nn_distance[is.finite(nn_distance) & nn_distance > 0]
    if (length(nn_distance) == 0L) {
        stop("Could not estimate nearest-neighbor distance.")
    }
    median_nn <- median(nn_distance)
    radius <- radius_multiplier * median_nn
    nb <- spdep::dnearneigh(coords, d1 = 0, d2 = radius, row.names = rownames(coords), longlat = FALSE)
    list(nb = nb, median_nn_distance = median_nn, radius = radius)
}

neighbor_mean_from_nb <- function(x, nb) {
    vapply(seq_along(nb), function(i) {
        ids <- nb[[i]]
        if (length(ids) == 0L) {
            return(NA_real_)
        }
        vals <- x[ids]
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0L) {
            return(NA_real_)
        }
        mean(vals)
    }, FUN.VALUE = numeric(1))
}

get_spatial_data_matrix <- function(obj) {
    if (!("Spatial" %in% names(obj@assays))) {
        stop("Spatial assay is absent from the Step64 object.")
    }
    mat <- tryCatch(SeuratObject::LayerData(object = obj, assay = "Spatial", layer = "data"), error = function(e) NULL)
    if (is.null(mat)) {
        mat <- tryCatch(Seurat::GetAssayData(object = obj, assay = "Spatial", slot = "data"), error = function(e) NULL)
    }
    if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
        stop("Could not retrieve log-normalized Spatial assay data.")
    }
    mat
}

match_gene_set <- function(requested_genes, available_genes) {
    req <- unique(requested_genes)
    avail_low <- tolower(available_genes)
    req_low <- tolower(req)
    idx <- match(req_low, avail_low)
    matched <- available_genes[idx[!is.na(idx)]]
    requested_matched <- req[!is.na(idx)]
    data.frame(requested_gene = req, matched_gene = ifelse(is.na(idx), NA_character_, available_genes[idx]), 
        detected = !is.na(idx), stringsAsFactors = FALSE)
}

lipid_gene_sets <- list(FA_uptake_transport = c("Cd36", "Slc27a1", "Slc27a2", "Slc27a3", "Slc27a4", "Slc27a5", 
    "Fabp3", "Fabp4", "Fabp5", "Lpl", "Lrp1", "Vldlr"), Cholesterol_efflux_transport = c("Abca1", "Abcg1", 
    "Apoe", "Pltp", "Npc1", "Npc2", "Lipa", "Soat1", "Lrp1", "Osbpl1a", "Nr1h3", "Nr1h2"), Lipoprotein_processing = c("Apoe", 
    "Apoc1", "Apoc2", "Apoc4", "Lpl", "Lrp1", "Vldlr", "Ldlr", "Abca1", "Abcg1", "Pltp", "Lipa", "Soat1"), 
    Lysosomal_lipid_catabolism = c("Lipa", "Npc1", "Npc2", "Gba1", "Galc", "Hexa", "Hexb", "Smpd1", "Asah1", 
        "Ctsb", "Ctsd", "Lamp1", "Lamp2", "Lgmn"))

lipid_gene_sets <- lapply(lipid_gene_sets, function(g) {
    setdiff(unique(g), "Clec7a")
})

core_lipid_genes <- unique(unlist(lipid_gene_sets, use.names = FALSE))

core_lipid_genes <- setdiff(core_lipid_genes, "Clec7a")

if (!file.exists(SPATIAL_OBJECT)) {
    stop("Step64 spatial object not found:\n", SPATIAL_OBJECT)
}

if (!file.exists(RCTD_FILE)) {
    stop("Step65 v06 normalized spot-weight file not found:\n", RCTD_FILE)
}

obj <- readRDS(SPATIAL_OBJECT)

md <- obj@meta.data

expr <- get_spatial_data_matrix(obj)

if (!identical(colnames(expr), rownames(md))) {
    if (!all(colnames(expr) %in% rownames(md))) {
        stop("Spatial expression columns cannot be aligned to Step64 metadata rownames.")
    }
    md <- md[colnames(expr), , drop = FALSE]
}

MD_SAMPLE_COL <- pick_first_existing_col(md, c("sample_label", "sample"), "Step64 sample_label")

MD_GSM_COL <- pick_first_existing_col(md, c("GSM", "gsm"), "Step64 GSM")

if ("barcode" %in% colnames(md)) {
    md$barcode_step73 <- as.character(md$barcode)
} else {
    md$barcode_step73 <- rownames(md)
}

md$sample_label_step73 <- as.character(md[[MD_SAMPLE_COL]])

md$GSM_step73 <- as.character(md[[MD_GSM_COL]])

if (!all(c("array_x", "array_y") %in% colnames(md))) {
    stop("Step64 metadata must contain array_x and array_y for formal spatial joining.")
}

md$array_x_step73 <- as.numeric(md$array_x)

md$array_y_step73 <- as.numeric(md$array_y)

if (any(!is.finite(md$array_x_step73)) || any(!is.finite(md$array_y_step73))) {
    stop("Step64 metadata contains non-finite array coordinates.")
}

gene_match_tables <- list()

for (nm in names(lipid_gene_sets)) {
    gene_match_tables[[nm]] <- match_gene_set(lipid_gene_sets[[nm]], rownames(expr)) %>% mutate(module = nm)
}

core_match <- match_gene_set(core_lipid_genes, rownames(expr)) %>% mutate(module = "Core_lipid_handling_union")

gene_match_table <- do.call(dplyr::bind_rows, c(unname(gene_match_tables), list(core_match))) %>% select(module, 
    requested_gene, matched_gene, detected)

GENE_MATCH_FILE <- file.path(OUT_DIR, paste0("01_73_lipid_program_gene_matching_", VERSION, ".csv"))

write_csv_utf8(gene_match_table, GENE_MATCH_FILE)

module_detection_summary <- gene_match_table %>% group_by(module) %>% summarise(n_requested = n(), n_detected = sum(detected), 
    detected_fraction = mean(detected), .groups = "drop")

print(module_detection_summary)

bad_module <- module_detection_summary %>% filter(module != "Core_lipid_handling_union", n_detected < 
    5L | detected_fraction < 0.5)

if (nrow(bad_module) > 0L) {
    stop("Insufficient spatial gene detection for at least one lipid module:\n", paste(paste0(bad_module$module, 
        " (", bad_module$n_detected, "/", bad_module$n_requested, ")"), collapse = "\n"))
}

score_df <- data.frame(object_cell = colnames(expr), barcode = md$barcode_step73, sample_label = md$sample_label_step73, 
    GSM = md$GSM_step73, array_x = md$array_x_step73, array_y = md$array_y_step73, stringsAsFactors = FALSE)

if (anyDuplicated(score_df[, c("GSM", "array_x", "array_y")]) > 0L) {
    stop("Step64 score table has duplicated GSM + array_x + array_y keys.")
}

score_names <- c(names(lipid_gene_sets), "Core_lipid_handling")

for (sn in score_names) {
    score_df[[paste0("score_", sn)]] <- NA_real_
}

sample_values <- unique(score_df$sample_label)

sample_values <- SAMPLE_ORDER[SAMPLE_ORDER %in% sample_values]

if (length(sample_values) != EXPECTED_CAPTURE_AREAS) {
    stop("Expected ", EXPECTED_CAPTURE_AREAS, " capture areas in Step64 object, found ", length(sample_values), 
        ".")
}

for (s in sample_values) {

  spot_idx <- which(
    score_df$sample_label ==
      s
  )

  if (length(spot_idx) < 50L) {
    stop(
      "Too few spots in capture area: ",
      s
    )
  }

  module_loop <- c(
    lipid_gene_sets,
    list(
      Core_lipid_handling =
        core_lipid_genes
    )
  )

  for (nm in names(module_loop)) {

    match_tab <- match_gene_set(
      module_loop[[nm]],
      rownames(expr)
    )

    genes_use <- unique(
      match_tab$matched_gene[
        match_tab$detected
      ]
    )

    if (length(genes_use) < 5L) {

      stop(
        "Fewer than 5 detected genes for ",
        nm,
        " in spatial expression matrix."
      )
    }

    m <- as.matrix(
      expr[
        genes_use,
        spot_idx,
        drop = FALSE
      ]
    )

    gene_mean <- rowMeans(
      m,
      na.rm = TRUE
    )

    gene_sd <- apply(
      m,
      1,
      sd,
      na.rm = TRUE
    )

    keep_gene <- is.finite(gene_sd) &
      gene_sd >
      1e-8

    if (sum(keep_gene) < 5L) {

      stop(
        "Fewer than 5 variable genes remain for ",
        nm,
        " in ",
        s,
        "."
      )
    }

    m <- m[
      keep_gene,
      ,
      drop = FALSE
    ]

    gene_mean <- gene_mean[
      keep_gene
    ]

    gene_sd <- gene_sd[
      keep_gene
    ]

    zmat <- sweep(
      m,
      1,
      gene_mean,
      "-"
    )

    zmat <- sweep(
      zmat,
      1,
      gene_sd,
      "/"
    )

    spot_score <- colMeans(
      zmat,
      na.rm = TRUE
    )

    score_df[
      spot_idx,
      paste0(
        "score_",
        nm
      )
    ] <- spot_score
  }
}

rctd0 <- read_csv_gz_base(RCTD_FILE)

RCTD_BARCODE_COL <- pick_first_existing_col(rctd0, c("barcode", "cell", "spot"), "Step65 barcode")

RCTD_SAMPLE_COL <- pick_first_existing_col(rctd0, c("sample_label", "sample"), "Step65 sample_label")

RCTD_GSM_COL <- pick_first_existing_col(rctd0, c("GSM", "gsm"), "Step65 GSM")

RCTD_TIME_COL <- pick_first_existing_col(rctd0, c("timepoint", "time"), "Step65 timepoint")

RCTD_X_COL <- pick_first_existing_col(rctd0, c("array_x", "x"), "Step65 array_x")

RCTD_Y_COL <- pick_first_existing_col(rctd0, c("array_y", "y"), "Step65 array_y")

RCTD_SECUID_COL <- pick_first_existing_col(rctd0, c("section_uid_internal", "section_uid"), "Step65 section_uid_internal")

if (!(TARGET_STATE %in% colnames(rctd0))) {
    target_hits <- grep("Clec7a.*lipid.*Microglia|Clec7a.*lipid", colnames(rctd0), ignore.case = TRUE, 
        value = TRUE)
    if (length(target_hits) != 1L) {
        stop("Could not uniquely identify Clec7a target-state column.\n", "Hits: ", paste(target_hits, 
            collapse = " | "))
    }
    TARGET_COL <- target_hits[[1]]
} else {
    TARGET_COL <- TARGET_STATE
}

SECTION_ELIGIBLE_COL <- if ("section_eligible" %in% colnames(rctd0)) {
    "section_eligible"
} else {
    NA_character_
}

rctd <- data.frame(barcode = as.character(rctd0[[RCTD_BARCODE_COL]]), sample_label = as.character(rctd0[[RCTD_SAMPLE_COL]]), 
    GSM = as.character(rctd0[[RCTD_GSM_COL]]), timepoint = as.character(rctd0[[RCTD_TIME_COL]]), array_x = as.numeric(rctd0[[RCTD_X_COL]]), 
    array_y = as.numeric(rctd0[[RCTD_Y_COL]]), section_uid_internal = as.character(rctd0[[RCTD_SECUID_COL]]), 
    Clec7a_RCTDnorm = as.numeric(rctd0[[TARGET_COL]]), stringsAsFactors = FALSE)

if (!is.na(SECTION_ELIGIBLE_COL)) {
    rctd$section_eligible <- as.logical(rctd0[[SECTION_ELIGIBLE_COL]])
} else {
    rctd$section_eligible <- !is.na(rctd$section_uid_internal) & nzchar(rctd$section_uid_internal)
}

score_key <- score_df %>% select(GSM, array_x, array_y, starts_with("score_"))

if (anyDuplicated(score_key[, c("GSM", "array_x", "array_y")]) > 0L) {
    stop("Step64 lipid-score table has duplicated GSM + array_x + array_y keys.")
}

if (anyDuplicated(rctd[, c("GSM", "array_x", "array_y")]) > 0L) {
    stop("Step65 v06 RCTD table has duplicated GSM + array_x + array_y keys.")
}

joined <- rctd %>% left_join(score_key, by = c("GSM", "array_x", "array_y"))

if (nrow(joined) != nrow(rctd)) {
    stop("Unexpected row-count change after formal GSM + array_x + array_y join.")
}

score_cols <- grep("^score_", colnames(joined), value = TRUE)

join_complete_fraction <- mean(complete.cases(joined[, score_cols, drop = FALSE]))

if (join_complete_fraction <
    0.95) {

  join_diag <- joined %>%
    mutate(
      lipid_score_joined =
        complete.cases(
          joined[
            ,
            score_cols,
            drop = FALSE
          ]
        )
    ) %>%
    count(
      sample_label,
      GSM,
      lipid_score_joined,
      name = "n_spots"
    )

  write_csv_utf8(
    join_diag,
    file.path(
      OUT_DIR,
      paste0(
        "00_73_coordinate_join_diagnostic_",
        VERSION,
        ".csv"
      )
    )
  )

  stop(
    "Lipid-program coordinate join completeness <95%: ",
    sprintf(
      "%.2f%%",
      100 *
        join_complete_fraction
    ),
    "\nDiagnostic table exported to Step73 output directory."
  )
}

dat <- joined %>% filter(section_eligible, timepoint %in% TIME_LEVELS, is.finite(Clec7a_RCTDnorm), is.finite(score_Core_lipid_handling))

n_sections <- length(unique(dat$section_uid_internal))

n_capture_areas <- length(unique(dat$sample_label))

if (n_sections != EXPECTED_MAJOR_SECTIONS) {
    stop("Expected ", EXPECTED_MAJOR_SECTIONS, " major coordinate-defined sections; found ", n_sections, 
        ".")
}

if (n_capture_areas != EXPECTED_CAPTURE_AREAS) {
    stop("Expected ", EXPECTED_CAPTURE_AREAS, " capture areas; found ", n_capture_areas, ".")
}

dat$lipid_neighbor_mean <- NA_real_

section_order <- dat %>% distinct(section_uid_internal, sample_label, GSM, timepoint) %>% arrange(factor(timepoint, 
    levels = TIME_LEVELS), factor(sample_label, levels = SAMPLE_ORDER), section_uid_internal) %>% pull(section_uid_internal)

section_rows <- list()
graph_rows <- list()
counter <- 1L

for (sec in section_order) {

  d <- dat %>%
    filter(
      section_uid_internal ==
        sec
    )

  if (nrow(d) <
      MIN_SECTION_SPOTS) {

    next
  }

  coords <- as.matrix(
    d[
      ,
      c(
        "array_x",
        "array_y"
      )
    ]
  )

  rownames(coords) <- as.character(
    seq_len(
      nrow(d)
    )
  )

  g <- build_visium_graph(
    coords
  )

  nb <- g$nb

  degree <- spdep::card(
    nb
  )

  lipid_neighbor <- neighbor_mean_from_nb(
    d$score_Core_lipid_handling,
    nb
  )

  rho_direct <- safe_spearman(
    d$Clec7a_RCTDnorm,
    d$score_Core_lipid_handling
  )

  rho_neighbor <- safe_spearman(
    d$Clec7a_RCTDnorm,
    lipid_neighbor
  )

  rho_FA <- safe_spearman(
    d$Clec7a_RCTDnorm,
    d$score_FA_uptake_transport
  )

  rho_chol <- safe_spearman(
    d$Clec7a_RCTDnorm,
    d$score_Cholesterol_efflux_transport
  )

  rho_lipo <- safe_spearman(
    d$Clec7a_RCTDnorm,
    d$score_Lipoprotein_processing
  )

  rho_lyso <- safe_spearman(
    d$Clec7a_RCTDnorm,
    d$score_Lysosomal_lipid_catabolism
  )

  section_rows[[counter]] <- data.frame(
    section_uid_internal =
      sec,
    sample_label =
      as.character(
        d$sample_label[[1]]
      ),
    GSM =
      as.character(
        d$GSM[[1]]
      ),
    timepoint =
      as.character(
        d$timepoint[[1]]
      ),
    n_spots =
      nrow(d),
    direct_spearman_rho =
      rho_direct,
    local_neighbor_spearman_rho =
      rho_neighbor,
    FA_uptake_transport_rho =
      rho_FA,
    Cholesterol_efflux_transport_rho =
      rho_chol,
    Lipoprotein_processing_rho =
      rho_lipo,
    Lysosomal_lipid_catabolism_rho =
      rho_lyso,
    direct_fisher_z =
      fisher_z(
        rho_direct
      ),
    local_neighbor_fisher_z =
      fisher_z(
        rho_neighbor
      ),
    stringsAsFactors = FALSE
  )

  graph_rows[[counter]] <- data.frame(
    section_uid_internal =
      sec,
    sample_label =
      as.character(
        d$sample_label[[1]]
      ),
    timepoint =
      as.character(
        d$timepoint[[1]]
      ),
    n_spots =
      nrow(d),
    median_nearest_neighbor_distance =
      g$median_nn_distance,
    graph_radius =
      g$radius,
    mean_degree =
      mean(degree),
    median_degree =
      median(degree),
    min_degree =
      min(degree),
    max_degree =
      max(degree),
    n_isolated =
      sum(
        degree ==
          0L
      ),
    isolated_fraction =
      mean(
        degree ==
          0L
      ),
    stringsAsFactors = FALSE
  )

  idx_global <- match(
    paste(
      d$GSM,
      d$barcode,
      sep = "||"
    ),
    paste(
      dat$GSM,
      dat$barcode,
      sep = "||"
    )
  )

  dat$lipid_neighbor_mean[
    idx_global
  ] <- lipid_neighbor

  counter <- counter + 1L
}

section_results <- bind_rows(section_rows)

graph_qc <- bind_rows(graph_rows)

if (nrow(section_results) != EXPECTED_MAJOR_SECTIONS) {
    stop("Section coupling result count is ", nrow(section_results), "; expected ", EXPECTED_MAJOR_SECTIONS, 
        ".")
}

time_summary <- section_results %>% group_by(timepoint) %>% summarise(n_sections = n(), n_capture_areas = n_distinct(sample_label), 
    median_direct_rho = median(direct_spearman_rho, na.rm = TRUE), q1_direct_rho = quantile(direct_spearman_rho, 
        0.25, na.rm = TRUE), q3_direct_rho = quantile(direct_spearman_rho, 0.75, na.rm = TRUE), median_local_neighbor_rho = median(local_neighbor_spearman_rho, 
        na.rm = TRUE), q1_local_neighbor_rho = quantile(local_neighbor_spearman_rho, 0.25, na.rm = TRUE), 
    q3_local_neighbor_rho = quantile(local_neighbor_spearman_rho, 0.75, na.rm = TRUE), .groups = "drop") %>% 
    arrange(match(timepoint, TIME_LEVELS))

blocked <- section_results %>% group_by(sample_label, GSM, timepoint) %>% summarise(n_sections_contributing_direct = sum(is.finite(direct_fisher_z)), 
    blocked_direct_fisher_z = median(direct_fisher_z[is.finite(direct_fisher_z)], na.rm = TRUE), blocked_direct_rho = inv_fisher_z(blocked_direct_fisher_z), 
    n_sections_contributing_local = sum(is.finite(local_neighbor_fisher_z)), blocked_local_fisher_z = median(local_neighbor_fisher_z[is.finite(local_neighbor_fisher_z)], 
        na.rm = TRUE), blocked_local_rho = inv_fisher_z(blocked_local_fisher_z), .groups = "drop")

formal_summary <- data.frame(analysis = FORMAL_ANALYSES, n_capture_areas = c(sum(is.finite(blocked$blocked_direct_fisher_z)), 
    sum(is.finite(blocked$blocked_local_fisher_z))), median_blocked_fisher_z = c(median(blocked$blocked_direct_fisher_z, 
    na.rm = TRUE), median(blocked$blocked_local_fisher_z, na.rm = TRUE)), median_blocked_rho = c(inv_fisher_z(median(blocked$blocked_direct_fisher_z, 
    na.rm = TRUE)), inv_fisher_z(median(blocked$blocked_local_fisher_z, na.rm = TRUE))), fraction_capture_areas_positive = c(mean(blocked$blocked_direct_fisher_z > 
    0, na.rm = TRUE), mean(blocked$blocked_local_fisher_z > 0, na.rm = TRUE)), one_sided_wilcox_p = c(wilcox_greater_zero(blocked$blocked_direct_fisher_z), 
    wilcox_greater_zero(blocked$blocked_local_fisher_z)), stringsAsFactors = FALSE)

formal_summary$FDR_BH <- p.adjust(formal_summary$one_sided_wilcox_p, method = "BH")

component_long <- section_results %>% select(section_uid_internal, sample_label, timepoint, FA_uptake_transport_rho, 
    Cholesterol_efflux_transport_rho, Lipoprotein_processing_rho, Lysosomal_lipid_catabolism_rho) %>% 
    pivot_longer(cols = ends_with("_rho"), names_to = "module", values_to = "spearman_rho") %>% mutate(module = gsub("_rho$", 
    "", module))

component_summary <- component_long %>% group_by(module) %>% summarise(n_sections = sum(is.finite(spearman_rho)), 
    median_section_rho = median(spearman_rho, na.rm = TRUE), fraction_sections_positive = mean(spearman_rho > 
        0, na.rm = TRUE), .groups = "drop") %>% arrange(desc(median_section_rho))

SPOT_OUTPUT_FILE <- file.path(OUT_DIR, paste0("02_73_Clec7a_lipid_spot_scores_section_aware_", VERSION, 
    ".csv.gz"))

write_csv_gz_base(dat, SPOT_OUTPUT_FILE)

SECTION_FILE <- file.path(OUT_DIR, paste0("03_73_Clec7a_lipid_section_coupling_", VERSION, ".csv"))

write_csv_utf8(section_results, SECTION_FILE)

TIME_FILE <- file.path(OUT_DIR, paste0("04_73_Clec7a_lipid_timepoint_descriptive_", VERSION, ".csv"))

write_csv_utf8(time_summary, TIME_FILE)

BLOCKED_FILE <- file.path(OUT_DIR, paste0("05_73_Clec7a_lipid_capture_area_blocked_effects_", VERSION, 
    ".csv"))

write_csv_utf8(blocked, BLOCKED_FILE)

FORMAL_FILE <- file.path(OUT_DIR, paste0("06_73_Clec7a_lipid_formal_blocked_tests_", VERSION, ".csv"))

write_csv_utf8(formal_summary, FORMAL_FILE)

COMPONENT_FILE <- file.path(OUT_DIR, paste0("07_73_lipid_component_coupling_summary_", VERSION, ".csv"))

write_csv_utf8(component_summary, COMPONENT_FILE)

GRAPH_FILE <- file.path(OUT_DIR, paste0("08_73_section_spatial_graph_QC_", VERSION, ".csv"))

write_csv_utf8(graph_qc, GRAPH_FILE)

plot_sec <- section_results %>% mutate(timepoint = factor(timepoint, levels = TIME_LEVELS, labels = c("3 dpi", 
    "7 dpi", "14 dpi", "28 dpi")))

p1 <- ggplot(plot_sec, aes(x = timepoint, y = direct_spearman_rho)) + geom_hline(yintercept = 0, linetype = "dashed", 
    linewidth = 0.45, color = "grey65") + geom_jitter(width = 0.13, height = 0, size = 1.8, alpha = 0.7, 
    color = "#8B4A78") + stat_summary(fun = median, geom = "point", shape = 21, size = 3.4, stroke = 0.8, 
    fill = "white", color = "#5D214F") + labs(title = "Same-spot coupling", x = NULL, y = "Spearman rho") + 
    theme_classic(base_family = "Arial", base_size = 9) + theme(plot.title = element_text(size = 10, 
    face = "bold"), axis.text.x = element_text(size = 8.5))

p2 <- ggplot(plot_sec, aes(x = timepoint, y = local_neighbor_spearman_rho)) + geom_hline(yintercept = 0, 
    linetype = "dashed", linewidth = 0.45, color = "grey65") + geom_jitter(width = 0.13, height = 0, 
    size = 1.8, alpha = 0.7, color = "#C56B45") + stat_summary(fun = median, geom = "point", shape = 21, 
    size = 3.4, stroke = 0.8, fill = "white", color = "#8D3D20") + labs(title = "Local spatial coupling", 
    x = NULL, y = "Spearman rho") + theme_classic(base_family = "Arial", base_size = 9) + theme(plot.title = element_text(size = 10, 
    face = "bold"), axis.text.x = element_text(size = 8.5))

p3 <- ggplot(component_summary, aes(x = reorder(module, median_section_rho), y = median_section_rho)) + 
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, color = "grey65") + geom_segment(aes(xend = module, 
    y = 0, yend = median_section_rho), linewidth = 1, color = "#D7A05C") + geom_point(size = 3, color = "#B95B3C") + 
    coord_flip() + labs(title = "Lipid-program components", x = NULL, y = "Median section rho") + theme_classic(base_family = "Arial", 
    base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 7.5))

diagnostic_plot <- (p1 | p2)/p3 + patchwork::plot_layout(heights = c(1, 0.85))

DIAG_PDF <- file.path(OUT_DIR, paste0("09_73_Clec7a_lipid_spatial_coupling_diagnostic_", VERSION, ".pdf"))

DIAG_PNG <- file.path(OUT_DIR, paste0("09_73_Clec7a_lipid_spatial_coupling_diagnostic_", VERSION, ".png"))

ggsave(DIAG_PDF, diagnostic_plot, width = 8.2, height = 6.2, units = "in", device = grDevices::cairo_pdf, 
    bg = "white")

ggsave(DIAG_PNG, diagnostic_plot, width = 8.2, height = 6.2, units = "in", dpi = 500, bg = "white")

all_section_size_ok <- all(section_results$n_spots >= MIN_SECTION_SPOTS)

graph_isolated_ok <- all(graph_qc$isolated_fraction < 0.05)

all_core_scores_finite <- all(is.finite(dat$score_Core_lipid_handling))

target_finite <- all(is.finite(dat$Clec7a_RCTDnorm))

clec7a_not_in_program <- !("Clec7a" %in% core_lipid_genes)

final_checks <- data.frame(check = c("Step64 v03 spatial object was used", "Step65 v06 section-aware normalized RCTD file was used", 
    "Spatial log-normalized assay data were available", "Clec7a was excluded from the lipid-handling program", 
    "All four curated lipid-handling component modules passed gene-detection QC", "Lipid program scores joined to at least 95% of Step65 spots by GSM + array_x + array_y", 
    "All 8 Visium capture areas were retained", "All 32 coordinate-defined major tissue sections were analyzed", 
    "All analyzed sections contain at least 100 spots", "Clec7a normalized RCTD weights are finite", 
    "Core lipid-handling scores are finite", "Spatial graph radius matches Step66 rule: 1.55 x section median nearest-neighbor distance", 
    "Graph isolated-spot fraction is <5% in every section", "No spot-level P values are used for manuscript inference", 
    "Formal inference uses capture-area-blocked effects", "No timepoint-specific hypothesis testing was performed", 
    "Diagnostic PDF generated", "Diagnostic PNG generated"), pass = c(identical(basename(SPATIAL_OBJECT), 
    "12_64_GSE195783_spatial_region_annotated_v03.rds"), identical(basename(RCTD_FILE), "03_65_RCTD_normalized_spot_weights_section_aware_v06.csv.gz"), 
    nrow(expr) > 0L && ncol(expr) > 0L, clec7a_not_in_program, nrow(bad_module) == 0L, join_complete_fraction >= 
        0.95, n_capture_areas == EXPECTED_CAPTURE_AREAS, n_sections == EXPECTED_MAJOR_SECTIONS, all_section_size_ok, 
    target_finite, all_core_scores_finite, isTRUE(all.equal(GRAPH_RADIUS_MULTIPLIER, 1.55)), graph_isolated_ok, 
    TRUE, nrow(formal_summary) == length(FORMAL_ANALYSES), TRUE, file.exists(DIAG_PDF), file.exists(DIAG_PNG)), 
    stringsAsFactors = FALSE)

QC_FILE <- file.path(OUT_DIR, paste0("10_73_GSE195783_Clec7a_lipid_coupling_FINAL_QC_", VERSION, ".csv"))

write_csv_utf8(final_checks, QC_FILE)

REPORT_FILE <- file.path(OUT_DIR, paste0("11_73_GSE195783_Clec7a_lipid_coupling_REPORT_", VERSION, ".txt"))

time_lines <- apply(time_summary, 1, function(z) {
    paste0("  ", z[["timepoint"]], ": n_sections=", z[["n_sections"]], "; median direct rho=", sprintf("%.3f", 
        as.numeric(z[["median_direct_rho"]])), "; median local-neighbor rho=", sprintf("%.3f", as.numeric(z[["median_local_neighbor_rho"]])))
})

formal_lines <- apply(formal_summary, 1, function(z) {
    paste0("  ", z[["analysis"]], ": n_capture_areas=", z[["n_capture_areas"]], "; median blocked rho=", 
        sprintf("%.3f", as.numeric(z[["median_blocked_rho"]])), "; fraction positive=", sprintf("%.3f", 
            as.numeric(z[["fraction_capture_areas_positive"]])), "; one-sided Wilcoxon P=", ifelse(is.finite(as.numeric(z[["one_sided_wilcox_p"]])), 
            sprintf("%.4g", as.numeric(z[["one_sided_wilcox_p"]])), "NA"), "; BH FDR=", ifelse(is.finite(as.numeric(z[["FDR_BH"]])), 
            sprintf("%.4g", as.numeric(z[["FDR_BH"]])), "NA"))
})

component_lines <- apply(component_summary, 1, function(z) {
    paste0("  ", z[["module"]], ": median section rho=", sprintf("%.3f", as.numeric(z[["median_section_rho"]])), 
        "; fraction positive=", sprintf("%.3f", as.numeric(z[["fraction_sections_positive"]])))
})

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    con = REPORT_FILE, useBytes = TRUE)

print(time_summary)

print(formal_summary)

print(component_summary)

print(final_checks)

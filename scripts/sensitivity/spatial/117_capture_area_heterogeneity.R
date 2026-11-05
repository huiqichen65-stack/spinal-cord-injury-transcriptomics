source("config.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1, scipen = 5)

ROOT <- PROJECT_ROOT
STEP <- 117L
VERSION <- "v02"
JOB <- "capture_area_heterogeneity"
OUT <- file.path(ROOT, "results", paste0(STEP, "_", JOB))
DIR_TABLES <- file.path(OUT, "tables")
DIR_FIGURES <- file.path(OUT, "figures")

INPUT <- list(qc_sample = file.path(ROOT, "results", "61_GSE195783_spatial_QC_normalization", "03_61_GSE195783_QC_retention_summary_v01.csv"),
  qc_spot = file.path(ROOT, "results", "61_GSE195783_spatial_QC_normalization", "02_61_GSE195783_spot_QC_decisions_v01.csv.gz"),
  pca_structure = file.path(ROOT, "results", "62_GSE195783_joint_embedding_batch_assessment_clustering",
    "02_62_PC_timepoint_vs_sample_structure_v01.csv"), harmony_structure = file.path(ROOT, "results",
    "63_GSE195783_Harmony_sensitivity_analysis", "01_63_global_sensitivity_metrics_v01.csv"), clec_section = file.path(ROOT,
    "results", "65_GSE195783_RCTD_myeloid_state_mapping", "04_65_Clec7a_section_summary_v06.csv"),
  section_manifest = file.path(ROOT, "results", "68_GSE195783_spatial_section_map_finalize", "01_68_GSE195783_spatial_section_map_section_manifest_v01.csv"),
  spot_section = file.path(ROOT, "results", "68_GSE195783_spatial_section_map_finalize", "02_68_GSE195783_spatial_section_map_spot_to_section_v01.csv.gz"),
  rep_sections = file.path(ROOT, "results", "69_Figure4A_spatial_transcriptional_regions", "02_69_Figure4A_representative_sections_v12.csv"),
  figure4B_spots = file.path(ROOT, "results", "70_Figure4B_Clec7a_spatial_abundance", "05_70_Figure4B_plot_spot_table_v08.csv.gz"),
  section_regions = file.path(ROOT, "results", "71_Figure4C_Clec7a_spatial_niche_heatmap", "01_71_Figure4C_section_region_spatial_burden_v04.csv"),
  signature_114 = file.path(ROOT, "results", "114_GSE195783_uninjured_section_validation", "03_114_section_signature_summary_v01.csv"))

MANDATORY <- c("qc_sample", "qc_spot", "clec_section", "section_manifest", "spot_section", "rep_sections",
  "section_regions", "signature_114")
STAGES <- c("3dpi", "7dpi", "14dpi", "28dpi")
SCENARIOS <- c("GEO_label_all32", "injury_only_inferred27")
EXPECTED_AREAS <- c("3dpi_rep1", "3dpi_rep2", "7dpi_rep1", "7dpi_rep2", "14dpi_rep1", "14dpi_rep2", "28dpi_rep1",
  "28dpi_rep2")
INFERRED_CONTROL_UIDS <- c("7dpi_rep2__upper_right", "14dpi_rep1__upper_right", "14dpi_rep2__lower_right",
  "28dpi_rep1__upper_right", "28dpi_rep2__lower_right")
EXPECTED_CORRECTED <- c(`3dpi` = 8L, `7dpi` = 7L, `14dpi` = 6L, `28dpi` = 6L)

stop_if_missing_cols <- function(x, cols, name) {
  absent <- setdiff(cols, names(x))
  if (length(absent))
    stop(name, " missing required columns: ", paste(absent, collapse = ", "))
}

read_data <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(connection), add = TRUE)
  read.csv(connection, check.names = FALSE, stringsAsFactors = FALSE)
}

save_csv <- function(x, prefix, name) {
  path <- file.path(DIR_TABLES, paste0(sprintf("%02d", prefix), "_", STEP, "_", name, "_", VERSION,
    ".csv"))
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  path
}

mean_ok <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x))
    return(NA_real_)
  mean(x)
}

median_ok <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x))
    return(NA_real_)
  median(x)
}

sd_ok <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L)
    return(NA_real_)
  stats::sd(x)
}

range_low <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x))
    return(NA_real_)
  min(x)
}

range_high <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x))
    return(NA_real_)
  max(x)
}

weighted_ok <- function(x, w) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep))
    return(NA_real_)
  sum(x[keep] * w[keep])/sum(w[keep])
}

safe_ratio <- function(a, b) {
  if (!is.finite(a) || !is.finite(b) || abs(b) < .Machine$double.eps)
    return(NA_real_)
  a/b
}

jsd <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (!length(a) || length(a) != length(b) || any(!is.finite(a)) || any(!is.finite(b)) || any(a < 0) ||
    any(b < 0) || sum(a) <= 0 || sum(b) <= 0)
    return(NA_real_)
  a <- a/sum(a)
  b <- b/sum(b)
  m <- (a + b)/2
  kl <- function(p, q) {
    use <- p > 0
    sum(p[use] * log2(p[use]/q[use]))
  }
  sqrt(max(0, 0.5 * kl(a, m) + 0.5 * kl(b, m)))
}

tv_distance <- function(a, b) {
  if (any(!is.finite(c(a, b))))
    return(NA_real_)
  a <- a/sum(a)
  b <- b/sum(b)
  0.5 * sum(abs(a - b))
}

if (!dir.exists(ROOT)) stop("Project root not found: ", ROOT)
missing_input <- MANDATORY[!vapply(INPUT[MANDATORY], file.exists, logical(1))]
if (length(missing_input)) {
  stop("Missing required CSV/CSV.GZ files (no ZIP fallback):\n", paste(sprintf("  %s: %s", missing_input,
    unlist(INPUT[missing_input])), collapse = "\n"))
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("R package ggplot2 is required. Install using install.packages('ggplot2').")
}
suppressPackageStartupMessages(library(ggplot2))
dir.create(DIR_TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

manifest <- read_data(INPUT$section_manifest)
spot_map <- read_data(INPUT$spot_section)
raw_regions <- read_data(INPUT$section_regions)
clec <- read_data(INPUT$clec_section)
s114 <- read_data(INPUT$signature_114)
qc_sample <- read_data(INPUT$qc_sample)
qc_spot <- read_data(INPUT$qc_spot)
representatives <- read_data(INPUT$rep_sections)

input_audit <- data.frame(input_name = names(INPUT), source_path = unlist(INPUT, use.names = FALSE),
  mandatory = names(INPUT) %in% MANDATORY, exists = vapply(INPUT, file.exists, logical(1)), bytes = vapply(INPUT,
    function(path) if (file.exists(path)) as.numeric(file.info(path)$size) else NA_real_, numeric(1)),
  input_type = "CSV or CSV.GZ (already generated local files)", stringsAsFactors = FALSE)
save_csv(input_audit, 0, "input_file_inventory")

figure4B_matching <- NA
if (file.exists(INPUT$figure4B_spots)) {
  fig4B <- read_data(INPUT$figure4B_spots)
  stop_if_missing_cols(fig4B, c("timepoint", "section_uid_internal"), "Step70 Figure4B spot table")
  fig4B_set <- unique(paste(fig4B$timepoint, fig4B$section_uid_internal, sep = "||"))
  fig4A_set <- unique(paste(representatives$timepoint, representatives$section_uid_internal, sep = "||"))
  figure4B_matching <- setequal(fig4B_set, fig4A_set)
  if (!isTRUE(figure4B_matching)) {
    stop("Figure4A and Figure4B use different section IDs. Inspect Step69/70 before continuing.")
  }
}

stop_if_missing_cols(manifest, c("sample_label", "GSM", "timepoint", "section_uid_internal", "section_position_internal",
  "n_spots"), "Step68 section manifest")
stop_if_missing_cols(spot_map, c("GSM", "barcode", "section_uid_internal", "is_major_component"), "Step68 spot map")
stop_if_missing_cols(raw_regions, c("GSM", "sample_label", "timepoint", "section_uid_internal", "spatial_region_v03",
  "n_spots", "total_section_spots", "region_target_weight_sum", "spatial_niche_burden"), "Step71 section-region")
stop_if_missing_cols(clec, c("summary_level", "section_uid_internal", "n_spots", "mean_target_weight",
  "fraction_weight_gt_0_10"), "Step65 section target")
stop_if_missing_cols(s114, c("section_uid_internal", "median_UMI_per_spot", "median_detected_features_per_spot",
  "normal_minus_injury_median_z", "injury_combined_median_z", "mean_RCTD_normal_weight", "mean_RCTD_injury_myeloid_weight"),
  "Step114 summary")
stop_if_missing_cols(qc_sample, c("sample_label", "n_spots_before", "n_pass", "retention_rate"), "Step61 QC summary")
stop_if_missing_cols(qc_spot, c("GSM", "barcode", "nCount_Spatial", "nFeature_Spatial", "percent_mt",
  "pass_qc"), "Step61 spot QC")
stop_if_missing_cols(representatives, c("timepoint", "sample_label", "section_uid_internal"), "Step69 representative manifest")

if (nrow(manifest) != 32L || anyDuplicated(manifest$section_uid_internal)) stop("Step68 manifest must have 32 unique internal section identifiers.")
if (!setequal(unique(manifest$sample_label), EXPECTED_AREAS)) stop("Step68 manifest does not contain the expected eight capture areas.")
if (!all(table(manifest$sample_label)[EXPECTED_AREAS] == 4L)) stop("Expected four major coordinate components per capture area.")
if (!all(INFERRED_CONTROL_UIDS %in% manifest$section_uid_internal)) stop("One or more inferred-uninjured section IDs are absent; do not guess replacements.")
if (anyDuplicated(s114$section_uid_internal)) stop("Step114 section summary has duplicate piece IDs.")
if (!setequal(s114$section_uid_internal, manifest$section_uid_internal)) stop("Step114 and Step68 section-ID sets differ.")
if (anyDuplicated(paste(raw_regions$section_uid_internal, raw_regions$spatial_region_v03, sep = "|"))) stop("Step71 section-by-region table contains duplicated section/region pairs.")

region_names <- unique(as.character(raw_regions$spatial_region_v03))
if (length(region_names) != 10L) stop("Expected ten Step71 transcriptional regions; found ", length(region_names))
TARGET_REGION <- grep("^Myeloid/immune", region_names, value = TRUE)
INJURY_REGIONS <- region_names[grepl("^Myeloid/immune|^Fibrotic/ECM|^Acute proliferative|^IFN-responsive",
  region_names)]
if (length(TARGET_REGION) != 1L || length(INJURY_REGIONS) != 4L) stop("Cannot identify target region and four injury-related region classes.")

manifest <- manifest[order(match(manifest$sample_label, EXPECTED_AREAS), match(manifest$section_position_internal,
  c("upper_left", "upper_right", "lower_left", "lower_right"))), ]
manifest$GEO_stage_original <- manifest$timepoint
manifest$inferred_uninjured <- manifest$section_uid_internal %in% INFERRED_CONTROL_UIDS
manifest$inferred_tissue_group <- ifelse(manifest$inferred_uninjured, "putative_uninjured", "putative_injured")
manifest$stage_after_inference <- ifelse(manifest$inferred_uninjured, "Uninjured_inferred", manifest$GEO_stage_original)
manifest$animal_id <- NA_character_
manifest$tissue_block_id <- NA_character_
manifest$original_section_id <- NA_character_
manifest$label_evidence <- ifelse(manifest$inferred_uninjured, "Step114 expression phenotype + original-paper count inference; not original animal metadata",
  "GEO original stage label; inferred putative injured after excluding five candidates")

counts_clean <- table(factor(manifest$stage_after_inference[!manifest$inferred_uninjured], levels = STAGES))
if (!identical(as.integer(counts_clean), as.integer(EXPECTED_CORRECTED[STAGES]))) stop("Unexpected corrected section totals; expected 8/7/6/6.")

sec_clec <- clec[clec$summary_level == "section", c("section_uid_internal", "n_spots", "mean_target_weight",
  "median_target_weight", "fraction_weight_gt_0_10")]
if (nrow(sec_clec) != 32L || anyDuplicated(sec_clec$section_uid_internal)) stop("Step65 section-level Clec7a summary must have exactly 32 unique pieces.")
names(sec_clec)[names(sec_clec) == "n_spots"] <- "n_spots_RCTD"
names(sec_clec)[names(sec_clec) == "mean_target_weight"] <- "Clec7a_mean_RCTD"
names(sec_clec)[names(sec_clec) == "median_target_weight"] <- "Clec7a_median_RCTD"

sec_sig <- s114[, c("section_uid_internal", "n_spots", "median_UMI_per_spot", "median_detected_features_per_spot",
  "normal_minus_injury_median_z", "injury_combined_median_z", "mean_RCTD_normal_weight", "mean_RCTD_injury_myeloid_weight")]
names(sec_sig)[names(sec_sig) == "n_spots"] <- "n_spots_114"

section <- merge(manifest, sec_clec, by = "section_uid_internal", all.x = TRUE, sort = FALSE)
section <- merge(section, sec_sig, by = "section_uid_internal", all.x = TRUE, sort = FALSE)
section <- section[match(manifest$section_uid_internal, section$section_uid_internal), ]
if (anyNA(section$Clec7a_mean_RCTD)) stop("Missing Step65 Clec7a section means.")

spot_map <- spot_map[!is.na(spot_map$section_uid_internal) & (spot_map$is_major_component %in% TRUE),
  c("GSM", "barcode", "section_uid_internal")]
if (anyDuplicated(paste(spot_map$GSM, spot_map$barcode))) stop("Duplicate GSM+barcode keys in Step68 spot map.")
qc_spot <- qc_spot[, c("GSM", "barcode", "nCount_Spatial", "nFeature_Spatial", "percent_mt", "pass_qc")]
if (anyDuplicated(paste(qc_spot$GSM, qc_spot$barcode))) stop("Duplicate GSM+barcode keys in Step61 spot QC.")
spot_join <- merge(spot_map, qc_spot, by = c("GSM", "barcode"), all.x = TRUE, sort = FALSE)
if (anyNA(spot_join$nCount_Spatial)) stop("Some Step68 spots are missing from Step61 spot QC.")
spot_join$pass_qc <- spot_join$pass_qc %in% TRUE

section_qc <- do.call(rbind, lapply(unique(section$section_uid_internal), function(uid) {
  z <- spot_join[spot_join$section_uid_internal == uid, ]
  data.frame(section_uid_internal = uid, n_spots_before_QC = nrow(z), n_spots_after_QC = sum(z$pass_qc),
    QC_retention_fraction = mean(z$pass_qc), median_UMI_QCpass = median_ok(z$nCount_Spatial[z$pass_qc]),
    median_genes_QCpass = median_ok(z$nFeature_Spatial[z$pass_qc]), median_percent_mt_QCpass = median_ok(z$percent_mt[z$pass_qc]),
    stringsAsFactors = FALSE)
}))
section <- merge(section, section_qc, by = "section_uid_internal", all.x = TRUE, sort = FALSE)
section <- section[match(manifest$section_uid_internal, section$section_uid_internal), ]

sec_total <- aggregate(raw_regions$n_spots, by = list(section_uid_internal = raw_regions$section_uid_internal),
  sum)
names(sec_total)[2] <- "n_spots_region_sum"
section <- merge(section, sec_total, by = "section_uid_internal", all.x = TRUE, sort = FALSE)
section <- section[match(manifest$section_uid_internal, section$section_uid_internal), ]
if (any(abs(section$n_spots_region_sum - section$n_spots_RCTD) > 0)) stop("Step71 region spot counts do not sum to Step65 section RCTD spot counts.")

raw_regions$region_fraction <- ifelse(raw_regions$total_section_spots > 0, raw_regions$n_spots/raw_regions$total_section_spots,
  NA_real_)
raw_regions$region_class <- ifelse(raw_regions$spatial_region_v03 %in% INJURY_REGIONS, "injury_associated_transcriptional",
  "other_transcriptional")
raw_regions$inferred_uninjured <- raw_regions$section_uid_internal %in% INFERRED_CONTROL_UIDS
raw_regions$analysis_stage <- ifelse(raw_regions$inferred_uninjured, "Uninjured_inferred", raw_regions$timepoint)

for (uid in section$section_uid_internal) {
  r <- raw_regions[raw_regions$section_uid_internal == uid, ]
  if (nrow(r) != 10L || length(unique(r$spatial_region_v03)) != 10L)
    stop("Section does not have ten explicit region rows: ", uid)
  if (length(unique(r$total_section_spots)) != 1L || sum(r$n_spots) != r$total_section_spots[1] ||
    abs(sum(r$region_fraction) - 1) > 1e-09)
    stop("Step71 region counts/denominators inconsistent: ", uid)
}

section$injury_region_fraction <- vapply(section$section_uid_internal, function(uid) {
  z <- raw_regions[raw_regions$section_uid_internal == uid & raw_regions$spatial_region_v03 %in% INJURY_REGIONS,
    ]
  sum(z$region_fraction)
}, numeric(1))
section$target_region_fraction <- vapply(section$section_uid_internal, function(uid) {
  z <- raw_regions[raw_regions$section_uid_internal == uid & raw_regions$spatial_region_v03 == TARGET_REGION,
    ]
  sum(z$region_fraction)
}, numeric(1))
section$target_region_spatial_burden <- vapply(section$section_uid_internal, function(uid) {
  z <- raw_regions[raw_regions$section_uid_internal == uid & raw_regions$spatial_region_v03 == TARGET_REGION,
    ]
  sum(z$spatial_niche_burden)
}, numeric(1))

section$RNA_depth_proxy <- section$median_UMI_QCpass
section$detected_genes_proxy <- section$median_genes_QCpass
section$mitochondrial_proxy <- section$median_percent_mt_QCpass

save_csv(section[, c("section_uid_internal", "GSM", "sample_label", "GEO_stage_original", "stage_after_inference",
  "section_position_internal", "inferred_tissue_group", "inferred_uninjured", "animal_id", "tissue_block_id",
  "original_section_id", "label_evidence")], 1, "frozen_32piece_status_and_provenance")
save_csv(raw_regions, 2, "section_region_10class_composition_all32")
save_csv(section, 3, "section_indicators_QC_all32")

capture_rows <- list()
capture_region_rows <- list()
kk <- 0L
rr <- 0L

for (scenario in SCENARIOS) {
    chosen <- if (scenario == "injury_only_inferred27") {
        section[!section$inferred_uninjured, ]
    }
    else section
    for (area in EXPECTED_AREAS) {
        d <- chosen[chosen$sample_label == area, ]
        if (!nrow(d)) 
            stop("No pieces retained in area: ", area, " scenario ", scenario)
        st <- unique(d$GEO_stage_original)
        if (length(st) != 1L) 
            stop("Inconsistent GEO timepoint within capture area: ", area)
        rr1 <- raw_regions[raw_regions$section_uid_internal %in% d$section_uid_internal, ]
        denominator <- sum(d$n_spots_RCTD)
        if (sum(rr1$n_spots)/denominator != 1) 
            stop("Capture-region counts mismatch: ", area)
        kk <- kk + 1L
        capture_rows[[kk]] <- data.frame(scenario = scenario, GEO_stage = st, capture_area = area, GSM = unique(d$GSM)[1], 
            n_sections = nrow(d), n_uninjured_included = sum(d$inferred_uninjured), n_RCTD_spots = denominator, 
            Clec7a_equal_section_mean = mean_ok(d$Clec7a_mean_RCTD), Clec7a_pooled_spot_mean = weighted_ok(d$Clec7a_mean_RCTD, 
                d$n_spots_RCTD), Clec7a_section_sd = sd_ok(d$Clec7a_mean_RCTD), Clec7a_section_min = range_low(d$Clec7a_mean_RCTD), 
            Clec7a_section_max = range_high(d$Clec7a_mean_RCTD), injury_region_fraction_pooled = weighted_ok(d$injury_region_fraction, 
                d$n_spots_RCTD), injury_region_fraction_equal_section = mean_ok(d$injury_region_fraction), 
            target_region_fraction_pooled = weighted_ok(d$target_region_fraction, d$n_spots_RCTD), target_region_fraction_equal_section = mean_ok(d$target_region_fraction), 
            target_region_burden_pooled = weighted_ok(d$target_region_spatial_burden, d$n_spots_RCTD), 
            target_region_burden_equal_section = mean_ok(d$target_region_spatial_burden), median_section_UMI = median_ok(d$RNA_depth_proxy), 
            median_section_genes = median_ok(d$detected_genes_proxy), median_section_percent_mt = median_ok(d$mitochondrial_proxy), 
            mean_section_retention = mean_ok(d$QC_retention_fraction), mean_RCTD_normal_weight = mean_ok(d$mean_RCTD_normal_weight), 
            mean_RCTD_injury_myeloid_weight = mean_ok(d$mean_RCTD_injury_myeloid_weight), stringsAsFactors = FALSE)
        for (region in region_names) {
            r <- rr1[rr1$spatial_region_v03 == region, ]
            rr <- rr + 1L
            capture_region_rows[[rr]] <- data.frame(scenario = scenario, GEO_stage = st, capture_area = area, 
                region = region, region_class = unique(r$region_class), n_sections = nrow(d), n_region_spots = sum(r$n_spots), 
                n_total_spots = denominator, fraction_pooled = sum(r$n_spots)/denominator, fraction_equal_section = mean_ok(r$region_fraction), 
                n_sections_with_region = sum(r$n_spots > 0), target_weight_sum_region = sum(r$region_target_weight_sum), 
                region_niche_burden_pooled = sum(r$region_target_weight_sum)/denominator, region_niche_burden_equal_section = mean_ok(r$spatial_niche_burden), 
                stringsAsFactors = FALSE)
        }
    }
}

capture <- do.call(rbind, capture_rows)
capture_region <- do.call(rbind, capture_region_rows)
rownames(capture) <- NULL
rownames(capture_region) <- NULL
if (nrow(capture) != 16L || nrow(capture_region) != 160L) stop("Capture area summary dimensions are incorrect (expected 16 areas and 160 region rows).")
save_csv(capture, 4, "capture_area_metrics_before_after")
save_csv(capture_region, 5, "capture_area_10class_regions_before_after")

get_capture_region <- function(scenario, area, column = "fraction_pooled") {
  z <- capture_region[capture_region$scenario == scenario & capture_region$capture_area == area, ]
  if (!setequal(z$region, region_names))
    stop("Incomplete capture region: ", area)
  unname(z[[column]][match(region_names, z$region)])
}
get_section_region <- function(uid) {
  z <- raw_regions[raw_regions$section_uid_internal == uid, ]
  unname(z$region_fraction[match(region_names, z$spatial_region_v03)])
}

heterogeneity_rows <- list()
hh <- 0L
region_delta_rows <- list()
dd <- 0L

for (scenario in SCENARIOS) for (stage in STAGES) {
  ca <- capture[capture$scenario == scenario & capture$GEO_stage == stage, ]
  ca <- ca[order(ca$capture_area), ]
  if (nrow(ca) != 2L)
    stop("Expected exactly two capture areas for stage ", stage)
  a <- ca$capture_area[1]
  b <- ca$capture_area[2]
  pa <- get_capture_region(scenario, a)
  pb <- get_capture_region(scenario, b)
  pa_equal <- get_capture_region(scenario, a, "fraction_equal_section")
  pb_equal <- get_capture_region(scenario, b, "fraction_equal_section")
  deltas <- pb - pa
  hh <- hh + 1L
  heterogeneity_rows[[hh]] <- data.frame(scenario = scenario, stage = stage, area_1 = a, area_2 = b,
    n_sections_area1 = ca$n_sections[1], n_sections_area2 = ca$n_sections[2], JS_distance_pooled = jsd(pa,
      pb), JS_distance_equal_section = jsd(pa_equal, pb_equal), TV_distance_pooled = tv_distance(pa,
      pb), max_region_fraction_absolute_difference = max(abs(deltas)), max_difference_region = region_names[which.max(abs(deltas))],
    Clec7a_mean_area1 = ca$Clec7a_equal_section_mean[1], Clec7a_mean_area2 = ca$Clec7a_equal_section_mean[2],
    Clec7a_abs_difference = abs(diff(ca$Clec7a_equal_section_mean)), Clec7a_area1_div_area2 = safe_ratio(ca$Clec7a_equal_section_mean[1],
      ca$Clec7a_equal_section_mean[2]), injury_fraction_area1 = ca$injury_region_fraction_pooled[1],
    injury_fraction_area2 = ca$injury_region_fraction_pooled[2], injury_fraction_abs_difference = abs(diff(ca$injury_region_fraction_pooled)),
    target_region_burden_area1 = ca$target_region_burden_pooled[1], target_region_burden_area2 = ca$target_region_burden_pooled[2],
    target_region_burden_abs_difference = abs(diff(ca$target_region_burden_pooled)), stringsAsFactors = FALSE)
  for (i in seq_along(region_names)) {
    dd <- dd + 1L
    region_delta_rows[[dd]] <- data.frame(scenario = scenario, stage = stage, region = region_names[i],
      area1 = a, area2 = b, proportion_area1 = pa[i], proportion_area2 = pb[i], delta_area2_minus_area1 = deltas[i],
      absolute_difference = abs(deltas[i]), stringsAsFactors = FALSE)
  }
}
heterogeneity <- do.call(rbind, heterogeneity_rows)
region_deltas <- do.call(rbind, region_delta_rows)
save_csv(heterogeneity, 6, "within_stage_capture_heterogeneity")
save_csv(region_deltas, 7, "region_resolved_between_capture_differences")

change_rows <- lapply(STAGES, function(stage) {
  old <- heterogeneity[heterogeneity$stage == stage & heterogeneity$scenario == "GEO_label_all32",
    ]
  new <- heterogeneity[heterogeneity$stage == stage & heterogeneity$scenario == "injury_only_inferred27",
    ]
  data.frame(stage = stage, JS_before = old$JS_distance_pooled, JS_after = new$JS_distance_pooled,
    JS_absolute_change = new$JS_distance_pooled - old$JS_distance_pooled, JS_percent_reduction = safe_ratio(old$JS_distance_pooled -
      new$JS_distance_pooled, old$JS_distance_pooled) * 100, TV_before = old$TV_distance_pooled,
    TV_after = new$TV_distance_pooled, Clec7a_gap_before = old$Clec7a_abs_difference, Clec7a_gap_after = new$Clec7a_abs_difference,
    Clec7a_gap_percent_reduction = safe_ratio(old$Clec7a_abs_difference - new$Clec7a_abs_difference,
      old$Clec7a_abs_difference) * 100, injury_fraction_gap_before = old$injury_fraction_abs_difference,
    injury_fraction_gap_after = new$injury_fraction_abs_difference, target_burden_gap_before = old$target_region_burden_abs_difference,
    target_burden_gap_after = new$target_region_burden_abs_difference, stringsAsFactors = FALSE)
})
change <- do.call(rbind, change_rows)
save_csv(change, 8, "effect_of_five_inferred_controls_on_capture_heterogeneity")

section_pairs <- list()
pp <- 0L
for (scenario in SCENARIOS) for (stage in STAGES) {
  d <- section[section$GEO_stage_original == stage & (scenario == "GEO_label_all32" | !section$inferred_uninjured),
    ]
  pairs <- utils::combn(seq_len(nrow(d)), 2)
  for (j in seq_len(ncol(pairs))) {
    x <- d[pairs[1, j], ]
    y <- d[pairs[2, j], ]
    pp <- pp + 1L
    section_pairs[[pp]] <- data.frame(scenario = scenario, stage = stage, section_1 = x$section_uid_internal,
      section_2 = y$section_uid_internal, capture_1 = x$sample_label, capture_2 = y$sample_label,
      same_capture = x$sample_label == y$sample_label, JS_distance = jsd(get_section_region(x$section_uid_internal),
        get_section_region(y$section_uid_internal)), absolute_Clec7a_difference = abs(x$Clec7a_mean_RCTD -
        y$Clec7a_mean_RCTD), stringsAsFactors = FALSE)
  }
}
section_pairs <- do.call(rbind, section_pairs)
save_csv(section_pairs, 9, "within_stage_section_pairwise_heterogeneity")

pair_groups <- split(section_pairs, interaction(section_pairs$scenario, section_pairs$stage, section_pairs$same_capture,
  drop = TRUE))
pair_summary <- do.call(rbind, lapply(pair_groups, function(x) {
  data.frame(scenario = x$scenario[1], stage = x$stage[1], same_capture = x$same_capture[1], n_pairs = nrow(x),
    JS_median = median_ok(x$JS_distance), JS_min = range_low(x$JS_distance), JS_max = range_high(x$JS_distance),
    Clec7a_absolute_difference_median = median_ok(x$absolute_Clec7a_difference), stringsAsFactors = FALSE)
}))
save_csv(pair_summary, 10, "within_vs_between_capture_section_variation")

qc_area <- merge(capture[capture$scenario == "injury_only_inferred27", ], qc_sample[, c("sample_label",
  "n_spots_before", "n_pass", "retention_rate")], by.x = "capture_area", by.y = "sample_label", all.x = TRUE,
  sort = FALSE)
save_csv(qc_area, 11, "corrected_capture_area_QC_technical_audit")

tech_predictors <- c("RNA_depth_proxy", "detected_genes_proxy", "mitochondrial_proxy", "QC_retention_fraction",
  "n_spots_RCTD")
tech_outcomes <- c("Clec7a_mean_RCTD", "injury_region_fraction", "target_region_fraction", "target_region_spatial_burden")
tech_rows <- list()
tt <- 0L
for (subset_name in c("all_inferred_injury27", STAGES)) {
  z <- section[!section$inferred_uninjured, ]
  if (subset_name != "all_inferred_injury27")
    z <- z[z$GEO_stage_original == subset_name, ]
  for (predictor in tech_predictors) for (outcome in tech_outcomes) {
    x <- as.numeric(z[[predictor]])
    y <- as.numeric(z[[outcome]])
    ok <- is.finite(x) & is.finite(y)
    rho <- if (sum(ok) >= 4L && length(unique(x[ok])) > 1L && length(unique(y[ok])) > 1L) {
      suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
    } else NA_real_
    tt <- tt + 1L
    tech_rows[[tt]] <- data.frame(subset = subset_name, predictor = predictor, outcome = outcome,
      n_pieces = sum(ok), n_capture_areas = length(unique(z$sample_label[ok])), spearman_rho_descriptive = rho,
      inference = "none: nested pieces and very few capture areas", stringsAsFactors = FALSE)
  }
}
tech <- do.call(rbind, tech_rows)
save_csv(tech, 12, "exploratory_QC_vs_spatial_associations")

batch_audit <- data.frame()
if (file.exists(INPUT$pca_structure)) {
  pca <- read_data(INPUT$pca_structure)
  batch_audit <- rbind(batch_audit, data.frame(source = "Step62_preexisting_PCA_sample_structure",
    metric = "sum_weighted_extra_sample_structure", value = if ("weighted_extra_sample_structure" %in%
      names(pca))
      sum(pca$weighted_extra_sample_structure, na.rm = TRUE) else NA_real_, detail = "PCA analysis before five-piece classification; global aggregate", stringsAsFactors = FALSE))
}
if (file.exists(INPUT$harmony_structure)) {
  hs <- read_data(INPUT$harmony_structure)
  if (all(c("method", "weighted_extra_sample_R2") %in% names(hs))) {
    batch_audit <- rbind(batch_audit, data.frame(source = "Step63_preexisting_Harmony_sensitivity",
      metric = paste0("weighted_extra_sample_R2_", hs$method), value = hs$weighted_extra_sample_R2,
      detail = "Preexisting, before inferred control exclusion; no causal attribution", stringsAsFactors = FALSE))
  }
}
if (nrow(batch_audit) > 0L) save_csv(batch_audit, 13, "preexisting_batch_structure_context")

if (nrow(representatives) != 4L || anyDuplicated(representatives$timepoint) || !setequal(representatives$timepoint,
  STAGES)) stop("Step69 Figure4A representative manifest must contain one piece for each stage.")

rep_all <- list()
rep_rank_rows <- list()
uu <- 0L
vv <- 0L
for (stage in STAGES) {
  d <- section[!section$inferred_uninjured & section$GEO_stage_original == stage, ]
  ca <- capture[capture$scenario == "injury_only_inferred27" & capture$GEO_stage == stage, ]
  selected_uid <- representatives$section_uid_internal[representatives$timepoint == stage]
  if (!selected_uid %in% d$section_uid_internal)
    stop("Original Step69 representative is not among corrected injury pieces: ", selected_uid)

  stage_center <- (get_capture_region("injury_only_inferred27", ca$capture_area[1]) + get_capture_region("injury_only_inferred27",
    ca$capture_area[2]))/2
  stage_mean_target <- mean_ok(ca$Clec7a_equal_section_mean)
  sample_medians <- median(d$Clec7a_mean_RCTD)

  ranks <- data.frame(stage = stage, section_uid_internal = d$section_uid_internal, sample_label = d$sample_label,
    originally_selected = d$section_uid_internal == selected_uid, JS_to_stage_center = vapply(d$section_uid_internal,
      function(uid) jsd(get_section_region(uid), stage_center), numeric(1)), Clec7a_abs_diff_to_stage_capture_equal_mean = abs(d$Clec7a_mean_RCTD -
      stage_mean_target), Clec7a_abs_diff_to_stage_section_median = abs(d$Clec7a_mean_RCTD - sample_medians),
    stringsAsFactors = FALSE)
  ranks$rank_JS_closer_is_better <- rank(ranks$JS_to_stage_center, ties.method = "min")
  ranks$rank_Clec7a_closer_is_better <- rank(ranks$Clec7a_abs_diff_to_stage_capture_equal_mean, ties.method = "min")
  ranks$stage_n_sections <- nrow(d)
  vv <- vv + 1L
  rep_rank_rows[[vv]] <- ranks

  r <- ranks[ranks$originally_selected, ]
  uu <- uu + 1L
  rep_all[[uu]] <- data.frame(stage = stage, original_Figure4A_section = selected_uid, n_corrected_injury_sections = nrow(d),
    n_capture_areas = nrow(ca), JS_to_stage_center = r$JS_to_stage_center, JS_rank_among_stage_sections = r$rank_JS_closer_is_better,
    JS_percentile_rank = (r$rank_JS_closer_is_better - 1)/max(1, nrow(d) - 1), median_JS_all_stage_sections = median_ok(ranks$JS_to_stage_center),
    Clec7a_representative = section$Clec7a_mean_RCTD[match(selected_uid, section$section_uid_internal)],
    Clec7a_stage_capture_equal_mean = stage_mean_target, Clec7a_abs_difference = r$Clec7a_abs_diff_to_stage_capture_equal_mean,
    Clec7a_closeness_rank = r$rank_Clec7a_closer_is_better, stage_n_capture_areas = 2L, interpretation = "Illustrative section only; formal representativeness not established",
    stringsAsFactors = FALSE)
}
rep_summary <- do.call(rbind, rep_all)
rep_section_ranks <- do.call(rbind, rep_rank_rows)
save_csv(rep_summary, 14, "Figure4A_representative_section_stage_consistency")
save_csv(rep_section_ranks, 15, "all_injury_sections_distance_to_stage_center")

qc_checks <- data.frame(check = c("32 unique original coordinate-defined pieces", "8 distinct GEO capture areas",
  "4 coordinate components per GEO capture area", "10 region classes per piece", "5 corrected Step114-derived inferred uninjured pieces",
  "27 retained putative injured pieces: 8/7/6/6", "All 32 Step65 section Clec7a means joined", "Region counts match the QC-passed Step65 piece denominators",
  "All 32 pieces have per-piece Step61 QC results", "Two capture areas per stage in both scenarios",
  "Four original Figure4A representatives checked", "Figure4B plotted sections agree with Figure4A (if Step70 is available)",
  "No original animal IDs or OCT-block IDs are manufactured"), pass = c(nrow(section) == 32L && !anyDuplicated(section$section_uid_internal),
  length(unique(section$sample_label)) == 8L, all(table(section$sample_label)[EXPECTED_AREAS] == 4L),
  all(table(raw_regions$section_uid_internal) == 10L), sum(section$inferred_uninjured) == 5L, sum(!section$inferred_uninjured) ==
    27L && all(as.integer(counts_clean) == as.integer(EXPECTED_CORRECTED[STAGES])), all(is.finite(section$Clec7a_mean_RCTD)),
  all(section$n_spots_region_sum == section$n_spots_RCTD), all(is.finite(section$QC_retention_fraction)),
  all(table(capture$scenario, capture$GEO_stage) == 2L), nrow(rep_summary) == 4L, is.na(figure4B_matching) ||
    isTRUE(figure4B_matching), all(is.na(section$animal_id)) && all(is.na(section$tissue_block_id))),
  stringsAsFactors = FALSE)
save_csv(qc_checks, 16, "FINAL_QC")
if (!all(qc_checks$pass)) stop("Step117 final QC failed; inspect 16_117_FINAL_QC_v01.csv")

plot_save <- function(plot_object, stem, width, height) {
  pdf <- file.path(DIR_FIGURES, paste0(stem, "_", VERSION, ".pdf"))
  png <- file.path(DIR_FIGURES, paste0(stem, "_", VERSION, ".png"))
  ggplot2::ggsave(pdf, plot = plot_object, width = width, height = height, units = "in", device = "pdf",
    bg = "white")
  ggplot2::ggsave(png, plot = plot_object, width = width, height = height, units = "in", dpi = 220,
    bg = "white")
}

fig20 <- ggplot(capture_region, aes(x = capture_area, y = fraction_pooled, fill = region)) + geom_col(width = 0.82) +
  facet_grid(scenario ~ GEO_stage, scales = "free_x", space = "free_x") + scale_y_continuous(labels = function(z) paste0(round(100 *
  z), "%")) + labs(title = "Capture-area transcriptional-region composition: 32 vs inferred 27 pieces",
  subtitle = "Spot-count pooled proportions; two capture areas per stage", x = "GEO capture area",
  y = "Regional spot fraction", fill = "Transcriptional region") + theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 50,
  hjust = 1, size = 7), legend.position = "bottom", legend.text = element_text(size = 7), strip.text = element_text(size = 8),
  plot.title = element_text(size = 11)) + guides(fill = guide_legend(ncol = 2))
plot_save(fig20, "20_117_capture_region_composition_before_after", 14, 8)

fig21 <- ggplot(capture, aes(x = scenario, y = Clec7a_equal_section_mean, group = capture_area, color = capture_area)) +
  geom_line(linewidth = 0.65) + geom_point(size = 2.4) + facet_wrap(~GEO_stage, ncol = 4) + labs(title = "Change in Clec7a RCTD weight after inferred-control exclusion",
  subtitle = "Equal section mean within capture area; no animal-level inference", x = NULL, y = "Mean of piece-level Clec7a RCTD means",
  color = "Capture area") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35,
  hjust = 1), legend.position = "bottom", plot.title = element_text(size = 11))
plot_save(fig21, "21_117_Clec7a_capture_area_sensitivity", 13, 5.5)

fig22 <- ggplot(section, aes(x = sample_label, y = Clec7a_mean_RCTD, shape = inferred_uninjured, color = inferred_uninjured)) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 117), size = 2.5) + facet_wrap(~GEO_stage_original,
  scales = "free_x", ncol = 4) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "Putative injured",
  `TRUE` = "Inferred uninjured")) + scale_color_manual(values = c(`FALSE` = "#2864A4", `TRUE` = "#D95F02"),
  labels = c(`FALSE` = "Putative injured", `TRUE` = "Inferred uninjured")) + labs(title = "Clec7a RCTD weights across the 32 coordinate-defined pieces",
  subtitle = "Each dot is a piece, not a capture area or verified animal-ID record", x = "Capture area",
  y = "Mean piece-level Clec7a weight", shape = NULL, color = NULL) + theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")
plot_save(fig22, "22_117_section_variability_inferred_control_flag", 12, 5)

fig23 <- ggplot(heterogeneity, aes(x = scenario, y = JS_distance_pooled, group = stage, color = stage)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.8) + labs(title = "Between-capture-area compositional heterogeneity by stage",
  subtitle = "Jensen-Shannon distance (base 2); 0=same composition; 1=maximally different", x = NULL,
  y = "Jensen-Shannon distance", color = "Stage") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 25,
  hjust = 1), legend.position = "bottom")
plot_save(fig23, "23_117_stage_Jensen_Shannon_before_after", 9, 5.5)

fig24 <- ggplot(rep_section_ranks, aes(x = reorder(section_uid_internal, JS_to_stage_center), y = JS_to_stage_center,
  color = originally_selected)) + geom_point(size = 2.3) + geom_segment(aes(xend = section_uid_internal,
  y = 0, yend = JS_to_stage_center), linewidth = 0.35) + coord_flip() + facet_wrap(~stage, scales = "free_y",
  ncol = 2) + scale_color_manual(values = c(`FALSE` = "#777777", `TRUE` = "#C04C36"), labels = c(`FALSE` = "Other section",
  `TRUE` = "Original Figure4A")) + labs(title = "How representative are the four original Figure 4A sections?",
  subtitle = "Section-to-stage equal-capture-center JS distance; smaller = more similar", x = NULL,
  y = "JS distance to corrected stage composition", color = NULL) + theme_classic(base_size = 9) +
  theme(legend.position = "bottom")
plot_save(fig24, "24_117_Figure4A_representative_consistency", 12, 8)

qc_plot <- section[!section$inferred_uninjured, ]
fig25 <- ggplot(qc_plot, aes(x = RNA_depth_proxy, y = Clec7a_mean_RCTD, color = sample_label)) + geom_point(size = 2) +
  facet_wrap(~GEO_stage_original, ncol = 4, scales = "free_x") + labs(title = "Clec7a weight and sample RNA complexity across putative injured pieces",
  subtitle = "Exploratory, nested pieces; no independent-animal effect estimation", x = "Median UMI of QC-passed spots per piece",
  y = "Mean Clec7a RCTD weight", color = "Capture area") + theme_classic(base_size = 9) + theme(legend.position = "bottom")
plot_save(fig25, "25_117_RNA_depth_vs_Clec7a_section_QC", 12, 5.5)

writeLines("Interpretation notes: ANALYSIS_NOTES.md. Numerical results: accompanying output tables.", 
    file.path(OUT, paste0("17_", STEP, "_REPORT_", VERSION, ".txt")), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, paste0("18_", STEP, "_sessionInfo_", VERSION,
  ".txt")))

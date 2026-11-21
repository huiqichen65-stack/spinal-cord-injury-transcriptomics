source("config.R", local = TRUE)
source("scripts/helpers/common.R", local = TRUE)

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- PROJECT_ROOT
SCRIPT_NAME <- "146_Reviewer1_RemainingMajor3_high_confidence_microglia_sensitivity_v01.R"
V <- "v01"

INPUT_IDENTITY <- file.path( ROOT, "results", "106_Reviewer1_Major2_identity_programs", "02_106_per_cell_identity_scores_v01.csv.gz" )

INPUT_AUCELL <- file.path( ROOT, "results", "20_AUCell_lipid_metabolism_scoring", "20_AUCell_cell_scores_v03.csv" )

OUT <- file.path( ROOT, "results", "146_Reviewer1_RemainingMajor3_high_confidence_microglia_sensitivity" )
FIG_DIR <- file.path(OUT, "figures")
LOG_DIR <- file.path(OUT, "logs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(LOG_DIR, "00_146_run_log_v01.log")
if (file.exists(LOG)) unlink(LOG)

SOURCE_ORDER <- c("17702045", "GSE162610", "GSE196928")
STAGE_ORDER <- c("Control", "4 h", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi", "38 dpi")

MAIN_FIGURE_STAGES <- c("Control", "1 dpi", "3 dpi", "7 dpi", "14 dpi", "28 dpi")

MG_STATES <- c( "Homeostatic Microglia", "Clec7a+ lipid-associated Microglia", "Activated inflammatory Microglia", "Ccl12+ acute Microglia", "Tnf+ inflammatory Microglia",
  "Proliferating Microglia" )

MAC_STATES <- c( "MHC-II+ antigen-presenting Macrophage", "Gpnmb+ tissue-repair Macrophage", "Acute-response Macrophage" )

FOCUS <- "Clec7a+ lipid-associated Microglia"
HOME <- "Homeostatic Microglia"

PRIMARY_LIPID_SCORE <- "AUC_FA_uptake_transport"

MG_MIN_SCORE <- 0

MIN_CELLS_PER_STATE_SOURCE <- 20L

pkgs <- c("dplyr", "ggplot2")
missing_pkg <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) {
  stop("Required R packages not installed: ", paste(missing_pkg, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

logit <- function(...) {
  x <- paste0(..., collapse = "")
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", x)
  cat(line, "\n")
  cat(line, "\n", file = LOG, append = TRUE)
}

write_csv <- function(x, filename) {
  f <- file.path(OUT, filename)
  utils::write.csv( x, f, row.names = FALSE, na = "", fileEncoding = "UTF-8" )
  logit("Wrote ", basename(f), " | rows=", nrow(x))
  invisible(f)
}

write_csv_gz <- function(x, filename) {
  f <- file.path(OUT, filename)
  con <- gzfile(f, open = "wt", compression = 8)
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "")
  logit("Wrote ", basename(f), " | rows=", nrow(x))
  invisible(f)
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  z <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(z))
  out[z %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[z %in% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  out
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

q25 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, 0.25, names = FALSE))
}

q75 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, 0.75, names = FALSE))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.finite(x), paste0(formatC(x, digits = digits, format = "f"), "%"), "NA")
}

file_md5 <- function(f) {
  if (!file.exists(f)) return(NA_character_)
  unname(tools::md5sum(f))
}

file_info_row <- function(role, f) {
  fi <- file.info(f)
  data.frame( input_role = role, path = normalizePath(f, winslash = "/", mustWork = FALSE), exists = file.exists(f), size_bytes = if (file.exists(f)) fi$size else NA_real_,
    modified_time = if (file.exists(f)) format(fi$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_, md5 = file_md5(f), stringsAsFactors = FALSE )
}

forbidden_input_re <- "\\.(rds|rdata|rda|rdb|rdx|mtx|loom|h5|h5ad|feather)(\\.gz)?$"
if (any(grepl(forbidden_input_re,
              tolower(c(INPUT_IDENTITY, INPUT_AUCELL)),
              perl = TRUE))) {
  stop("Forbidden RDS/matrix-like input configured. Step146 must remain CSV-only.")
}

for (f in c(INPUT_IDENTITY, INPUT_AUCELL)) {
  if (!file.exists(f)) stop("Required direct input missing: ", f)
}

logit("Step 146 started")
logit("SCRIPT_NAME: ", SCRIPT_NAME)
logit("INPUT_IDENTITY: ", INPUT_IDENTITY)
logit("INPUT_AUCELL: ", INPUT_AUCELL)
logit("OUTPUT: ", OUT)
logit("No ZIP, RDS, expression matrix, Seurat, or trajectory input is used.")

logit("Reading Step106 per-cell identity scores...")
identity <- utils::read.csv( gzfile(INPUT_IDENTITY), check.names = FALSE, stringsAsFactors = FALSE )

need_identity <- c( "cell_id", "source_dataset", "timepoint_main", "sample_id", "sample_id_source", "myeloid_state", "eligible_identity_scoring", "MG_score", "MDM_score", "BAM_score" )
miss_identity <- setdiff(need_identity, names(identity))
if (length(miss_identity)) {
  stop("Step106 missing required columns: ", paste(miss_identity, collapse = ", "))
}
if (anyDuplicated(identity$cell_id)) stop("Duplicated cell_id in Step106 input.")

identity$source_dataset <- as.character(identity$source_dataset)
identity$timepoint_main <- as.character(identity$timepoint_main)
identity$myeloid_state <- as.character(identity$myeloid_state)
identity$eligible_identity_scoring <- as_bool(identity$eligible_identity_scoring)
identity$MG_score <- safe_num(identity$MG_score)
identity$MDM_score <- safe_num(identity$MDM_score)
identity$BAM_score <- safe_num(identity$BAM_score)

if (any(!identity$source_dataset %in% SOURCE_ORDER)) {
  stop( "Unexpected source_dataset values: ", paste(setdiff(unique(identity$source_dataset), SOURCE_ORDER), collapse = ", ") )
}
if (any(!identity$timepoint_main %in% STAGE_ORDER)) {
  stop( "Unexpected timepoint_main values: ", paste(setdiff(unique(identity$timepoint_main), STAGE_ORDER), collapse = ", ") )
}

expected_states <- c(MG_STATES, MAC_STATES)
unexpected_states <- setdiff(unique(identity$myeloid_state), expected_states)
if (length(unexpected_states)) {
  stop( "Unexpected frozen myeloid-state labels: ", paste(unexpected_states, collapse = " | ") )
}

logit("Reading Step20 precomputed AUCell scores...")
auc <- utils::read.csv( INPUT_AUCELL, check.names = FALSE, stringsAsFactors = FALSE )

need_auc <- c("cell", "myeloid_state", "source_dataset", "timepoint_main", PRIMARY_LIPID_SCORE)
miss_auc <- setdiff(need_auc, names(auc))
if (length(miss_auc)) {
  stop("Step20 missing required columns: ", paste(miss_auc, collapse = ", "))
}
if (anyDuplicated(auc$cell)) stop("Duplicated cell IDs in Step20 input.")

score_cols <- grep("^AUC_", names(auc), value = TRUE)
if (!length(score_cols)) stop("No AUC_ score columns found in Step20.")
if (!PRIMARY_LIPID_SCORE %in% score_cols) {
  stop("Primary lipid score missing: ", PRIMARY_LIPID_SCORE)
}

for (z in score_cols) auc[[z]] <- safe_num(auc[[z]])
auc$source_dataset <- as.character(auc$source_dataset)
auc$timepoint_main <- as.character(auc$timepoint_main)
auc$myeloid_state <- as.character(auc$myeloid_state)

if (!setequal(identity$cell_id, auc$cell)) {
  stop( "Step106 and Step20 cell sets differ. ", "Step146 will not silently subset inconsistent inputs." )
}

auc <- auc[match(identity$cell_id, auc$cell), , drop = FALSE]
if (!identical(identity$cell_id, auc$cell)) stop("Cell-order alignment failed.")

for (z in c("source_dataset", "timepoint_main", "myeloid_state")) {
  if (!identical(as.character(identity[[z]]), as.character(auc[[z]]))) {
    stop("Step106 vs Step20 metadata mismatch after cell matching: ", z)
  }
}

input_audit <- rbind( file_info_row("Step106_identity_scores", INPUT_IDENTITY), file_info_row("Step20_AUCell_scores", INPUT_AUCELL) )
input_audit$n_rows <- c(nrow(identity), nrow(auc))
input_audit$n_unique_cells <- c(length(unique(identity$cell_id)), length(unique(auc$cell)))
write_csv(input_audit, "01_146_input_audit_v01.csv")

identity$score_available <-
  !is.na(identity$eligible_identity_scoring) &
  identity$eligible_identity_scoring &
  is.finite(identity$MG_score) &
  is.finite(identity$MDM_score) &
  is.finite(identity$BAM_score)

identity$rule_MG_nonnegative <-
  identity$score_available & identity$MG_score >= MG_MIN_SCORE

identity$rule_MG_gt_MDM <-
  identity$score_available & identity$MG_score > identity$MDM_score

identity$rule_MG_gt_BAM <-
  identity$score_available & identity$MG_score > identity$BAM_score

identity$HC_microglia_consistent <-
  identity$score_available &
  identity$rule_MG_nonnegative &
  identity$rule_MG_gt_MDM &
  identity$rule_MG_gt_BAM

identity$original_MG_labeled <- identity$myeloid_state %in% MG_STATES
identity$original_macrophage_labeled <- identity$myeloid_state %in% MAC_STATES
identity$original_focus_state <- identity$myeloid_state == FOCUS

identity$confidence_status <- ifelse( identity$HC_microglia_consistent, "expression_based_HC_microglia_consistent", "did_not_meet_HC_criterion" )

identity$fail_score_unavailable <- !identity$score_available
identity$fail_MG_below_zero <- identity$score_available & identity$MG_score < MG_MIN_SCORE
identity$fail_MG_not_gt_MDM <- identity$score_available & !(identity$MG_score > identity$MDM_score)
identity$fail_MG_not_gt_BAM <- identity$score_available & !(identity$MG_score > identity$BAM_score)

rule_spec <- data.frame( rule_name = "expression_based_HC_microglia_consistent", component = c( "input_eligibility", "microglia_score_floor", "microglia_vs_MDM", "microglia_vs_BAM",
  "original_label_requirement", "interpretation" ), criterion = c( "eligible_identity_scoring == TRUE and all three program scores finite", "MG_score >= 0", "MG_score > MDM_score",
  "MG_score > BAM_score", "none; rule is applied to all frozen myeloid cells", paste0( "expression-based high-confidence microglia-consistent; ", "not proof of resident microglial lineage" )
  ), fixed_before_outcome_comparison = TRUE, stringsAsFactors = FALSE )
write_csv(rule_spec, "02_146_confidence_rule_specification_v01.csv")

per_cell_out <- identity[, c( "cell_id", "source_dataset", "timepoint_main", "sample_id", "sample_id_source", "myeloid_state", "eligible_identity_scoring", "score_available",
  "MG_score", "MDM_score", "BAM_score", "rule_MG_nonnegative", "rule_MG_gt_MDM", "rule_MG_gt_BAM", "HC_microglia_consistent", "confidence_status",
  "original_MG_labeled", "original_macrophage_labeled", "original_focus_state", "fail_score_unavailable", "fail_MG_below_zero", "fail_MG_not_gt_MDM", "fail_MG_not_gt_BAM" )]
write_csv_gz( per_cell_out, "03_146_per_cell_confidence_assignment_v01.csv.gz" )

retention_by_state <- identity %>%
  group_by(myeloid_state) %>%
  summarise( n_total = n(), n_score_available = sum(score_available), n_HC = sum(HC_microglia_consistent), n_not_HC = n_total - n_HC, pct_HC_of_total = 100 * n_HC / n_total,
    pct_HC_of_score_available = ifelse( n_score_available > 0, 100 * n_HC / n_score_available, NA_real_ ), .groups = "drop" ) %>%
  mutate( state_class = ifelse( myeloid_state %in% MG_STATES, "original_microglia_labeled_state", "original_macrophage_labeled_state" ) )

state_order <- c(MG_STATES, MAC_STATES)
retention_by_state$myeloid_state <- factor( retention_by_state$myeloid_state, levels = state_order )
retention_by_state <- retention_by_state %>%
  arrange(myeloid_state) %>%
  mutate(myeloid_state = as.character(myeloid_state))

write_csv( retention_by_state, "04_146_confidence_retention_by_state_v01.csv" )

retention_source_stage_state <- identity %>%
  group_by(source_dataset, timepoint_main, myeloid_state) %>%
  summarise( n_total = n(), n_score_available = sum(score_available), n_HC = sum(HC_microglia_consistent), n_not_HC = n_total - n_HC, pct_HC_of_total = 100 * n_HC / n_total,
    pct_HC_of_score_available = ifelse( n_score_available > 0, 100 * n_HC / n_score_available, NA_real_ ), .groups = "drop" )

write_csv( retention_source_stage_state, "05_146_confidence_retention_by_source_stage_state_v01.csv" )

identity$original_label_class <- ifelse( identity$original_MG_labeled, "original_microglia_labeled", "original_macrophage_labeled" )

label_concordance <- identity %>%
  group_by(original_label_class, confidence_status) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(original_label_class) %>%
  mutate( class_total = sum(n_cells), pct_within_original_label_class = 100 * n_cells / class_total ) %>%
  ungroup()

write_csv( label_concordance, "06_146_original_label_concordance_v01.csv" )

plot_ret <- retention_by_state
plot_ret$myeloid_state <- factor(plot_ret$myeloid_state, levels = rev(state_order))

pA <- ggplot(plot_ret, aes(x = myeloid_state, y = pct_HC_of_total)) +
  geom_col(width = 0.72) +
  coord_flip() +
  labs( x = NULL, y = "Cells meeting expression-based HC criterion (%)", title = "Reviewer 1 Remaining Major Comment 3",
    subtitle = "High-confidence microglia-consistent assignment by original myeloid state" ) +
  theme_classic(base_size = 11)

ggsave( file.path(FIG_DIR, "Fig_146A_high_confidence_retention_by_state_v01.pdf"), pA, width = 8.3, height = 5.5 )
ggsave( file.path(FIG_DIR, "Fig_146A_high_confidence_retention_by_state_v01.png"), pA, width = 8.3, height = 5.5, dpi = 300 )

make_temporal_row <- function(df, cohort_name, denominator_definition) {
  out <- df %>%
    filter(timepoint_main %in% MAIN_FIGURE_STAGES) %>%
    group_by(timepoint_main) %>%
    summarise( denominator_n = n(), Clec7a_focal_n = sum(myeloid_state == FOCUS), Clec7a_focal_fraction = Clec7a_focal_n / denominator_n, Clec7a_focal_percent = 100 * Clec7a_focal_fraction,
      .groups = "drop" )
  out$cohort <- cohort_name
  out$denominator_definition <- denominator_definition
  out
}

temporal_all <- make_temporal_row( identity, "Manuscript reference: all myeloid cells", "all frozen myeloid cells at each displayed stage" )

temporal_mg_label <- make_temporal_row( identity[identity$original_MG_labeled, , drop = FALSE], "Original microglia-labeled states",
  "six original microglia-labeled states at each displayed stage" )

temporal_hc <- make_temporal_row( identity[identity$HC_microglia_consistent, , drop = FALSE], "Expression-based HC microglia-consistent",
  "all myeloid cells satisfying the expression-based HC rule at each displayed stage" )

temporal <- bind_rows(temporal_all, temporal_mg_label, temporal_hc)
temporal$timepoint_main <- factor( temporal$timepoint_main, levels = MAIN_FIGURE_STAGES )
temporal <- temporal %>%
  arrange(cohort, timepoint_main) %>%
  mutate(timepoint_main = as.character(timepoint_main))

write_csv( temporal, "07_146_Clec7a_temporal_enrichment_original_vs_HC_v01.csv" )

get_stage_pct <- function(dat, stage) {
  z <- dat$Clec7a_focal_percent[dat$timepoint_main == stage]
  if (length(z) != 1L) return(NA_real_)
  z
}

temporal_metric <- function(dat, cohort_name) {
  p3 <- get_stage_pct(dat, "3 dpi")
  p7 <- get_stage_pct(dat, "7 dpi")
  p14 <- get_stage_pct(dat, "14 dpi")
  p28 <- get_stage_pct(dat, "28 dpi")

  early_max <- max(c(p3, p7), na.rm = TRUE)
  late_min <- min(c(p14, p28), na.rm = TRUE)

  data.frame( cohort = cohort_name, pct_3dpi = p3, pct_7dpi = p7, pct_14dpi = p14, pct_28dpi = p28, early_3_7_max_percent = early_max, late_14_28_min_percent = late_min,
    late_minus_early_separation_pp = late_min - early_max, late_enrichment_pattern = is.finite(late_min) && is.finite(early_max) && late_min > early_max, stringsAsFactors = FALSE )
}

t_all <- temporal_metric( temporal[temporal$cohort == "Manuscript reference: all myeloid cells", ], "Manuscript reference: all myeloid cells" )
t_mg <- temporal_metric( temporal[temporal$cohort == "Original microglia-labeled states", ], "Original microglia-labeled states" )
t_hc <- temporal_metric( temporal[temporal$cohort == "Expression-based HC microglia-consistent", ], "Expression-based HC microglia-consistent" )

ref6 <- temporal %>%
  filter(cohort == "Manuscript reference: all myeloid cells") %>%
  arrange(match(timepoint_main, MAIN_FIGURE_STAGES))
hc6 <- temporal %>%
  filter(cohort == "Expression-based HC microglia-consistent") %>%
  arrange(match(timepoint_main, MAIN_FIGURE_STAGES))

rank_concordance <- if (
  nrow(ref6) == length(MAIN_FIGURE_STAGES) &&
  nrow(hc6) == length(MAIN_FIGURE_STAGES)
) {
  suppressWarnings(stats::cor( ref6$Clec7a_focal_percent, hc6$Clec7a_focal_percent, method = "spearman" ))
} else {
  NA_real_
}

temporal_pattern_summary <- bind_rows(t_all, t_mg, t_hc)
temporal_pattern_summary$rank_concordance_with_manuscript_reference <- c( 1, NA_real_, rank_concordance )
temporal_pattern_summary$HC_sensitivity_conclusion <- c( NA_character_, NA_character_, ifelse( t_hc$late_enrichment_pattern, "principal_late_enrichment_pattern_preserved",
  "principal_late_enrichment_pattern_not_preserved" ) )

write_csv( temporal_pattern_summary, "08_146_temporal_pattern_summary_v01.csv" )

plot_temporal <- temporal %>%
  filter(cohort %in% c( "Manuscript reference: all myeloid cells", "Expression-based HC microglia-consistent" ))

plot_temporal$timepoint_main <- factor( plot_temporal$timepoint_main, levels = MAIN_FIGURE_STAGES )

pB <- ggplot( plot_temporal, aes( x = timepoint_main, y = Clec7a_focal_percent, group = cohort, linetype = cohort, shape = cohort ) ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  labs( x = "Injury stage", y = "Clec7a+ focal-state representation (%)", title = "Temporal pattern before and after expression-based HC restriction",
    subtitle = "Percentages use the cohort-specific denominators stated in the audit table", linetype = NULL, shape = NULL ) +
  theme_classic(base_size = 11) +
  theme( legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1) )

ggsave( file.path(FIG_DIR, "Fig_146B_Clec7a_temporal_enrichment_original_vs_HC_v01.pdf"), pB, width = 7.6, height = 5.2 )
ggsave( file.path(FIG_DIR, "Fig_146B_Clec7a_temporal_enrichment_original_vs_HC_v01.png"), pB, width = 7.6, height = 5.2, dpi = 300 )

auc$HC_microglia_consistent <- identity$HC_microglia_consistent

lipid_base <- auc %>%
  filter(myeloid_state %in% c(FOCUS, HOME))

summarise_lipid_long <- function(df, cohort_name) {
  out <- lapply(score_cols, function(sc) {
    do.call(rbind, lapply(c(FOCUS, HOME), function(st) {
      z <- df[df$myeloid_state == st, sc]
      z <- z[is.finite(z)]
      data.frame( cohort = cohort_name, score = sc, myeloid_state = st, n_cells = length(z), mean = safe_mean(z), median = safe_median(z), sd = safe_sd(z), q25 = q25(z), q75 = q75(z),
        stringsAsFactors = FALSE )
    }))
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

lipid_summary_original <- summarise_lipid_long( lipid_base, "Original operational states" )

lipid_summary_hc <- summarise_lipid_long( lipid_base[lipid_base$HC_microglia_consistent, , drop = FALSE], "Expression-based HC restriction" )

lipid_summary <- bind_rows(lipid_summary_original, lipid_summary_hc)

write_csv( lipid_summary, "09_146_lipid_program_summary_original_vs_HC_v01.csv" )

summarise_lipid_source <- function(df, cohort_name) {
  ans <- list()
  k <- 0L

  for (src in SOURCE_ORDER) {
    ds <- df[df$source_dataset == src, , drop = FALSE]

    for (sc in score_cols) {
      f <- ds[ds$myeloid_state == FOCUS, sc]
      h <- ds[ds$myeloid_state == HOME, sc]
      f <- f[is.finite(f)]
      h <- h[is.finite(h)]

      k <- k + 1L
      ans[[k]] <- data.frame( cohort = cohort_name, source_dataset = src, score = sc, n_focus = length(f), n_homeostatic = length(h), median_focus = safe_median(f),
        median_homeostatic = safe_median(h), median_difference_focus_minus_homeostatic = safe_median(f) - safe_median(h), enough_cells_for_descriptive_direction =
        length(f) >= MIN_CELLS_PER_STATE_SOURCE && length(h) >= MIN_CELLS_PER_STATE_SOURCE, stringsAsFactors = FALSE )
    }
  }

  do.call(rbind, ans)
}

lipid_source <- bind_rows( summarise_lipid_source(lipid_base, "Original operational states"), summarise_lipid_source( lipid_base[lipid_base$HC_microglia_consistent, , drop = FALSE],
  "Expression-based HC restriction" ) )

write_csv( lipid_source, "10_146_lipid_program_by_source_original_vs_HC_v01.csv" )

make_lipid_comparison <- function(summary_df, cohort_name) {
  z <- summary_df[summary_df$cohort == cohort_name, , drop = FALSE]

  out <- lapply(score_cols, function(sc) {
    zz <- z[z$score == sc, , drop = FALSE]
    f <- zz[zz$myeloid_state == FOCUS, , drop = FALSE]
    h <- zz[zz$myeloid_state == HOME, , drop = FALSE]

    if (nrow(f) != 1L || nrow(h) != 1L) {
      return(data.frame( score = sc, cohort = cohort_name, n_focus = NA_integer_, n_homeostatic = NA_integer_, median_focus = NA_real_, median_homeostatic = NA_real_,
        median_difference_focus_minus_homeostatic = NA_real_, mean_difference_focus_minus_homeostatic = NA_real_, focus_median_gt_homeostatic = NA, stringsAsFactors = FALSE ))
    }

    data.frame( score = sc, cohort = cohort_name, n_focus = f$n_cells, n_homeostatic = h$n_cells, median_focus = f$median, median_homeostatic = h$median,
      median_difference_focus_minus_homeostatic = f$median - h$median, mean_difference_focus_minus_homeostatic = f$mean - h$mean, focus_median_gt_homeostatic =
      is.finite(f$median) && is.finite(h$median) && f$median > h$median, stringsAsFactors = FALSE )
  })

  do.call(rbind, out)
}

lipid_cmp_original <- make_lipid_comparison( lipid_summary, "Original operational states" )
lipid_cmp_hc <- make_lipid_comparison( lipid_summary, "Expression-based HC restriction" )

lipid_robustness <- merge( lipid_cmp_original, lipid_cmp_hc, by = "score", suffixes = c("_original", "_HC"), all = TRUE, sort = FALSE )

lipid_robustness$direction_preserved <- with( lipid_robustness, is.finite(median_difference_focus_minus_homeostatic_original) & is.finite(median_difference_focus_minus_homeostatic_HC) &
  sign(median_difference_focus_minus_homeostatic_original) == sign(median_difference_focus_minus_homeostatic_HC) )

lipid_robustness$primary_Figure2G_score <-
  lipid_robustness$score == PRIMARY_LIPID_SCORE

lipid_robustness$primary_result_preserved <- ifelse( lipid_robustness$primary_Figure2G_score, lipid_robustness$direction_preserved &
  lipid_robustness$median_difference_focus_minus_homeostatic_HC > 0, NA )

hc_src <- lipid_source[ lipid_source$cohort == "Expression-based HC restriction" & lipid_source$enough_cells_for_descriptive_direction, , drop = FALSE ]
source_dir <- hc_src %>%
  group_by(score) %>%
  summarise( n_informative_sources = n(), n_sources_focus_gt_homeostatic = sum(median_difference_focus_minus_homeostatic > 0, na.rm = TRUE), .groups = "drop" )

lipid_robustness <- merge( lipid_robustness, source_dir, by = "score", all.x = TRUE, sort = FALSE )

write_csv( lipid_robustness, "11_146_lipid_robustness_summary_v01.csv" )

primary_plot_df <- lipid_summary %>%
  filter(score == PRIMARY_LIPID_SCORE) %>%
  mutate( state_display = ifelse( myeloid_state == FOCUS, "Clec7a+ focal state", "Homeostatic microglia" ), cohort = factor( cohort, levels = c( "Original operational states",
    "Expression-based HC restriction" ) ) )

pC <- ggplot( primary_plot_df, aes( x = state_display, y = median, ymin = q25, ymax = q75, shape = state_display ) ) +
  geom_linerange(linewidth = 0.7) +
  geom_point(size = 2.7) +
  facet_wrap(~cohort) +
  labs( x = NULL, y = "FA uptake/transport AUCell score (median, IQR)", title = "Primary lipid-handling result before and after HC restriction",
    subtitle = "Precomputed Step20 AUCell scores; no expression matrix is reloaded", shape = NULL ) +
  theme_classic(base_size = 11) +
  theme( legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1) )

ggsave( file.path(FIG_DIR, "Fig_146C_FA_uptake_transport_original_vs_HC_v01.pdf"), pC, width = 7.6, height = 4.8 )
ggsave( file.path(FIG_DIR, "Fig_146C_FA_uptake_transport_original_vs_HC_v01.png"), pC, width = 7.6, height = 4.8, dpi = 300 )

n_all <- nrow(identity)
n_score_available <- sum(identity$score_available)
n_hc <- sum(identity$HC_microglia_consistent)
n_not_hc <- n_all - n_hc

n_original_mg <- sum(identity$original_MG_labeled)
n_original_mg_hc <- sum( identity$original_MG_labeled & identity$HC_microglia_consistent )
n_original_focus <- sum(identity$original_focus_state)
n_original_focus_hc <- sum( identity$original_focus_state & identity$HC_microglia_consistent )

temporal_hc_ok <- isTRUE(t_hc$late_enrichment_pattern)

primary_row <- lipid_robustness[ lipid_robustness$score == PRIMARY_LIPID_SCORE, , drop = FALSE ]
if (nrow(primary_row) != 1L) {
  stop("Primary lipid robustness row is not uniquely defined.")
}
lipid_primary_ok <- isTRUE(primary_row$primary_result_preserved)

final_qc <- data.frame( check = c( "Step106_and_Step20_cell_sets_identical", "Expected_total_myeloid_cells_58533", "All_three_identity_scores_available_for_at_least_one_cell",
  "HC_subset_nonempty", "HC_rule_applied_without_original_label_requirement", "Main_display_stages_present_in_HC_subset", "Clec7a_late_enrichment_pattern_preserved_in_HC_subset",
  "Primary_FA_uptake_transport_direction_preserved_in_HC_subset", "Overall_focused_sensitivity_supports_principal_results" ), value = c( TRUE, n_all == 58533L, n_score_available > 0,
  n_hc > 0, TRUE, all(MAIN_FIGURE_STAGES %in% unique(identity$timepoint_main[identity$HC_microglia_consistent])), temporal_hc_ok, lipid_primary_ok, temporal_hc_ok && lipid_primary_ok ),
  stringsAsFactors = FALSE )

write_csv(final_qc, "12_146_FINAL_QC_v01.csv")

report_lines <- c( "Reviewer 1 Remaining Major Comment 3 - focused sensitivity analysis", paste0("Script: ", SCRIPT_NAME), "", "Expression-based high-confidence criterion:",
  "  eligible_identity_scoring == TRUE;", "  finite MG_score, MDM_score, BAM_score;", "  MG_score >= 0;", "  MG_score > MDM_score;", "  MG_score > BAM_score.",
  "  The rule was applied to all frozen myeloid cells and did not require the original microglia label.", "", "Cell accounting:", paste0("  Total frozen myeloid cells: ", n_all),
  paste0("  Cells with available identity-program scores: ", n_score_available), paste0( "  Expression-based HC microglia-consistent cells: ",
  n_hc, " / ", n_all, " (", fmt_pct(100 * n_hc / n_all), ")" ), paste0( "  Cells not meeting the HC criterion: ", n_not_hc, " / ", n_all, " (", fmt_pct(100 * n_not_hc / n_all), ")" ),
  paste0( "  Original microglia-labeled cells retained by HC rule: ", n_original_mg_hc, " / ", n_original_mg, " (", fmt_pct(100 * n_original_mg_hc / n_original_mg), ")" ), paste0(
  "  Existing Clec7a+ focal-state cells retained by HC rule: ", n_original_focus_hc, " / ", n_original_focus, " (", fmt_pct(100 * n_original_focus_hc / n_original_focus), ")" ), "",
  "Temporal robustness:", paste0( "  HC Clec7a+ representation at 3, 7, 14, 28 dpi: ", fmt_pct(t_hc$pct_3dpi), ", ", fmt_pct(t_hc$pct_7dpi), ", ", fmt_pct(t_hc$pct_14dpi), ", ",
  fmt_pct(t_hc$pct_28dpi), "." ), paste0( "  Prespecified late-versus-early pattern preserved: ", temporal_hc_ok ), paste0(
  "  Descriptive six-stage Spearman rank concordance with manuscript reference: ", fmt_num(rank_concordance, 3) ), "", "Primary lipid-handling robustness:", paste0(
  "  Score: ", PRIMARY_LIPID_SCORE ), paste0( "  Original median difference (focal - homeostatic): ", fmt_num(primary_row$median_difference_focus_minus_homeostatic_original, 4) ), paste0(
  "  HC median difference (focal - homeostatic): ", fmt_num(primary_row$median_difference_focus_minus_homeostatic_HC, 4) ), paste0( "  Direction preserved and remains focal > homeostatic: ",
  lipid_primary_ok ), "", paste0( "Focused sensitivity conclusion: ", ifelse( temporal_hc_ok && lipid_primary_ok, paste0( "the principal late-stage enrichment and primary lipid-handling ",
  "transcriptional findings persist after restriction to the ", "expression-based high-confidence microglia-consistent subset." ), paste0(
  "one or more prespecified principal patterns were not preserved; ", "inspect the audit tables before drafting the response." ) ) ), "", "Interpretation limit:", paste0(
  "This is an expression-based robustness analysis. It does not establish ", "resident microglial lineage, and cells failing the criterion should not ", "be described as non-microglia." ) )

writeLines( report_lines, con = file.path(OUT, "13_146_REPORT_v01.txt"), useBytes = TRUE )
logit("Wrote 13_146_REPORT_v01.txt")

methods_lines <- c( "Step 146 methods and limitations", "", "Purpose", paste0( "This analysis addresses Reviewer 1 Remaining Major Comment 3 by testing ",
  "whether the principal temporal and lipid-handling findings persist after ", "restriction to an expression-based high-confidence microglia-consistent subset." ), "", "Inputs", paste0(
  "Only two existing project outputs are read: Step106 per-cell identity-program ", "scores (CSV.GZ) and Step20 per-cell AUCell scores (CSV). No ZIP, RDS, ",
  "expression matrix, Seurat object, pseudotime object, or differential-expression ", "result is loaded or regenerated." ), "", "High-confidence rule", paste0(
  "The rule was fixed before outcome comparison and was applied to every frozen ", "myeloid cell with available Step106 scores, irrespective of its original state label. ",
  "A cell was classified as expression-based high-confidence microglia-consistent ", "when MG_score >= 0, MG_score > MDM_score, and MG_score > BAM_score." ), "",
  "Why this is not lineage validation", paste0( "Step106 identity-program scores are transcriptional summaries. The high-confidence ",
  "classification therefore indicates relative expression consistency with the ", "microglial program under the stated rule; it does not establish resident origin, ",
  "exclude infiltrating-cell ancestry, or substitute for lineage tracing." ), "", "Temporal sensitivity", paste0(
  "The existing Clec7a+ focal-state label was retained unchanged. Its stage-wise ", "representation was recalculated within the expression-based HC subset and compared ",
  "descriptively with the manuscript-reference all-myeloid proportions. Because the ", "denominators differ by design, absolute percentages are not treated as directly ",
  "interchangeable. Robustness is defined qualitatively by preservation of the late ", "(14/28 dpi) enrichment relative to earlier injury stages (3/7 dpi)." ), "",
  "Lipid-handling sensitivity", paste0( "Precomputed Step20 AUCell values were matched by exact cell identifier. The existing ",
  "Clec7a+ focal state was compared descriptively with homeostatic microglia before and ", "after HC restriction. AUC_FA_uptake_transport was prespecified as the primary score ",
  "because it is the program displayed in revised Figure 2G. Other AUC programs are ", "reported as secondary descriptive robustness summaries." ), "", "Statistical interpretation", paste0(
  "No cell-level hypothesis-testing P values are generated in Step146 because individual ", "cells are not independent biological replicates. Medians, IQRs, source-specific ",
  "directions, temporal patterns, and effect-direction concordance are used for ", "sensitivity assessment." ), "", "Critical wording", paste0(
  "Cells that fail the expression-based high-confidence criterion should be described ", "only as not meeting the criterion. They should not be called non-microglia. ",
  "Likewise, persistence of the findings within the HC subset supports robustness to ", "an alternative expression-based confidence definition but does not independently ",
  "confirm lineage identity." ) )

writeLines( methods_lines, con = file.path(OUT, "14_146_methods_and_limitations_v01.txt"), useBytes = TRUE )
logit("Wrote 14_146_methods_and_limitations_v01.txt")

sink(file.path(OUT, "15_146_sessionInfo_v01.txt"))
cat("Step146 sessionInfo\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Script: ", SCRIPT_NAME, "\n", sep = "")
cat("No ZIP/RDS/expression matrix was read.\n\n")
print(sessionInfo())
sink()
logit("Wrote 15_146_sessionInfo_v01.txt")

all_outputs <- list.files( OUT, recursive = TRUE, full.names = TRUE, all.files = FALSE, include.dirs = FALSE )

all_outputs <- all_outputs[ basename(all_outputs) != "16_146_OUTPUT_INDEX_v01.csv" ]

rel_path <- function(x) {
  out_norm <- gsub("\\\\", "/", normalizePath(OUT, winslash = "/", mustWork = FALSE))
  x_norm <- gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
  prefix <- paste0(out_norm, "/")
  if (startsWith(tolower(x_norm), tolower(prefix))) {
    return(substr(x_norm, nchar(prefix) + 1L, nchar(x_norm)))
  }
  basename(x_norm)
}

logit("Step 146 completed.")
logit("Total cells: ", n_all)
logit("HC cells: ", n_hc, " (", fmt_pct(100 * n_hc / n_all), ")")
logit("Temporal late-enrichment pattern preserved: ", temporal_hc_ok)
logit("Primary lipid direction preserved: ", lipid_primary_ok)
logit( "Overall focused sensitivity support: ", temporal_hc_ok && lipid_primary_ok )
logit("Output directory: ", OUT)

all_outputs <- list.files( OUT, recursive = TRUE, full.names = TRUE, all.files = FALSE, include.dirs = FALSE )
all_outputs <- all_outputs[ basename(all_outputs) != "16_146_OUTPUT_INDEX_v01.csv" ]

output_index <- data.frame( relative_path = vapply(all_outputs, rel_path, character(1)), file_name = basename(all_outputs), extension = tolower(tools::file_ext(all_outputs)),
  size_bytes = file.info(all_outputs)$size, md5 = unname(tools::md5sum(all_outputs)), stringsAsFactors = FALSE )

output_index <- output_index[order(output_index$relative_path), , drop = FALSE]
utils::write.csv( output_index, file.path(OUT, "16_146_OUTPUT_INDEX_v01.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8" )

cat("\n============================================================\n")
cat("Step 146 completed\n")
cat("============================================================\n")
cat("Output: ", OUT, "\n", sep = "")
cat("No ZIP, RDS, expression matrix, Seurat object, or trajectory object was read.\n")
cat("See 13_146_REPORT_v01.txt and 12_146_FINAL_QC_v01.csv first.\n")

source("scripts/helpers/common.R", local = TRUE)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

ROOT <- "D:/zz"

DATA_DIR <- file.path(ROOT, "data", "GSE195783")

OUT_DIR <- file.path(ROOT, "results", "59_GSE195783_spatial_input_audit")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VERSION <- "v02"

MINIML_TGZ <- file.path(DATA_DIR, "GSE195783_family.xml.tgz")

RAW_TAR <- file.path(DATA_DIR, "GSE195783_RAW.tar")

EXPECTED_GSE <- "GSE195783"
EXPECTED_TIMEPOINTS <- c("3dpi", "7dpi", "14dpi", "28dpi")

TMP_DIR <- file.path(OUT_DIR, "_tmp_GSE195783_audit")

if (dir.exists(TMP_DIR)) {
    unlink(TMP_DIR, recursive = TRUE, force = TRUE)
}

dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

on.exit({
    if (dir.exists(TMP_DIR)) {
        unlink(TMP_DIR, recursive = TRUE, force = TRUE)
    }
}, add = TRUE)

required_pkgs <- c("xml2", "dplyr", "stringr", "tibble")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall them first, e.g.:\n", 
        "install.packages(c(", paste(sprintf("\"%s\"", missing_pkgs), collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(xml2)
  library(dplyr)
  library(stringr)
  library(tibble)
})

write_csv_utf8 <- write_csv_utf8

collapse_unique <- function(x, sep = " | ") {
    x <- x[!is.na(x) & nzchar(trimws(x))]
    x <- unique(trimws(x))
    paste(x, collapse = sep)
}

infer_timepoint <- function(text) {

  z <- tolower(
    paste(text, collapse = " ")
  )

  patterns <- list(
    "3dpi" = c(
      "time\\s*point\\s*:\\s*day\\s*3([^0-9]|$)",
      "(^|[^0-9])3\\s*dpi([^0-9]|$)",
      "(^|[^0-9])3\\s*days?\\s*(post[- ]?injury|after\\s*injury)([^0-9]|$)"
    ),
    "7dpi" = c(
      "time\\s*point\\s*:\\s*day\\s*7([^0-9]|$)",
      "(^|[^0-9])7\\s*dpi([^0-9]|$)",
      "(^|[^0-9])7\\s*days?\\s*(post[- ]?injury|after\\s*injury)([^0-9]|$)"
    ),
    "14dpi" = c(
      "time\\s*point\\s*:\\s*day\\s*14([^0-9]|$)",
      "(^|[^0-9])14\\s*dpi([^0-9]|$)",
      "(^|[^0-9])14\\s*days?\\s*(post[- ]?injury|after\\s*injury)([^0-9]|$)"
    ),
    "28dpi" = c(
      "time\\s*point\\s*:\\s*day\\s*28([^0-9]|$)",
      "(^|[^0-9])28\\s*dpi([^0-9]|$)",
      "(^|[^0-9])28\\s*days?\\s*(post[- ]?injury|after\\s*injury)([^0-9]|$)"
    )
  )

  candidates <- names(patterns)[
    vapply(
      patterns,
      function(pp) {
        any(
          vapply(
            pp,
            function(pat) {
              stringr::str_detect(
                z,
                stringr::regex(
                  pat,
                  ignore_case = TRUE
                )
              )
            },
            logical(1)
          )
        )
      },
      logical(1)
    )
  ]

  if (length(candidates) == 1L) {
    return(candidates)
  }

  NA_character_
}

infer_replicate <- function(text) {
    z <- paste(text, collapse = " ")
    x <- str_match(z, regex("(?:rep(?:licate)?|biological\\s*rep(?:licate)?)\\s*[_-]?\\s*([0-9]+)", ignore_case = TRUE))
    if (!is.na(x[1, 2])) {
        return(as.integer(x[1, 2]))
    }
    NA_integer_
}

detect_component <- function(files, pattern) {
    any(str_detect(tolower(files), regex(pattern, ignore_case = TRUE)))
}

input_audit <- data.frame(file = c(MINIML_TGZ, RAW_TAR), exists = file.exists(c(MINIML_TGZ, RAW_TAR)), 
    size_bytes = ifelse(file.exists(c(MINIML_TGZ, RAW_TAR)), file.info(c(MINIML_TGZ, RAW_TAR))$size, 
        NA_real_), stringsAsFactors = FALSE)

input_audit$size_MB <- input_audit$size_bytes / 1024^2

write_csv_utf8(input_audit, file.path(OUT_DIR, paste0("00_59_GSE195783_input_files_", VERSION, ".csv")))

if (!all(input_audit$exists)) {
    stop("Required input file(s) are missing.\n", paste(input_audit$file[!input_audit$exists], collapse = "\n"))
}

MINIML_DIR <- file.path(TMP_DIR, "miniml")

dir.create(MINIML_DIR, recursive = TRUE, showWarnings = FALSE)

utils::untar(MINIML_TGZ, exdir = MINIML_DIR)

xml_files <- list.files(MINIML_DIR, pattern = "\\.xml$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

if (length(xml_files) == 0L) {
    stop("No XML file was found after extracting the MINiML archive.")
}

if (length(xml_files) > 1L) {
    warning("More than one XML file found. Using the first one:\n", xml_files[1])
}

xml_file <- xml_files[1]

doc <- xml2::read_xml(xml_file)

sample_nodes <- xml2::xml_find_all(doc, ".//*[local-name()='Sample']")

if (length(sample_nodes) == 0L) {
    stop("No <Sample> node was detected in the GEO MINiML XML.")
}

parse_geo_sample <- function(node) {
    gsm_attr <- xml2::xml_attr(node, "iid")
    accession_nodes <- xml2::xml_find_all(node, ".//*[local-name()='Accession']")
    accessions <- xml2::xml_text(accession_nodes)
    gsm_text <- accessions[str_detect(accessions, "^GSM[0-9]+$")]
    gsm <- if (!is.na(gsm_attr) && str_detect(gsm_attr, "^GSM[0-9]+$")) {
        gsm_attr
    }
    else if (length(gsm_text) > 0L) {
        gsm_text[1]
    }
    else {
        NA_character_
    }
    title <- xml2::xml_text(xml2::xml_find_first(node, ".//*[local-name()='Title']"))
    source <- collapse_unique(xml2::xml_text(xml2::xml_find_all(node, ".//*[local-name()='Source']")))
    organism <- collapse_unique(xml2::xml_text(xml2::xml_find_all(node, ".//*[local-name()='Organism']")))
    characteristics_nodes <- xml2::xml_find_all(node, ".//*[local-name()='Characteristics']")
    characteristics <- character(0)
    if (length(characteristics_nodes) > 0L) {
        characteristics <- vapply(characteristics_nodes, function(x) {
            tag <- xml2::xml_attr(x, "tag")
            val <- trimws(xml2::xml_text(x))
            if (!is.na(tag) && nzchar(tag)) {
                paste0(tag, ": ", val)
            }
            else {
                val
            }
        }, character(1))
    }
    all_text <- collapse_unique(c(title, source, organism, characteristics, xml2::xml_text(xml2::xml_find_all(node, 
        ".//*[local-name()='Description']"))))
    data.frame(GSM = gsm, title = title, source = source, organism = organism, characteristics = collapse_unique(characteristics), 
        metadata_text = all_text, inferred_timepoint = infer_timepoint(all_text), inferred_replicate = infer_replicate(all_text), 
        stringsAsFactors = FALSE)
}

geo_samples <- do.call(rbind, lapply(sample_nodes, parse_geo_sample))

geo_samples <- geo_samples[!is.na(geo_samples$GSM), , drop = FALSE]

geo_samples <- geo_samples[!duplicated(geo_samples$GSM), , drop = FALSE]

geo_samples <- geo_samples %>% arrange(GSM)

write_csv_utf8(geo_samples, file.path(OUT_DIR, paste0("01_59_GSE195783_MINiML_samples_", VERSION, ".csv")))

tar_members <- utils::untar(RAW_TAR, list = TRUE)

if (length(tar_members) == 0L) {
    stop("GSE195783_RAW.tar contains no listed members.")
}

tar_df <- data.frame(tar_member = tar_members, basename = basename(tar_members), extension = tolower(tools::file_ext(tar_members)), 
    GSM = str_extract(basename(tar_members), "GSM[0-9]+"), stringsAsFactors = FALSE)

write_csv_utf8(tar_df, file.path(OUT_DIR, paste0("02_59_GSE195783_RAW_tar_members_", VERSION, ".csv")))

RAW_TMP_DIR <- file.path(TMP_DIR, "raw_tar")

dir.create(RAW_TMP_DIR, recursive = TRUE, showWarnings = FALSE)

utils::untar(RAW_TAR, exdir = RAW_TMP_DIR)

zip_files <- list.files(RAW_TMP_DIR, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

if (length(zip_files) == 0L) {
    warning(paste0("No ZIP files were found after extracting RAW.tar. ", "The script will inspect extracted files directly, but ", 
        "the expected GEO packaging is TAR-of-ZIP."))
}

zip_manifest_list <- list()
zip_contents_list <- list()

if (length(zip_files) > 0L) {

  for (i in seq_along(zip_files)) {

    zf <- zip_files[i]
    zbase <- basename(zf)
    gsm <- str_extract(
      zbase,
      "GSM[0-9]+"
    )

    zlist <- utils::unzip(
      zf,
      list = TRUE
    )

    zcontents <- as.character(
      zlist$Name
    )

    zcontents_norm <- gsub(
      "\\\\",
      "/",
      zcontents
    )

    zip_contents_list[[i]] <- data.frame(
      GSM = gsm,
      zip_file = zbase,
      zip_member = zcontents_norm,
      compressed_length = if (
        "CompressedLength" %in% colnames(zlist)
      ) {
        zlist$CompressedLength
      } else {
        NA_real_
      },
      uncompressed_length = if (
        "Length" %in% colnames(zlist)
      ) {
        zlist$Length
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )

    zip_manifest_list[[i]] <- data.frame(
      GSM = gsm,
      zip_file = zbase,
      n_members = length(zcontents_norm),

      has_filtered_feature_bc_matrix_h5 =
        detect_component(
          zcontents_norm,
          "filtered_feature_bc_matrix\\.h5$"
        ),

      has_matrix_mtx =
        detect_component(
          zcontents_norm,
          "(^|/)matrix\\.mtx(\\.gz)?$"
        ),

      has_barcodes =
        detect_component(
          zcontents_norm,
          "(^|/)barcodes\\.tsv(\\.gz)?$"
        ),

      has_features =
        detect_component(
          zcontents_norm,
          "(^|/)(features|genes)\\.tsv(\\.gz)?$"
        ),

      has_tissue_positions =
        detect_component(
          zcontents_norm,
          "tissue_positions(_list)?\\.csv$|tissue_positions\\.csv$"
        ),

      has_scalefactors =
        detect_component(
          zcontents_norm,
          "scalefactors_json\\.json$"
        ),

      has_hires_image =
        detect_component(
          zcontents_norm,
          "tissue_hires_image\\.(png|jpg|jpeg)$"
        ),

      has_lowres_image =
        detect_component(
          zcontents_norm,
          "tissue_lowres_image\\.(png|jpg|jpeg)$"
        ),

      has_spatial_directory =
        detect_component(
          zcontents_norm,
          "(^|/)spatial/"
        ),

      has_any_coordinate_like_file =
        detect_component(
          zcontents_norm,
          "coord|coordinate|position|spatial|tissue_pos|spot"
        ),

      has_any_image_like_file =
        detect_component(
          zcontents_norm,
          "\\.(png|jpg|jpeg|tif|tiff)$"
        ),

      stringsAsFactors = FALSE
    )
  }

  zip_manifest <- bind_rows(
    zip_manifest_list
  )

  zip_contents <- bind_rows(
    zip_contents_list
  )

} else {

  zip_manifest <- data.frame(
    GSM = character(0),
    zip_file = character(0),
    n_members = integer(0),
    stringsAsFactors = FALSE
  )

  zip_contents <- data.frame(
    GSM = character(0),
    zip_file = character(0),
    zip_member = character(0),
    compressed_length = numeric(0),
    uncompressed_length = numeric(0),
    stringsAsFactors = FALSE
  )
}

write_csv_utf8(zip_manifest, file.path(OUT_DIR, paste0("03_59_GSE195783_processed_zip_manifest_", VERSION, 
    ".csv")))

write_csv_utf8(zip_contents, file.path(OUT_DIR, paste0("04_59_GSE195783_processed_zip_contents_", VERSION, 
    ".csv")))

sample_manifest <- geo_samples %>% select(GSM, title, source, organism, characteristics, metadata_text, 
    inferred_timepoint, inferred_replicate)

if (nrow(zip_manifest) > 0L) {
    sample_manifest <- sample_manifest %>% left_join(zip_manifest, by = "GSM")
} else {
    sample_manifest$zip_file <- NA_character_
}

sample_manifest <- sample_manifest %>% mutate(timepoint_order = match(inferred_timepoint, EXPECTED_TIMEPOINTS)) %>% 
    arrange(timepoint_order, inferred_replicate, GSM) %>% select(-timepoint_order)

write_csv_utf8(sample_manifest, file.path(OUT_DIR, paste0("05_59_GSE195783_sample_manifest_", VERSION, 
    ".csv")))

timepoint_summary <- sample_manifest %>% count(inferred_timepoint, name = "n_samples")

expected_summary <- data.frame(inferred_timepoint = EXPECTED_TIMEPOINTS, stringsAsFactors = FALSE)

timepoint_summary <- expected_summary %>% left_join(timepoint_summary, by = "inferred_timepoint") %>% 
    mutate(n_samples = ifelse(is.na(n_samples), 0L, n_samples))

n_unresolved <- sum(is.na(sample_manifest$inferred_timepoint))

timepoint_summary <- bind_rows(timepoint_summary, data.frame(inferred_timepoint = "UNRESOLVED", n_samples = n_unresolved))

write_csv_utf8(timepoint_summary, file.path(OUT_DIR, paste0("06_59_GSE195783_timepoint_summary_", VERSION, 
    ".csv")))

qc <- data.frame(check = c("MINiML file exists", "RAW TAR exists", "GEO samples parsed", "Expected total GEO samples = 8", 
    "RAW TAR contains files", "Processed ZIP files detected", "Every GEO sample linked to a ZIP", "All 3/7/14/28 dpi resolved from GEO metadata", 
    "No unresolved GEO sample timepoint"), pass = c(file.exists(MINIML_TGZ), file.exists(RAW_TAR), nrow(geo_samples) > 
    0L, nrow(geo_samples) == 8L, length(tar_members) > 0L, length(zip_files) > 0L, if (nrow(zip_manifest) > 
    0L) {
    all(geo_samples$GSM %in% zip_manifest$GSM)
} else {
    FALSE
}, all(EXPECTED_TIMEPOINTS %in% sample_manifest$inferred_timepoint), n_unresolved == 0L), stringsAsFactors = FALSE)

write_csv_utf8(qc, file.path(OUT_DIR, paste0("07_59_GSE195783_audit_QC_", VERSION, ".csv")))

report_file <- file.path(OUT_DIR, paste0("08_59_GSE195783_audit_REPORT_", VERSION, ".txt"))

con <- file(report_file, open = "wt", encoding = "UTF-8")

writeLines(c("============================================================", "GSE195783 spatial input audit", 
    paste0("Version: ", VERSION), "============================================================", "", 
    paste0("Input directory: ", DATA_DIR), paste0("Output directory: ", OUT_DIR), "", "INPUT FILES", 
    paste0("  MINiML: ", MINIML_TGZ, " | ", round(file.info(MINIML_TGZ)$size/1024^2, 3), " MB"), paste0("  RAW TAR: ", 
        RAW_TAR, " | ", round(file.info(RAW_TAR)$size/1024^2, 3), " MB"), "", "GEO SAMPLE STRUCTURE", 
    paste0("  Parsed GEO samples: ", nrow(geo_samples)), paste0("  TAR members: ", length(tar_members)), 
    paste0("  ZIP files detected: ", length(zip_files)), paste0("  Unresolved timepoints: ", n_unresolved), 
    "", "TIMEPOINT COUNTS", paste(paste0(timepoint_summary$inferred_timepoint, ": ", timepoint_summary$n_samples), 
        collapse = "\n"), "", "QC", paste(paste0(ifelse(qc$pass, "[PASS] ", "[CHECK] "), qc$check), collapse = "\n"), 
    "", "IMPORTANT INTERPRETATION RULE", paste0("Timepoint labels in this audit are inferred only from explicit ", 
        "GEO metadata text. Ambiguous samples are left unresolved rather ", "than guessed."), "", "NEXT STEP", 
    paste0("Proceed to the next spatial step only after checking ", "05_59_GSE195783_sample_manifest_", 
        VERSION, ".csv and confirming the sample/timepoint mapping."), "", "============================================================"), 
    con = con)

close(con)

print(timepoint_summary)
print(qc)

if (n_unresolved > 0L) {
    warning(n_unresolved, " sample(s) have unresolved timepoints. ", "Sample manifest/report written; ", 
        "do not manually guess the mapping.")
}

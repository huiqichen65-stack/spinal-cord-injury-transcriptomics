# Project paths
PROJECT_ROOT <- Sys.getenv(
  "SCI_PROJECT_ROOT",
  unset = normalizePath(".", winslash = "/", mustWork = FALSE)
)

RESULTS_DIR <- file.path(PROJECT_ROOT, "results")

REFERENCE_DIR <- Sys.getenv(
  "SCI_REFERENCE_DIR",
  unset = file.path(PROJECT_ROOT, "reference")
)

PYSCENIC_PYTHON <- Sys.getenv(
  "PYSCENIC_PYTHON",
  unset = Sys.which("python")
)

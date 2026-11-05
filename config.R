# Central path configuration for reproducible execution.
# Run scripts from the repository root, or set the environment variables below.

PROJECT_ROOT <- Sys.getenv(
  "SCI_PROJECT_ROOT",
  unset = normalizePath(".", winslash = "/", mustWork = FALSE)
)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = FALSE)

REFERENCE_DIR <- Sys.getenv(
  "SCI_REFERENCE_DIR",
  unset = file.path(PROJECT_ROOT, "reference")
)
REFERENCE_DIR <- normalizePath(REFERENCE_DIR, winslash = "/", mustWork = FALSE)

PYSCENIC_PYTHON <- Sys.getenv("PYSCENIC_PYTHON", unset = Sys.which("python"))

FONT_REGULAR <- Sys.getenv("SCI_FONT_REGULAR", unset = "")
FONT_BOLD <- Sys.getenv("SCI_FONT_BOLD", unset = "")

# Preserve trailing slashes supplied by legacy paste0() calls while making
# all paths relative to the configurable project/reference roots.
project_path <- function(path = "") {
  if (!nzchar(path)) return(PROJECT_ROOT)
  paste0(sub("/+$", "", PROJECT_ROOT), "/", sub("^/+", "", path))
}

reference_path <- function(path = "") {
  if (!nzchar(path)) return(REFERENCE_DIR)
  paste0(sub("/+$", "", REFERENCE_DIR), "/", sub("^/+", "", path))
}

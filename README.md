# Clec7a+ myeloid states after spinal cord injury

Analysis code for **Integrated transcriptomic analyses characterize a late-stage Clec7a+ lipid-associated myeloid state after spinal cord injury**.

## Data and environment

Public datasets: GSE162610, GSE196928, Figshare 17702045, GSE195783 and GSE281205. Accessions and study roles are listed in `public_datasets.csv`; analysis settings, random seeds and software versions are recorded in the corresponding root-level CSV files.

Raw data and large intermediate objects are not included. Supply the inputs required by each script, including reference gene sets where used, and adapt the local data, Python environment and font paths to your machine. Manuscript Supplementary Tables contain the released tabular results.

## Path configuration

Scripts are intended to be run from the repository root. Local machine-specific paths are not hard-coded. By default, the repository root is used as the project root; alternative locations can be supplied through the following environment variables:

- `SCI_PROJECT_ROOT`: project/data root
- `SCI_REFERENCE_DIR`: external reference-data directory (defaults to `<project root>/reference`)
- `PYSCENIC_PYTHON`: Python executable used for pySCENIC (defaults to `python` on `PATH`)
- `SCI_FONT_REGULAR` / `SCI_FONT_BOLD`: optional font-file paths for figure scripts; if unset, the scripts fall back to the default sans-serif font

Large input and intermediate files are not included in the repository. Reference resources should be downloaded separately and placed under the configured reference directory where required by the relevant script.

## Analysis scripts

| Directory | Analysis |
| --- | --- |
| `scripts/atlas` | Input checks, QC, integration and major cell types |
| `scripts/myeloid` | Myeloid states, differential expression and trajectories |
| `scripts/bulk_wgcna` | Bulk RNA-seq and WGCNA |
| `scripts/spatial` | Visium, RCTD and spatial associations |
| `scripts/cellchat` | Cell communication |
| `scripts/regulatory` | pySCENIC and network perturbation |
| `scripts/external_cohort` | Independent GSE281205 validation |
| `scripts/sensitivity` | Sensitivity analyses |
| `scripts/figures` | Main and supplementary figures |
| `scripts/helpers` | Shared R utilities and Python programs |

Run scripts from the repository root, starting with the input check:

```sh
Rscript scripts/atlas/01_dataset_check.R
```

Step numbers identify the study workflow. Each script lists its input objects; complete those upstream steps before running the corresponding analysis or figure. Interpretation notes are collected in [ANALYSIS_NOTES.md](ANALYSIS_NOTES.md). Legacy prose-report files contain a pointer to these notes; QC tables, numerical outputs, session information and file paths remain available.

## Checks

```sh
Rscript scripts/helpers/verify_release.R
```

This checks R syntax, shared-function environments and synthetic spatial input alignment. It does not rerun the biological analyses.

The script organization draws on the published SCI companion repositories of [Milich et al. (2021)](https://github.com/JaeLeeLab/sci_scRNAseq) and [Brennan et al. (2022)](https://github.com/OSU-BMBL/Spinal-cord-scRNAseq).

# Environment

This document describes how to set up the computational environment for the analyses in this repository.

## Tested platform

- macOS on Apple Silicon (arm64)
- Python 3.12 via conda-forge (Miniforge)
- R 4.6.0, installed as a macOS framework from CRAN
- The scVI training code uses `accelerator="mps"` (Apple GPU). In its current form it requires an Apple Silicon Mac.

## Installation

### 1. Install R 4.6.0

Install R 4.6.0 for macOS (arm64) from CRAN. rpy2 is compiled against this R installation, so R must be installed before the Python environment is created.

Install the R packages used by the analysis:
```r
install.packages("BiocManager")
BiocManager::install(c("miQC", "scDblFinder", "scry"))
```

ALRA is installed from source. Record the commit hash used.

Record the installed versions with `packageVersion("<package>")`. They are documented in `environment_R.json`.

### 2. Set R_HOME and create the Python environment

rpy2 is built from source, so `R_HOME` must point to the active R installation before the environment is created:
```bash
cd /path/to/repository
export R_HOME="$(R RHOME)"
conda env create -f environment.yml -n cap
conda activate cap
```

### 3. Verify the installation
```bash
python -m pip check
python -c "import cap, torch, scanpy, scvi, rpy2.robjects as ro; print(cap.__file__, torch.__version__, scanpy.__version__, scvi.__version__, ro.r('R.version.string')[0])"
```

Expected: `pip check` reports no broken requirements, and the import line shows R 4.6.0.

To check which R framework rpy2 is linked against:
```bash
otool -L "$(python -c 'import site; print(site.getsitepackages()[0])')/_rinterface_cffi_api.abi3.so" | grep -i "R.framework"
```

The path must match the R version from step 1.

## Environment files

| File | Purpose |
| --- | --- |
| `environment.yml` | Installation specification (conda and pip packages) |
| `environment_full.yml` | Complete export of the analysis environment `conda_py312` (no build strings), for exact reproduction |
| `environment_R.json` | R versions of the packages used (miQC, scDblFinder, scry, ALRA) |

`requirements-lock.txt` is not published because it contains machine-specific paths.

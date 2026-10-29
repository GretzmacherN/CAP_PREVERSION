# Breast cancer single-cell and proteomic analysis

Code for the analysis of epithelial cell states, antigen expression and protein associations in breast cancer, with a normal-tissue protein reference analysis.

The single-cell workflow reanalyses filtered count matrices from GSE161529. Tumour-derived epithelium is defined by specimen origin; CNV-high status is a separate sensitivity definition, not a complete malignancy classifier. Antigen detection is defined from raw UMI counts. ALRA values are used only for exploratory diagnostics.

## Installation

Create the Python environment and register its notebook kernel:

```bash
# Run from the repository root.
conda env create --file environment.yml
conda activate cap
python -m pip check
python -m ipykernel install --user --name cap --display-name "Python (cap)"
```

The notebooks also require an external R installation with miQC, scDblFinder, scry and ALRA, accessed through rpy2. R and package versions, including the ALRA source commit, are recorded in [environment_R.json](environment_R.json). The Python YAML does not install these R dependencies; the JSON records their provenance and is not an installation lockfile.

The recorded environment is Apple Silicon macOS. scVI training uses `accelerator="mps"`; select a supported accelerator when running on other hardware. The full platform-specific environment snapshot is optional. A fresh installation has not yet been validated.

## Input data

Large input datasets are stored outside the repository. Set their locations before starting Jupyter or the notebook editor:

```bash
# Replace these example paths with the local data directories.
export CAP_PROJECT="/path/to/project-data"
export CAP_RAW_DATA="/path/to/filtered-count-matrices"
export CAP_PROTEOME="/path/to/protein-data"
```

`CAP_PROJECT` contains `Data/`, `metadata/`, `Adata_Saves/` and `Dataframe_Saves/`. The shared `features.tsv` must be one directory above `CAP_RAW_DATA`. `CAP_PROTEOME` contains the protein matrices and associated metadata.

The [input manifest](docs/INPUTS.md) lists the required files, source studies and expected paths. Frozen signatures and supporting barcode inputs are included in [data/](data/README.md). Keep [scripts/](scripts/) alongside the notebooks; it contains required analysis and plotting functions.

## Notebooks

Run 000 through 003c in the order below, from the repository root, using the `cap` kernel and a fresh kernel for each notebook.

| Notebook | Analysis | Main output |
|---|---|---|
| [000_Raw_Data_And_Cohort.ipynb](000_Raw_Data_And_Cohort.ipynb) | Library inventory and cohort selection | `metadata/samples.csv` |
| [001_Preprocessing.ipynb](001_Preprocessing.ipynb) | Per-library QC, doublet detection and preprocessing | `Adata_Saves/head/<sample_id>.h5ad` |
| [002_Celltype_ID.ipynb](002_Celltype_ID.ipynb) | Atlas integration and cell-type annotation | `Adata_Saves/epithelial.h5ad` |
| [003a_Epithelial_Integration_Malignancy_Publication.ipynb](003a_Epithelial_Integration_Malignancy_Publication.ipynb) | Epithelial integration and inferred CNV | `Adata_Saves/epithelial_tumour_scvi.h5ad` |
| [003b_Meta_Programmes_Publication.ipynb](003b_Meta_Programmes_Publication.ipynb) | Meta-programme scoring and cell-state comparisons | `Adata_Saves/epithelial_meta_programmes.h5ad` |
| [003c_Antigen_Coexpression_Publication.ipynb](003c_Antigen_Coexpression_Publication.ipynb) | Antigen detection, donor comparisons and sensitivity analyses | `Adata_Saves/epithelial_antigen_coexpression.h5ad` |
| [003d_HPA_Normal_Tissue_Liability.ipynb](003d_HPA_Normal_Tissue_Liability.ipynb) | Normal-tissue protein co-representation | `Dataframe_Saves/hpa_dvp_*.csv` |
| [Three_Proteome.ipynb](Three_Proteome.ipynb) | Cell-line and biopsy protein associations | `$CAP_PROTEOME/output/` |

Notebook 001 is a Papermill template run once per included library. Its `patient` parameter takes the full library `sample_id`, not the donor identifier. Sample-specific executed notebooks are saved in `papermill_runs/`.

003d and `Three_Proteome` are independent of the single-cell chain. The primary Anurag protein analysis uses baseline specimens only.

## Results

Single-cell and HPA tables are written to `$CAP_PROJECT/Dataframe_Saves/` and figures to `$CAP_PROJECT/Figures/`. Protein results are written to `$CAP_PROTEOME/output/`. Derived H5ADs and trained models remain in the external project-data directory.

The [figure notebooks](figure_notebooks/README.md) export figures from saved results without retraining models. Their plotting-only objects must not be used to overwrite analytical H5ADs.

Population definitions, expression layers and statistical conventions are documented in [docs/ANALYSIS.md](docs/ANALYSIS.md). Validation commands and execution records are described in [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

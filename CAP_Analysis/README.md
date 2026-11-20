# Breast cancer single-cell and proteomic analysis

This repository contains the analysis notebooks and shared Python helpers for breast-cancer single-cell and proteomic analyses. The workflow reanalyses filtered count matrices from GSE161529 and includes separate analyses of epithelial cell states, antigen expression, normal-tissue protein reference data, and published proteomics cohorts.

Tumour-derived epithelium is defined by specimen origin. CNV-high status is a sensitivity definition, not a complete malignancy classifier. Antigen detection is based on raw UMI counts; ALRA values are used only for exploratory diagnostics. Protein abundance associations do not demonstrate molecular interaction or same-cell co-expression.

## Repository layout and data

The repository root is this directory, identified by `pyproject.toml`. Keep large inputs and run outputs outside the Git repository under `CAP_DATA_ROOT`. The repository contains selected frozen reference files under `data/frozen/` and their checksum manifest under `data/manifests/`.

The path helper supports these configuration values:

```text
CAP_DATA_ROOT=C:/path/to/CAP_Data
CAP_RUN_ID=your_run_id
```

The same `CAP_RUN_ID` identifies all stages of one analysis run. Use a new ID for a new run. The path helper derives the input and output locations from `CAP_DATA_ROOT`; separate environment variables for individual data and output subdirectories are not needed.

Expected external inputs include:

- Filtered GSE161529 matrices: `$CAP_DATA_ROOT/raw/gse161529/raw_counts/`
- Shared feature file: `$CAP_DATA_ROOT/raw/gse161529/features.tsv`
- GENCODE annotation: `$CAP_DATA_ROOT/references/gencode/gencode.v50.annotation.gtf.gz`
- CellTypist model: `$CAP_DATA_ROOT/references/celltypist/model_from_pal.pkl`
- HPA reference archive: `$CAP_DATA_ROOT/references/hpa/dvp_cell_type_group_data.tsv.zip`
- Three-proteome supplements: `$CAP_DATA_ROOT/references/proteomics/`

The proteomics directory must contain `Lee_MOESM2_ESM.csv`, `Lehmann.csv`, `mmc2.xlsx`, `cd-22-0200_table_s3_suppst3.xlsx`, `cd-22-0200_table_s1_suppst1.csv`, `h.all.v2026.1.Hs.symbols.gmt`, and `lee_uniprot_gene_mapping.csv`.

Do not commit local `.env` files, raw data, or run outputs. Check that `.env` is ignored by Git before creating or committing local configuration; the standalone `.gitignore` must include that rule. The checked-in `.env.example` should contain placeholders only and should be kept aligned with the two supported settings above.

## Environment and setup

The reported Windows environment is named `cap-rerun`. From PowerShell, activate it, install this repository in editable mode, and register the notebook kernel:

```powershell
conda activate cap-rerun
Set-Location C:\path\to\CAP_Analysis
python -m pip install -e .
python -m ipykernel install --user --name cap-rerun --display-name "Python (cap-rerun)"
```

Start Jupyter or the IDE from inside the repository directory so the path helper can locate `pyproject.toml`. Restart the notebook kernel after changing the package or environment.

The checked-in `environment.yml` and `environment_R.json` record an Apple Silicon/macOS environment; they are not validated Windows lockfiles. The Windows R setup and current validation limits are described in [ENVIRONMENT.md](ENVIRONMENT.md). In particular, a fresh Windows installation and the complete R-dependent single-cell workflow have not yet been confirmed.

## Notebook workflow

The core single-cell chain is `000` → `001` → `002` → `003a` → `003b` → `003c` → `004`. Notebook `001_Papermill.ipynb` is the driver that runs `001_Preprocessing.ipynb` for included samples; the preprocessing notebook uses the full library `sample_id` parameter. Notebooks `003d_Normal_Tissue_Protein.ipynb` and `005_Three_Proteome.ipynb` are independent of that single-cell chain.

| Notebook | Stage | Purpose |
|---|---|---|
| [000_Raw_Data_And_Cohort.ipynb](notebooks/000_Raw_Data_And_Cohort.ipynb) | `00_cohort` | Inventory raw libraries and write the cohort table. |
| [001_Papermill.ipynb](notebooks/001_Papermill.ipynb) | `01_preprocessing` | Run per-sample preprocessing through Papermill. |
| [001_Preprocessing.ipynb](notebooks/001_Preprocessing.ipynb) | `01_preprocessing` | Preprocess one library; normally invoked by the Papermill driver. |
| [002_Celltype_ID.ipynb](notebooks/002_Celltype_ID.ipynb) | `02_cell_annotation` | Integrate the atlas and annotate cell types. |
| [003a_Epithelial_Integration.ipynb](notebooks/003a_Epithelial_Integration.ipynb) | `03a_epithelial_integration` | Integrate epithelial cells and perform the CNV sensitivity analysis. |
| [003b_Meta_Programmes.ipynb](notebooks/003b_Meta_Programmes.ipynb) | `03b_meta_programmes` | Score meta-programmes and compare cell states. |
| [003c_Antigen_Coexpression.ipynb](notebooks/003c_Antigen_Coexpression.ipynb) | `03c_antigen_coexpression` | Analyse antigen detection, donor comparisons, and sensitivities. |
| [003d_Normal_Tissue_Protein.ipynb](notebooks/003d_Normal_Tissue_Protein.ipynb) | `03d_normal_tissue_protein` | Analyse normal-tissue protein reference data. |
| [004_Programme_Summary.ipynb](notebooks/004_Programme_Summary.ipynb) | `04_programme_summary` | Summarise programme scores across epithelial clusters. |
| [005_Three_Proteome.ipynb](notebooks/005_Three_Proteome.ipynb) | `05_three_proteome` | Analyse prespecified protein-abundance associations across three proteomics cohorts. |

Use a fresh kernel for each notebook. Run the core chain in order; run `003d` and `005` separately when their inputs are available. See [VALIDATION.md](VALIDATION.md) for what has and has not been validated.

## Run outputs

Stage outputs are stored outside the repository under:

```text
$CAP_DATA_ROOT/runs/$CAP_RUN_ID/<stage>/
```

The stage directories are `00_cohort`, `01_preprocessing`, `02_cell_annotation`, `03a_epithelial_integration`, `03b_meta_programmes`, `03c_antigen_coexpression`, `03d_normal_tissue_protein`, `04_programme_summary`, and `05_three_proteome`. Depending on the stage, outputs include H5AD objects, models, tables, figures, logs, and executed notebooks. The three-proteome stage writes its tables, figures, `run_provenance.json`, and `input_checksums.json` in its stage directory.

**Rerun warning:** the current cohort notebook and full Papermill preprocessing run contain stage-reset/deletion logic. Do not rerun them against an existing run ID or treat an existing run as immutable until those code paths are corrected. Use a new `CAP_RUN_ID` for any test or rerun.

Figures should be saved through the shared helpers in [`src/cap/figures.py`](src/cap/figures.py); path and run configuration is implemented in [`src/cap/paths.py`](src/cap/paths.py), with provenance helpers in [`src/cap/provenance.py`](src/cap/provenance.py).

# Reproducibility and validation

## Execution records

Use the analytical notebooks at the release root. Notebook 001 is a Papermill template; its library-specific execution records belong in `papermill_runs/`. Retain the original executed notebooks and their session reports.

Each analytical stage writes a named output consumed by the next stage. Model directories, parameter dictionaries, cohort/eligibility tables and source-data files accompany the H5ADs. Modification timestamps alone are insufficient to identify a run.

The environment snapshots record the software installed when captured. They do not replace the version and parameter records attached to individual results. See [environment configuration](../ENVIRONMENT.md).

## Sensitivity analyses

The analytical notebooks contain the study-specific sensitivity cells. They can also be executed from the release directory:

```bash
python -m scripts.publication_targeted_checks
```

The command requires the current final epithelial object, programme tables, protein scores and relevant source inputs. It writes to `results/publication_checks/`. It does not overwrite the primary objects or result tables. Doublet-score filtering retains the original programme calls; alternative class ranking is a separate sensitivity.

The helper supports section-specific execution, for example:

```bash
python -m scripts.publication_targeted_checks --sections doublet protein
```

## Export and parity checks

```bash
python -m scripts.validate_figure_repairs
```

This validator checks notebook schema and Python syntax, selected nonblank image exports, HPA formula parity and selected primary-statistic parity. It requires the saved analytical tables and the generated validation artifacts under `results/repaired_figures/`, `results/hpa_updated/` and `results/publication_checks/`. For an existing development checkout, the validator also recognises those artifacts in the parent `results/` directory. These artifacts can be supplied separately; they need not be committed to Git. The validator is not a data-free test or a validation of every notebook result.

Numerical checks complement visual inspection of final assembled figures. Inspect scales, legends, population labels, denominators, missingness symbols and panel correspondence with the manuscript.

## Public release records

A reproducible release should include:

- A versioned code snapshot and documented notebook/working-directory configuration.
- Source accessions, input checksums, signature provenance and reference-model metadata.
- Python and R package records, plus a fresh-install setup test.
- Exact specimen, donor and cell denominators with the selected source-data tables.
- Versioned links/checksums for processed objects, models and executed artifacts stored outside Git.
- Author-approved license and citation metadata, with source-data reuse conditions recorded separately.

Outstanding provenance fields are identified in the [input manifest](INPUTS.md). Historical audit reports describe development checks; the current notebooks, input records and exported tables define the release.

## Repository scope

The computational publication workflow comprises 000–003d and `Three_Proteome`, with the required scoring, figure and validation helpers. Functional assay source records and biological replicate verification are separate from these computational notebooks.

Legacy notebooks, migration scripts, local manuscript drafts, download caches and bulk figure directories are not part of the required execution path. Excluding them from Git does not require deleting the local files. Choose the tracked representation of notebook outputs explicitly and retain executed copies outside any output-cleaning operation.

# Validation record

## Current CAP_Analysis status

This section describes the current project layout and the latest user-reported Windows result. It supersedes the older path guidance below, but does not replace or rewrite the historical numerical findings.

- The canonical repository is `CAP_Analysis/`; the external data root is configured as `CAP_DATA_ROOT`, and stages share one `CAP_RUN_ID` under `CAP_DATA_ROOT/runs/<CAP_RUN_ID>/`.
- Stage `05_three_proteome` is registered in `cap.paths` and was included in the same run as the existing stages. Notebook 005 completed successfully on the Windows `cap-rerun` environment, according to the user.
- The confirmed proteomics inputs are under `CAP_DATA_ROOT/references/proteomics/`. Stage 005 writes `tables/`, `figures/`, `run_provenance.json` and `input_checksums.json` in its stage directory.
- `.env.example` has been added. The local `.env` is not a repository input and must remain ignored; submitted documentation must not contain a machine-specific absolute data path.
- A fresh Windows installation and full end-to-end execution of the R-dependent 000–004 workflow under the new paths have not yet been confirmed. The duplicate-directory inventory and cleanup are also still pending.

## Historical validation record: 13 November 2026, previous layout

The sections below preserve the earlier validation record for notebooks 002–004 and its then-current data, results and review directories. It documents the prior macOS/MPS run and does **not** validate the current Windows environment or the new `CAP_Analysis` path migration. Paths such as `results/publication_final/` and `../outputs/` below are historical references, not current output instructions.

Updated on **13 November 2026** after execution of 002, 003a–d and 004. The
preparation checks below describe the 12 November source preparation. Earlier
saved executions were archived before those copies were prepared.

## Packaged result copies

On 13 November 2026, all 44 current SVG exports were archived in
`figures_svg.zip`. All seven H5AD objects, four scVI model files and matching
numerical outputs were copied into `saved_objects/`, retaining their stage
structure. The 192 copied result files and every archived SVG were checked
against their original SHA256 values. The ZIP integrity check passed.

The result copies are independent regular files, with no symbolic or hard
links. Source file identities, sizes and modification times were unchanged.
APFS may share storage until a copy is edited. The source notebooks and their
input/output paths were not changed. Packaging does not execute any analysis.
The full copy and figure inventory is `saved_objects/manifest.json`.

## Completed execution on 13 November 2026

| Notebook | Executed nonempty code cells | Saved errors | Completion time, Europe/Berlin |
| --- | --- | --- | --- |
| 002 | 247/247 | 0 | 01:19, including final exports and diagnostics |
| 003a | 386/386 | 0 | 03:04 |
| 003b | 553/553 | 0 | 05:08 |
| 003c | 397/397 | 0 | 05:33 |
| 003d | 45/45 | 0 | 09:28 |
| 004 | 16/16 | 0 | 09:49 |

The 003a–c notebooks ran sequentially in fresh `cap` kernels through the
installed Jupyter executor. Their actual counters and outputs were saved in
place after successful completion. The run used explicit MPS acceleration
for the atlas model and all three epithelial models. All four model-weight
files are present. 000/001 remain source copies for the selected existing
preprocessing outputs. 003d and 004 also completed in fresh `cap` kernels.

003d reproduces all six earlier CTLI scores exactly from the same HPA archive,
whose SHA256 is recorded in its output provenance. The 24-group coverage and
exploratory total-protein interpretation remain documented in that notebook.

004 generates Figure 1d from the completed 003b class scores and the current
epithelial cluster metadata. All 97,506 scored cells, 29 libraries and 28
donors match the primary eligibility table. Technical libraries and donors
receive equal weights at their respective aggregation steps. Its thresholds
follow the existing cluster-summary convention of ten cells inside and
outside a cluster per library and three contributing donors. Eleven clusters
meet the donor threshold. Five are shown with coverage markers.

The new panel exports all 112 cluster-by-class summaries, contributing
library and donor scores, and exact input hashes. The SVG retains editable
text and vector colour rendering. Its standalone layout was visually checked.

The new atlas contains 229,043 cells from 43 libraries and 41 individuals.
Its manually selected epithelial compartment contains 123,639 cells. The
primary programme analysis contains 97,506 tumour-derived cells from
29 libraries and 28 donor-condition units. SH0106 contributes 59 cells and
falls below the existing 100-cell library eligibility threshold. The full
epithelial object retains those cells. The new resolution-0.5 epithelial
clustering has 16 clusters.

Read-only checks verify identical ordered cells, source cluster labels,
manual annotations, primary inclusion flags and sample/patient/condition
values across the four new epithelial objects. A streaming comparison verifies
every raw count on all 27,270 shared genes, with original gene IDs aligned.
Counts are finite nonnegative integers. All five saved embeddings agree from
003a through 003c. All 100 programme fields agree between 003b and 003c.
Every SVG under the four completed stages has a PDF counterpart. This checks
export completeness. Final panel assembly still requires visual review.

The global programme test gives R² 0.167016 and P 0.0133. Omitting MH0135
gives P 0.0686 in the leave-one-donor-out check. All six primary Lim/MP13
contrasts retain positive TN differences and Holm-adjusted P below 0.05.
The primary ROR1–EGFR TN–normal exact donor comparison gives P 0.085618.
The conditional TN overlap test gives P 0.120131. The corresponding
manuscript values require updating for this new run.

The detailed comparison and check record is
`../outputs/publication_final_review_2026-11-13/`. These checks establish
continuity of the current completed chain from the explicitly selected
preprocessing objects. Figure 1d is now generated. Final figure/manuscript
reconciliation, data/model distribution and the final release manifest remain
pending. The review records the remaining collaborator checklist items.

## Preparation checks

- All seven notebook schemas and Python/IPython syntax are checked.
- A static symbol-table scan checks for unbound global names; this does not
  replace execution.
- The fresh 002-to-003a-to-003b-to-003c output/input paths agree exactly.
  003c also points to the matching 002 epithelial-inclusion table.
- All analysis stages use distinct folders under `results/publication_final/`.
  The selected new root is checked to be absent during preparation.
- Papermill parameter injection is tested using an isolated preprocessing
  destination. 001 derives output paths after the injected `run_dir` value.
- Actual stage guards are exercised against occupied scratch directories.
  The files in those directories must remain unchanged.
- An empty 002 manual-annotation mapping must stop execution before labels are
  applied. Earlier numerical cluster IDs are not treated as new identities.
- The fixed gene-policy, meta-programme workbook, Lim signatures and their
  source-mapping files retain their original checksums.
- Retained computational function bodies are compared with the archived
  source. Removed functions are presentation functions; scoring, statistical
  inference and sensitivity calculations are retained.
- Import/setup checks use separate fresh processes in the existing `cap`
  environment. Changed figure labels/layouts are reviewed in PNG previews
  using explicitly selected completed objects. Those previews are not final
  results from the forthcoming run.

Detailed checks and their actual outcomes are recorded in
`../outputs/publication_rerun_preparation_2026-11-12/`.

## Final checks before rerunning

The additional record is `../outputs/final_preflight_2026-11-12/`.

- All 43 selected post-QC H5ADs were checked. Every nonzero raw count is
  finite, positive and integer-valued; cell totals and detected-gene counts
  match their QC metadata. Library identities, gene order, original Ensembl
  IDs, singlet status and the gene-policy record agree. The total remains
  229,043 cells from 41 patients. No downstream labels or scVI embeddings
  occur in these inputs.
- Regenerating 000 reproduces every existing cohort column; its additional
  `files_present` audit column is true for every library.
- A complete 001 execution for `GSM4909292_HER2-MH0069` reproduces the same
  558 retained cells and raw counts. An explicit PCA seed of zero reproduces
  identical PCA coordinates and matches the installed Scanpy default.
- Seven isolated stage executions pass in fresh kernels. Downstream tests use
  25,758 cells from all 43 libraries, arbitrary fixture labels, two training
  epochs, 999 null signatures, reduced Monte Carlo budgets and at most 1,200
  silhouette evaluation cells. The 003a retry caps normal references at 50
  cells per library while processing all query libraries. These are execution
  tests, not estimates for the manuscript.
- All four fresh test model files reload and reproduce saved latent means
  on 128 cells each. Four saved epithelial handoffs retain matching ordered
  cells, identity metadata, gene IDs and every raw count on 24,577 common
  features.
- Focused checks cover Holm correction, enumerated overlap probabilities,
  exact donor permutations, equal weighting of technical libraries, and
  generic/named annotation cases. A separate check verifies that discarding
  inherited per-library PCA coordinates leaves counts and metadata unchanged.
  The corrected marker dendrogram uses scVI and exports SVG/PDF successfully.

The initial 003a test exposed an assumption that a generic `Epithelial` label
must exist. That diagnostic now handles its absence and the absence of named
subtypes. Annotation silhouettes are undefined when only one category exists.
An obsolete fixed cell count was removed from a comment, and the diagnostic
was split into smaller cells. No manual identity is inferred by these changes.

The primary 003c permutation test now includes tied absolute differences
within `1e-12` percentage points, consistently with its other exact tests.
A rational-arithmetic example validates the correction. The archived primary
P-value remains 0.07058086392451718. QC, signature and matched-null settings
remain unchanged in the publication sources.

The publication run uses explicit MPS acceleration in 002 and for all three
003a models. The 002 setting was changed from automatic selection after the
reduced execution checks, at the user's request. The full atlas run and the
three epithelial MPS models subsequently completed on 13 November, as recorded
above.

The 12 November prepared-source hashes are in `source_after.json` in the
final-check record. The current executed notebooks and main artifacts are
identified in `../outputs/publication_final_review_2026-11-13/artifact_manifest.json`.
Earlier preparation manifests describe their dated snapshots, with matching
source copies retained in that record's `source_snapshot/` directory.

## Scope of the edits

Path configuration, early overwrite guards, R discovery, kernel metadata,
figure labels/layout and duplicate presentation code were changed. 003a also
saves protected-feature and second-seed model weights before those model
objects are deleted. Its existing training settings are retained.

The misleading `miqc_mt_cutoff` field name is replaced with
`max_retained_pct_mt`, with conversion of the older selected input CSV headers
in 002. The underlying values and QC rules are unchanged. Marker colourbars
identify the actual scaling. Donor-test annotations are placed outside the
data region.

All quantitative programme and lineage analyses remain, including the
luminal-progenitor and Lim et al. stem-cell-associated signatures, coverage
audits, donor-level contrasts and sensitivity families. The primary
matched-null procedure still uses 10,000 null signatures. Removing duplicate
figures does not remove those calculations.

## Selected starting data and earlier evidence

The explicitly selected 000/001 outputs contain 43 libraries, 41 patients and
229,043 post-QC cells, with the intended MH0114 exclusion and Tum0554 retention.
The preprocessing object/QC checks and prior per-library executions remain
available. The downstream rerun starts from these documented upstream inputs.

The earlier four epithelial handoff files had matching ordered cell IDs,
sample/patient/condition metadata and every raw count on their common feature
set. These checks describe the archived completed chain, not a newly executed
chain under `publication_final/`.

The CellTypist gene-ID correction was independently checked before this
preparation: original predictions reproduced exactly, all 1,396 model genes
matched after correction, and 591/229,043 supporting labels changed. Cluster
dominant labels, manual annotations and cell inclusion were unchanged. The
completed new 002 used original `gene_id` values and saved its reference
predictions directly into its own annotated H5ADs.

## Required final-run checks

The completed-run checks above verify the selected cohort and the saved manual
labels and inclusion flags. The comparison record identifies the changed cell
counts and primary statistics. Final figure and manuscript reconciliation
should use the selected release's actual values. The current manual mapping
was reviewed against this run's markers and reference predictions.

The archived analysis had donor-sensitive global programme results and
non-significant primary TN-normal/conditional-overlap tests. Preserve those
interpretive limits when assessing the rerun. Transcript detection does not
establish accessible surface-protein co-display, and the Lim-associated
signatures do not demonstrate functional stemness or a state transition.
HPA-DVP includes no breast tissue and does not establish clinical safety.

Freeze the final executed notebooks and numerical outputs together. Record
input/output checksums, provide accessible data/model locations, and match
each manuscript panel and value to its source table. The CellTypist training
source and historical HPA archive metadata remain documentation limitations.

The 12 November preparation checks used the existing environment and reduced
downstream fixtures. The completed 13 November execution began from the
explicitly selected post-QC objects. The execution and handoff validation are
recorded above. A fresh environment installation and complete raw-data
reproduction remain separate verification tasks.

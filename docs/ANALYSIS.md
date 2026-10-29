# Analysis definitions

## Samples and statistical units

A library is identified by `sample`/`sample_id`; a donor is identified by `patient`. Multiple libraries from one donor are not assumed to be technical replicates. Notebook 000 records cohort membership, and later notebooks record analysis-specific eligibility.

Programme composition uses per-library fractions averaged within donor. Raw antigen prevalence in 003c pools eligible cells within donor and specimen-origin group. Donor-level comparisons therefore differ from pooled cell proportions; both require their own denominators.

## Quality control and integration

Notebook 001 applies miQC per library, with an override retaining cells below 10% mitochondrial counts subject to other QC rules and a fixed ≤20% fallback if the mixture fails. Additional limits are 500–40,000 UMI and ≤3% haemoglobin counts. scDblFinder-predicted doublets are removed. No primary minimum detected-gene cutoff is imposed; ≥500 genes is a downstream sensitivity.

Whole-atlas integration uses 2,000 deviance-ranked features and a 30-dimensional scVI representation. Epithelial integration uses 4,000 deviance-ranked features and 20 latent dimensions. Both use library identity as the batch covariate, raw counts as model input, two hidden layers and a negative-binomial likelihood. Primary epithelial Leiden resolution is 0.5. Parameter and sensitivity settings remain explicit in the notebooks.

## Epithelial populations

| Field | Definition |
|---|---|
| `tumour_derived` | Epithelial cells from TN, HER2+ or ER+ specimens |
| `cnv_high` | CNV burden strictly above the normal epithelial q95 |
| `tumour_derived_cnv_high` | Intersection of tumour-derived and CNV-high cells |
| `g1_cell` | Cells assigned G1 by the expression-based cell-cycle classifier |

Tumour-derived epithelium can contain non-malignant cells. CNV-high is a conservative burden-based subgroup, not a complete malignancy classifier. A G1 expression assignment is not direct measurement of quiescence or DNA copy number.

CNV inference is performed per library using normal-tissue epithelium as reference, excluding a query normal library from its own reference. infercnvpy uses a 100-gene window, 10-gene step and dynamic threshold 1.5. Burden averages chromosome RMS signals with weights proportional to inference-window count. The infercnvpy dynamic threshold and the normal-reference q95 selection are distinct parameters.

## Meta-programmes

Notebook 003b projects 41 published signatures using an adapted matched-null framework. Activity requires within-cell Holm-adjusted empirical P ≤0.05 against 10,000 expression-matched null signatures and detection of at least 20% of retained signature genes. Signature coverage and analysis eligibility are exported.

Thirteen biological classes are accompanied by `low_activity` and `mixed_active_but_ambiguous`. These are uncertainty categories, not cell types; they remain in composition denominators. Cells outside the scoring population are distinct from low-activity cells.

Cluster labels describe donor-centred relative enrichment. In the profile bubble plot, colour is mean donor-centred MP deviation and circle area is the fraction of contributing donors with positive deviation. Neither encoding is a fraction of individual cells expressing a programme.

## Antigen detection and imputation

Primary antigen detection is ≥1 raw UMI; double positivity requires both genes to meet the threshold. Thresholds of ≥2 and ≥3 UMI are sensitivities. Primary donor prevalence comparisons retain donors with at least 100 cells and include zero-event donors.

The exact donor permutation test compares mean donor prevalences. The conditional overlap test fixes gene-positive margins within donor, library and depth strata; it addresses a different question from the prevalence comparison.

ALRA is run per library for exploratory visualisation. Imputed values are not counts, positivity calls, measured surface abundance or input to the raw-count statistical tests. Empty-droplet-based ambient-RNA estimation is not possible from the provided filtered matrices.

## Protein scores

ROR1 and EGFR abundances are standardised separately within each cohort; their mean z-score forms the primary pair-abundance score. EMT uses the frozen Hallmark signature. The BCSC-associated composite uses CD44, ALDH1A1, ALDH1A3, ITGA6 and PROM1. These are specimen-level protein summaries, not direct plasticity or functional stemness measurements.

Primary correlations use two-sided Spearman permutation tests with 99,999 permutations, 10,000 specimen-bootstrap resamples and BH correction across nine cohort–association tests. Missing values are not imputed. Basal-only, subtype-adjusted and tumour-content-adjusted results are supporting sensitivities.

## HPA co-representation

For measured intensity `I(g,c)`, scaled abundance is `A(g,c) = sqrt(I(g,c)) / max(sqrt(I))`, using one maximum across all four targets and 24 groups. Local pair score is `min(A_i, A_j)`; CTLI sums this score across groups.

Missing values remain marked as unquantified and contribute zero only to the declared primary score. Coverage, transformations, group exclusions and missing-value scenarios are exported separately. CTLI describes exploratory total-protein co-representation, not same-cell surface co-display, calibrated AND-gate activation or clinical safety.

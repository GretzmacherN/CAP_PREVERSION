# Frozen analytical inputs

These files are read directly by the publication notebooks. Checksums identify the exact local inputs.

| File | Use | SHA-256 |
|---|---|---|
| `cancer_meta_programme_signatures.xlsx` | Forty-one published cancer meta-programme signatures used by 003a and 003b | `a980e87edf6b216f591e29452f85d16bdc87b510c9b5989f524096e1862ac377` |
| `source_tumour_cluster_barcodes.csv.gz` | Supporting source-cluster concordance in 003a | `26af400ae38a5e878873c1cb30e83b8447440518a8727fdf7a8808d7780ede5c` |
| `source_transfer_exclusions.csv` | Explicit sample exclusions from barcode transfer | Versioned with the code |

The signature workbook derives from the cancer meta-programmes published by Gavish et al. (2023), *Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours*.

The source barcode table was extracted from the released breast-atlas objects:

| Source object | Figshare file identifier |
|---|---|
| `SeuratObject_TNBCTum.rds` | 31546544 |
| `SeuratObject_HER2Tum.rds` | 31544936 |
| `SeuratObject_ERTotalTum.rds` | 31545167 |

These objects cover the source study's jointly analysed TN, HER2+ and ER+ tumour groups. Transfer requires an audited sample identity; recurring 10x barcode sequences alone do not establish a match. The ER-MH0001 mapping is excluded. This supporting annotation does not define the primary population in 003b or 003c.

For external matrices, annotations, reference models and protein inputs, see the [input manifest](../docs/INPUTS.md). Refer to each source resource's access and reuse terms.

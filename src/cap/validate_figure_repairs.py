"""Check repaired exports, notebook syntax and parity with the frozen analysis."""

from pathlib import Path
import json

import nbformat
import numpy as np
import pandas as pd
from PIL import Image, ImageStat

from .publication_figure_state import project_paths
from .publication_targeted_checks import condition_tests


def validate():
    project, repository = project_paths()
    for directory in (repository, repository/"figure_notebooks"):
        for path in directory.glob("*.ipynb"):
            nb = nbformat.read(path, as_version=4)
            nbformat.validate(nb)
            for index, cell in enumerate(nb.cells):
                if cell.cell_type == "code" and not any(line.lstrip().startswith(("%", "!")) for line in cell.source.splitlines()):
                    compile(cell.source, f"{path.name}:cell{index}", "exec")
    reference = repository/"results"
    if not (reference/"repaired_figures").is_dir():
        reference = repository.parent/"results"
    folder = reference/"repaired_figures"
    images = sorted(folder.glob("*.png"))
    required = {"epithelial_umap_condition.png", "epithelial_umap_sample.png",
                "epithelial_cluster_markers_top5.png", "infercnv_score_by_condition.png",
                "infercnv_cnv_high_threshold.png", "tumour_derived_cnv_high_fraction_by_sample.png",
                "antigen_expression_umap.png", "supplement_2_ror1_raw_and_imputed.png",
                "epithelial_umap_cnv_burden.png", "epithelial_umap_cnv_high.png"}
    if not required.issubset({path.name for path in images}):
        raise AssertionError("A required figure export is missing.")
    for path in images + [reference/"hpa_updated/figure5_hpa_liability.png"]:
        with Image.open(path) as image:
            image.thumbnail((128, 128))
            if max(ImageStat.Stat(image.convert("RGB")).stddev) < 1:
                raise AssertionError(f"Blank figure: {path.name}")
    hpa = pd.read_csv(reference/"hpa_updated/hpa_dvp_liability_score.csv", index_col=0).ctli
    expected = {"ROR1+EGFR": .2588060163522138, "ROR1+ERBB2": .23164866482420873,
                "ROR1+EPCAM": .23164866482420873, "EGFR+ERBB2": 2.544719600266159,
                "EGFR+EPCAM": 3.013639904740072, "ERBB2+EPCAM": 1.685624127444305}
    assert np.allclose(hpa.reindex(expected).values, list(expected.values()), atol=1e-12)
    primary = pd.read_csv(project/"Dataframe_Saves/meta_programme_condition_tests.csv").set_index("programme")
    profiles = pd.read_csv(project/"Dataframe_Saves/meta_programme_donor_profiles.csv")
    recomputed = condition_tests(profiles, list(primary.index), permutations=9999).set_index("programme")
    assert np.allclose(primary.p_permutation, recomputed.p_permutation.reindex(primary.index), atol=1e-12)
    assert np.allclose(primary.p_holm, recomputed.p_holm.reindex(primary.index), atol=1e-12)
    baseline = pd.read_csv(reference/"publication_checks/doublet_score_codetection_sensitivity.csv")
    baseline = baseline[(baseline.retained_score_quantile == 1) & (baseline.min_umi == 1)].iloc[0]
    assert np.isclose(baseline.exact_permutation_p, .18021155830753355)
    report = {"notebook_syntax_and_schema": "passed", "nonblank_pngs": len(images)+1,
              "hpa_formula_parity": "passed", "primary_programme_test_parity": "passed",
              "primary_donor_test_parity": "passed"}
    (reference/"figure_repair_validation.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    validate()

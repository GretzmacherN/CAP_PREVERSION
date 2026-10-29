"""Load saved plotting state without training models or reading full expression layers."""

from __future__ import annotations

import ast
import json
import os
import re
from pathlib import Path

import anndata as ad
from anndata.io import read_elem
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns


def project_paths(project=None):
    """CAP_PROJECT is the directory containing Adata_Saves and Dataframe_Saves."""
    if project is None:
        project = os.environ.get("CAP_PROJECT")
    if project is None:
        for parent in (Path.cwd().resolve(), *Path.cwd().resolve().parents):
            if (parent / "Adata_Saves").is_dir() and (parent / "Dataframe_Saves").is_dir():
                project = parent
                break
    if project is None:
        raise FileNotFoundError("Set CAP_PROJECT to the directory containing Adata_Saves.")
    project = Path(project).expanduser().resolve()
    repository = Path(__file__).resolve().parents[1]
    return project, repository


def protein_paths(project=None):
    """Locate external protein inputs without depending on the code-folder name."""
    project, repository = project_paths(project)
    configured = os.environ.get("CAP_PROTEOME")
    candidates = ([Path(configured).expanduser()] if configured else [
        project / "Data" / "proteomics",
        repository.parent / "bulk",
    ])
    for candidate in candidates:
        if (candidate / "data" / "mmc2.xlsx").is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "Set CAP_PROTEOME to the protein-data directory containing data/mmc2.xlsx."
    )


def read_plot_object(path, genes=None):
    """Read metadata and selected genes in CSR row chunks to limit memory use.

    This object is explicitly for plots. Its expression matrix contains only
    requested genes, and must never replace a complete analysis object.
    """
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        uns = read_elem(handle["uns"])
        if genes is None:
            genes = ["ROR1", "ROR2", "EGFR", "ERBB2", "EPCAM", "MKI67", "ACTB",
                     "KRT8", "KRT18", "KRT19", "KRT7", "KRT5", "KRT14", "KRT17", "MUC1"]
            ranked = uns.get("rank_genes_groups", {}).get("names")
            if ranked is not None:
                for group in ranked.dtype.names:
                    genes.extend(str(value) for value in ranked[group][:5])
        positions = np.flatnonzero(var.index.isin(genes))
        if len(positions) == 0:
            raise ValueError("No requested plotting genes are present.")
        layers = {}
        for name in ("counts", "lognorm"):
            group = handle[f"layers/{name}"]
            if group.attrs["encoding-type"] != "csr_matrix":
                raise ValueError(f"Expected a CSR layer: {name}")
            pointers = group["indptr"][:]
            blocks = []
            for start in range(0, len(obs), 2048):
                stop = min(start + 2048, len(obs))
                lo, hi = int(pointers[start]), int(pointers[stop])
                block = sp.csr_matrix((group["data"][lo:hi], group["indices"][lo:hi],
                                       pointers[start:stop + 1] - lo),
                                      shape=(stop - start, len(var)))
                blocks.append(block[:, positions])
            layers[name] = sp.vstack(blocks, format="csr")
        epi = ad.AnnData(X=layers["lognorm"].copy(), obs=obs, var=var.iloc[positions].copy())
        epi.layers.update(layers)
        epi.layers["log1p"] = epi.layers["lognorm"]
        epi.uns.update(uns)
        epi.uns["figure_only"] = True
        for name in ("X_umap",):
            epi.obsm[name] = handle[f"obsm/{name}"][:]
    exclusions = Path(__file__).resolve().parents[1] / "data/source_transfer_exclusions.csv"
    if exclusions.exists() and "source_tumour_cluster" in epi.obs:
        excluded_samples = pd.read_csv(exclusions)["sample"]
        mask = epi.obs["sample"].astype(str).isin(excluded_samples)
        epi.obs.loc[mask, "source_tumour_cluster"] = False
        epi.obs.loc[mask, "source_tumour_cluster_eligible"] = False
        epi.obs.loc[mask, "source_tumour_cluster_call"] = "not evaluated"
        epi.uns["figure_state_transfer_exclusions"] = excluded_samples.to_list()
    return epi


def _run_tagged_setup(notebook, namespace):
    """Only explicitly tagged, non-training setup cells are executed."""
    for cell in json.loads(Path(notebook).read_text())["cells"]:
        if "figure_setup" in cell.get("metadata", {}).get("tags", []):
            exec(compile("".join(cell["source"]), str(notebook), "exec"), namespace)


def load_figure_state(stage="003a", project=None, output_dir=None):
    """Return the namespace needed by the existing tagged figure blocks."""
    project, repository = project_paths(project)
    filenames = {"003a": "epithelial_tumour_scvi.h5ad",
                 "003b": "epithelial_meta_programmes.h5ad",
                 "003c": "epithelial_antigen_coexpression.h5ad"}
    if stage not in filenames:
        raise ValueError(f"Unknown stage {stage!r}; use 003a, 003b, or 003c.")
    tables = project / "Dataframe_Saves"
    figures = Path(output_dir or os.environ.get("CAP_FIGURE_OUTPUT", project / "Figures"))
    figures.mkdir(parents=True, exist_ok=True)
    epi = read_plot_object(project / "Adata_Saves" / filenames[stage])
    ns = dict(np=np, pd=pd, plt=plt, sc=sc, sp=sp, sns=sns, re=re,
              Line2D=Line2D, TwoSlopeNorm=TwoSlopeNorm, Path=Path,
              epi=epi, project=project, notebook_dir=repository,
              df_dir=tables, input_df_dir=tables, fig_dir=figures, save_dir=project / "Adata_Saves",
              random_seed=0, find_publication=lambda: repository,
              display=__import__("IPython.display", fromlist=["display"]).display)
    sns.set_theme(style="whitegrid", context="paper", palette="colorblind")
    sc.settings.figdir = figures
    np.random.seed(0)
    if stage == "003a":
        ns.update(cnv_adata=epi, cnv_high_threshold=float(
            epi.uns["epithelial_population_definitions"]["cnv_high_threshold"]))
        ns["sample_population_summary"] = pd.read_csv(tables / "epithelial_population_summary_by_sample.csv")
        ns["plot_data"] = ns["sample_population_summary"].query("condition != 'N'").copy()
    elif stage == "003b":
        mapping = {"unit_composition": "donor_profiles", "cluster_stability": "cluster_donor_stability",
                   "group_mp_summary": "cluster_enrichment_bubble", "margin_sensitivity": "top2_margin_sensitivity"}
        for variable, suffix in mapping.items():
            ns[variable] = pd.read_csv(tables / f"meta_programme_{suffix}.csv", dtype={"epi_leiden": str})
        ns["class_scores"] = pd.read_csv(tables / "meta_programme_class_scores.csv", index_col=0)
        ns["scored_cells"] = ns["class_scores"].index
        ns["scored_obs"] = epi.obs.loc[ns["scored_cells"]].copy()
        ns["cluster_unit_deviation"] = ns["class_scores"].join(
            ns["scored_obs"][["epi_leiden", "analysis_unit"]]
        ).groupby(["epi_leiden", "analysis_unit"], observed=True).mean()
        means = ns["class_scores"].join(ns["scored_obs"][["analysis_unit"]]).groupby("analysis_unit", observed=True).mean()
        ns["cluster_unit_deviation"] -= means.reindex(
            ns["cluster_unit_deviation"].index.get_level_values("analysis_unit")
        ).set_axis(ns["cluster_unit_deviation"].index)
        ns["meta_programme_states"] = {name: [] for name in epi.uns["meta_programme_mp_names"].values()}
        ns["mp_id_to_name"] = dict(epi.uns["meta_programme_mp_names"])
        ns["mp_global_z"] = pd.DataFrame(columns=list(ns["meta_programme_states"]))
        ns["composition_categories"] = list(epi.uns["meta_programme_hallmark_order"]) + [
            "low_activity", "mixed_active_but_ambiguous"]
    else:
        ns["donor_levels"] = pd.read_csv(tables / "coexpression_donor_prevalence.csv")
        ns["donor_inference"] = dict(epi.uns["coexpression_window"])
        ns["group_order"] = ["TN tumour-derived epithelium", "ER+ tumour-derived epithelium",
                             "HER2+ tumour-derived epithelium", "Normal-tissue epithelium"]
        primary = ns["donor_levels"].query("pair == 'ROR1+EGFR' and min_umi == 1")
        for label in ("TN tumour-derived epithelium", "Normal-tissue epithelium"):
            values = primary.loc[primary.group.eq(label), "prevalence"]
            if values.empty or not np.isfinite(values).all():
                raise ValueError(f"Missing or non-finite donor group: {label}")
        p_value = ns["donor_inference"]["exact_permutation_p"]
        if not np.isfinite(p_value) or not 0 <= p_value <= 1:
            raise ValueError("The saved donor-level result is invalid.")
    notebook = next(repository.glob(f"{stage}_*Publication.ipynb"))
    _run_tagged_setup(notebook, ns)
    # Presentation-only objects carry an explicit guard against accidental saving.
    return ns


def run_figure_blocks(stage, project=None, output_dir=None):
    """Run only tagged plot cells. No model, null scoring, or ALRA cell is run."""
    project, repository = project_paths(project)
    namespace = load_figure_state(stage, project, output_dir)
    source_tables = Path(namespace["fig_dir"]) / "source_tables"
    source_tables.mkdir(parents=True, exist_ok=True)
    notebook = next(repository.glob(f"{stage}_*Publication.ipynb"))
    for index, cell in enumerate(json.loads(notebook.read_text())["cells"]):
        if "publication_figure" in cell.get("metadata", {}).get("tags", []):
            namespace["df_dir"] = namespace["input_df_dir"] if "scvi_epithelial_history.csv" in "".join(cell["source"]) else source_tables
            exec(compile("".join(cell["source"]), f"{notebook.name}:cell{index}", "exec"), namespace)
            plt.close("all")
    return namespace


def save_cnv_panels(epi, output_dir):
    """Standalone CNV panels for figure assembly, using the saved embedding."""
    output_dir = Path(output_dir)
    fig = sc.pl.umap(epi, color="cnv_score", cmap="viridis", vmax="p99.9",
                     frameon=False, title="CNV burden", return_fig=True, show=False)
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"epithelial_umap_cnv_burden.{suffix}", dpi=300, bbox_inches="tight")
    plt.show()
    coordinates = epi.obsm["X_umap"]
    high = epi.obs["tumour_derived_cnv_high"].to_numpy(bool)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(*coordinates.T, s=1, color="#D0D0D0", linewidths=0, rasterized=True)
    ax.scatter(*coordinates[high].T, s=1, color="#D55E00", linewidths=0, rasterized=True)
    ax.set_title("CNV-high tumour-derived epithelium", loc="left")
    ax.set_axis_off()
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="none", color="#D55E00",
                               label="tumour-derived, above normal epithelial q95"),
                       Line2D([], [], marker="o", linestyle="none", color="#D0D0D0", label="other epithelial cells")],
              frameon=False, fontsize=8, loc="lower left")
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"epithelial_umap_cnv_high.{suffix}", dpi=300, bbox_inches="tight")
    plt.show()


def save_supplementary_ror1(epi, fig_dir, output_tables=None):
    """Export six raw/imputed diagnostics without rerunning imputation."""
    counts = {gene: np.asarray(epi[:, gene].layers["counts"].toarray()).ravel()
              for gene in ("ROR1", "EPCAM", "EGFR", "ERBB2")}
    condition_order = ["N", "TN", "HER2+", "ER+"]
    palette = {"N": "#031A6B", "ER+": "#F46036", "HER2+": "#2CA02C", "TN": "#976391"}
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    rates = pd.Series({gene: 100 * (values > 0).mean() for gene, values in counts.items()})
    axes[0, 0].bar(rates.index, rates.values, color="#4477AA")
    axes[0, 0].set(ylabel="cells with ≥1 raw UMI (%)", title="a  Raw transcript detection")
    raw = counts["ROR1"]
    comparison = epi.copy()
    comparison.obs["ROR1_log_normalised"] = np.asarray(epi[:, "ROR1"].layers["lognorm"].toarray()).ravel()
    sc.pl.umap(comparison, color="ROR1_log_normalised", ax=axes[0, 1], cmap="viridis",
               vmax="p99.9", frameon=False, show=False, title="b  ROR1, non-imputed")
    sc.pl.umap(epi, color="ROR1_alra", ax=axes[0, 2], cmap="viridis", vmax="p99",
               frameon=False, show=False, title="c  ROR1, ALRA (exploratory)")
    bins = np.arange(1, max(int(raw.max()), 3) + 2) - 0.5
    for condition in condition_order:
        mask = epi.obs["condition"].astype(str).eq(condition).to_numpy()
        axes[1, 0].hist(raw[mask & (raw > 0)], bins=bins, histtype="step",
                        label=condition, color=palette[condition])
    axes[1, 0].set(xlabel="ROR1 UMI among detecting cells", ylabel="cells", title="d  Raw UMI distribution")
    axes[1, 0].legend(frameon=False)
    rows = []
    for condition in condition_order:
        mask = epi.obs["condition"].astype(str).eq(condition).to_numpy()
        rows.append({"condition": condition, "n_cells": int(mask.sum()),
                     "zero": int((raw[mask] == 0).sum()), "one": int((raw[mask] == 1).sum()),
                     "two_plus": int((raw[mask] >= 2).sum())})
    distribution = pd.DataFrame(rows).set_index("condition")
    distribution[["zero", "one", "two_plus"]].div(distribution.n_cells, axis=0).plot.bar(
        stacked=True, ax=axes[1, 1], color=["#D9D9D9", "#56B4E9", "#0072B2"])
    axes[1, 1].set(ylabel="fraction of all epithelial cells", xlabel="", title="e  Raw detection categories")
    diagnostic = pd.DataFrame({"condition": epi.obs.condition.astype(str),
                               "ALRA value": epi.obs["ROR1_alra"].to_numpy()})
    sns.boxplot(data=diagnostic, x="condition", y="ALRA value", order=condition_order,
                showfliers=False, color="#BDD7E7", ax=axes[1, 2])
    axes[1, 2].set(xlabel="", title="f  ALRA distribution (diagnostic only)")
    for suffix in ("png", "pdf"):
        fig.savefig(Path(fig_dir) / f"supplement_2_ror1_raw_and_imputed.{suffix}", dpi=300, bbox_inches="tight")
    if output_tables is not None:
        distribution.to_csv(Path(output_tables) / "supplement_2_ror1_raw_detection.csv")
    plt.show()
    return distribution

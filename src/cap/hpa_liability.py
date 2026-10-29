"""Exploratory total-protein co-representation from a frozen HPA-DVP snapshot."""

from __future__ import annotations

import hashlib
import itertools
import json
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


gene_order = ["ROR1", "EGFR", "ERBB2", "EPCAM"]
gene_colors = {"ROR1": "#d62728", "EGFR": "#9467bd", "ERBB2": "#2ca02c", "EPCAM": "#1f77b4"}
organ_groups = {
    "Blood / Immune": ["b-cells", "t-cells", "macrophages", "neutrophils"],
    "Brain": ["neurons", "astrocytes and neuropil", "microglia and neuropil"],
    "Heart": ["cardiomyocytes", "capillaries"],
    "Kidney": ["glomeruli", "proximal tubules", "distal tubules", "collecting ducts"],
    "Liver": ["hepatocytes"], "Lung": ["alveolar cell types"],
    "Ovary": ["granulosa cells", "oocytes"],
    "Pancreas": ["pancreatic islets", "pancreas epithelial cells"],
    "Salivary gland": ["large ducts"], "Skin": ["keratinocytes"],
    "Skeletal muscle": ["skeletal myofibers"], "Fallopian tube": ["fallopian tube cells"],
    "Vasculature": ["smooth muscle cells"],
}


def load_hpa(path):
    """Keep non-detection as NaN; the zero convention is applied only in scoring."""
    path = Path(path)
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.endswith(".tsv")]
        if len(members) != 1:
            raise ValueError("Expected exactly one TSV in the cached HPA archive.")
        with archive.open(members[0]) as handle:
            table = pd.read_csv(handle, sep="\t")
    targets = table.loc[table["Gene name"].isin(gene_order)].copy()
    if targets.duplicated(["Cell type name", "Gene name"]).any():
        raise ValueError("Duplicate target/cell-type measurements need an explicit aggregation rule.")
    abundance = targets.pivot(index="Cell type name", columns="Gene name", values="Intensity").reindex(columns=gene_order)
    expected = [cell for members in organ_groups.values() for cell in members]
    if set(expected) != set(abundance.index):
        raise ValueError(f"HPA cell-type mapping changed: {set(expected) ^ set(abundance.index)}")
    abundance = abundance.reindex(expected)
    if abundance.lt(0).any().any() or abundance.notna().sum().eq(0).any():
        raise ValueError("Invalid intensities or a wholly unquantified target.")
    provenance = {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "formula": "A=sqrt(I)/max(sqrt(I)); L=min(A_i,A_j); CTLI=sum_c L",
                  "missing_policy": "zero contribution to the primary score; biological abundance remains unknown",
                  "n_groups": len(abundance), "source_url": "https://www.proteinatlas.org/download/tsv/dvp_cell_type_group_data.tsv.zip",
                  "source_version": "not encoded in cached archive; verify before deposit",
                  "interpretation": "exploratory total-protein co-representation; not same-cell surface co-display or safety"}
    return abundance, provenance


def score_hpa(abundance, transform="sqrt", missing_value=0.0):
    """Use one global denominator for all targets and retain the measured mask."""
    transforms = {"sqrt": np.sqrt, "linear": lambda x: x, "log1p": np.log1p}
    transformed = transforms[transform](abundance)
    denominator = float(np.nanmax(transformed.to_numpy()))
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("The common intensity scale must be positive and finite.")
    filled = transforms[transform](abundance.fillna(missing_value)) / denominator
    local = pd.DataFrame({f"{a}+{b}": filled[[a, b]].min(axis=1)
                          for a, b in itertools.combinations(gene_order, 2)})
    observed = pd.DataFrame({f"{a}+{b}": abundance[[a, b]].notna().all(axis=1)
                             for a, b in itertools.combinations(gene_order, 2)})
    return transformed / denominator, local, observed


def hpa_summaries(abundance):
    scaled, local, observed = score_hpa(abundance)
    summary = pd.DataFrame({"ctli_sum": local.sum(), "ctli_max": local.max(),
                            "top_cell_type": local.idxmax(), "n_joint_quantified": observed.sum(),
                            "n_groups": len(local),
                            "joint_quantified_mean": local.where(observed).mean()}).sort_values("ctli_sum")
    summary["top_share_of_sum"] = summary.ctli_max.div(summary.ctli_sum.replace(0, np.nan))
    summary["rank"] = summary.ctli_sum.rank(method="min").astype(int)
    transforms = []
    for transform in ("linear", "sqrt", "log1p"):
        _, values, _ = score_hpa(abundance, transform)
        for pair, total in values.sum().items():
            transforms.append({"transform": transform, "pair": pair, "ctli_sum": total,
                               "rank": values.sum().rank(method="min")[pair]})
    omitted = []
    for cell_type in local.index:
        totals = local.drop(index=cell_type).sum()
        for pair, total in totals.items():
            omitted.append({"dropped": cell_type, "pair": pair, "ctli_sum": total,
                            "rank": totals.rank(method="min")[pair]})
    # These are declared perturbation scenarios, not estimated detection limits.
    positive = abundance.to_numpy()[np.isfinite(abundance.to_numpy()) & (abundance.to_numpy() > 0)]
    missing_rows = []
    for name, fill in (("zero contribution", 0.0), ("pooled detected q10", float(np.quantile(positive, .10))),
                       ("pooled detected median", float(np.median(positive)))):
        _, values, _ = score_hpa(abundance, missing_value=fill)
        totals = values.sum()
        for pair, total in totals.items():
            missing_rows.append({"scenario": name, "assumed_missing_intensity": fill, "pair": pair,
                                 "ctli_sum": total, "rank": totals.rank(method="min")[pair]})
    coverage = pd.DataFrame({"n_quantified": abundance.notna().sum(),
                             "n_unquantified": abundance.isna().sum(), "n_groups": len(abundance)})
    return dict(scaled=scaled, local=local, observed=observed, summary=summary,
                coverage=coverage, transforms=pd.DataFrame(transforms),
                leave_one_out=pd.DataFrame(omitted), missing_scenarios=pd.DataFrame(missing_rows))


def export_hpa(abundance, results, provenance, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    abundance.to_csv(output_dir / "hpa_dvp_raw_intensities.csv")
    names = {"local": "hpa_dvp_coexpression_by_celltype", "summary": "hpa_dvp_liability_summary",
             "coverage": "hpa_dvp_target_coverage", "transforms": "hpa_dvp_transform_sensitivity",
             "leave_one_out": "hpa_dvp_leave_one_group_out", "missing_scenarios": "hpa_dvp_missingness_scenarios"}
    for key, name in names.items():
        results[key].to_csv(output_dir / f"{name}.csv", index=key in ("local", "summary", "coverage"))
    results["summary"].ctli_sum.rename("ctli").to_csv(output_dir / "hpa_dvp_liability_score.csv")
    (output_dir / "hpa_dvp_input_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def draw_hpa(abundance, results, output_dir):
    """One cohesive figure; missing measurements are visible as crosses."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scaled, local, observed = (results[key] for key in ("scaled", "local", "observed"))
    fig, axes = plt.subplots(1, 3, figsize=(17, 12), gridspec_kw={"width_ratios": [1.5, 1.25, 1]}, layout="constrained")
    positions, labels, headers = {}, [], []
    cursor = 0
    for organ, members in organ_groups.items():
        headers.append((cursor, organ))
        labels.append("")
        cursor += 1
        for cell in members:
            positions[cell] = cursor
            labels.append(cell)
            cursor += 1
    for column, gene in enumerate(gene_order):
        for cell, y in positions.items():
            value = scaled.loc[cell, gene]
            if pd.isna(value):
                axes[0].scatter(column, y, marker="x", s=14, color="#AAAAAA", linewidths=.7)
            elif value > 0:
                axes[0].scatter(column, y, s=400 * value, color=gene_colors[gene], edgecolors="none")
    axes[0].set(xticks=range(4), xticklabels=gene_order, yticks=range(cursor), yticklabels=labels,
                ylim=(cursor-.3, -.8), title="a  Measured total-protein intensity",
                xlabel="square-root intensity / common global maximum")
    for y, organ in headers:
        axes[0].text(-.57, y, organ, transform=axes[0].get_yaxis_transform(), fontweight="bold", va="center")
    points = None
    for column, pair in enumerate(local):
        for cell, y in positions.items():
            if not observed.loc[cell, pair]:
                axes[1].scatter(column, y, marker="x", s=12, color="#AAAAAA", linewidths=.6)
            else:
                value = local.loc[cell, pair]
                points = axes[1].scatter(column, y, s=400 * value, c=[value], cmap="viridis", vmin=0, vmax=1)
    axes[1].set(xticks=range(len(local.columns)), xticklabels=[x.replace("+", "–") for x in local],
                yticks=[], ylim=(cursor-.3, -.8), title="b  Limiting-antigen score")
    axes[1].tick_params(axis="x", rotation=55)
    fig.colorbar(points, ax=axes[1], orientation="horizontal", pad=.09, fraction=.035, label="minimum scaled intensity")
    totals = results["summary"].ctli_sum.sort_values(ascending=False)
    axes[2].barh(range(len(totals)), totals, color=["#1f77b4" if name=="ROR1+EGFR" else "#CCCCCC" for name in totals.index])
    axes[2].set(yticks=range(len(totals)), yticklabels=[x.replace("+", "–") for x in totals.index],
                xlabel="CTLI (zero contribution for unquantified values)", title="c  Exploratory co-representation")
    axes[2].invert_yaxis()
    axes[2].set_xlim(0, totals.max()*1.20)
    for y, value in enumerate(totals):
        axes[2].text(value + totals.max()*.025, y, f"{value:.2f}", va="center")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(False)
    size_handles = [Line2D([], [], marker="o", linestyle="none", color="#777777",
                          markersize=np.sqrt(400*value), label=f"{value:.2f}")
                    for value in (.25, .50, .75, 1.0)]
    size_handles.append(Line2D([], [], marker="x", linestyle="none", color="#AAAAAA", label="not quantified"))
    axes[0].legend(handles=size_handles, title="scaled intensity (bubble area)", ncol=3,
                   loc="upper left", bbox_to_anchor=(0, -.06), frameon=False, fontsize=8)
    fig.suptitle("HPA-DVP: 24 pooled cell-type groups; no breast tissue; one donor\nTotal protein, not surface antigen density or a validated safety estimate", fontsize=11)
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"figure5_hpa_liability.{suffix}", dpi=300, bbox_inches="tight")
    plt.show()
    return fig

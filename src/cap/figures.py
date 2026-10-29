# src/cap/figures.py
"""Shared figure style, palettes and output for all notebooks."""
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

from cap.paths import figures_dir

CONDITION_ORDER = ["N", "TN", "HER2+", "ER+"]
CONDITION_PALETTE = {"N": "#031A6B", "TN": "#976391", "HER2+": "#2CA02C", "ER+": "#F46036"}
CELLTYPE_PALETTE = {
    # immune
    "T Cells":                        "#2166AC",
    "T/NK Cells":                       "#67A9CF",
    "B Cells":                        "#FFD92F",
    "Plasma Cells":                   "#FF7F00",
    "Macrophages":                    "#33A02C",
    "Dendritic Cells":                "#6A3D9A",
    "Mast Cells":                     "#8C510A",
    # stromal / vascular
    "Fibroblasts":                    "#D73027",
    "Perivascular":                   "#CAB2D6",
    "Endothelial":                    "#F781BF",
    "Endothelial: Lymphatic":         "#FCC5C0",
    # epithelial
    "Epithelial":                     "#00BFC4",
    "Epithelial: Basal":              "#00441B",
    "Epithelial: Luminal Progenitor": "#A6CEE3",
    "Epithelial: Mature Luminal":     "#1F9E89",
    "Epithelial: Ciliated":           "#7FCDBB",
    "Epithelial: Cycling":            "#B2182B",
    # descriptive
    "ribo-high":                 "#CCCCCC",
}
# CELLTYPE_PALETTE: copy here the single version that the final figures use (see below)


RASTER_MIN_POINTS = 1000   # scatter collections with at least this many points become images in SVG/PDF
def set_style(font_size=10, dpi=300):
    """Call once per notebook, in the setup cell."""
    sns.set_theme(style="whitegrid", context="paper", palette="colorblind")
    plt.rcParams.update({
        "font.size": font_size, "axes.titlesize": font_size,
        "axes.labelsize": font_size, "xtick.labelsize": font_size, "ytick.labelsize": font_size,
        "figure.dpi": 150, "savefig.dpi": dpi,
        # fonts: first available is used; Liberation Sans is metric-compatible with Arial
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        # text stays as editable text in SVG, not as outlines
        "svg.fonttype": "none",
        # TrueType fonts in PDF, editable in Illustrator
        "pdf.fonttype": 42,
        # reproducible SVG ids, so re-running produces identical files
        "svg.hashsalt": "cap",
    })
    sc.settings.figdir = figures_dir()


def _rasterize_dense_artists(fig, min_points=RASTER_MIN_POINTS):
    """Mark large scatter collections as rasterized, so vector formats embed them as images."""
    for ax in fig.axes:
        for collection in ax.collections:
            if hasattr(collection, "get_offsets") and len(collection.get_offsets()) >= min_points:
                collection.set_rasterized(True)


def save_figure(fig, name, formats=("png", "svg"), dpi=300, raster_dpi=200):
    """Write one figure to results/figures/ in every requested format.

    name: file name without extension.
    dpi: resolution for PNG; raster_dpi: resolution of embedded point clouds in PDF/SVG.
    """
    out = figures_dir()
    _rasterize_dense_artists(fig)
    paths = []
    for ext in formats:
        path = out / f"{name}.{ext}"
        fig.savefig(path, dpi=raster_dpi if ext in ("svg", "pdf") else dpi,
                    bbox_inches="tight")
        if ext == "svg" and path.stat().st_size > 20e6:
            print(f"warning: {path.name} is {path.stat().st_size / 1e6:.0f} MB; "
                  "consider fewer points or a lower raster_dpi")
        paths.append(path)
    return paths

# add to src/cap/figures.py (seaborn import needed at the top: import seaborn as sns)
def categorical_palette(categories):
    """One distinct colour per category. Call with the same categories to get the same colours."""
    cats = list(categories)
    return dict(zip(cats, sns.color_palette("husl", n_colors=len(cats))))


def sample_order_by_condition(obs):
    """Sample IDs grouped by CONDITION_ORDER, then sorted by name."""
    rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
    table = obs[["sample", "condition"]].drop_duplicates()
    return (table.assign(rank=table["condition"].map(rank))
            .sort_values(["rank", "sample"])["sample"].tolist())
"""Shared figure style, palettes, and portable figure exports for CAP."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Literal
from warnings import warn

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.figure import Figure

FigureFormat = Literal["png", "svg"]

CONDITION_ORDER = ["N", "TN", "HER2+", "ER+"]
CONDITION_PALETTE = {
    "N": "#031A6B",
    "TN": "#976391",
    "HER2+": "#2CA02C",
    "ER+": "#F46036",
}

CELLTYPE_PALETTE = {
    # Immune
    "T Cells": "#2166AC",
    "T/NK Cells": "#67A9CF",
    "B Cells": "#FFD92F",
    "Plasma Cells": "#FF7F00",
    "Macrophages": "#33A02C",
    "Dendritic Cells": "#6A3D9A",
    "Mast Cells": "#8C510A",
    # Stromal / vascular
    "Fibroblasts": "#D73027",
    "Perivascular": "#CAB2D6",
    "Endothelial": "#F781BF",
    "Endothelial: Lymphatic": "#FCC5C0",
    # Epithelial
    "Epithelial": "#00BFC4",
    "Epithelial: Basal": "#00441B",
    "Epithelial: Luminal Progenitor": "#A6CEE3",
    "Epithelial: Mature Luminal": "#1F9E89",
    "Epithelial: Ciliated": "#7FCDBB",
    "Epithelial: Cycling": "#B2182B",
    # Descriptive
    "ribo-high": "#CCCCCC",
}

VALID_FORMATS = frozenset({"png", "svg"})
RASTER_MIN_POINTS = 1_000
SVG_SIZE_WARNING_MB = 20


def configure_figure_style(
    font_size: float = 10,
    *,
    figure_dpi: int = 150,
    save_dpi: int = 300,
) -> None:
    """Apply the publication-wide Matplotlib and Seaborn settings once per notebook.

    This function intentionally does not set ``scanpy.settings.figdir``. Each
    notebook must provide its run-, stage-, and sample-specific output path.
    """
    sns.set_theme(style="whitegrid", context="paper", palette="colorblind")

    plt.rcParams.update(
        {
            "figure.dpi": figure_dpi,
            "savefig.dpi": save_dpi,
            "font.size": font_size,
            "axes.titlesize": font_size,
            "axes.labelsize": font_size,
            "xtick.labelsize": font_size,
            "ytick.labelsize": font_size,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            # Retain text as editable text, rather than converting it to paths.
            "svg.fonttype": "none",
            # Deterministic element identifiers for reproducible SVG output.
            "svg.hashsalt": "cap",
        }
    )


def rasterize_dense_artists(
    figure: Figure,
    *,
    min_points: int = RASTER_MIN_POINTS,
) -> int:
    """Rasterize dense scatter collections while retaining vector axes and text.

    Matplotlib embeds each rasterized collection as an image inside SVG
    exports. The resulting files retain editable vector labels, axes, legends,
    and non-rasterized panels, while avoiding one SVG element per cell.

    Returns
    -------
    int
        Number of collections marked for rasterization.
    """
    if min_points < 1:
        raise ValueError("min_points must be at least one.")

    rasterized = 0
    for axis in figure.axes:
        for collection in axis.collections:
            get_offsets = getattr(collection, "get_offsets", None)
            if get_offsets is None:
                continue

            offsets = get_offsets()
            if offsets is None or len(offsets) < min_points:
                continue

            collection.set_rasterized(True)
            rasterized += 1

    return rasterized


def _normalise_destination(destination: Path | str) -> Path:
    path = Path(destination)
    if not path.name:
        raise ValueError("destination must include a figure filename.")
    return path.with_suffix("")


def _validate_formats(formats: Sequence[FigureFormat]) -> tuple[FigureFormat, ...]:
    requested = tuple(formats)
    if not requested:
        raise ValueError("At least one output format is required.")

    unknown = set(requested).difference(VALID_FORMATS)
    if unknown:
        raise ValueError(f"Unsupported figure format(s): {sorted(unknown)}")

    if len(set(requested)) != len(requested):
        raise ValueError("Output formats must not contain duplicates.")

    return requested


def save_figure(
    figure: Figure,
    destination: Path | str,
    *,
    formats: Sequence[FigureFormat] = ("svg", "png"),
    rasterize_points: bool = False,
    min_points: int = RASTER_MIN_POINTS,
    png_dpi: int = 300,
    embedded_raster_dpi: int = 300,
    close: bool = True,
    overwrite: bool = False,
) -> list[Path]:
    """Save a figure without relying on a global output directory.

    Parameters
    ----------
    figure
        Matplotlib figure to export.
    destination
        Target filename without an extension. Its parent must be the explicit
        run-, stage-, or sample-specific figure directory selected by a notebook.
    formats
        Requested file formats. The default is SVG plus PNG. SVG is vector,
        while PNG is a raster image.
    rasterize_points
        If true, dense scatter collections are embedded as raster images in SVG
        output. This is useful for UMAPs: points become a manageable image
        object, while axes, labels, legends, and annotations remain vector
        elements. PNG output is a standalone raster image.
    min_points
        Minimum scatter-collection size for rasterization.
    png_dpi
        Resolution used only for standalone PNG output.
    embedded_raster_dpi
        Resolution of rasterized scatter collections embedded in SVG output.
    close
        Close the figure after a successful save. Keep this true for repeated
        per-library QC exports to avoid accumulating open Matplotlib figures.
    overwrite
        If false, fail before writing any output when a requested target exists.
    """
    output = _normalise_destination(destination)
    requested_formats = _validate_formats(formats)
    output.parent.mkdir(parents=True, exist_ok=True)

    output_files = [output.with_suffix(f".{extension}") for extension in requested_formats]
    existing = [path for path in output_files if path.exists()]
    if existing and not overwrite:
        existing_names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Figure output already exists: {existing_names}")

    if rasterize_points:
        rasterize_dense_artists(figure, min_points=min_points)

    temporary_files: list[Path] = []
    try:
        for output_file, extension in zip(output_files, requested_formats, strict=True):
            temporary_file = output_file.with_name(
                f"{output_file.stem}.tmp{output_file.suffix}"
            )
            temporary_files.append(temporary_file)

            dpi = png_dpi if extension == "png" else embedded_raster_dpi
            figure.savefig(
                temporary_file,
                dpi=dpi,
                bbox_inches="tight",
                facecolor="white",
                transparent=False,
            )

        for temporary_file, output_file in zip(temporary_files, output_files, strict=True):
            temporary_file.replace(output_file)

    except Exception:
        for temporary_file in temporary_files:
            temporary_file.unlink(missing_ok=True)
        raise
    finally:
        if close:
            plt.close(figure)

    for output_file in output_files:
        if (
            output_file.suffix == ".svg"
            and output_file.stat().st_size > SVG_SIZE_WARNING_MB * 1_000_000
        ):
            size_mb = output_file.stat().st_size / 1_000_000
            warn(
                f"{output_file.name} is {size_mb:.1f} MB. For dense scatter plots, "
                "call save_figure(..., rasterize_points=True).",
                stacklevel=2,
            )

    return output_files


def save_umap_figure(
    figure: Figure,
    destination: Path | str,
    *,
    embedded_raster_dpi: int = 300,
    close: bool = True,
    overwrite: bool = False,
) -> list[Path]:
    """Export a UMAP as SVG plus PNG.

    In the SVG, dense points are embedded as one raster image while axes, labels,
    legends, and annotations remain vector. The PNG is exported as a standalone
    raster image.
    """
    return save_figure(
        figure,
        destination,
        formats=("svg", "png"),
        rasterize_points=True,
        embedded_raster_dpi=embedded_raster_dpi,
        close=close,
        overwrite=overwrite,
    )


def categorical_palette(categories: Iterable[str]) -> dict[str, tuple[float, float, float]]:
    """Return deterministic colors in the supplied, unique category order."""
    category_list = list(categories)
    if len(category_list) != len(set(category_list)):
        raise ValueError("categories must not contain duplicates.")

    colors = sns.color_palette("husl", n_colors=len(category_list))
    return dict(zip(category_list, colors, strict=True))


def sample_order_by_condition(obs: pd.DataFrame) -> list[str]:
    """Return sample IDs ordered by clinical condition and then sample identifier."""
    required_columns = {"sample", "condition"}
    missing = required_columns.difference(obs.columns)
    if missing:
        raise KeyError(f"obs is missing required column(s): {sorted(missing)}")

    table = obs[["sample", "condition"]].drop_duplicates()
    unknown_conditions = set(table["condition"]).difference(CONDITION_ORDER)
    if unknown_conditions:
        raise ValueError(
            f"Unexpected condition label(s): {sorted(unknown_conditions)}. "
            f"Expected only: {CONDITION_ORDER}."
        )

    rank = {condition: index for index, condition in enumerate(CONDITION_ORDER)}
    return (
        table.assign(rank=table["condition"].map(rank))
        .sort_values(["rank", "sample"])["sample"]
        .tolist()
    )

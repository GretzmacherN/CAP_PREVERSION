"""Central, CWD-independent paths for the CAP analysis."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

_ROOT_MARKER = "pyproject.toml"
_STAGE_NAMES = {
    "00_cohort",
    "01_preprocessing",
    "02_cell_annotation",
    "03a_epithelial_integration",
    "03b_meta_programmes",
    "03c_antigen_coexpression",
    "03d_normal_tissue_protein",
    "04_programme_summary",
    "05_three_proteome"
}


def _read_dotenv(dotenv_path: Path) -> dict[str, str]:
    """Read simple KEY=VALUE pairs without adding a runtime dependency."""
    if not dotenv_path.is_file():
        return {}

    values: dict[str, str] = {}
    for raw_line in dotenv_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", maxsplit=1)
        values[key.strip()] = value.strip().strip("'\"")
    return values


def repository_root(start: Path | None = None) -> Path:
    """Find the repository root from the current location or a supplied path."""
    location = (start or Path.cwd()).resolve()
    for candidate in (location, *location.parents):
        if (candidate / _ROOT_MARKER).is_file():
            return candidate
    raise FileNotFoundError(
        f"No {_ROOT_MARKER} found above {location}. "
        "Open a directory inside the CAP repository."
    )


def _setting(name: str, root: Path) -> str | None:
    """Environment variables take precedence over a local ignored .env file."""
    return os.environ.get(name) or _read_dotenv(root / ".env").get(name)

@dataclass(frozen=True)
class ProjectPaths:
    """Resolved repository, external-data, and run locations."""

    root: Path
    data_root: Path
    run_id: str

    @property
    def frozen_dir(self) -> Path:
        return self.root / "data" / "frozen"

    @property
    def manifest_dir(self) -> Path:
        return self.root / "data" / "manifests"

    @property
    def raw_counts_dir(self) -> Path:
        return self.data_root / "raw" / "gse161529" / "raw_counts"

    @property
    def features_file(self) -> Path:
        return self.raw_counts_dir.parent / "features.tsv"

    @property
    def gencode_file(self) -> Path:
        return self.data_root / "references" / "gencode" / "gencode.v50.annotation.gtf.gz"

    @property
    def celltypist_model_file(self) -> Path:
        return self.data_root / "references" / "celltypist" / "model_from_pal.pkl"

    @property
    def hpa_dvp_file(self) -> Path:
        return self.data_root / "references" / "hpa" / "dvp_cell_type_group_data.tsv.zip"

    @property
    def supplied_release_dir(self) -> Path:
        return self.data_root / "processed" / "supplied_submission_release"

    @property
    def run_root(self) -> Path:
        return self.data_root / "runs" / self.run_id

    def stage_dir(self, stage: str) -> Path:
        if stage not in _STAGE_NAMES:
            allowed = ", ".join(sorted(_STAGE_NAMES))
            raise ValueError(f"Unknown stage {stage!r}; expected one of: {allowed}")
        return self.run_root / stage

    def ensure_stage(self, stage: str) -> Path:
        """Create the requested run-stage directory if it does not yet exist."""
        output_dir = self.stage_dir(stage)
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def prepare_stage(self, stage: str) -> Path:
        """Create an empty output directory; refuse accidental overwrites."""
        output_dir = self.stage_dir(stage)
        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError(
                f"Stage output already exists and is non-empty: {output_dir}. "
                "Set a new CAP_RUN_ID rather than overwriting it."
            )
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def require_file(self, path: Path, description: str) -> Path:
        if not path.is_file():
            raise FileNotFoundError(f"Missing {description}: {path}")
        return path

    def require_dir(self, path: Path, description: str) -> Path:
        if not path.is_dir():
            raise FileNotFoundError(f"Missing {description}: {path}")
        return path

    def validate_base_inputs(self) -> None:
        """Fail before expensive computation if required single-cell inputs are absent."""
        self.require_dir(self.raw_counts_dir, "GSE161529 filtered count directory")
        self.require_file(self.features_file, "shared features.tsv")
        self.require_file(self.gencode_file, "GENCODE v50 annotation")
        self.require_file(self.frozen_dir / "gene_policy.json", "frozen gene policy")
        self.require_file(
            self.frozen_dir / "mammary_lineage_signatures.json",
            "frozen lineage signatures",
            )
        self.require_file(
            self.frozen_dir / "cancer_meta_programme_signatures.xlsx",
            "frozen cancer meta-programme workbook",
            )


def get_paths(start: Path | None = None) -> ProjectPaths:
    """Return validated settings without creating input directories."""
    root = repository_root(start)
    data_root_raw = _setting("CAP_DATA_ROOT", root)
    run_id = _setting("CAP_RUN_ID", root)

    if not data_root_raw:
        raise EnvironmentError("CAP_DATA_ROOT is not set in the environment or .env.")
    if not run_id:
        raise EnvironmentError("CAP_RUN_ID is not set in the environment or .env.")

    data_root = Path(data_root_raw).expanduser().resolve()
    if not data_root.is_dir():
        raise FileNotFoundError(f"CAP_DATA_ROOT does not exist: {data_root}")

    return ProjectPaths(root=root, data_root=data_root, run_id=run_id)
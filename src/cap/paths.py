"""Locate the project root and the external raw data without hard-coded paths."""
import os
from pathlib import Path

ROOT_MARKER = "pyproject.toml"


def project_root(start=None) -> Path:
    env = os.environ.get("CAP_PROJECT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(start or Path.cwd()).resolve()
    for parent in [here, *here.parents]:
        if (parent / ROOT_MARKER).is_file():
            return parent
    raise FileNotFoundError(
        f"No {ROOT_MARKER} found above {here}. Run from inside the repository "
        "or set CAP_PROJECT to the repository root."
    )

def _ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def figures_dir() -> Path:
    return _ensure(project_root() / "results" / "figures")


def dataframe_dir() -> Path:
    return _ensure(project_root() / "results" / "dataframe_saves")

def results_dir() -> Path:
    return _ensure(project_root() / "results")


def runs_dir(step: str) -> Path:
    """Papermill outputs (executed notebooks and logs) for one pipeline step."""
    return _ensure(project_root() / "runs" / step)

def data_dir() -> Path:
    return _ensure(project_root() / "data")

def log_dir() -> Path:
    return _ensure(project_root() / "runs" / "001_Preprocessing" / "logs")

def interim_head_dir() -> Path:
    return _ensure(project_root() / "data" / "Adata_Saves" / "Interim_head")

def save_dir() -> Path:
    return _ensure(project_root() / "data" / "Adata_Saves")

def model_dir() -> Path:
    return _ensure(project_root() / "data" / "Scvi_Models")

def notebook_dir() -> Path:
    return _ensure(project_root() / "Notebooks")

def metadata_dir() -> Path:
    return _ensure(project_root() / "Metadata")



def sample_sheet() -> Path:
    path = metadata_dir() / "samples.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Sample sheet missing: {display_path(path)}")
    return path


def features_file() -> Path:
    """features.tsv is expected next to the raw_counts folder (shared by all samples)."""
    path = raw_data_dir().parent / "features.tsv"
    if not path.is_file():
        raise FileNotFoundError("features.tsv not found next to the raw_counts folder")
    return path

def _read_env_file(root: Path) -> dict:
    env_file = root / ".env"
    if not env_file.is_file():
        return {}
    pairs = (line.split("=", 1) for line in env_file.read_text().splitlines()
             if "=" in line and not line.lstrip().startswith("#"))
    return {k.strip(): v.strip().strip("'\"") for k, v in pairs}


def raw_data_dir() -> Path:
    env = os.environ.get("CAP_RAW_DATA") or _read_env_file(project_root()).get("CAP_RAW_DATA")
    if not env:
        raise EnvironmentError(
            "CAP_RAW_DATA is not set. Either export it or put a line "
            "CAP_RAW_DATA=/path/to/raw_counts in a .env file in the project root. "
            "See data/README.md for the GEO download."
        )
    path = Path(env).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"CAP_RAW_DATA does not exist: {path.name}")
    return path


def display_path(path, root=None) -> Path:
    """Project-relative where possible, so external absolute paths are never printed."""
    path, root = Path(path).resolve(), root or project_root()
    try:
        return path.relative_to(root)
    except ValueError:
        return Path("<external>") / path.name
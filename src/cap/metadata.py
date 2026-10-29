# src/cap/metadata.py
"""Sample-sheet loading and validation shared by all notebooks."""
from pathlib import Path
import pandas as pd
from cap.paths import sample_sheet

ALLOWED_CONDITIONS = {"N", "TN", "HER2+", "ER+"}
REQUIRED_COLUMNS = {"sample_id", "condition", "patient", "include", "exclusion_reason"}


def validate_sample_sheet(df: pd.DataFrame) -> pd.DataFrame:
    """Check a sample sheet before it is written or used. Returns df unchanged."""
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"sample sheet lacks columns: {sorted(missing)}")
    if df["sample_id"].duplicated().any():
        dup = df.loc[df["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"duplicate sample_id: {dup}")
    if df["include"].dtype != bool:
        raise ValueError(f"'include' must be True/False, found dtype {df['include'].dtype}")
    included = df[df["include"]]
    bad = set(included["condition"]) - ALLOWED_CONDITIONS
    if bad:
        raise ValueError(f"unexpected condition labels among included samples: {sorted(bad)}")
    if included["patient"].isna().any():
        raise ValueError("included sample without patient")
    return df


def load_sample_sheet(path=None) -> pd.DataFrame:
    """Read Metadata/samples.csv (or the given path), validate it, index by sample_id."""
    path = Path(path) if path is not None else sample_sheet()
    df = validate_sample_sheet(pd.read_csv(path))
    return df.set_index("sample_id")
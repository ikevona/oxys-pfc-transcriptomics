"""Small shared helpers."""
from pathlib import Path

import pandas as pd


def read_table(path, sheet=None):
    """Read .csv/.tsv/.xls/.xlsx into a DataFrame (first sheet unless `sheet` is given)."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_excel(path, sheet_name=sheet if sheet is not None else 0)


def norm_symbols(series):
    """Strip whitespace and upper-case, so 'Slc7a9 ' and 'SLC7A9' match."""
    return series.dropna().astype(str).str.strip().str.upper()

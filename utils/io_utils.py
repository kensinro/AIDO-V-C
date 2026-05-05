from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable
import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def detect_sep(path: str | Path) -> str:
    suffix = Path(path).suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return "\t"
    return ","


def read_table(path: str | Path, index_col: int | str | None = 0) -> pd.DataFrame:
    path = Path(path)
    sep = detect_sep(path)
    return pd.read_csv(path, sep=sep, index_col=index_col)


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".tsv":
        df.to_csv(path, sep="\t", index=True)
    else:
        df.to_csv(path, index=True)


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_gmt(path: str | Path) -> Dict[str, list[str]]:
    gene_sets: Dict[str, list[str]] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            set_name = parts[0]
            genes = [g for g in parts[2:] if g]
            gene_sets[set_name] = genes
    return gene_sets


def save_text(text: str, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")

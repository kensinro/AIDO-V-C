from __future__ import annotations

from typing import Iterable
import numpy as np
import pandas as pd


def tcga_patient_id(sample_id: str) -> str:
    return str(sample_id)[:12]


def is_primary_tumor(sample_id: str) -> bool:
    s = str(sample_id)
    parts = s.split("-")
    if len(parts) >= 4 and len(parts[3]) >= 2:
        return parts[3][:2] == "01"
    return True


def keep_primary_tumor_columns(df: pd.DataFrame) -> pd.DataFrame:
    keep_cols = [c for c in df.columns if is_primary_tumor(c)]
    return df.loc[:, keep_cols]


def patient_collapse(df: pd.DataFrame) -> pd.DataFrame:
    patient_map = {}
    for c in df.columns:
        pid = tcga_patient_id(c)
        patient_map.setdefault(pid, []).append(c)
    collapsed = {}
    for pid, cols in patient_map.items():
        if len(cols) == 1:
            collapsed[pid] = df[cols[0]]
        else:
            collapsed[pid] = df[cols].mean(axis=1)
    out = pd.DataFrame(collapsed)
    out.index = df.index
    return out


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    means = df.mean(axis=1)
    stds = df.std(axis=1).replace(0, np.nan)
    z = df.sub(means, axis=0).div(stds, axis=0)
    return z.fillna(0.0)


def align_modalities_by_patient(expr: pd.DataFrame, clin: pd.DataFrame, cn: pd.DataFrame | None = None,
                                patient_id_col: str = "patient_id"):
    common = set(expr.columns).intersection(set(clin[patient_id_col].astype(str)))
    if cn is not None:
        common = common.intersection(set(cn.columns))
    common = sorted(common)
    expr2 = expr.loc[:, common]
    clin2 = clin.set_index(patient_id_col).loc[common].reset_index()
    if cn is None:
        return expr2, clin2
    cn2 = cn.loc[:, common]
    return expr2, clin2, cn2

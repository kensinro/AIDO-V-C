from __future__ import annotations

from pathlib import Path
from typing import Literal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


Endpoint = Literal["OS", "PFI"]


def endpoint_cols(endpoint: str) -> tuple[str, str]:
    endpoint = endpoint.upper()
    if endpoint == "OS":
        return "OS_time", "OS_event"
    if endpoint == "PFI":
        return "PFI_time", "PFI_event"
    raise ValueError(f"Unsupported endpoint: {endpoint}")


def split_by_median(scores: pd.Series) -> pd.Series:
    threshold = scores.median()
    return pd.Series(np.where(scores >= threshold, "High", "Low"), index=scores.index)


def split_by_quartiles(scores: pd.Series) -> pd.Series:
    return pd.qcut(scores.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])


def logrank_p_from_groups(clin: pd.DataFrame, groups: pd.Series, endpoint: str = "OS") -> float:
    tcol, ecol = endpoint_cols(endpoint)
    df = clin.copy()
    df.index = df["patient_id"].astype(str)
    idx = groups.index.intersection(df.index)
    df = df.loc[idx]
    groups = groups.loc[idx]
    unique = pd.unique(groups)
    if len(unique) != 2:
        raise ValueError("logrank_p_from_groups expects exactly 2 groups.")
    g1, g2 = unique[0], unique[1]
    mask1 = groups == g1
    mask2 = groups == g2
    res = logrank_test(
        df.loc[mask1, tcol], df.loc[mask2, tcol],
        event_observed_A=df.loc[mask1, ecol],
        event_observed_B=df.loc[mask2, ecol],
    )
    return float(res.p_value)


def discriminability_from_p(p: float) -> float:
    p = max(float(p), 1e-300)
    return -np.log10(p)


def plot_km_by_groups(clin: pd.DataFrame, groups: pd.Series, endpoint: str, title: str, outpath: str | Path) -> None:
    tcol, ecol = endpoint_cols(endpoint)
    df = clin.copy()
    df.index = df["patient_id"].astype(str)
    idx = groups.index.intersection(df.index)
    df = df.loc[idx]
    groups = groups.loc[idx]

    plt.figure(figsize=(7, 5))
    kmf = KaplanMeierFitter()
    for label in pd.unique(groups):
        mask = groups == label
        kmf.fit(df.loc[mask, tcol], df.loc[mask, ecol], label=str(label))
        kmf.plot(ci_show=False)
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300)
    plt.close()

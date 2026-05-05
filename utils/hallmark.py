from __future__ import annotations

from typing import Dict
import numpy as np
import pandas as pd


def filter_gene_sets(gene_sets: Dict[str, list[str]], available_genes: set[str], min_genes: int = 5) -> Dict[str, list[str]]:
    out = {}
    for name, genes in gene_sets.items():
        keep = [g for g in genes if g in available_genes]
        if len(keep) >= min_genes:
            out[name] = keep
    return out


def compute_gene_set_scores(matrix: pd.DataFrame, gene_sets: Dict[str, list[str]]) -> pd.DataFrame:
    scores = {}
    for name, genes in gene_sets.items():
        scores[name] = matrix.loc[genes].mean(axis=0)
    return pd.DataFrame(scores)


def fuse_modalities(expr: pd.DataFrame, cn: pd.DataFrame, alpha: float = 1.0, beta: float = 1.0) -> pd.DataFrame:
    common_genes = expr.index.intersection(cn.index)
    common_samples = expr.columns.intersection(cn.columns)
    expr2 = expr.loc[common_genes, common_samples]
    cn2 = cn.loc[common_genes, common_samples]
    return alpha * expr2 + beta * cn2

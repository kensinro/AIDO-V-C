from __future__ import annotations

import random
from typing import Iterable
import pandas as pd

from .survival_utils import split_by_median, logrank_p_from_groups, discriminability_from_p


def random_gene_sets(all_genes: list[str], set_size: int, repeats: int, seed: int = 42) -> list[list[str]]:
    rng = random.Random(seed)
    out = []
    for _ in range(repeats):
        out.append(rng.sample(all_genes, set_size))
    return out


def random_baseline_distribution(matrix: pd.DataFrame, clin: pd.DataFrame, set_size: int, repeats: int,
                                 endpoint: str = "OS", seed: int = 42) -> list[float]:
    all_genes = list(matrix.index)
    dvals = []
    for i, genes in enumerate(random_gene_sets(all_genes, set_size, repeats, seed=seed)):
        scores = matrix.loc[genes].mean(axis=0)
        groups = split_by_median(scores)
        p = logrank_p_from_groups(clin, groups, endpoint=endpoint)
        dvals.append(discriminability_from_p(p))
    return dvals

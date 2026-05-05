from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from utils.io_utils import load_config, read_table, read_gmt, write_table, ensure_dir
from utils.preprocess import keep_primary_tumor_columns, patient_collapse, zscore_rows, align_modalities_by_patient
from utils.hallmark import filter_gene_sets, compute_gene_set_scores
from utils.survival_utils import split_by_median, logrank_p_from_groups, discriminability_from_p


def direction_from_group_means(scores: pd.Series, clin: pd.DataFrame, endpoint: str) -> int:
    # Simple pragmatic rule: high-score group with worse event-weighted median time is risk direction.
    groups = split_by_median(scores)
    df = clin.copy().set_index("patient_id")
    idx = groups.index.intersection(df.index)
    groups = groups.loc[idx]
    df = df.loc[idx]
    tcol = "OS_time" if endpoint.upper() == "OS" else "PFI_time"
    high_mean = df.loc[groups == "High", tcol].mean()
    low_mean = df.loc[groups == "Low", tcol].mean()
    return -1 if high_mean > low_mean else 1


def main(config_path: str) -> None:
    cfg = load_config(config_path)
    outdir = ensure_dir(Path(cfg["outputs"]["base_dir"]) / "exp3_prhs")

    expr = read_table(cfg["inputs"]["expression_path"], index_col=0)
    clin = read_table(cfg["inputs"]["clinical_path"], index_col=None)
    gene_sets = read_gmt(cfg["inputs"]["hallmark_gmt_path"])

    if cfg["preprocessing"]["keep_primary_tumor_only"]:
        expr = keep_primary_tumor_columns(expr)
    expr = patient_collapse(expr)
    if cfg["preprocessing"]["zscore_expression"]:
        expr = zscore_rows(expr)

    expr, clin = align_modalities_by_patient(expr, clin, patient_id_col=cfg["clinical"]["patient_id_col"])
    clin = clin.rename(columns={cfg["clinical"]["patient_id_col"]: "patient_id"})

    gene_sets = filter_gene_sets(gene_sets, set(expr.index), cfg["analysis"]["min_genes_per_set"])
    scores = compute_gene_set_scores(expr, gene_sets)

    rng = np.random.default_rng(int(cfg["random_seed"]))
    patients = np.array(scores.index)
    rng.shuffle(patients)
    split = int(len(patients) * float(cfg["analysis"]["train_fraction"]))
    train_ids = patients[:split]
    test_ids = patients[split:]

    endpoint = cfg["analysis"]["endpoint"]
    clin_idx = clin.set_index("patient_id")
    train_clin = clin_idx.loc[train_ids].reset_index()
    test_clin = clin_idx.loc[test_ids].reset_index()
    train_scores = scores.loc[train_ids]
    test_scores = scores.loc[test_ids]

    ranked = []
    for hallmark in train_scores.columns:
        s = train_scores[hallmark]
        p = logrank_p_from_groups(train_clin, split_by_median(s), endpoint=endpoint)
        d = discriminability_from_p(p)
        direction = direction_from_group_means(s, train_clin, endpoint)
        ranked.append((hallmark, p, d, direction))
    ranked.sort(key=lambda x: x[2], reverse=True)
    top_m = int(cfg["analysis"]["top_m_hallmarks"])
    selected = ranked[:top_m]

    weights = np.array([x[2] for x in selected], dtype=float)
    weights = weights / weights.sum() if weights.sum() > 0 else np.ones_like(weights) / len(weights)

    comp = pd.Series(0.0, index=test_scores.index)
    for w, (hallmark, _, _, direction) in zip(weights, selected):
        comp = comp + w * direction * test_scores[hallmark]

    p_test = logrank_p_from_groups(test_clin, split_by_median(comp), endpoint=endpoint)
    d_test = discriminability_from_p(p_test)

    summary = pd.DataFrame(selected, columns=["hallmark", "p_train", "D_train", "direction"])
    summary["weight"] = weights
    write_table(summary.set_index("hallmark"), outdir / "selected_hallmarks.csv")
    write_table(pd.DataFrame({"patient_id": comp.index, "PRHS": comp.values}).set_index("patient_id"), outdir / "prhs_test_scores.csv")
    write_table(pd.DataFrame({"metric": ["p_test", "D_test"], "value": [p_test, d_test]}).set_index("metric"), outdir / "prhs_test_summary.csv")
    print("EXP3 completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from utils.io_utils import load_config, read_table, read_gmt, write_table, ensure_dir
from utils.preprocess import keep_primary_tumor_columns, patient_collapse, zscore_rows, align_modalities_by_patient
from utils.hallmark import filter_gene_sets, compute_gene_set_scores
from utils.survival_utils import split_by_median, logrank_p_from_groups, discriminability_from_p
from utils.random_baseline import random_baseline_distribution


def main(config_path: str) -> None:
    cfg = load_config(config_path)
    outdir = ensure_dir(Path(cfg["outputs"]["base_dir"]) / "exp2_baseline")

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

    rows = []
    repeats = int(cfg["analysis"]["random_repeats"])
    endpoint = cfg["analysis"]["endpoint"]
    seed = int(cfg["random_seed"])

    for hallmark, genes in gene_sets.items():
        obs_scores = scores[hallmark]
        obs_groups = split_by_median(obs_scores)
        p_obs = logrank_p_from_groups(clin, obs_groups, endpoint=endpoint)
        d_obs = discriminability_from_p(p_obs)

        d_rand = random_baseline_distribution(expr, clin, set_size=len(genes), repeats=repeats, endpoint=endpoint, seed=seed)
        mu_rand = float(np.mean(d_rand))
        sd_rand = float(np.std(d_rand))
        d_norm = d_obs - mu_rand
        z = d_norm / sd_rand if sd_rand > 0 else np.nan
        percentile = float(np.mean(np.array(d_rand) <= d_obs))

        rows.append({
            "hallmark": hallmark,
            "n_genes": len(genes),
            "p_obs": p_obs,
            "D_obs": d_obs,
            "mu_rand": mu_rand,
            "sd_rand": sd_rand,
            "D_norm": d_norm,
            "Z": z,
            "percentile": percentile,
        })

    res = pd.DataFrame(rows).sort_values("D_norm", ascending=False).set_index("hallmark")
    write_table(res, outdir / "baseline_calibrated_discriminability.csv")
    print("EXP2 completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)

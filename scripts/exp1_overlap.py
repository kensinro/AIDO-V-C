from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from utils.io_utils import load_config, read_table, read_gmt, write_table, ensure_dir
from utils.preprocess import keep_primary_tumor_columns, patient_collapse, zscore_rows, align_modalities_by_patient
from utils.hallmark import filter_gene_sets, compute_gene_set_scores
from utils.survival_utils import split_by_quartiles, plot_km_by_groups


def main(config_path: str) -> None:
    cfg = load_config(config_path)
    outdir = ensure_dir(Path(cfg["outputs"]["base_dir"]) / "exp1_overlap")

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

    # Use the first hallmark as placeholder, or change manually later.
    target = scores.columns[0]
    groups = split_by_quartiles(scores[target])
    plot_km_by_groups(clin, groups, cfg["analysis"]["endpoint"], f"EXP1 Overlap: {target}", outdir / f"km_{target}.png")

    out = pd.DataFrame({"patient_id": scores.index, "score": scores[target].values, "quartile": groups.values}).set_index("patient_id")
    write_table(out, outdir / f"quartiles_{target}.csv")
    print(f"EXP1 completed. Target hallmark: {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)

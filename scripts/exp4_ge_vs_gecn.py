from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from utils.io_utils import load_config, read_table, read_gmt, write_table, ensure_dir
from utils.preprocess import keep_primary_tumor_columns, patient_collapse, zscore_rows, align_modalities_by_patient
from utils.hallmark import filter_gene_sets, compute_gene_set_scores, fuse_modalities
from utils.survival_utils import split_by_median, logrank_p_from_groups, discriminability_from_p


def main(config_path: str) -> None:
    cfg = load_config(config_path)
    outdir = ensure_dir(Path(cfg["outputs"]["base_dir"]) / "exp4_ge_vs_gecn")

    expr = read_table(cfg["inputs"]["expression_path"], index_col=0)
    cn = read_table(cfg["inputs"]["copy_number_path"], index_col=0)
    clin = read_table(cfg["inputs"]["clinical_path"], index_col=None)
    gene_sets = read_gmt(cfg["inputs"]["hallmark_gmt_path"])

    if cfg["preprocessing"]["keep_primary_tumor_only"]:
        expr = keep_primary_tumor_columns(expr)
        cn = keep_primary_tumor_columns(cn)
    expr = patient_collapse(expr)
    cn = patient_collapse(cn)
    if cfg["preprocessing"]["zscore_expression"]:
        expr = zscore_rows(expr)
    if cfg["preprocessing"]["zscore_copy_number"]:
        cn = zscore_rows(cn)

    expr, clin, cn = align_modalities_by_patient(expr, clin, cn=cn, patient_id_col=cfg["clinical"]["patient_id_col"])
    clin = clin.rename(columns={cfg["clinical"]["patient_id_col"]: "patient_id"})

    gene_sets_expr = filter_gene_sets(gene_sets, set(expr.index), cfg["analysis"]["min_genes_per_set"])
    ge_scores = compute_gene_set_scores(expr, gene_sets_expr)

    fusion = fuse_modalities(expr, cn, alpha=float(cfg["analysis"]["alpha_ge"]), beta=float(cfg["analysis"]["beta_cn"]))
    gene_sets_fusion = filter_gene_sets(gene_sets, set(fusion.index), cfg["analysis"]["min_genes_per_set"])
    fused_scores = compute_gene_set_scores(fusion, gene_sets_fusion)

    common_hallmarks = sorted(set(ge_scores.columns).intersection(fused_scores.columns))
    rows = []
    endpoint = cfg["analysis"]["endpoint"]
    for hallmark in common_hallmarks:
        p_ge = logrank_p_from_groups(clin, split_by_median(ge_scores[hallmark]), endpoint=endpoint)
        p_f = logrank_p_from_groups(clin, split_by_median(fused_scores[hallmark]), endpoint=endpoint)
        d_ge = discriminability_from_p(p_ge)
        d_f = discriminability_from_p(p_f)
        rows.append({
            "hallmark": hallmark,
            "p_GE": p_ge,
            "D_GE": d_ge,
            "p_GECN": p_f,
            "D_GECN": d_f,
            "Delta_D": d_f - d_ge,
        })

    out = pd.DataFrame(rows).sort_values("Delta_D", ascending=False).set_index("hallmark")
    write_table(out, outdir / "ge_vs_gecn.csv")
    print("EXP4 completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)

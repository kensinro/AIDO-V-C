# AIDO-VALIDATION-CORE Python Package

This package provides a modular starting codebase for running the four core experiments in **AIDO-VALIDATION-CORE**:

- **EXP1**: Observable overlap analysis (quartiles / Kaplan–Meier / Q2–Q3 overlap)
- **EXP2**: Baseline-calibrated discriminability using random gene-set ensembles
- **EXP3**: PRHS (train/test composite observable)
- **EXP4**: GE vs GE+CN multi-omics comparison

## Folder structure

- `scripts/exp1_overlap.py` — EXP1
- `scripts/exp2_baseline_calibration.py` — EXP2
- `scripts/exp3_prhs.py` — EXP3
- `scripts/exp4_ge_vs_gecn.py` — EXP4
- `utils/io_utils.py` — file loading / writing helpers
- `utils/preprocess.py` — sample alignment and preprocessing
- `utils/hallmark.py` — gene set parsing and score computation
- `utils/survival_utils.py` — Kaplan–Meier, log-rank, D(O)
- `utils/random_baseline.py` — random gene-set ensemble generation
- `config/config_template.yaml` — configurable parameters

## Expected inputs

You should provide:

1. Gene expression matrix (`csv` / `tsv`)
   - rows = genes
   - columns = samples
2. Copy number matrix (`csv` / `tsv`) for EXP4
3. Clinical table (`csv` / `tsv`)
4. Hallmark gene sets (`.gmt` preferred; csv/txt can be adapted)

## Required clinical columns

By default the code expects the clinical table to contain:

- `patient_id`
- `OS_time`
- `OS_event`
- `PFI_time`
- `PFI_event`

These names can be changed in `config/config_template.yaml`.

## Sample ID logic

The code assumes TCGA-like sample barcodes and truncates sample IDs to the first 12 characters for patient-level matching.
Only primary tumor samples are kept by default (`sample type = 01`) when this can be inferred from the barcode.

## Installation

Recommended Python version: 3.10+

Suggested packages:

```bash
pip install pandas numpy scipy matplotlib lifelines pyyaml
```

## Quick start

Edit `config/config_template.yaml`, then run for example:

```bash
python scripts/exp1_overlap.py --config config/config_template.yaml
python scripts/exp2_baseline_calibration.py --config config/config_template.yaml
python scripts/exp3_prhs.py --config config/config_template.yaml
python scripts/exp4_ge_vs_gecn.py --config config/config_template.yaml
```

## Notes

- This package is designed to be **transparent and editable**, not a black box.
- Default settings are intentionally conservative.
- You should inspect outputs before using them in a manuscript.
- If your data format differs slightly, most changes will be confined to `io_utils.py` and `preprocess.py`.

# Observable–Endpoint Survival Audit

Code and reproducibility materials for the manuscript:

> **Auditing observable–endpoint alignment before survival-model escalation across cancer systems**  
> Sin Guan Kong

This repository implements a transparent computational workflow for auditing whether a molecular representation contains reproducible information about a specified survival endpoint **before** additional model complexity or omics layers are introduced.

The repository supports two linked analyses:

1. a **cross-cancer transcriptomic audit** based on Hallmark observables, size-matched random gene-set calibration, and held-out Personalized Relative Hallmark Score (PRHS) evaluation; and
2. a **selected multi-omics extension** in BRCA, GBM, and KIRC that evaluates whether copy-number, mutation, and RPPA layers provide incremental survival-relevant information beyond gene expression after accounting for matched-sample attrition and metric dependence.

The central purpose is diagnostic rather than algorithmic. A weak result may reflect:

- inadequate extraction of information that is present;
- mismatch between the measured representation and the clinical endpoint;
- insufficient endpoint information, such as sparse events or limited follow-up; or
- loss of effective information after multi-omics matching and fusion.

Accordingly, the workflow asks whether further model escalation is justified, or whether the representation, modality set, or endpoint should be reconsidered first.

---

## 1. Scope

### Primary transcriptomic audit

The primary analysis evaluates Hallmark-based gene-expression observables across selected TCGA cohorts:

- BRCA
- HNSC
- LIHC
- LUAD
- GBM
- COAD/READ
- PRAD

The workflow includes:

- patient-level barcode harmonization;
- primary-tumor filtering;
- gene-wise standardization;
- Hallmark mean-score construction;
- median-split Kaplan–Meier and two-sided log-rank analysis;
- endpoint-conditioned evidence display;
- size-matched random gene-set calibration; and
- held-out PRHS evaluation.

### Extended multi-omics audit

The extended analysis evaluates the following matched configurations in BRCA, GBM, and KIRC:

- GE
- GE+CN
- GE+MU
- GE+RPPA
- GE+CN+MU
- GE+CN+RPPA
- GE+CN+MU+RPPA

where:

- **GE** = gene expression;
- **CN** = copy-number alteration;
- **MU** = somatic mutation;
- **RPPA** = reverse-phase protein array abundance.

The multi-omics analysis is intentionally simple and transparent. It is not presented as a benchmark of all integration algorithms and does not imply that untested modalities are uninformative.

---

## 2. Audit questions

The computational workflow is organized around five questions:

1. Does a structured Hallmark observable exceed an equally sized random gene set?
2. Does its survival separation generalize to held-out patients?
3. Does an added modality improve on the simpler gene-expression reference?
4. Is any apparent gain retained after matched-sample loss and across complementary survival summaries?
5. Does weak performance point primarily to predictor inadequacy, representation mismatch, or endpoint limitation?

---

## 3. Key quantities

### 3.1 Hallmark observable

For patient \(i\) and Hallmark gene set \(G_h\), the Hallmark score is

\[
H_h(i)=\frac{1}{m_h}\sum_{g=1}^{m_h}z_{g,i},
\]

where \(z_{g,i}\) is the standardized expression of gene \(g\) in patient \(i\), and \(m_h\) is the number of measured genes in Hallmark set \(h\).

### 3.2 Discriminability evidence scale

The display quantity is

\[
D(o,e)=-\log_{10}\left[P_{\mathrm{LR}}(o,e)\right],
\]

where \(P_{\mathrm{LR}}\) is the two-sided log-rank \(P\) value for observable \(o\) and endpoint \(e\).

\(D\) is a monotonic evidence scale for visualization. It is **not** a replacement for hazard ratios, confidence intervals, concordance measures, calibration, or predictive accuracy metrics.

### 3.3 Personalized Relative Hallmark Score

For held-out patient \(i\),

\[
\mathrm{PRHS}_{i}
=
\frac{\sum_{h\in\mathcal{H}^{*}}w_hs_{h,i}}
{\sum_{h\in\mathcal{H}^{*}}|w_h|},
\]

where:

- \(\mathcal{H}^{*}\) is the Hallmark set selected only in the training data;
- \(s_{h,i}\) is the standardized, direction-aligned score for Hallmark \(h\) in patient \(i\); and
- \(w_h\) is the corresponding training-derived weight.

### 3.4 Multi-omics retention and increment

For configuration \(C\), sample retention relative to GE is

\[
R_n(C)=\frac{n(C)}{n_{\mathrm{GE}}},
\]

and the change in displayed log-rank evidence is

\[
\Delta D(C)=D(C)-D(\mathrm{GE}).
\]

Neither quantity is interpreted alone. Incremental value is judged together with event count, continuous Cox association, matched-sample loss, and consistency across metrics.

---

## 4. Repository structure

```text
AIDO-V-C/
├── README.md
├── requirements.txt
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── config/
│   ├── config_template.yaml
│   ├── config_cbm.yaml
│   └── config_demo.yaml
├── main_pipeline/
│   └── run_aido_v_all_cancers.py
├── legacy_modular_pipeline/
│   ├── README.md
│   ├── run_all.py
│   ├── scripts/
│   │   ├── exp1_overlap.py
│   │   ├── exp2_baseline_calibration.py
│   │   ├── exp3_prhs.py
│   │   └── exp4_ge_vs_gecn.py
│   └── utils/
│       ├── __init__.py
│       ├── hallmark.py
│       ├── io_utils.py
│       ├── preprocess.py
│       ├── random_baseline.py
│       └── survival_utils.py
├── scripts/
│   ├── 01_prepare_expression.py
│   ├── 02_prepare_survival.py
│   ├── 03_align_patients.py
│   ├── 04_run_hallmark_audit.py
│   ├── 05_run_random_baseline.py
│   ├── 06_run_prhs.py
│   ├── 07_run_multiomics_audit.py
│   ├── 08_run_pca_comparison.py
│   └── 09_build_manuscript_outputs.py
├── utils/
│   ├── data_io.py
│   ├── preprocessing.py
│   ├── hallmark_scoring.py
│   ├── survival_analysis.py
│   ├── random_sets.py
│   ├── prhs.py
│   ├── multiomics.py
│   ├── audit_logging.py
│   └── plotting.py
├── data_preparation_helpers/
│   ├── STEP0.2_ExpressionCSV.py
│   ├── STEP0.3_SurvivalCSV.py
│   └── STEP0.4_Alignment.py
├── data/
│   ├── README.md
│   ├── manifest.tsv
│   ├── raw/
│   ├── processed/
│   ├── hallmark/
│   └── demo/
├── outputs/
│   ├── tables/
│   ├── figures/
│   ├── audit/
│   └── logs/
├── demo/
│   ├── demo_expression.csv
│   ├── demo_clinical.csv
│   ├── run_demo.py
│   └── expected_outputs/
└── docs/
    ├── CODE_PACKAGE_MANIFEST.txt
    ├── MANUSCRIPT_CODE_MAP.md
    ├── DATA_DICTIONARY.md
    └── REPRODUCIBILITY_NOTES.md
```

The existing `legacy_modular_pipeline/` is retained for provenance. The manuscript-facing workflow should use the current scripts and configuration files described below.

---

## 5. Software requirements

### Recommended environment

- Python 3.10 or later
- 16 GB RAM or more for full TCGA-scale analyses
- Windows, Linux, or macOS

### Installation

```bash
git clone https://github.com/kensinro/AIDO-V-C.git
cd AIDO-V-C
python -m venv .venv
```

Activate the environment.

**Windows**

```bash
.venv\Scripts\activate
```

**Linux/macOS**

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

The exact package versions used for the manuscript release should be recorded in `requirements.txt` and in the versioned release notes.

---

## 6. Data availability and acquisition

Raw TCGA molecular and clinical data are not redistributed in this repository.

The manuscript analyses use publicly available TCGA data obtained through UCSC Xena or the corresponding TCGA public resources. Users must download the required molecular and clinical files from the original source.

### Required inputs

For the primary transcriptomic audit:

- gene-expression matrix;
- survival annotations;
- TCGA sample or patient identifiers; and
- Hallmark gene-set definitions.

For the extended multi-omics audit:

- gene expression;
- copy number;
- mutation calls;
- RPPA abundance; and
- overall-survival annotations.

### Recommended local layout

```text
data/
├── raw/
│   ├── BRCA/
│   │   ├── expression.tsv
│   │   ├── copy_number.tsv
│   │   ├── mutation.tsv
│   │   ├── rppa.tsv
│   │   └── survival.tsv
│   ├── GBM/
│   ├── HNSC/
│   ├── KIRC/
│   ├── LIHC/
│   ├── LUAD/
│   ├── COADREAD/
│   └── PRAD/
├── hallmark/
│   └── hallmark_gene_sets.gmt
└── manifest.tsv
```

### Data manifest

`data/manifest.tsv` should contain at least:

```text
cohort	modality	source_url	original_filename	local_filename	access_date	checksum
```

This file provides an auditable link between the public source and each local input.

---

## 7. Input format

### Expression, copy-number, and RPPA matrices

Expected orientation:

- rows: features;
- columns: TCGA samples or patients;
- first column: gene or feature identifier.

### Mutation data

Mutation calls are converted to a binary gene-by-patient matrix:

- 1 = at least one qualifying mutation event observed;
- 0 = no qualifying event observed.

### Survival file

Minimum required columns:

```text
patient_id	time	event	endpoint
```

where:

- `patient_id` is the first 12 characters of the TCGA barcode;
- `time` is the follow-up or survival time;
- `event` is coded as 1 for event and 0 for censored; and
- `endpoint` identifies OS or another prespecified endpoint.

### Patient harmonization

TCGA barcodes are truncated to the first 12 characters for patient-level matching. Only primary tumor samples are retained for the transcriptomic analysis. Where multiple molecular samples map to one patient, continuous measurements are averaged unless otherwise specified in the configuration file.

---

## 8. Configuration

Copy the template:

```bash
cp config/config_template.yaml config/config_cbm.yaml
```

Edit paths and analysis settings in `config/config_cbm.yaml`.

Illustrative structure:

```yaml
project:
  name: observable_endpoint_survival_audit
  seed: 20260715

paths:
  data_root: data/raw
  hallmark_gmt: data/hallmark/hallmark_gene_sets.gmt
  output_root: outputs

cohorts:
  primary:
    - BRCA
    - HNSC
    - LIHC
    - LUAD
    - GBM
    - COADREAD
    - PRAD
  multiomics:
    - BRCA
    - GBM
    - KIRC

preprocessing:
  primary_tumor_code: "01"
  tcga_patient_barcode_length: 12
  standardize_features: true

hallmark:
  minimum_measured_genes: 10
  score_method: mean_z

survival:
  endpoint: OS
  split_rule: median
  logrank_two_sided: true

random_baseline:
  n_random_sets: 1000
  size_matched: true

prhs:
  train_fraction: 0.70
  n_selected_hallmarks: 5
  use_training_direction: true
  use_training_weights: true

multiomics:
  continuous_top_variance_features: 300
  mutation_top_variance_features: 300
  continuous_summary: mean_absolute_z
  mutation_summary: log1p_burden
  fusion: unweighted_mean
  configurations:
    - GE
    - GE_CN
    - GE_MU
    - GE_RPPA
    - GE_CN_MU
    - GE_CN_RPPA
    - GE_CN_MU_RPPA
```

Values shown above are examples and must match the final manuscript analysis before release.

---

## 9. Running the analyses

### 9.1 Full manuscript workflow

```bash
python main_pipeline/run_aido_v_all_cancers.py \
  --config config/config_cbm.yaml
```

or, if using the manuscript-facing modular scripts:

```bash
python scripts/04_run_hallmark_audit.py --config config/config_cbm.yaml
python scripts/05_run_random_baseline.py --config config/config_cbm.yaml
python scripts/06_run_prhs.py --config config/config_cbm.yaml
python scripts/07_run_multiomics_audit.py --config config/config_cbm.yaml
python scripts/08_run_pca_comparison.py --config config/config_cbm.yaml
python scripts/09_build_manuscript_outputs.py --config config/config_cbm.yaml
```

### 9.2 Legacy workflow

The original modular pipeline is retained in `legacy_modular_pipeline/`:

```bash
cd legacy_modular_pipeline
python scripts/run_all.py --config config/config_template.yaml
```

The legacy workflow is provided for provenance and may not reproduce every output of the upgraded CBM manuscript without synchronization.

### 9.3 Demonstration run

A lightweight demonstration should be runnable without downloading TCGA-scale data:

```bash
python demo/run_demo.py --config config/config_demo.yaml
```

The demonstration verifies software execution and output structure only. It is not intended to reproduce the manuscript’s numerical results.

---

## 10. Manuscript-to-code map

| Manuscript item | Analysis | Script or entry point | Principal output |
|---|---|---|---|
| Figure 1 | Study workflow | documentation figure | `outputs/figures/figure1_workflow.png` |
| Table 1 | Cross-cancer audit summary | `04_run_hallmark_audit.py`, `05_run_random_baseline.py`, `06_run_prhs.py` | `outputs/tables/table1_cross_cancer_audit.csv` |
| Figure 2 | Endpoint-conditioned discriminability landscape | `09_build_manuscript_outputs.py` | `outputs/figures/figure2_discriminability_landscape.png` |
| Figure 3 | Diagnostic regimes schematic | documentation figure | `outputs/figures/figure3_diagnostic_regimes.png` |
| Table 2 | Main multi-omics audit findings | `07_run_multiomics_audit.py` | `outputs/tables/table2_multiomics_summary.csv` |
| Supplementary Table S1 | Complete multi-omics configurations | `07_run_multiomics_audit.py` | `outputs/tables/table_s1_complete_multiomics.csv` |
| Supplementary Table S2 | Cross-cancer multi-omics synthesis | `09_build_manuscript_outputs.py` | `outputs/tables/table_s2_multiomics_audit_decisions.csv` |
| Supplementary Figure S1 | Multi-omics survival-separation evidence | `07_run_multiomics_audit.py` | `outputs/figures/figure_s1_multiomics_discriminability.png` |
| Supplementary Figure S2 | GE versus integrated PC1 correspondence | `08_run_pca_comparison.py` | `outputs/figures/figure_s2_pca_correspondence.png` |

The exact filenames in this table should be synchronized with the released repository before manuscript submission.

---

## 11. Expected outputs

### Primary transcriptomic audit

```text
outputs/
├── tables/
│   ├── hallmark_results_by_cohort.csv
│   ├── random_baseline_results.csv
│   ├── prhs_heldout_results.csv
│   └── table1_cross_cancer_audit.csv
├── figures/
│   ├── km_curves/
│   ├── random_baseline/
│   └── figure2_discriminability_landscape.png
└── audit/
    ├── patient_alignment_log.csv
    ├── hallmark_gene_coverage.csv
    └── run_manifest.json
```

### Extended multi-omics audit

```text
outputs/
├── tables/
│   ├── table2_multiomics_summary.csv
│   ├── table_s1_complete_multiomics.csv
│   ├── sample_retention_by_configuration.csv
│   └── metric_agreement_by_configuration.csv
├── figures/
│   ├── figure_s1_multiomics_discriminability.png
│   └── figure_s2_pca_correspondence.png
└── audit/
    ├── multiomics_patient_matching.csv
    ├── modality_feature_counts.csv
    └── multiomics_run_manifest.json
```

---

## 12. Audit trail and reproducibility outputs

Each run should record:

- analysis timestamp;
- software versions;
- random seed;
- configuration file;
- input filenames;
- optional input checksums;
- cohort and modality;
- patient counts before and after matching;
- event counts;
- feature counts before and after filtering;
- Hallmark gene coverage;
- selected Hallmarks and training-derived directions;
- random-set parameters;
- survival statistics;
- output filenames; and
- warnings, exclusions, and failed configurations.

Recommended audit files:

```text
outputs/audit/run_manifest.json
outputs/audit/patient_alignment_log.csv
outputs/audit/sample_retention_by_configuration.csv
outputs/audit/feature_selection_log.csv
outputs/audit/prhs_training_selection.csv
outputs/audit/warnings.log
```

Sample attrition is treated as an analysis result rather than a hidden preprocessing detail.

---

## 13. Interpretation boundaries

The following constraints apply to the released workflow:

1. The primary score is a transparent Hallmark mean-z construction and may not capture nonlinear or distribution-sensitive pathway behavior.
2. The displayed quantity \(D\) is derived from a log-rank \(P\) value and is not a complete measure of predictive performance.
3. Median splitting can attenuate continuous risk information; Cox results are therefore reported separately where applicable.
4. Multi-omics fusion is an unweighted mean of standardized layer summaries and is not optimized using survival outcomes.
5. Matched sample size and event count can differ across configurations.
6. Hazard ratios are interpreted within configurations and are not ranked across differently scaled fused scores.
7. The analysis does not evaluate all cancers, endpoints, omics modalities, or integration algorithms.
8. miRNA, methylation, metabolomics, spatial, single-cell, longitudinal, and treatment-specific data are outside the present scope.
9. PCA correspondence describes unsupervised variance structure and does not establish survival relevance.
10. The results support an observable–endpoint audit; they do not prove intrinsic non-identifiability of a cancer system.

---

## 14. Reproducing the manuscript release

For a formal manuscript release:

1. Confirm that all scripts run from a clean environment.
2. Lock package versions in `requirements.txt`.
3. Confirm that `config/config_cbm.yaml` matches the submitted manuscript.
4. Regenerate all tables and figures.
5. Compare regenerated outputs with the manuscript values.
6. Update `docs/MANUSCRIPT_CODE_MAP.md`.
7. Update `CHANGELOG.md`.
8. Create a versioned GitHub release.
9. Archive the release in Zenodo.
10. Add the Zenodo DOI to this README and the manuscript Code Availability statement.

Recommended release name:

```text
v1.0.0-cbm
```

Recommended release title:

```text
Code release accompanying the observable–endpoint survival audit manuscript
```

---

## 15. Code availability statement for the manuscript

After the repository and Zenodo release are synchronized, the following statement may be used:

> **Code Availability**  
> The analysis code and reproducibility materials are publicly available at the GitHub repository listed below and archived in a versioned Zenodo release. The repository includes scripts for patient-level harmonization, Hallmark score construction, size-matched random gene-set calibration, held-out PRHS evaluation, matched multi-omics integration, median-split and continuous survival analyses, sample-retention auditing, and generation of the principal manuscript and supplementary outputs. A manuscript-to-code map, configuration templates, software dependencies, data-acquisition instructions, and a lightweight demonstration dataset are provided. Raw TCGA data are not redistributed and must be obtained from the original public sources described in the Data Availability statement.

Repository:

```text
https://github.com/kensinro/AIDO-V-C
```

Archived release:

```text
https://doi.org/[ADD-ZENODO-DOI]
```

Do not add the Zenodo sentence until the archived release exists.

---

## 16. Citation

Please cite the accompanying manuscript when using this code.

```bibtex
@article{KongObservableEndpointAudit,
  author  = {Kong, Sin Guan},
  title   = {Auditing observable--endpoint alignment before survival-model escalation across cancer systems},
  journal = {Computers in Biology and Medicine},
  year    = {2026},
  note    = {Submitted manuscript}
}
```

Update the journal, year, volume, pages, and DOI after publication.

A machine-readable citation record should also be provided in `CITATION.cff`.

---

## 17. License

Add a clear software license before release. A permissive license such as MIT or BSD-3-Clause is commonly used for academic code, subject to institutional requirements.

The license applies to the code in this repository. It does not alter the terms imposed by the original TCGA, UCSC Xena, MSigDB, or other data providers.

---

## 18. Contact

**Sin Guan Kong**  
Department of Electrical and Electronic Engineering  
Lee Kong Chian Faculty of Engineering and Science  
Universiti Tunku Abdul Rahman  
Malaysia

Email: `kongsg@utar.edu.my`

Repository: `https://github.com/kensinro/AIDO-V-C`

---

## 19. Changelog summary

### v1.0.0-cbm

- Reframed the repository around an observable–endpoint survival audit.
- Added manuscript-to-code mapping.
- Clarified the role of the log-rank-derived evidence scale.
- Expanded the multi-omics audit to all reported GE/CN/MU/RPPA configurations.
- Added matched-sample retention and metric-agreement outputs.
- Distinguished median-split log-rank results from continuous Cox associations.
- Added reproducible environment, configuration, data-manifest, and demo guidance.
- Removed unsupported patient-level resampling claims from the evidentiary chain.
- Added bounded interpretation and reproducibility notes.

---

## 20. Acknowledgement of prior review

The repository should not reproduce reviewer comments or identify individual reviewers. Improvements made after peer review should be documented neutrally through the changelog and version history.

The manuscript may separately acknowledge anonymous reviewers for constructive feedback.

# DDS-Hayling Paper Pipeline

Fully reproducible analysis pipeline for:

> **Fronto-cingulate Information Dynamics During Semantic Inhibition: A Dual Damped Sinusoid Decomposition of the N450 ERP**

This repository contains every script needed to reproduce the figures, tables, and statistical results reported in the paper, starting from the pre-computed derivative files.

---

## Repository structure

```
DDS-Hayling-paper/
├── config/
│   └── analysis.yaml          # single source of truth for all parameters
├── pipeline/
│   ├── _utils.py              # shared helpers (paths, stats, style)
│   ├── 00_validate_environment.py
│   ├── 01_erp_grand_average.py
│   ├── 02_dds_fit.py          # derivative check / regeneration instructions
│   ├── 03_info_dynamics_ais.py
│   ├── 04_info_dynamics_te.py
│   ├── 05_stats_main.py       # primary paper statistics (Table 2)
│   ├── 06_stats_group_condition.py
│   ├── 07_stats_sensitivity.py    # Supp Tables S2/S3
│   ├── 08_stats_alt_alignment.py  # Supp Table S1
│   ├── 09_correlations_dds_clinical.py  # Table 3
│   └── run_all.py             # master orchestrator
├── figures/
│   ├── make_fig01_publication.py  # Fig 1: ERP + topomap + bar chart
│   ├── make_fig02_publication.py  # Fig 2: DDS fit + boxplots + PCA scree
│   ├── make_fig03_AB_publication.py  # Fig 3: AIS (Frontal, cACC)
│   ├── make_fig04_TE_publication.py  # Fig 4: TE frontal→cACC
│   └── make_fig05_concept_publication.py  # Fig 5: conceptual diagram
├── outputs/
│   ├── figures/               # .png and .svg outputs (git-ignored)
│   └── tables/                # .csv and .txt outputs (git-ignored)
├── environment.yml
├── CITATION.cff
├── LICENSE
└── .github/workflows/ci.yml
```

---

## Quick start

### 1. Set up the environment

```bash
conda env create -f environment.yml
conda activate dds-hayling
```

### 2. Verify the environment

```bash
python pipeline/00_validate_environment.py
```

This checks Python version, required packages, and all derivative input files.

### 3. Run the full pipeline

```bash
python pipeline/run_all.py
```

Or for statistics and figures only (derivatives already present):

```bash
python pipeline/run_all.py --stats-only
```

For figures only:

```bash
python pipeline/run_all.py --figs-only
```

Dry run (print commands without executing):

```bash
python pipeline/run_all.py --dry-run
```

---

## Data requirements

This pipeline reads pre-computed derivatives from the companion `DDS-Hayling/` data directory. The following files must exist relative to `../DDS-Hayling/` (one level above this repo):

| File | Description |
|------|-------------|
| `derivatives/dds_peak_aligned_n450/dds_n450_results.csv` | DDS fit parameters per trial |
| `derivatives/ais_n450_residual_r2pos/ais_n450_residual_r2pos_results.csv` | AIS values (R²>0 filter) |
| `derivatives/te_n450_residual_r2pos/te_n450_residual_r2pos_results.csv` | TE values (R²>0 filter) |
| `derivatives/hayling_N450_subject_table.csv` | Subject → group mapping |
| `derivatives/results_n450_clinical/n450_master_clinical_results_with_group.csv` | DDS params + clinical scores |
| `EEG/P*/HYL_*_90_Hz-2-ar-den-epo.fif` | Preprocessed epoch files (for Figs 1–2 only) |

Raw EEG and large derivative files are **not** included in this repository (see `.gitignore`).

---

## Figure → script mapping

| Figure | Script | Key data source |
|--------|--------|-----------------|
| Fig 1: ERP waveform, topomap, bar chart | `figures/make_fig01_publication.py` | Epoch .fif files |
| Fig 2: DDS fit, boxplots, PCA scree | `figures/make_fig02_publication.py` | `dds_n450_results.csv` |
| Fig 3: AIS Frontal + cACC | `figures/make_fig03_AB_publication.py` | `ais_n450_residual_r2pos_results.csv` |
| Fig 4: TE Frontal→cACC | `figures/make_fig04_TE_publication.py` | `te_n450_residual_r2pos_results.csv` |
| Fig 5: Conceptual network diagram | `figures/make_fig05_concept_publication.py` | Hardcoded paper values |

---

## Statistical results

### Primary results (Table 2)

| Measure | ROI | t | df | p | direction |
|---------|-----|---|----|----|-----------|
| AIS | Frontal | 2.13 | 46 | 0.038 | INIT > INHIB |
| AIS | cACC | −2.58 | 46 | 0.013 | INHIB > INIT |
| TE | Frontal→cACC | 2.46 | 46 | 0.018 | INIT > INHIB |

All p-values FDR-corrected (Benjamini-Hochberg). No DDS parameter × clinical correlations survive FDR correction (Table 3).

### Participants

- N = 47 total (after excluding P4, P5, P19)
- CTL: n = 22 (BDI-II ≤ 8)
- DEP: n = 25 (BDI-II ≥ 14)

---

## Configuration

All analysis parameters are stored in `config/analysis.yaml`. This includes ROI channel definitions, N450 time windows, DDS model parameters, statistical thresholds, and figure style settings. Edit this file to modify any parameter; all pipeline scripts read from it automatically.

---

## Citation

If you use this code, please cite:

```bibtex
@software{dds_hayling_paper,
  title  = {DDS-Hayling: Analysis Pipeline},
  year   = {2026},
  url    = {https://github.com/<owner>/DDS-Hayling-paper}
}
```

See also `CITATION.cff` for structured citation metadata.

---

## License

MIT — see `LICENSE`.

#!/usr/bin/env python3
"""
07_stats_sensitivity.py
=======================
Supplementary sensitivity analyses.

  Table S2 — TE lag sensitivity
    Repeat frontal→cACC TE test at lag=1 and lag=4 samples using the
    unfiltered te_n450_results.csv (bins=4, raw signal — not residuals).
    Note: only one binning (k=4) was used in the original computation;
    multi-bin comparison is not available from archived derivatives.

  Table S3 — R² threshold sensitivity
    Repeat primary AIS/TE tests at R²>0.0, 0.1, 0.2, 0.3, 0.5.

Outputs (outputs/tables/):
  stats_sensitivity_bins.csv     — S2: lag sensitivity (proxy for TE estimation)
  stats_sensitivity_r2.csv       — S3: R² threshold sweep

Usage:
  python pipeline/07_stats_sensitivity.py
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline._utils import (
    CFG, OUT_TABLES, AIS_CSV, TE_CSV, DERIV_ROOT,
    load_subject_table,
    cohens_d_paired, sem, p_to_stars,
)

OUT_TABLES.mkdir(parents=True, exist_ok=True)
df_grp = load_subject_table()

print("=" * 60)
print("07_stats_sensitivity.py — Sensitivity Analyses")
print("=" * 60)

EXCLUDE = set(CFG["participants"]["exclude"])


def paired_test(pivot: pd.DataFrame) -> dict:
    n = len(pivot)
    if n < 5:
        return dict(n=n, t=np.nan, p=np.nan, d=np.nan)
    t, p = scipy_stats.ttest_rel(pivot["INIT"], pivot["INHIB"])
    d    = cohens_d_paired(pivot["INIT"].values, pivot["INHIB"].values)
    return dict(
        n=n,
        init_mean=round(pivot["INIT"].mean(), 4),
        inhib_mean=round(pivot["INHIB"].mean(), 4),
        delta=round(pivot["INIT"].mean() - pivot["INHIB"].mean(), 4),
        t=round(t, 4), p=round(p, 5), d=round(d, 3),
    )


# ---------------------------------------------------------------------------
# Table S2: TE lag sensitivity  (lag=1 vs lag=4, bins=4 throughout)
# ---------------------------------------------------------------------------
print("\n── Table S2: TE lag sensitivity (lag=1 vs lag=4, bins=4) ───")

# Use te_n450_results.csv which contains both lags and a 'bins' column
te_raw_path = DERIV_ROOT / "te_n450" / "te_n450_results.csv"

if not te_raw_path.exists():
    print(f"  [SKIP] {te_raw_path} not found — Table S2 not produced.")
    s2_rows = []
else:
    df_te_raw = pd.read_csv(te_raw_path)
    # Exclude bad subjects
    df_te_raw = df_te_raw[~df_te_raw["subject"].isin(EXCLUDE)]
    bins_val  = df_te_raw["bins"].unique()
    lags_avail = sorted(df_te_raw["lag_samples"].unique())
    print(f"  bins in data: {bins_val}  |  lags (samples): {lags_avail}")
    if len(bins_val) == 1:
        print(f"  Note: only bins={bins_val[0]} archived — multi-bin comparison N/A.")

    s2_rows = []
    for lag in lags_avail:
        df_lag = (df_te_raw[
            (df_te_raw["source_roi"] == "frontal") &
            (df_te_raw["target_roi"] == "cacc") &
            (df_te_raw["lag_samples"] == lag)
        ].groupby(["subject", "cond"])["te_bits"].mean().reset_index()
         .merge(df_grp, on="subject", how="left"))

        pivot = (df_lag.pivot_table(index="subject", columns="cond",
                                     values="te_bits")
                 .dropna(subset=["INIT", "INHIB"]))
        res = paired_test(pivot)
        s2_rows.append(dict(
            metric="TE", direction="frontal→cACC",
            bins=int(bins_val[0]), lag_samples=int(lag),
            **res,
        ))

    df_s2 = pd.DataFrame(s2_rows)
    df_s2["stars"] = df_s2["p"].apply(p_to_stars)
    print(df_s2.to_string(index=False))

s2_path = OUT_TABLES / "stats_sensitivity_bins.csv"
if s2_rows:
    pd.DataFrame(s2_rows).assign(
        stars=lambda d: d["p"].apply(p_to_stars)
    ).to_csv(s2_path, index=False)
    print(f"\n[OK] saved: {s2_path}")
else:
    # Write an empty placeholder with a note
    pd.DataFrame(
        [{"note": "te_n450_results.csv not found; only bins=4 was computed"}]
    ).to_csv(s2_path, index=False)
    print(f"[OK] saved (placeholder): {s2_path}")


# ---------------------------------------------------------------------------
# Table S3: R² threshold sensitivity
# ---------------------------------------------------------------------------
print("\n── Table S3: R² threshold sensitivity ───────────────────────")

df_dds_full = pd.read_csv(
    DERIV_ROOT / "dds_peak_aligned_n450" / "dds_n450_results.csv"
)
df_dds_full = df_dds_full[~df_dds_full["subject"].isin(EXCLUDE)]

df_ais_full = pd.read_csv(AIS_CSV)
df_te_full  = pd.read_csv(TE_CSV)

R2_THRESHOLDS = [0.0, 0.1, 0.2, 0.3, 0.5]
ROIS_FOCUS    = ["frontal", "cacc"]

s3_rows = []
for thr in R2_THRESHOLDS:
    valid_uids = set(df_dds_full.loc[df_dds_full["r2"] > thr, "trial_uid"])

    for roi in ROIS_FOCUS:
        sub = (df_ais_full[
                   (df_ais_full["trial_uid"].isin(valid_uids)) &
                   (df_ais_full["roi"] == roi)
               ].groupby(["subject", "cond"])["ais_bits"].mean().reset_index()
               .merge(df_grp, on="subject", how="left"))
        pivot = (sub.pivot_table(index="subject", columns="cond", values="ais_bits")
                 .dropna(subset=["INIT", "INHIB"]))
        res = paired_test(pivot)
        s3_rows.append(dict(
            metric="AIS", target=roi, r2_threshold=thr, **res
        ))

    df_te_t = (df_te_full[
        (df_te_full["trial_uid"].isin(valid_uids)) &
        (df_te_full["source_roi"] == "frontal") &
        (df_te_full["target_roi"] == "cacc") &
        (df_te_full["lag"] == 1)
    ].groupby(["subject", "cond"])["te_bits"].mean().reset_index()
     .merge(df_grp, on="subject", how="left"))
    pivot_te = (df_te_t.pivot_table(index="subject", columns="cond",
                                     values="te_bits")
                .dropna(subset=["INIT", "INHIB"]))
    res_te = paired_test(pivot_te)
    s3_rows.append(dict(metric="TE", target="frontal→cACC",
                        r2_threshold=thr, **res_te))

df_s3 = pd.DataFrame(s3_rows)
df_s3["stars"] = df_s3["p"].apply(p_to_stars)
print(df_s3.to_string(index=False))

s3_path = OUT_TABLES / "stats_sensitivity_r2.csv"
df_s3.to_csv(s3_path, index=False)
print(f"\n[OK] saved: {s3_path}")
print("\nDone.")

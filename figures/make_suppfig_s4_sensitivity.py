#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_suppfig_s4_sensitivity.py
===============================
Supplementary Figure S4: Sensitivity analysis visualizations

Graphical complement to Tables S2 and S3.

  (A) R² threshold sensitivity — t-statistic and p-value for
      Frontal AIS, cACC AIS, and TE frontal→cACC as R² threshold
      increases from 0.0 to 0.5. Demonstrates result stability.

  (B) TE lag sensitivity — mean TE frontal→cACC for INIT and INHIB
      at lag=1 and lag=4 (bins=4 throughout; from te_n450_results.csv).

  (C) DDS R² distribution by condition — histogram of R² values
      for INIT vs INHIB, showing fit quality is comparable between
      conditions.

  (D) N valid trials per subject (R²>0) — bar chart per subject,
      INIT and INHIB side by side; confirms sufficient trial counts.

Output: outputs/figures/suppfig_s4_sensitivity.{png,svg}
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline._utils import (
    CFG, OUT_TABLES, DERIV_ROOT, DDS_CSV,
    load_subject_table, load_dds,
    ROI_ORDER, sem, add_panel_label, apply_style, save_figure,
)

apply_style()

C_INIT    = CFG["figures"]["colors"]["init"]
C_INHIB   = CFG["figures"]["colors"]["inhib"]
ALPHA_FDR = CFG["statistics"]["alpha"]

print("Loading data for Figure S4…")
df_grp = load_subject_table()
EXCLUDE = set(CFG["participants"]["exclude"])

# ── Load sensitivity tables ───────────────────────────────────────────────────
s3_path = OUT_TABLES / "stats_sensitivity_r2.csv"
s2_path = OUT_TABLES / "stats_sensitivity_bins.csv"

if not s3_path.exists():
    print(f"  [warn] {s3_path} not found — run pipeline/07_stats_sensitivity.py first")
    sys.exit(1)

df_s3 = pd.read_csv(s3_path)
df_s2 = pd.read_csv(s2_path) if s2_path.exists() else pd.DataFrame()

# ── Load DDS for panels C, D ──────────────────────────────────────────────────
df_dds_all = pd.read_csv(DDS_CSV)
df_dds_all = df_dds_all[~df_dds_all["subject"].isin(EXCLUDE)].copy()

# ── Figure layout ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(15, 10), facecolor="white",
                          gridspec_kw={"hspace": 0.42, "wspace": 0.38})
ax_a, ax_b = axes[0]
ax_c, ax_d = axes[1]

# ── Panel A: R² threshold sensitivity (t-statistic) ──────────────────────────
TARGETS_A = [
    ("AIS",    "frontal",     "Frontal AIS",   C_INIT,  "o-"),
    ("AIS",    "cacc",        "cACC AIS",      C_INHIB, "s-"),
    ("TE",     "frontal→cACC","TE front.→cACC","#2a9d8f","^-"),
]

ax_a2 = ax_a.twinx()   # right y-axis for p-values

for metric, target, label, col, style in TARGETS_A:
    sub = df_s3[(df_s3["metric"] == metric) & (df_s3["target"] == target)]
    if sub.empty:
        continue
    ax_a.plot(sub["r2_threshold"], sub["t"],
              style, color=col, lw=1.8, ms=6, label=label, zorder=4)
    ax_a2.plot(sub["r2_threshold"], sub["p"],
               style, color=col, lw=1.0, ms=4, alpha=0.4, zorder=2)

ax_a2.axhline(ALPHA_FDR, color="0.5", ls="--", lw=1.0)
ax_a2.set_ylabel("p-value (faint lines)", fontsize=10, color="0.5")
ax_a2.tick_params(axis="y", labelsize=9, labelcolor="0.5")
ax_a2.set_ylim(0, 0.12)

ax_a.axhline(0, color="black", lw=0.5, alpha=0.3)
ax_a.set_xlabel("R² threshold", fontsize=13)
ax_a.set_ylabel("t-statistic", fontsize=13)
ax_a.tick_params(labelsize=11)
ax_a.spines["top"].set_visible(False)
ax_a.legend(frameon=False, fontsize=10, loc="upper left")
add_panel_label(ax_a, "(A)")

# ── Panel B: TE lag sensitivity ───────────────────────────────────────────────
if not df_s2.empty and "lag_samples" in df_s2.columns:
    lag_labels = {1: "lag=1 (~4 ms)", 4: "lag=4 (~16 ms)"}
    x_b = np.arange(2)
    bar_w_b = 0.35
    for j, (cond_col, color, label) in enumerate([
        ("init_mean",  "#bdc3c7", "Initiation"),
        ("inhib_mean", C_INHIB,   "Inhibition"),
    ]):
        means_b = df_s2[cond_col].values if cond_col in df_s2.columns else []
        if len(means_b) == len(x_b):
            ax_b.bar(x_b + (j-0.5)*bar_w_b, means_b, bar_w_b,
                     color=color, label=label)

    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels([lag_labels.get(int(l), str(l))
                           for l in df_s2["lag_samples"].values], fontsize=11)
    ax_b.set_ylabel("TE frontal→cACC (bits)", fontsize=13)
    ax_b.set_xlabel("TE estimation lag", fontsize=13)
    ax_b.tick_params(axis="y", labelsize=11)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.legend(frameon=False, fontsize=11)
    note = f"bins = {df_s2['bins'].iloc[0]}" if "bins" in df_s2.columns else ""
    ax_b.text(0.97, 0.97, note, transform=ax_b.transAxes, fontsize=9,
              ha="right", va="top", color="0.5")
else:
    ax_b.text(0.5, 0.5, "Lag sensitivity data not available",
              transform=ax_b.transAxes, ha="center", va="center",
              fontsize=11, color="0.5")
    ax_b.axis("off")
add_panel_label(ax_b, "(B)", x=-0.12)

# ── Panel C: R² distribution by condition ─────────────────────────────────────
r2_init  = df_dds_all.loc[df_dds_all["cond"] == "INIT",  "r2"].dropna()
r2_inhib = df_dds_all.loc[df_dds_all["cond"] == "INHIB", "r2"].dropna()

bins_hist = np.linspace(-0.5, 1.0, 40)
ax_c.hist(r2_init,  bins=bins_hist, color=C_INIT,  alpha=0.6, label="Initiation",
          density=True, histtype="stepfilled")
ax_c.hist(r2_inhib, bins=bins_hist, color=C_INHIB, alpha=0.6, label="Inhibition",
          density=True, histtype="stepfilled")
ax_c.axvline(0, color="black", lw=1.0, ls="--", alpha=0.6)
ax_c.set_xlabel("DDS goodness-of-fit (R²)", fontsize=13)
ax_c.set_ylabel("Density", fontsize=13)
ax_c.tick_params(labelsize=11)
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.legend(frameon=False, fontsize=11)

pct_init  = 100 * (r2_init  > 0).mean()
pct_inhib = 100 * (r2_inhib > 0).mean()
ax_c.text(0.97, 0.95,
          f"INIT  R²>0: {pct_init:.1f}%\nINHIB R²>0: {pct_inhib:.1f}%",
          transform=ax_c.transAxes, fontsize=9, ha="right", va="top",
          bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc",
                    alpha=0.85, lw=0.8))
add_panel_label(ax_c, "(C)")

# ── Panel D: N valid trials per subject ───────────────────────────────────────
subj_counts = (df_dds_all[df_dds_all["r2"] > 0]
               .groupby(["subject", "cond"])["trial"].count()
               .unstack("cond").fillna(0))
subj_sorted = subj_counts.mean(axis=1).sort_values().index.tolist()
subj_counts = subj_counts.loc[subj_sorted]

x_d = np.arange(len(subj_sorted))
bw  = 0.38
ax_d.bar(x_d - bw/2, subj_counts.get("INIT",  pd.Series(dtype=float)).values,
         bw, color=C_INIT,  label="Initiation", alpha=0.85)
ax_d.bar(x_d + bw/2, subj_counts.get("INHIB", pd.Series(dtype=float)).values,
         bw, color=C_INHIB, label="Inhibition", alpha=0.85)

ax_d.set_xticks(x_d[::5])
ax_d.set_xticklabels(subj_sorted[::5], fontsize=8, rotation=45, ha="right")
ax_d.set_xlabel("Subject", fontsize=13)
ax_d.set_ylabel("N valid DDS fits (R²>0)", fontsize=13)
ax_d.tick_params(axis="y", labelsize=11)
ax_d.spines["top"].set_visible(False)
ax_d.spines["right"].set_visible(False)
ax_d.legend(frameon=False, fontsize=11, loc="upper left")
add_panel_label(ax_d, "(D)")

save_figure(fig, "suppfig_s4_sensitivity")
plt.close(fig)
print("Done.")

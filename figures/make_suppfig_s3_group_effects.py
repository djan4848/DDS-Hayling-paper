#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_suppfig_s3_group_effects.py
=================================
Supplementary Figure S3: Group (CTL vs DEP) effects on AIS and TE

Four panels showing CTL/DEP breakdown across conditions:

  (A) Frontal AIS — CTL vs DEP, INIT and INHIB separately
  (B) cACC AIS   — CTL vs DEP, INIT and INHIB separately
  (C) TE frontal→cACC — CTL vs DEP, INIT and INHIB
  (D) AIS INIT−INHIB difference score — CTL vs DEP, all ROIs

Individual data points shown with group means ± SEM.

Output: outputs/figures/suppfig_s3_group_effects.{png,svg}
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats as scipy_stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline._utils import (
    CFG, load_ais, load_te, load_subject_table,
    ROI_ORDER, ROI_DISPLAY_NAMES, sem,
    add_panel_label, apply_style, save_figure,
)

apply_style()

C_CTL      = CFG["figures"]["colors"]["ctl"]
C_DEP      = CFG["figures"]["colors"]["dep"]
ALPHA_PT   = 0.45
ALPHA_LINE = CFG["figures"]["alpha_subject_line"]
COND_ORDER = ["INIT", "INHIB"]
COND_LABELS = {"INIT": "Initiation", "INHIB": "Inhibition"}
X_POS       = {c: i for i, c in enumerate(COND_ORDER)}

print("Loading data for Figure S3…")
df_ais = load_ais(r2_positive_only=True)
df_te  = load_te(lag=1)
df_grp = load_subject_table()

ais_subj = (df_ais.groupby(["subject", "cond", "roi"])["ais_bits"]
            .mean().reset_index().merge(df_grp, on="subject", how="left"))
te_fc    = (df_te[(df_te["source_roi"] == "frontal") &
                   (df_te["target_roi"] == "cacc")]
            .groupby(["subject", "cond"])["te_bits"]
            .mean().reset_index().merge(df_grp, on="subject", how="left"))

fig, axes = plt.subplots(2, 2, figsize=(12, 9.5), constrained_layout=True,
                          facecolor="white")
ax_a, ax_b = axes[0]
ax_c, ax_d = axes[1]


def panel_cond_group(ax, data, value_col, roi_or_dir, panel_label, ylabel):
    """INIT / INHIB on x-axis, CTL/DEP colour-coded, with means ± SEM."""
    pivot = (data.pivot_table(index=["subject", "group"],
                               columns="cond", values=value_col)
             .reset_index().dropna(subset=COND_ORDER))

    # Individual subject lines
    for _, row in pivot.iterrows():
        col = C_CTL if row["group"] == "CTL" else C_DEP
        ax.plot([0, 1], [row["INIT"], row["INHIB"]],
                color=col, lw=0.7, alpha=ALPHA_LINE, zorder=2)

    # Individual dots
    for grp, col in [("CTL", C_CTL), ("DEP", C_DEP)]:
        gp = pivot[pivot["group"] == grp]
        for cond in COND_ORDER:
            ax.scatter([X_POS[cond]] * len(gp), gp[cond].values,
                       c=col, s=18, alpha=ALPHA_PT, zorder=3)

    # Group means ± SEM
    for grp, col in [("CTL", C_CTL), ("DEP", C_DEP)]:
        gp = pivot[pivot["group"] == grp]
        for cond in COND_ORDER:
            vals = gp[cond].values
            ax.errorbar(X_POS[cond], vals.mean(), yerr=sem(vals),
                        fmt="o", color=col, ms=10, lw=2.2,
                        capsize=5, capthick=2.0, zorder=6,
                        markeredgecolor="white", markeredgewidth=1.0)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([COND_LABELS["INIT"], COND_LABELS["INHIB"]])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.tick_params(axis="x", length=0)

    # Welch t-test CTL vs DEP within each condition
    for cond in COND_ORDER:
        ctl = pivot.loc[pivot["group"] == "CTL", cond].values
        dep = pivot.loc[pivot["group"] == "DEP", cond].values
        if len(ctl) >= 3 and len(dep) >= 3:
            _, p = scipy_stats.ttest_ind(ctl, dep, equal_var=False)
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            if star:
                ymax = pivot[COND_ORDER].values.max()
                ymin = pivot[COND_ORDER].values.min()
                rng  = ymax - ymin
                ax.text(X_POS[cond], ymax + rng * 0.07, star,
                        ha="center", va="bottom", fontsize=13, color="black")

    # Auto-range y
    all_vals = pivot[COND_ORDER].values.flatten()
    ypad = np.ptp(all_vals) * 0.15
    ax.set_ylim(all_vals.min() - ypad, all_vals.max() + ypad * 1.8)

    ax.text(-0.15, 1.04, panel_label,
            transform=ax.transAxes, fontsize=14, fontweight="bold",
            va="top", ha="left")
    ax.text(0.97, 0.97, roi_or_dir,
            transform=ax.transAxes, fontsize=9, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc",
                      alpha=0.85, lw=0.8))


# Panels A, B, C
panel_cond_group(ax_a, ais_subj[ais_subj["roi"] == "frontal"],
                 "ais_bits", "Frontal ROI AIS", "(A)", "AIS (bits)")
panel_cond_group(ax_b, ais_subj[ais_subj["roi"] == "cacc"],
                 "ais_bits", "cACC AIS", "(B)", "AIS (bits)")
panel_cond_group(ax_c, te_fc,
                 "te_bits", "TE Frontal→cACC", "(C)", "TE (bits)")

# ── Panel D: INIT−INHIB difference scores, all ROIs, CTL vs DEP ──────────────
rng_order = ROI_ORDER
x_d = np.arange(len(rng_order))
bar_w_d = 0.35

for offset, grp, col, label in [
    (-bar_w_d/2, "CTL", C_CTL, "CTL"),
    (+bar_w_d/2, "DEP", C_DEP, "DEP"),
]:
    means_d, sems_d = [], []
    for roi in rng_order:
        sub = ais_subj[ais_subj["roi"] == roi]
        pivot = (sub.pivot_table(index="subject", columns="cond", values="ais_bits")
                 .reset_index().dropna(subset=COND_ORDER))
        pivot = pivot.merge(df_grp, on="subject", how="left")
        diff = (pivot.loc[pivot["group"] == grp, "INIT"].values -
                pivot.loc[pivot["group"] == grp, "INHIB"].values)
        means_d.append(diff.mean() if len(diff) else np.nan)
        sems_d.append(sem(diff) if len(diff) > 1 else np.nan)

    ax_d.bar(x_d + offset, means_d, bar_w_d, color=col, label=label,
             yerr=sems_d, capsize=4, error_kw=dict(elinewidth=1.1))

ax_d.axhline(0, color="black", lw=0.7, alpha=0.5)
ax_d.set_xticks(x_d)
ax_d.set_xticklabels([ROI_DISPLAY_NAMES[r] for r in rng_order], fontsize=11)
ax_d.set_xlabel("ROI", fontsize=13)
ax_d.set_ylabel("AIS  INIT − INHIB (bits)", fontsize=13)
ax_d.tick_params(axis="y", labelsize=11)
ax_d.spines["top"].set_visible(False)
ax_d.spines["right"].set_visible(False)
ax_d.legend(frameon=False, fontsize=11, loc="upper right")
ax_d.text(-0.15, 1.04, "(D)",
          transform=ax_d.transAxes, fontsize=14, fontweight="bold",
          va="top", ha="left")

# Shared legend
patch_ctl = mpatches.Patch(facecolor=C_CTL, label=f"CTL  (n={CFG['participants']['n_ctl']})")
patch_dep = mpatches.Patch(facecolor=C_DEP, label=f"DEP  (n={CFG['participants']['n_dep']})")
fig.legend(handles=[patch_ctl, patch_dep],
           loc="upper center", ncol=2, frameon=False, fontsize=11,
           bbox_to_anchor=(0.5, 1.02))

save_figure(fig, "suppfig_s3_group_effects")
plt.close(fig)
print("Done.")

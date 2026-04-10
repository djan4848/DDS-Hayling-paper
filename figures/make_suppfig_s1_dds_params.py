#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_suppfig_s1_dds_params.py
==============================
Supplementary Figure S1: DDS parameter distributions

Violin + box plots showing the empirical distribution of all 8 DDS
parameters (A₁, γ₁, f₁, φ₁, A₂, γ₂, f₂, φ₂) across the 4 ROIs
(Frontal, cACC, LH, RH). Only valid fits (R²>0) included.

Layout: 2 rows × 4 columns.  One subplot per parameter.
        X-axis: 4 ROIs.  Y-axis: parameter value.
        Bold panel labels (A)–(H) top-left each subplot.
        Clipped at 1st–99th percentile per parameter for legibility.

Output: outputs/figures/suppfig_s1_dds_params.{png,svg}
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
    CFG, load_dds, ROI_ORDER, ROI_DISPLAY_NAMES,
    apply_style, save_figure,
)

apply_style()

PARAMS = CFG["dds"]["params"]
PARAM_LABELS = {
    "A1":     r"$A_1$ (V)",
    "gamma1": r"$\gamma_1$ (s⁻¹)",
    "f1":     r"$f_1$ (Hz)",
    "phi1":   r"$\varphi_1$ (rad)",
    "A2":     r"$A_2$ (V)",
    "gamma2": r"$\gamma_2$ (s⁻¹)",
    "f2":     r"$f_2$ (Hz)",
    "phi2":   r"$\varphi_2$ (rad)",
}
PANEL_LABELS = list("ABCDEFGH")
ROI_LABELS   = [ROI_DISPLAY_NAMES[r] for r in ROI_ORDER]

# Neutral colour per ROI for violin fill
ROI_COLORS = ["#6fa8c7", "#c47070", "#7bbf8e", "#c4a55a"]

print("Loading DDS data…")
df = load_dds(r2_positive_only=True)
df = df[df["roi"].isin(ROI_ORDER)].copy()
valid_subjects = {f for f in df["subject"].unique()}
n_subj  = len(valid_subjects)
n_fits  = len(df)
print(f"  N={n_subj} subjects, n={n_fits:,} trial–ROI fits (R²>0)")

fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), facecolor="white",
                          gridspec_kw={"hspace": 0.42, "wspace": 0.38})
axes_flat = axes.flatten()

for ax_idx, param in enumerate(PARAMS):
    ax     = axes_flat[ax_idx]
    label  = PARAM_LABELS[param]
    plabel = PANEL_LABELS[ax_idx]

    # Per-ROI data with 1%–99% clip for legibility
    q01 = df[param].quantile(0.01)
    q99 = df[param].quantile(0.99)
    data_per_roi = [
        df.loc[(df["roi"] == roi) & df[param].between(q01, q99), param].values
        for roi in ROI_ORDER
    ]

    # Violin
    parts = ax.violinplot(
        data_per_roi,
        positions=range(len(ROI_ORDER)),
        widths=0.7,
        showmedians=False,
        showextrema=False,
    )
    for pc, col in zip(parts["bodies"], ROI_COLORS):
        pc.set_facecolor(col)
        pc.set_alpha(0.55)
        pc.set_edgecolor("none")

    # Box overlay (thin, dark)
    bp = ax.boxplot(
        data_per_roi,
        positions=range(len(ROI_ORDER)),
        widths=0.18,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.8),
        boxprops=dict(facecolor="white", linewidth=1.0),
        whiskerprops=dict(linewidth=0.9, color="#444444"),
        capprops=dict(linewidth=0.9, color="#444444"),
    )

    ax.set_xticks(range(len(ROI_ORDER)))
    ax.set_xticklabels(ROI_LABELS, fontsize=10)
    ax.set_ylabel(label, fontsize=11)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # panel label
    ax.text(-0.14, 1.06, f"({plabel})",
            transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="top", ha="left")

fig.text(
    0.5, 0.01,
    f"DDS parameters (R²>0 valid fits only)  |  "
    f"N = {n_subj} subjects, n = {n_fits:,} trial–ROI fits  |  "
    f"Y-axes clipped at 1st–99th percentile",
    ha="center", fontsize=10, style="italic",
    bbox=dict(facecolor="white", alpha=0.85, edgecolor="#cccccc",
              boxstyle="round,pad=0.3"),
)

save_figure(fig, "suppfig_s1_dds_params")
plt.close(fig)
print("Done.")

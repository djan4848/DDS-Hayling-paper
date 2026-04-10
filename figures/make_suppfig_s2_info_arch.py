#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_suppfig_s2_info_arch.py
=============================
Supplementary Figure S2: Full residual information architecture

Four panels showing the complete AIS and TE landscape across all ROIs:

  (A) AIS — all 4 ROIs × 2 conditions, grouped bar chart ± SEM
  (B) TE heatmap — mean TE (bits) for all source→target pairs,
      separately for INIT and INHIB; difference INIT − INHIB
  (C) TE condition comparison — all directions, paired dot plot
      INIT vs INHIB (overall mean across subjects)
  (D) AIS − TE correlation — scatter of subject-mean AIS (frontal)
      vs subject-mean TE (frontal→cACC) for each condition

Output: outputs/figures/suppfig_s2_info_arch.{png,svg}
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline._utils import (
    CFG, load_ais, load_te, load_subject_table,
    ROI_ORDER, ROI_DISPLAY_NAMES, sem, p_to_stars,
    add_panel_label, apply_style, save_figure,
)

apply_style()

C_INIT  = CFG["figures"]["colors"]["init"]
C_INHIB = CFG["figures"]["colors"]["inhib"]
C_CTL   = CFG["figures"]["colors"]["ctl"]
C_DEP   = CFG["figures"]["colors"]["dep"]

print("Loading data for Figure S2…")
df_ais  = load_ais(r2_positive_only=True)
df_te   = load_te(lag=1)
df_grp  = load_subject_table()

# Subject-level aggregates
ais_subj = (df_ais.groupby(["subject", "cond", "roi"])["ais_bits"]
            .mean().reset_index().merge(df_grp, on="subject", how="left"))

te_subj  = (df_te.groupby(["subject", "cond", "source_roi", "target_roi"])["te_bits"]
            .mean().reset_index().merge(df_grp, on="subject", how="left"))

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 13), facecolor="white")
gs  = gridspec.GridSpec(2, 2, figure=fig,
                        hspace=0.42, wspace=0.36,
                        left=0.07, right=0.97,
                        bottom=0.07, top=0.94)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

# ── Panel A: AIS all ROIs × conditions ───────────────────────────────────────
bar_w  = 0.35
x_pos  = np.arange(len(ROI_ORDER))

for offset, cond, color, label in [
    (-bar_w/2, "INIT",  "#bdc3c7", "Initiation"),
    (+bar_w/2, "INHIB", C_INHIB,   "Inhibition"),
]:
    means, sems = [], []
    for roi in ROI_ORDER:
        vals = ais_subj.loc[(ais_subj["cond"] == cond) &
                             (ais_subj["roi"]  == roi), "ais_bits"].values
        means.append(vals.mean())
        sems.append(sem(vals))
    ax_a.bar(x_pos + offset, means, bar_w, color=color, label=label,
             yerr=sems, capsize=4, error_kw=dict(elinewidth=1.1))

ax_a.set_xticks(x_pos)
ax_a.set_xticklabels([ROI_DISPLAY_NAMES[r] for r in ROI_ORDER], fontsize=11)
ax_a.set_ylabel("AIS (bits)", fontsize=13)
ax_a.set_xlabel("ROI", fontsize=13)
ax_a.tick_params(axis="y", labelsize=11)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.legend(frameon=False, fontsize=11, loc="lower right")
add_panel_label(ax_a, "(A)")

# ── Panel B: TE heatmap (INIT − INHIB difference) ────────────────────────────
# Build matrix: rows = source ROIs, cols = target ROIs
all_sources = sorted(df_te["source_roi"].unique())
all_targets = sorted(df_te["target_roi"].unique())

def te_matrix(cond):
    m = np.full((len(all_sources), len(all_targets)), np.nan)
    for si, src in enumerate(all_sources):
        for ti, tgt in enumerate(all_targets):
            sub = te_subj.loc[(te_subj["cond"] == cond) &
                               (te_subj["source_roi"] == src) &
                               (te_subj["target_roi"] == tgt), "te_bits"]
            if len(sub) > 0:
                m[si, ti] = sub.mean()
    return m

mat_init  = te_matrix("INIT")
mat_inhib = te_matrix("INHIB")
mat_diff  = mat_init - mat_inhib

vlim = np.nanmax(np.abs(mat_diff)) * 1.1
im = ax_b.imshow(mat_diff, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")

source_labels = [ROI_DISPLAY_NAMES.get(s, s) for s in all_sources]
target_labels = [ROI_DISPLAY_NAMES.get(t, t) for t in all_targets]
ax_b.set_xticks(range(len(all_targets))); ax_b.set_xticklabels(target_labels, fontsize=11)
ax_b.set_yticks(range(len(all_sources))); ax_b.set_yticklabels(source_labels, fontsize=11)
ax_b.set_xlabel("Target ROI", fontsize=13)
ax_b.set_ylabel("Source ROI", fontsize=13)

# Annotate cells with numeric values
for si in range(len(all_sources)):
    for ti in range(len(all_targets)):
        val = mat_diff[si, ti]
        if not np.isnan(val):
            ax_b.text(ti, si, f"{val:.3f}", ha="center", va="center",
                      fontsize=9, color="black" if abs(val) < vlim*0.6 else "white")

cbar = plt.colorbar(im, ax=ax_b, shrink=0.85, pad=0.02)
cbar.set_label("TE INIT − INHIB (bits)", fontsize=10)
cbar.ax.tick_params(labelsize=9)
add_panel_label(ax_b, "(B)", x=-0.12)

# ── Panel C: TE all directions, paired INIT vs INHIB ─────────────────────────
dirs = (te_subj.assign(direction=te_subj["source_roi"] + "→" + te_subj["target_roi"])
        ["direction"].unique())
dirs_sorted = sorted(dirs)
x_c = np.arange(len(dirs_sorted))
bar_w_c = 0.35

for offset, cond, color, label in [
    (-bar_w_c/2, "INIT",  "#bdc3c7", "Initiation"),
    (+bar_w_c/2, "INHIB", C_INHIB,   "Inhibition"),
]:
    means_c, sems_c = [], []
    for d in dirs_sorted:
        src, tgt = d.split("→")
        vals = te_subj.loc[(te_subj["cond"] == cond) &
                            (te_subj["source_roi"] == src) &
                            (te_subj["target_roi"] == tgt), "te_bits"].values
        means_c.append(vals.mean() if len(vals) else np.nan)
        sems_c.append(sem(vals)    if len(vals) > 1 else np.nan)
    ax_c.bar(x_c + offset, means_c, bar_w_c, color=color, label=label,
             yerr=sems_c, capsize=3, error_kw=dict(elinewidth=1.0))

ax_c.set_xticks(x_c)
dir_display = [d.replace("frontal", "Frontal").replace("cacc", "cACC")
               .replace("lh", "LH").replace("rh", "RH")
               for d in dirs_sorted]
ax_c.set_xticklabels(dir_display, fontsize=9, rotation=30, ha="right")
ax_c.set_ylabel("TE (bits)", fontsize=13)
ax_c.set_xlabel("TE direction", fontsize=13)
ax_c.tick_params(axis="y", labelsize=11)
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.legend(frameon=False, fontsize=11, loc="upper right")
add_panel_label(ax_c, "(C)")

# ── Panel D: Subject-mean frontal AIS vs frontal→cACC TE ──────────────────────
ais_f = (ais_subj[ais_subj["roi"] == "frontal"]
         .pivot_table(index="subject", columns="cond", values="ais_bits")
         .rename(columns={"INIT": "ais_init", "INHIB": "ais_inhib"})
         .reset_index())

te_fc = (te_subj[(te_subj["source_roi"] == "frontal") &
                  (te_subj["target_roi"] == "cacc")]
         .pivot_table(index="subject", columns="cond", values="te_bits")
         .rename(columns={"INIT": "te_init", "INHIB": "te_inhib"})
         .reset_index())

scatter_df = ais_f.merge(te_fc, on="subject").merge(df_grp, on="subject")

for cond_a, cond_t, color, marker, label in [
    ("ais_init",  "te_init",  C_INIT,  "o", "Initiation"),
    ("ais_inhib", "te_inhib", C_INHIB, "s", "Inhibition"),
]:
    sub = scatter_df.dropna(subset=[cond_a, cond_t])
    ax_d.scatter(sub[cond_a], sub[cond_t], c=color, marker=marker,
                 s=30, alpha=0.65, label=label, zorder=3)
    # Trend line
    m, b = np.polyfit(sub[cond_a], sub[cond_t], 1)
    x_line = np.linspace(sub[cond_a].min(), sub[cond_a].max(), 50)
    ax_d.plot(x_line, m*x_line + b, color=color, lw=1.2, alpha=0.7)

ax_d.set_xlabel("AIS Frontal (bits)", fontsize=13)
ax_d.set_ylabel("TE Frontal→cACC (bits)", fontsize=13)
ax_d.tick_params(labelsize=11)
ax_d.spines["top"].set_visible(False)
ax_d.spines["right"].set_visible(False)
ax_d.legend(frameon=False, fontsize=11)
add_panel_label(ax_d, "(D)")

save_figure(fig, "suppfig_s2_info_arch")
plt.close(fig)
print("Done.")

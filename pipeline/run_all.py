#!/usr/bin/env python3
"""
run_all.py
==========
Master pipeline orchestrator for DDS-Hayling.

Runs all numbered pipeline steps in sequence and then generates all
publication figures.  Each step is run as a subprocess so that import
errors or crashes in individual steps are clearly isolated.

Flags:
  --dry-run       Print the command sequence without executing
  --stats-only    Run only steps 00 + 05–09 (skip ERP/DDS/AIS/TE checks)
  --figs-only     Run only step 00 and the figures directory scripts
  --no-figs       Run pipeline steps but skip figure generation
  --step N        Run a single step N (0–9) only
  --verbose       Show full stdout of each step

Usage:
  python pipeline/run_all.py
  python pipeline/run_all.py --dry-run
  python pipeline/run_all.py --stats-only
  python pipeline/run_all.py --figs-only
  python pipeline/run_all.py --step 5

Exit code:
  0 — all steps passed
  N — N steps failed
"""

import sys
import argparse
import subprocess
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_DIR = PROJECT_ROOT / "pipeline"
FIGURES_DIR  = PROJECT_ROOT / "figures"

# ── step definitions ──────────────────────────────────────────────────────────
PIPELINE_STEPS = [
    (0,  "pipeline/00_validate_environment.py", "Environment validation"),
    (1,  "pipeline/01_erp_grand_average.py",    "ERP grand averages"),
    (2,  "pipeline/02_dds_fit.py",              "DDS fit check"),
    (3,  "pipeline/03_info_dynamics_ais.py",    "AIS check"),
    (4,  "pipeline/04_info_dynamics_te.py",     "TE check"),
    (5,  "pipeline/05_stats_main.py",           "Primary statistics"),
    (6,  "pipeline/06_stats_group_condition.py","Group × condition stats"),
    (7,  "pipeline/07_stats_sensitivity.py",    "Sensitivity analyses"),
    (8,  "pipeline/08_stats_alt_alignment.py",  "Alt-alignment sensitivity"),
    (9,  "pipeline/09_correlations_dds_clinical.py", "DDS × clinical correlations"),
]

FIGURE_SCRIPTS = [
    # Main figures
    ("figures/make_fig01_publication.py",          "Figure 1 — ERP"),
    ("figures/make_fig02_publication.py",          "Figure 2 — DDS evidence"),
    ("figures/make_fig03_AB_publication.py",       "Figure 3 — AIS"),
    ("figures/make_fig04_TE_publication.py",       "Figure 4 — TE"),
    ("figures/make_fig05_concept_publication.py",  "Figure 5 — Concept diagram"),
    # Supplementary figures
    ("figures/make_suppfig_s1_dds_params.py",      "Supp Fig S1 — DDS param distributions"),
    ("figures/make_suppfig_s2_info_arch.py",       "Supp Fig S2 — Full info architecture"),
    ("figures/make_suppfig_s3_group_effects.py",   "Supp Fig S3 — Group effects"),
    ("figures/make_suppfig_s4_sensitivity.py",     "Supp Fig S4 — Sensitivity analysis"),
]


def run_script(rel_path: str, label: str, dry_run: bool, verbose: bool) -> bool:
    """Run a single script. Returns True on success."""
    full_path = PROJECT_ROOT / rel_path
    cmd = [sys.executable, str(full_path)]

    if dry_run:
        print(f"  [DRY-RUN] {label}")
        print(f"            {' '.join(cmd)}")
        return True

    print(f"\n{'─'*60}")
    print(f"  Running: {label}")
    print(f"  Script : {rel_path}")
    t0 = time.time()

    result = subprocess.run(
        cmd,
        cwd=str(PROJECT_ROOT),
        capture_output=not verbose,
        text=True,
    )

    elapsed = time.time() - t0
    if result.returncode == 0:
        print(f"  [PASS]  {label}  ({elapsed:.1f}s)")
        if verbose and result.stdout:
            print(result.stdout)
    else:
        print(f"  [FAIL]  {label}  (exit {result.returncode}, {elapsed:.1f}s)")
        if result.stdout:
            print("  STDOUT:")
            for line in result.stdout.splitlines()[-30:]:
                print(f"    {line}")
        if result.stderr:
            print("  STDERR:")
            for line in result.stderr.splitlines()[-20:]:
                print(f"    {line}")

    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(
        description="DDS-Hayling pipeline orchestrator",
    )
    parser.add_argument("--dry-run",    action="store_true")
    parser.add_argument("--stats-only", action="store_true")
    parser.add_argument("--figs-only",  action="store_true")
    parser.add_argument("--no-figs",    action="store_true")
    parser.add_argument("--step",       type=int, default=None,
                        help="Run a single step (0–9)")
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    # ── select which steps to run ─────────────────────────────────────────────
    if args.step is not None:
        steps_to_run = [(n, p, l) for n, p, l in PIPELINE_STEPS if n == args.step]
        if not steps_to_run:
            print(f"Error: no step {args.step} defined.")
            sys.exit(1)
        run_figs = False
    elif args.figs_only:
        steps_to_run = [(n, p, l) for n, p, l in PIPELINE_STEPS if n == 0]
        run_figs = True
    elif args.stats_only:
        steps_to_run = [(n, p, l) for n, p, l in PIPELINE_STEPS if n in (0, 5, 6, 7, 8, 9)]
        run_figs = not args.no_figs
    else:
        steps_to_run = PIPELINE_STEPS
        run_figs = not args.no_figs

    # ── run ───────────────────────────────────────────────────────────────────
    n_fail = 0
    t_start = time.time()

    if not args.dry_run:
        print("=" * 60)
        print("DDS-Hayling Analysis Pipeline")
        print("=" * 60)

    for n, path, label in steps_to_run:
        ok = run_script(path, f"Step {n:02d}: {label}", args.dry_run, args.verbose)
        if not ok:
            n_fail += 1

    if run_figs:
        if not args.dry_run:
            print(f"\n{'─'*60}")
            print("  Generating figures…")
        for path, label in FIGURE_SCRIPTS:
            ok = run_script(path, label, args.dry_run, args.verbose)
            if not ok:
                n_fail += 1

    # ── summary ───────────────────────────────────────────────────────────────
    elapsed = time.time() - t_start
    if not args.dry_run:
        print(f"\n{'='*60}")
        total = len(steps_to_run) + (len(FIGURE_SCRIPTS) if run_figs else 0)
        n_pass = total - n_fail
        status = "ALL PASSED" if n_fail == 0 else f"{n_fail} FAILED"
        print(f"  {status}  ({n_pass}/{total} steps, {elapsed:.0f}s total)")

    sys.exit(n_fail)


if __name__ == "__main__":
    main()

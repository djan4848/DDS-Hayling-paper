"""Frozen machinery, copied verbatim from the Hayling and MDD pipelines.
   - AIS/TE: 04_extract_ais_n450_residual.py / 05_extract_te_n450_residual.py (NBINS=8, LAG=1)
   - IAAFT : DEP_DECO/reproduce/nonlinearity/iaaft_analysis.py (N_SURR=100, N_ITER=20, SEED=42)
   - EPR   : DEP_DECO/scripts/utils/epr_estimator.py (n_bins=4, lag=1)
   No estimator logic is modified."""
import numpy as np, pandas as pd
from scipy.stats import entropy

NBINS = 8; LAG = 1; MIN_SAMPLES = 20
N_SURR = 100; N_ITER = 20; SEED = 42

# ---- frozen Hayling AIS ---------------------------------------------------
def safe_qcut(x, q):
    try: x_disc = pd.qcut(x, q=q, labels=False, duplicates="drop")
    except ValueError: return None
    if x_disc is None: return None
    x_disc = np.asarray(x_disc, dtype=float)
    if np.all(np.isnan(x_disc)): return None
    if (~np.isnan(x_disc)).sum() < 2: return None
    return x_disc.astype(int)

def entropy_from_counts(arr):
    _, counts = np.unique(arr, return_counts=True, axis=0)
    return entropy(counts, base=2)

def calculate_ais_shannon(series, bins=NBINS, lag=LAG):
    x = np.asarray(series, dtype=float)
    if len(x) <= lag + 1: return np.nan
    if np.allclose(np.std(x), 0): return np.nan
    x_disc = safe_qcut(x, q=bins)
    if x_disc is None: return np.nan
    past = x_disc[:-lag]; current = x_disc[lag:]
    if len(past) < 2 or len(current) < 2: return np.nan
    return max(0.0, float(entropy_from_counts(past) + entropy_from_counts(current)
                          - entropy_from_counts(np.stack((past, current), axis=1))))

def transfer_entropy_discrete(source, target, lag=LAG, bins=NBINS):
    x = np.asarray(source, dtype=float); y = np.asarray(target, dtype=float)
    if len(x) != len(y) or len(x) <= lag + 1: return np.nan
    if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0): return np.nan
    xd = safe_qcut(x, q=bins); yd = safe_qcut(y, q=bins)
    if xd is None or yd is None: return np.nan
    y_t = yd[lag:]; y_past = yd[:-lag]; x_past = xd[:-lag]
    if len(y_t) < 2: return np.nan
    te = (entropy_from_counts(np.stack((y_t, y_past), axis=1))
          + entropy_from_counts(np.stack((y_past, x_past), axis=1))
          - entropy_from_counts(y_past)
          - entropy_from_counts(np.stack((y_t, y_past, x_past), axis=1)))
    return max(0.0, float(te))

# ---- frozen EPR -----------------------------------------------------------
def compute_epr(signal, n_bins=4, lag=1):
    signal = np.asarray(signal, dtype=float); signal = signal[np.isfinite(signal)]
    if len(signal) < n_bins*n_bins*5 + 2*lag: return np.nan
    if np.std(signal) < 1e-12: return np.nan
    edges = np.unique(np.percentile(signal, np.linspace(0, 100, n_bins + 1)))
    if len(edges) < 3: return np.nan
    d = np.digitize(signal, edges[1:-1])
    counts = np.zeros((n_bins, n_bins))
    for i, j in zip(d[:-lag], d[lag:]):
        counts[min(int(i), n_bins-1), min(int(j), n_bins-1)] += 1.
    tot = counts.sum()
    if tot < 1: return np.nan
    P = counts/tot; epr = 0.
    for i in range(n_bins):
        for j in range(n_bins):
            if P[i,j] > 1e-14 and P[j,i] > 1e-14: epr += P[i,j]*np.log2(P[i,j]/P[j,i])
    return float(max(0., epr))

# ---- frozen IAAFT (univariate) -------------------------------------------
def iaaft_batch(signal, n_surr=N_SURR, n_iter=N_ITER, seed=SEED):
    rng = np.random.default_rng(seed); n = len(signal)
    sorted_amp = np.sort(signal); target_mags = np.abs(np.fft.rfft(signal))
    surrs = np.stack([rng.permutation(signal) for _ in range(n_surr)])
    for _ in range(n_iter):
        S = np.fft.rfft(surrs, axis=1)
        S_a = target_mags[np.newaxis,:]*np.exp(1j*np.angle(S))
        surrs = np.fft.irfft(S_a, n=n, axis=1)
        surrs = sorted_amp[np.argsort(np.argsort(surrs, axis=1), axis=1)]
    return surrs

# ---- bivariate IAAFT: SAME phase randomisation on both channels ----------
# Prichard & Theiler (1994) multivariate surrogates. Preserves both auto-spectra
# AND the cross-spectrum (hence all linear cross-predictability / Granger
# structure) while destroying higher-order (nonlinear) coupling.
def iaaft_batch_bivariate(x, y, n_surr=N_SURR, n_iter=N_ITER, seed=SEED):
    rng = np.random.default_rng(seed); n = len(x)
    sx, sy = np.sort(x), np.sort(y)
    mx, my = np.abs(np.fft.rfft(x)), np.abs(np.fft.rfft(y))
    perm = np.stack([rng.permutation(n) for _ in range(n_surr)])
    SX = x[perm].copy(); SY = y[perm].copy()
    for _ in range(n_iter):
        FX = np.fft.rfft(SX, axis=1); FY = np.fft.rfft(SY, axis=1)
        # common phase: average the two phase fields -> identical phase applied
        # to both channels, so the cross-spectrum phase is left at its original
        # value while the absolute phases are randomised.
        ph = np.angle(FX)
        SX = np.fft.irfft(mx[None,:]*np.exp(1j*ph), n=n, axis=1)
        SY = np.fft.irfft(my[None,:]*np.exp(1j*(ph + (np.angle(np.fft.rfft(y)) - np.angle(np.fft.rfft(x)))[None,:])), n=n, axis=1)
        SX = sx[np.argsort(np.argsort(SX, axis=1), axis=1)]
        SY = sy[np.argsort(np.argsort(SY, axis=1), axis=1)]
    return SX, SY

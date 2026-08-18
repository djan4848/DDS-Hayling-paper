# Ficheros a nivel de trial — NO versionados

Fuera de git por la política declarada en `.gitignore`: *«exclude trial-level files»*.

| Fichero | Tamaño | Contenido |
|---|---|---|
| `fix_verify.csv` | 1.3 MB | verificación por unidad de la reparación: R² recuperado vs almacenado y var(residuo)/var(ERP) para las 16 304 unidades |
| `fase2_unit.csv` | ~4 MB | (en scratchpad) AIS residuo/crudo, EPR, χ, ΔA_NL por unidad |
| `fase2_series.pkl` | ~30 MB | (en scratchpad) series temporales por unidad |

**Regeneración** (desde `code/`, con los derivados `*_FIXED` en su sitio):
`build2.py` → `fase2.py` → `ladder.py` → `sweeps.py` → `famB_robust.py` → `closing.py` → `chi_raw.py`.

Los agregados que sostienen **todas** las cifras del informe sí están versionados:
`fase2_families.csv`, `ladder.csv`, `ladder_diffs.csv`, `auc.csv`, `trialvar.csv`, `sweep_ais.csv`,
`sweep_te.csv`, `famB_robust.csv`, `closing_raw_iaaft.csv`, `chi_both.csv`, `chi_residuo.csv`,
`chi_EEG.csv`, `fase2_chi.csv`, `ladder_lmm_auc.txt`.

Los resúmenes de la reparación (mediana 2.1·10⁻¹⁰, máx 8.1·10⁻⁹, var ratio < 1 en el 100 %) están
en `../../hayling_fix/CHANGELOG.md` §3 y en `../report.md` §1.1.

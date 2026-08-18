# Ficheros a nivel de trial — NO versionados

Estos ficheros existen en disco pero quedan fuera de git, siguiendo la política declarada en
`.gitignore` del repositorio: *«Keep small summary CSVs, exclude trial-level files»*.

| Fichero | Tamaño | Contenido |
|---|---|---|
| `unit_battery.csv` | 4.3 MB | 16 304 unidades trial × ROI: AIS, EPR, r₁, χ, ΔA_NL (variantes `pub` y `cor`) |
| `unit_te.csv` | 926 KB | 3 836 pares frontal/cACC: TE real + nulo bivariado de cross-espectro |
| `keys.csv` | 268 KB | índice sujeto × condición × trial × ROI de las 16 304 unidades |

**Regeneración** (desde `code/`, con los derivados `*_FIXED` en su sitio):
`build_arrays.py` → `battery.py` → `infer.py`.
Semillas fijas (IAAFT `seed=42`, permutaciones `default_rng(2026)`), por lo que la regeneración es
determinista salvo el ruido de Monte Carlo ya documentado en el informe.

**Todas las cifras del informe son verificables sin estos ficheros**, a partir de
`battery_results.csv` (resumen por contraste), que sí está versionado.

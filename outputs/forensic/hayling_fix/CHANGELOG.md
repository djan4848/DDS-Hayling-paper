# CHANGELOG — reparación de la construcción del residuo DDS

**Fecha:** 2026-08-18 · **Alcance:** `DDS-Hayling/02_dds_modeling/03b_build_residual_trial_roi_timeseries.py`
y todos los derivados que dependen de él. **El manuscrito no se ha tocado.**

---

## 1. Qué estaba mal

`03_dds_peak_aligned_trial_by_tiral.py` ajusta el modelo DDS sobre la ventana **pico N450 ± 200 ms**
definiendo el tiempo del modelo **relativo al borde izquierdo de la ventana**:

```python
# 03_dds_peak_aligned_trial_by_tiral.py, línea ~266
t_win = times[mask] - tmin          # t = 0 en el borde izquierdo, t ∈ [0, 0.4] s
res   = fit_dds_blind(t_win, y_roi)
```

`03b_build_residual_trial_roi_timeseries.py` reconstruía el modelo con los parámetros ajustados
pero evaluándolo en **tiempo absoluto de época**, y sobre las 301 muestras de la época completa:

```python
# 03b (versión defectuosa)
t_ms = g[time_col_erp].to_numpy(dtype=float)
t    = t_ms / 1000.0                # t ∈ [-0.2, +1.0] s  ← BASE TEMPORAL DISTINTA
y_hat = dds_model(t, A1, gamma1, f1, phi1, A2, gamma2, f2, phi2)
residual = y - y_hat
```

Las dos bases temporales no coinciden. El término `A·exp(−γ·t)` con γ hasta 150 s⁻¹ y t = −0.2 s
diverge: `exp(30) ≈ 10¹³`. En consecuencia la columna `dds_fit` **no era el modelo ajustado** y la
columna `residual` **no era el residuo del ajuste**, ni dentro ni fuera de la ventana de ajuste.

### Evidencia (16 304 unidades trial × ROI, filtro R² > 0)

| Diagnóstico | Valor |
|---|---|
| R² del `dds_fit` publicado, **dentro de su propia ventana de ajuste** | mediana **−0.891** (IQR −1.740 … −0.371) |
| Unidades con R² < 0 dentro de su ventana | **92.2 %** |
| var(`residual` publicado) / var(ERP), dentro de la ventana | mediana **1.38**; > 1 en el **76.5 %** |
| Valores absolutos máximos del `residual` fuera de la ventana | hasta ~10³⁴ |
| Modelo reevaluado con la base temporal **correcta** vs. R² almacenado | mediana \|Δ\| = **2.1·10⁻⁹**, r = 0.9998 |

La última fila es la prueba diagnóstica: al reevaluar el modelo con `t` relativo a la ventana se
recupera el R² almacenado con precisión numérica. **Los parámetros DDS son correctos; el defecto
estaba únicamente en el paso de resta.**

### Defecto secundario (documentado, sin efecto en los datos vigentes)

El script de ajuste llama a `fit_dds_blind`, cuyo modelo (`dds_model_free`) incluye un **offset C**;
`C` no se guarda en `dds_n450_results.csv` y `03b` reconstruye el modelo **sin C**. En la práctica
es inocuo: el CSV vigente se reproduce exactamente con el modelo sin offset (es decir, procede de
`fit_dds_phi`, la línea comentada). Queda como divergencia código/datos a sellar en el repositorio.

---

## 2. Qué se ha cambiado

`03b_build_residual_trial_roi_timeseries.py` (respaldo del original en `*.pre_fix.bak`):

1. **Base temporal.** El modelo se evalúa en `t = times[mask] − (peak − 0.200)`, idéntica a la del
   ajuste.
2. **Ventana.** La salida se restringe a la ventana de ajuste (100–101 muestras), sin extrapolar.
3. **Máscara exacta.** La ventana se reconstruye sobre la **rejilla temporal float de MNE**
   (leída de los propios `-epo.fif`; idéntica en los 50 sujetos), no sobre `time_ms` redondeado.
   Motivo: el ajuste guardó `n_samples = 100` en el **12.2 %** de los trials por error de
   redondeo en `times >= peak − 0.200`; reconstruir desde el ms redondeado dejaba una muestra de
   borde de más en esos trials. Con la rejilla exacta el `n_samples` reconstruido coincide con el
   del ajuste en el **100 %** de los 17 496 trials.
4. **Dos variantes** en una sola pasada: con filtro `R² > 0` y sin filtro (necesaria para la
   sensibilidad de alineación/umbral del paper).
5. Columna nueva `t_win_s` (tiempo del modelo) para trazabilidad.

`04_extract_ais_n450_residual.py` y `05_extract_te_n450_residual.py` (respaldos `*.pre_fix.bak`):
**no se ha tocado la lógica de los estimadores.** Solo se han parametrizado entrada/salida y
`NBINS`/`LAG` mediante variables de entorno, con **valores por defecto idénticos a los originales**
(`NBINS = 8`, `LAG = 1`), para poder correr los barridos de sensibilidad sin duplicar código.

---

## 3. Verificación del arreglo

Sobre las 16 304 unidades del derivado corregido:

| Verificación | Criterio | Resultado |
|---|---|---|
| R² recuperado = R² almacenado | \|Δ\| ~ 2·10⁻⁹ | mediana **2.1·10⁻¹⁰**, p95 9.3·10⁻¹⁰, **máx 8.1·10⁻⁹**; 0 % por encima de 10⁻⁶ |
| var(residuo) / var(ERP) | < 1 | mediana **0.402**, p95 0.805, **máx 0.9992**; **< 1 en el 100 %** |
| Longitud de ventana | 100–101 muestras | 100 % coincide con `n_samples` del ajuste |

---

## 4. Qué cambia en los derivados

Directorios renombrados a `*_BUGGED` (conservados, **no borrados**, ver `derivatives/_README_BUGGED.txt`):

```
trial_roi_timeseries_residual_BUGGED         trial_roi_timeseries_residual_r2pos_BUGGED
ais_n450_residual_BUGGED                     ais_n450_residual_r2pos_BUGGED
te_n450_residual_BUGGED                      te_n450_residual_r2pos_BUGGED
```

Sustitutos correctos, regenerados con la misma maquinaria congelada:

```
trial_roi_timeseries_residual_FIXED          trial_roi_timeseries_residual_r2pos_FIXED
ais_n450_residual_FIXED                      ais_n450_residual_r2pos_FIXED
te_n450_residual_FIXED                       te_n450_residual_r2pos_FIXED
```

`DDS-Hayling-paper/config/analysis.yaml` apunta ahora a los derivados `_FIXED`
(respaldo del original en `config/analysis.yaml.pre_fix.bak`).

### Consecuencia estadística (detalle completo en `../hayling_incremental/report.md`)

Cambian todos los números de las Tablas 2 y 3 y de las Figuras 3 y 4. Adelanto:
el contraste **AIS frontal** sobrevive y se refuerza; **AIS cACC** y **TE frontal→cACC** dejan de
alcanzar significación; aparece un **TE cACC→frontal** en dirección opuesta que el manuscrito no
reporta.

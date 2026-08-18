# LÍNEA RETIRADA Y ARCHIVADA — NO REANUDAR. Ver §5 Justificaciones.

**Entorno:** `DDS-HAYLING/` (repos `DDS-Hayling` y `DDS-Hayling-paper`)
**Línea de análisis:** descomposición DDS del N450 → residuo → AIS / Transfer Entropy → «dinámica
fronto-cingulada de la información» en la tarea Hayling.
**Estado:** **ARCHIVADA**. Cerrada por falta de validez incremental, no por falta de trabajo.
**Fecha de cierre:** 2026-08-18 · **n = 47** sujetos · permutaciones 10 000 · bootstrap 10 000.
**El manuscrito NO se ha editado. No se ha borrado ningún dato ni derivado. No se ha hecho push.**

Documento autocontenido: puede leerse sin abrir nada más del repositorio. Las fuentes primarias
están en `outputs/forensic/` (ver §8).

---

## 1. Resumen en cinco líneas

1. El pipeline contenía un **error de alineación temporal** que hacía que el «residuo DDS» no fuera
   el residuo del ajuste. Se reparó y se verificó (§2).
2. Reparado, **solo un efecto sobrevive**: AIS frontal INIT > INHIB (q = 0.031). Los otros dos
   resultados primarios del manuscrito desaparecen y aparece un TE en **dirección opuesta** (§3).
3. Ese efecto superviviente **también está en el EEG crudo, sin DDS, con el mismo tamaño**
   → la residualización DDS **no aporta nada** (§4).
4. Y **no es un efecto de información**: surrogados IAAFT que solo preservan el espectro lo
   reproducen; el excedente no lineal es el 0.5 % del AIS; el AUC validado es 0.508 = azar (§4).
5. Lo único con vida propia es un **exponente aperiódico 1/f** que no pertenece a este entorno (§6).

---

## 2. El defecto: alineación temporal del residuo (`03b`)

### 2.1 Qué hacía

`03_dds_peak_aligned_trial_by_tiral.py` ajusta el modelo DDS en la ventana **pico N450 ± 200 ms**
con el tiempo del modelo **relativo al borde izquierdo de la ventana**:

```python
t_win = times[mask] - tmin        # t = 0 en el borde izquierdo,  t ∈ [0, 0.4] s
```

`03b_build_residual_trial_roi_timeseries.py` reconstruía ese modelo **en tiempo absoluto de época**,
y sobre las 301 muestras de la época completa:

```python
t = t_ms / 1000.0                 # t ∈ [-0.2, +1.0] s   ← BASE TEMPORAL DISTINTA
residual = y - dds_model(t, ...)
```

El término `A·exp(−γ·t)` con γ hasta 150 s⁻¹ y t = −0.2 s diverge (`exp(30) ≈ 10¹³`). La columna
`dds_fit` no era el modelo ajustado y `residual` no era su residuo — ni dentro ni fuera de la ventana.

### 2.2 Evidencia (16 304 unidades trial × ROI, filtro R² > 0)

| Diagnóstico | Valor medido |
|---|---|
| R² del `dds_fit` publicado, **dentro de su propia ventana de ajuste** | **mediana −0.891** (IQR −1.740 … −0.371) |
| Unidades con R² < 0 dentro de su ventana | **92.2 %** |
| var(`residual`) / var(ERP) dentro de la ventana | mediana 1.38; **> 1 en el 76.5 %** |
| Magnitud del `residual` fuera de la ventana | hasta ~10³⁴ |

### 2.3 Prueba de que los parámetros DDS son correctos

Reevaluando el mismo modelo con la base temporal **correcta** se recupera el R² almacenado con
precisión numérica:

| Verificación | Objetivo | Medido |
|---|---|---|
| \|R²_recuperado − R²_almacenado\| (primera reconstrucción, tolerancia ±0.5 ms) | ~2·10⁻⁹ | mediana **2.1·10⁻⁹**, r = 0.9998 |
| \|R²_recuperado − R²_almacenado\| (reconstrucción exacta sobre la rejilla float de MNE) | ~2·10⁻⁹ | mediana **2.1·10⁻¹⁰**, máx **8.1·10⁻⁹**, **0 %** por encima de 10⁻⁶ |
| var(residuo) / var(ERP) tras la reparación | < 1 | mediana **0.402**, máx 0.9992, **< 1 en el 100 %** |

**Conclusión: los parámetros DDS eran correctos; el defecto estaba únicamente en el paso de resta.**
Detalle completo, incluida la discrepancia secundaria del offset `C`, en
`outputs/forensic/hayling_fix/CHANGELOG.md`.

---

## 3. Efectos: publicado → tras reparar → veredicto

**Familia FDR pre-declarada (sellada antes de mirar resultados, `outputs/forensic/hayling_incremental/PREREG.md`):**
**familia P = AIS × 4 ROIs + TE × 2 direcciones = 6 contrastes**, corrección Benjamini-Hochberg
dentro de P. Frontal se reporta dentro de esa familia, no en otra.

| Contraste | Publicado (`*_BUGGED`) | Tras reparar (`*_FIXED`) | q en familia P (6) | Veredicto |
|---|---|---|---|---|
| **AIS frontal** INIT > INHIB | t = 2.135, p = 0.0382, **p_fdr = 0.0763 n.s.** | t = 2.682, p = 0.0101, dz = 0.391 [0.114, 0.705] | **q = 0.031** ✔ | **Único superviviente** |
| **AIS cACC** INHIB > INIT | t = −2.582, p = 0.0131, **p_fdr = 0.0523 n.s.** | t = −1.792, p = 0.0798, dz = −0.261 | q = 0.148 | **Desaparece** |
| **TE frontal→cACC** INIT > INHIB | t = 2.461, p = 0.0177 (sin corregir en el paper) | t = 1.622, p = 0.1116, dz = 0.237 | q = 0.163 | **Desaparece** |
| **TE cACC→frontal** (no reportado en el manuscrito) | — | dz = **−0.392** [−0.674, −0.125], p_perm = 0.0090 | **q = 0.031** ✔ | **Efecto INVERTIDO respecto al publicado** |
| AIS LH / AIS RH | n.s. | dz = −0.214 / −0.088 | q = 0.188 / 0.559 | n.s. |

**Sensibilidad del superviviente (todos sobre el residuo corregido):**
umbral R² 0.0–0.5 → p = 0.0044–0.0108 en los cinco; sin filtro R² → p = 0.0089, q = 0.035;
bins 4/6/8/10 → q = 0.022 / 0.016 / 0.042 / 0.013.
**Ventana:** ±200 ms (la del ajuste) q = 0.040; ±100 ms q = 0.139; ±50 ms q = 0.916 — el efecto se
disuelve al acortar, pero eso es **el estimador saturando** con 25 muestras y 8 bins, no evidencia
en contra. El TE invertido cACC→frontal es estable en lags 1/2/4 (bins 8): q = 0.015 / 0.001 / 0.003.

---

## 4. Batería de validez incremental — las tres rutas y su cifra de cierre

Criterios fijados en `PREREG.md` **antes** de ejecutar nada.

### (c) ¿Aporta algo la residualización DDS? — **NO**

AIS sobre **EEG crudo single-trial**, misma ventana peak-aligned, mismos ROIs, mismo estimador
congelado (8 bins, lag 1). Familia C = 4 contrastes.

| ROI | AIS residuo (familia P) | AIS **crudo** (familia C) | Veredicto |
|---|---|---|---|
| **frontal** | dz = **0.391** [0.114, 0.705], q = 0.031 | dz = **0.365** [0.089, 0.674], **q = 0.039** | mismo signo, tamaños casi idénticos, IC95 solapados → **(c1) EL DDS ES INERTE** |
| **cACC** | dz = −0.261, q = 0.148 (n.s.) | dz = **−0.403** [−0.739, −0.115], **q = 0.037** | crudo significativo, residual no → **(c3) EL DDS DESTRUYE SEÑAL** |
| LH / RH | n.s. | q = 0.757 / 0.757 | n.s. |

Los dos veredictos posibles se dan **a la vez**: donde hay efecto el DDS es prescindible; donde el
DDS actúa, elimina un efecto real.

### (a) Escalera de modelos — el efecto se absorbe por completo

Modelo mixto a nivel de trial (4 082 trials frontales, intercepto aleatorio por sujeto):

| Modelo | β(condición) | p |
|---|---|---|
| `AIS_resid ~ cond` | +0.03324 | 1.0·10⁻⁵ |
| `+ amplitud N450` | +0.03695 | 1.1·10⁻⁶ |
| `+ pendiente 1/f` | +0.01126 (**−66 %**) | 0.043 |
| `+ AIS crudo` | +0.01253 (**−62 %**) | 0.057 |
| **`+ amplitud + 1/f + AIS crudo`** | **+0.00531 (−84 %)** | **0.31** |

*(A nivel de sujeto, con 47 diferencias, ninguna covariable disuelve a ninguna otra — pero todas se
refuerzan mutuamente, con r = 0.76–0.82 entre ellas: supresión por colinealidad, modelo no
identificable con n = 47. El nivel de trial, con dos órdenes de magnitud más de observaciones, sí
resuelve.)*

### (b) Parámetros DDS como variable dependiente — no robusto

Familia B = {\|A1\|, γ1, f1, \|A2\|, γ2, f2} × 4 ROIs = 24 contrastes.

| | \|A1\| frontal | \|A2\| frontal |
|---|---|---|
| Agregación **media** (pre-declarada) | dz = 0.355, **q = 0.028** | dz = 0.353, **q = 0.028** |
| Agregación **mediana** | dz = 0.336, q = **0.383** | dz = 0.183, q = **0.676** |
| **Media recortada 10 %** | dz = 0.147, q = **0.379** | dz = 0.147, q = **0.379** |

Causa: `fit_dds_blind` ejecuta `curve_fit` **sin cotas** y produce amplitudes de hasta **0.246 V**
(mediana ≈ 1·10⁻⁵ V). El «efecto» lo llevan unos pocos ajustes desbocados. Los otros 22 contrastes:
q ≥ 0.28.

### Discriminación validada — azar

| Métrica | Valor | Criterio pre-declarado |
|---|---|---|
| AUC agrupado (AIS estandarizado within-subject) | 0.5322, p_perm = 0.0004 | — |
| **AUC LOSO** | **0.5075, IC95 [0.474, 0.539]** | IC95 inferior > 0.55 → **NO SE CUMPLE** |
| Fracción de sujetos con INIT > INHIB | 0.660 (31/47), binomial p = 0.040 | — |

Estadísticamente detectable, prácticamente inútil.

### Naturaleza del efecto — enteramente lineal-espectral

Surrogados IAAFT (100, 20 iteraciones, seed 42) sobre el **EEG crudo**:

| ROI | AIS real | AIS en surrogados IAAFT | Contraste del excedente no lineal ΔA_NL |
|---|---|---|---|
| frontal | dz = 0.365 | **dz = 0.359** | dz = 0.024, **p = 0.870** |
| cACC | dz = −0.403 | **dz = −0.394** | dz = 0.026, **p = 0.858** |

Los surrogados lineales reproducen el **98 %** del tamaño del efecto. El excedente no lineal es el
**0.53 %** del AIS (conjunto de datos; 1.21 % en frontal). Es decir, **≈ 99.5 % del AIS es
estructura lineal-espectral**.

**EPR** (irreversibilidad temporal, sobre el residuo): ningún contraste sobrevive FDR
(mínimos p = 0.036 y 0.049 → **q = 0.099** en su familia de 4).
**Dip test de Hartigan**: **unimodalidad en las 24 pruebas** (p mínima = 0.425). Sin subgrupos latentes.

---

## 5. Rutas descartadas y por qué NO reabrir ninguna

Cada línea está **resuelta y en negativo**. Ninguna es trabajo pendiente.

| Ruta | Falsador que la cierra | Estado |
|---|---|---|
| **Reencaminar como marcador informacional** («information storage / transfer», reorganización fronto-cingulada) | **(c1)**: el AIS del EEG crudo da el mismo efecto frontal que el del residuo (dz 0.365 vs 0.391, IC95 solapados, ambos q < 0.05). El DDS no aporta nada sobre el EEG crudo. Y los surrogados IAAFT reproducen el 98 % del efecto: no es información, es espectro. | **RESUELTO — NEGATIVO** |
| **Reencaminar como marcador lineal honesto estilo MDD** | **AUC LOSO = 0.5075, IC95 [0.474, 0.539]** — indistinguible del azar, por debajo del umbral pre-declarado de 0.55. Un dz ≈ 0.4 nominal que no clasifica nada no es un marcador. | **RESUELTO — NEGATIVO** |
| **Reencuadrar a paper de parámetros DDS (dinámico)** | **(b) no robusto**: \|A1\|/\|A2\| frontales pasan de q = 0.028 (media) a q = 0.38 (mediana o media recortada). Lo produce `curve_fit` sin cotas, con amplitudes de hasta 0.246 V. Es un artefacto de ajuste, no un efecto. | **RESUELTO — NEGATIVO** |
| **Corregir y reintentar el AIS/TE residual** | **Ya está corregido.** El defecto se reparó y se verificó (§2.3: \|Δ\| ≤ 8.1·10⁻⁹, var ratio < 1 en el 100 %). Con el dato corregido, el único efecto que queda es lineal (§4) y prescindible (§4c). **Reparar más no puede cambiar el veredicto**: el comparador —el EEG crudo sin DDS— ya iguala o supera al residuo. | **RESUELTO — NEGATIVO** |
| **Rescatar el TE / la asimetría direccional** | Tras reparar, el TE frontal→cACC del manuscrito es n.s. (q = 0.163) y el único TE robusto va en **dirección contraria** (cACC→frontal, q = 0.031, estable en lags 1/2/4). Sus surrogados de cross-espectro también lo reproducen. No hay asimetría direccional que defender. | **RESUELTO — NEGATIVO** |
| **Ruta aperiódica (1/f)** | **No pertenece a este entorno.** Es un hallazgo post-hoc, ajeno al pre-registro, que exige un estudio nuevo, pre-registrado y replicado. Ver §6. | **FUERA DE ALCANCE — no es trabajo pendiente aquí** |

**Regla operativa para cualquier sesión futura:** si una idea sobre este entorno consiste en
recalcular AIS o TE sobre residuos DDS del N450 de Hayling —con otros bins, otro lag, otra ventana,
otro umbral de R², otro alineamiento, otro estimador o más sujetos de la misma muestra— **ya está
ejecutada y documentada aquí, y el resultado es negativo**. No la reabras.

---

## 6. Hecho relevante: el resultado ya era n.s. antes de la reparación

En la versión **publicada** (sobre el residuo defectuoso, hoy `*_BUGGED`), la propia salida del
repositorio marcaba ambos efectos AIS como no significativos tras FDR:

```
frontal   t = 2.1346   p_raw = 0.03815   p_fdr = 0.0763   n.s.   d = 0.311
cacc      t = -2.5821  p_raw = 0.01307   p_fdr = 0.0523   n.s.   d = -0.377
```

*(Cifras verificadas recomputando desde `derivatives/ais_n450_residual_r2pos_BUGGED/`. El fichero
`outputs/tables/stats_main_ais.csv` ha sido regenerado desde el dato corregido y ahora contiene los
valores post-reparación, no estos.)*

La premisa de que «el efecto está establecido y sobrevive FDR» nunca fue cierta en el propio
repositorio. El cierre no destruye un resultado sólido: formaliza uno que nunca lo fue.

---

## 7. Único residuo con vida propia — y por qué NO es trabajo pendiente aquí

El **exponente aperiódico 1/f** del segmento N450 difiere entre condiciones:

| Señal | ROI | χ INIT | χ INHIB | dz | p_perm | q (familia de 4) | GOF R² mediana |
|---|---|---|---|---|---|---|---|
| Residuo corregido | frontal | 1.062 | 0.911 | **0.405** | **0.0071** | 0.028 | 0.949 |
| **EEG crudo** | frontal | 1.848 | 1.697 | **0.399** | **0.0090** | 0.018 | 0.981 |
| **EEG crudo** | cACC | 1.516 | 1.600 | **−0.430** | **0.0045** | 0.018 | 0.981 |

> **Corrección de etiqueta importante.** La cifra «dz = 0.405, p = 0.008» que circuló en el
> encargo de cierre corresponde al **residuo**, no al EEG crudo. Al detectarlo se calculó también
> sobre el crudo con el mismo procedimiento: **dz = 0.399, p = 0.0090**. Las conclusiones no
> cambian —y de hecho el hallazgo se refuerza, porque aparece igual en el crudo y además en cACC—
> pero la cifra correcta para «χ crudo frontal» es **0.399 / p = 0.0090**, no 0.405 / 0.008.
> Artefactos: `outputs/forensic/hayling_incremental/data/chi_both.csv`, `chi_raw.py`.

**Por qué no es trabajo pendiente de este entorno:**

1. Es **post-hoc**. No figura en el pre-registro; salió como covariable de la escalera de modelos.
2. **No necesita el DDS ni la teoría de la información**: es un exponente espectral del EEG crudo.
   Vive fuera de esta línea, que es precisamente la que se archiva.
3. Es **pequeño y no validado**: dz ≈ 0.4 con n = 47, sin discriminación demostrada, en la misma
   muestra que produjo todos los falsos positivos anteriores.
4. Requiere **estudio nuevo, pre-registrado y replicado en muestra independiente**, con decisión
   del grupo sobre si merece la pena. Reabrir este entorno para perseguirlo repetiría exactamente
   el error que lo cerró: exploración post-hoc en la misma muestra.

**Acción correcta:** anotarlo como hipótesis para un estudio futuro. **Acción incorrecta:**
tratarlo como continuación de este análisis.

---

## 8. Dónde está la evidencia

```
DDS-Hayling-paper/
├── CLOSURE_REPORT.md          ← este documento (canónico)
├── STATUS.md                  ← señal machine-readable: ARCHIVED / RESUME: NO
├── README.md                  ← banner de retirada al principio
└── outputs/forensic/
    ├── hayling_battery/       batería de falsación sobre el residuo publicado
    │   ├── report.md              IAAFT, EPR, 1/f, dip, nulos de cross-espectro
    │   └── code/ data/
    ├── hayling_fix/
    │   └── CHANGELOG.md           defecto, evidencia, arreglo, verificación
    └── hayling_incremental/
        ├── PREREG.md              familias FDR y criterios, sellados antes de FASE 2
        ├── report.md              batería de validez incremental completa
        ├── code/                  build2 fase2 ladder sweeps famB_robust closing chi_raw core
        └── data/                  fase2_families, ladder, ladder_diffs, auc, trialvar, sweep_ais,
                                   sweep_te, famB_robust, closing_raw_iaaft, chi_both, chi_residuo,
                                   chi_EEG, fase2_chi, ladder_lmm_auc.txt + tablas regeneradas
```

**Nota sobre lo versionado.** Siguiendo la política declarada en `.gitignore` del repositorio
(«keep small summary CSVs, exclude trial-level files»), cuatro ficheros a nivel de trial
(`unit_battery.csv` 4.3 MB, `unit_te.csv` 926 KB, `keys.csv` 268 KB, `fix_verify.csv` 1.3 MB)
**permanecen en disco pero fuera de git**. Cada directorio `data/` incluye un `_LOCAL_ONLY.md` con
su contenido y las instrucciones exactas de regeneración (determinista: semillas fijas).
**Todas las cifras de este informe son verificables a partir de los agregados sí versionados.**

**Derivados** (`DDS-Hayling/derivatives/`, ninguno borrado):
`*_BUGGED` = producidos por el pipeline defectuoso, **no usar** (ver `_README_BUGGED.txt`);
`*_FIXED` = correctos, generados con la misma maquinaria congelada.

**Scripts modificados** (respaldos `*.pre_fix.bak` junto a cada uno):
`DDS-Hayling/02_dds_modeling/03b_build_residual_trial_roi_timeseries.py` (la reparación);
`DDS-Hayling/04_info_dynamics/04_extract_ais_n450_residual.py` y `05_extract_te_n450_residual.py`
(**solo** parametrización de E/S y bins/lag — la lógica de los estimadores está intacta);
`DDS-Hayling-paper/pipeline/08_stats_alt_alignment.py` y `config/analysis.yaml` (rutas `_FIXED`).

**Semillas:** IAAFT `seed = 42`; permutaciones y bootstrap `numpy.random.default_rng(2026)`.

---

## 9. Estado del control de versiones

**No se ha hecho push.** El repositorio queda listo para revisión humana:

- Rama `main`, último commit `942289b` («Remove outputs/ from tracking, ensure gitignored»).
- La mayor parte del árbol (`pipeline/`, `figures/`, `config/`, `outputs/`) está **sin trackear**.
- `.gitignore` excluye `outputs/figures/*` y `outputs/tables/*`, pero **no** `outputs/forensic/`:
  los informes forenses **sí son commiteables**. `CLOSURE_REPORT.md` y `STATUS.md` están en la raíz
  y tampoco están ignorados.
- Antes de cualquier commit, revisar que los `.bak` y los derivados grandes queden fuera.

---

## 10. Barrido de ficheros que pueden leerse como «trabajo pendiente»

Encontrados y señalizados (ninguno borrado):

| Fichero | Por qué se lee como pendiente | Acción tomada |
|---|---|---|
| `DDS-Hayling-paper/README.md` | describe el pipeline como activo y reproducible | banner de retirada insertado al principio |
| `DDS-Hayling/README_STAGEPLAN.md` | cuatro etapas marcadas «(to be created)» | banner de retirada insertado al principio |
| `DDS-Hayling/README.md` | además de describir el proyecto como activo, **describe un experimento Stroop, no Hayling** (copia-pega de otro repositorio) | banner insertado, con la incoherencia señalada |
| `DDS-Hayling-paper/config/analysis.yaml` | bloque `primary_tests` con los valores publicados | **se conserva a propósito**: hace que `05_stats_main.py` emita `WARN` al correr sobre datos corregidos. Ya comentado en el fichero |
| `DDS-HAYLING/Mejoras_por_hacer.docx` | «mejoras por hacer» (marzo 2026) | **no modificado** (binario). Obsoleto: precede al cierre |
| `DDS-HAYLING/Sugerencias_deepseek.docx` | sugerencias de análisis | **no modificado** (binario). Obsoleto |
| `DDS-HAYLING/DDS-Hayling_informe_retirada.docx` | creado por el usuario el 2026-08-18 11:27 | **no modificado** — documento del propio usuario |
| `DDS-Hayling-paper/prueba.py`, `outputs/tables/prueba.py` | scripts de prueba sueltos | **no modificados**; sin relevancia analítica |

---

## 11. Firma del cierre

Esta línea se cierra **con resultado**, no por abandono: el pipeline se reparó, se verificó, se
regeneró por completo y se sometió a una batería de falsación pre-registrada. El resultado es
negativo y está documentado hasta el nivel de la unidad trial × ROI.

**No reanudar. Si alguien cree tener una idea nueva sobre este entorno, léase §5 primero.**

# Hayling / N450 — reparación del residuo y batería de validez incremental

**Postura:** falsación-primero, pre-registro estricto. **El manuscrito no se ha tocado.**
Fecha: 2026-08-18 · n = 47 sujetos · permutaciones 10 000 · bootstrap 10 000.

---

# BLOQUE FIJO — FAMILIAS FDR PRE-DECLARADAS

> Declaradas y selladas en `PREREG.md` **antes** de ejecutar ningún análisis de FASE 2.
> **Frontal se reporta dentro de la familia P**, no en ninguna otra.

| Familia | Contenido | nº contrastes | Corrección |
|---|---|---|---|
| **P** — primaria (información) | AIS × 4 ROIs (frontal, cACC, LH, RH) + TE × 2 direcciones (frontal→cACC, cACC→frontal) | **6** | BH dentro de P |
| **C** — AIS sin DDS | AIS sobre EEG crudo × 4 ROIs | **4** | BH dentro de C |
| **B** — parámetros DDS | {\|A1\|, γ1, f1, \|A2\|, γ2, f2} × 4 ROIs | **24** | BH dentro de B |
| **S** — sensibilidad | cada barrido (bins, ventana, umbral R², lag) se corrige dentro de sí mismo | 4 por barrido | BH dentro del barrido |

**Fuera de toda familia FDR, marcados como exploratorios:** la escalera de modelos (FASE 2a: tests
anidados del *mismo* efecto), las fases DDS (circulares), el AUC validado y la variabilidad
trial-a-trial.

**Transparencia sobre P:** el paper corrige AIS sobre 4 ROIs y deja el TE sin corregir. Aquí se
declaró la familia más exigente (6 = todas las medidas de teoría de la información calculadas).
La columna «q (solo AIS, 4)» se incluye únicamente como nota informativa; **la decisión se toma
con q sobre 6**.

Contraste primario en todo el documento: **INIT vs INHIB, within-subject**, medias por sujeto,
dz pareado + IC95 bootstrap + p por permutación de signo (10 000).
Estimadores congelados sin modificar: **cuantiles, 8 bins, lag = 1 muestra (4 ms)**.

---

# FASE 1 — Reparación y regeneración

## 1.1 El arreglo

Detalle completo en `../hayling_fix/CHANGELOG.md`. Resumen: el ajuste DDS define el tiempo del
modelo **relativo al borde izquierdo de la ventana** (pico N450 ± 200 ms), mientras que
`03b_build_residual_trial_roi_timeseries.py` lo evaluaba en **tiempo absoluto de época** sobre las
301 muestras. El `dds_fit` guardado no era el modelo ajustado y el `residual` no era su residuo.

El arreglo evalúa el modelo en la base temporal del ajuste, restringido a su ventana, y reconstruye
la máscara sobre la **rejilla temporal float de MNE** (necesario: el ajuste guardó `n_samples = 100`
en el 12.2 % de los trials por redondeo; con la rejilla exacta la coincidencia es del 100 %).

**Verificación exigida:**

| Criterio | Objetivo | Resultado (16 304 unidades) |
|---|---|---|
| \|R²_recuperado − R²_almacenado\| | ~ 2·10⁻⁹ | mediana **2.1·10⁻¹⁰**, p95 9.3·10⁻¹⁰, máx 8.1·10⁻⁹; **0 %** por encima de 10⁻⁶ |
| var(residuo) / var(ERP) | < 1 | mediana **0.402**, p95 0.805, máx 0.9992; **< 1 en el 100 %** |

Derivados antiguos conservados como `*_BUGGED` (no borrados); sustitutos `*_FIXED`;
`config/analysis.yaml` reapuntado (respaldo `.pre_fix.bak`). Los valores publicados quedan
codificados en `primary_tests` a propósito, de modo que `05_stats_main.py` emite **WARN** al correr
sobre el dato corregido:

```
[WARN]  Frontal AIS:    t=2.682 (publicado 2.130)   p=0.0101 (publicado 0.0380)
[WARN]  cACC AIS:       t=-1.792 (publicado -2.580) p=0.0798 (publicado 0.0131)
[WARN]  TE front→cACC:  t=1.622 (publicado 2.460)   p=0.1116 (publicado 0.0177)
```

## 1.2 Tablas 2/3 y Figuras 3/4 regeneradas

`outputs/tables/stats_main_{ais,te}.csv`, `stats_group_condition_*.csv`,
`stats_sensitivity_r2.csv`, `stats_alt_alignment_*.csv`;
`outputs/figures/fig03_*`, `fig04_*` (regeneradas desde el dato corregido).

Efectos del paper, antes y después:

| Resultado primario del manuscrito | publicado (`_BUGGED`) | corregido (`_FIXED`) |
|---|---|---|
| AIS frontal INIT > INHIB | t = 2.13, p = 0.038, p_fdr = 0.076 **n.s.** | t = **2.68**, p = **0.010**, p_fdr = **0.041** ✔ |
| AIS cACC INHIB > INIT | t = −2.58, p = 0.013, p_fdr = 0.052 **n.s.** | t = −1.79, p = 0.080, p_fdr = 0.160 **n.s.** |
| TE frontal→cACC INIT > INHIB | t = 2.46, p = 0.018 | t = 1.62, p = **0.112** **n.s.** |
| (no reportado) TE cACC→frontal | — | dz = **−0.39**, p_perm = 0.009, **q = 0.031** ✔ |

Es decir: al reparar el defecto, **el único efecto AIS que pasa a ser significativo tras FDR es el
frontal**; el de cACC y el TE frontal→cACC dejan de serlo; y emerge un TE en **dirección opuesta**
que el manuscrito no reporta.

## 1.3 Sensibilidad (obligatoria, sobre el residuo corregido)

- **Umbral R²** (0.0 → 0.5): frontal estable, p = 0.0044–0.0108 en los cinco umbrales; cACC n.s.
  en los cinco; TE frontal→cACC n.s. en los cinco.
- **Sin filtro R²** (alineación alternativa del paper): frontal p = 0.0089, **q = 0.035**; cACC
  q = 0.071 n.s.; TE frontal→cACC p = 0.123 n.s.
- **Bins** (4, 6, 8, 10) y **ventana** y **lag**: tablas completas en §3.5. Adelanto: frontal
  estable en bins; **se disuelve al acortar la ventana**; TE cACC→frontal estable en lags 1/2/4.

---

# FASE 2 — Validez incremental

## 2.0 Resultado primario tras la reparación

#### Familia P (6 contrastes) — residuo corregido

| contraste | media INIT | media INHIB | dz [IC95] | p_perm | **q (familia 6)** | q (solo AIS, 4) |
|---|---|---|---|---|---|---|
| AIS frontal | 0.8922 | 0.8543 | 0.391 [0.114, 0.705] | 0.0102 | **0.0306** | 0.0408 |
| AIS cacc | 0.9287 | 0.9450 | -0.261 [-0.507, 0.006] | 0.0742 | **0.1484** | 0.1484 |
| AIS lh | 0.7645 | 0.7796 | -0.214 [-0.456, 0.066] | 0.1567 | **0.1880** | 0.2089 |
| AIS rh | 0.7530 | 0.7586 | -0.088 [-0.343, 0.223] | 0.5594 | **0.5594** | 0.5594 |
| TE frontal->cacc | 1.1812 | 1.1692 | 0.237 [-0.037, 0.479] | 0.1087 | **0.1630** | — |
| TE cacc->frontal | 1.2082 | 1.2348 | -0.392 [-0.674, -0.125] | 0.0090 | **0.0306** | — |

## 2.1 (c) AIS SIN DDS — ¿aporta algo la residualización?

#### Familia C (4 contrastes) — AIS sobre EEG CRUDO, misma ventana/ROIs/estimador

| contraste | media INIT | media INHIB | dz [IC95] | p_perm | **q (familia 4)** |
|---|---|---|---|---|---|
| AIS crudo frontal | 1.2703 | 1.2256 | 0.365 [0.089, 0.673] | 0.0196 | **0.0392** |
| AIS crudo cacc | 1.2862 | 1.3132 | -0.403 [-0.738, -0.114] | 0.0092 | **0.0368** |
| AIS crudo lh | 1.0947 | 1.0899 | 0.049 [-0.249, 0.316] | 0.7568 | **0.7568** |
| AIS crudo rh | 1.0642 | 1.0580 | 0.060 [-0.236, 0.345] | 0.6930 | **0.7568** |

## 2.2 (b) Parámetros DDS como variable dependiente

#### Familia B (24 contrastes) — parámetros DDS. Solo los q < 0.20:

| contraste | dz [IC95] | p_perm | q (familia 24) |
|---|---|---|---|
| absA1 frontal | 0.355 [0.243, 0.525] | 0.0021 | 0.0276 |
| absA2 frontal | 0.353 [0.236, 0.524] | 0.0023 | 0.0276 |

Los 20 contrastes restantes: q ≥ 0.28 (todos n.s.).

#### Familia B — robustez de la amplitud frontal a la regla de agregación

| agregación por sujeto | parámetro | dz | p_perm | q (dentro de su barrido) |
|---|---|---|---|---|
| media (pre-declarado) | absA1 | 0.355 | 0.0017 | 0.0264 |
| media (pre-declarado) | absA2 | 0.353 | 0.0022 | 0.0264 |
| mediana | absA1 | 0.336 | 0.0238 | 0.3825 |
| mediana | absA2 | 0.183 | 0.2259 | 0.6758 |
| media recortada 10% | absA1 | 0.147 | 0.0794 | 0.3788 |
| media recortada 10% | absA2 | 0.147 | 0.0533 | 0.3788 |
| media de log10|A| | logA1 | 0.370 | 0.0135 | 0.0720 |
| media de log10|A| | logA2 | 0.359 | 0.0180 | 0.0720 |

## 2.3 (a) Escalera de modelos — frontal

#### Escalera de modelos — frontal (ANCOVA within-subject sobre diferencias, n = 47)

| modelo | dz del efecto de condición [IC95] | p_perm |
|---|---|---|
| M0  AIS_resid ~ 1 | 0.391 [0.115, 0.711] | 0.0095 |
| M-i  + amplitud N450 | 0.407 [0.129, 0.725] | 0.0073 |
| M-ii + exponente 1/f (specparam) | 0.471 [0.153, 0.905] | 0.0026 |
| M-ii' + exponente 1/f (OLS trial) | 0.675 [0.185, 1.301] | 0.0001 |
| M-iii + AIS crudo | 0.600 [0.189, 1.055] | 0.0003 |
| M-full + amp + 1/f + AIS crudo | 0.638 [0.220, 1.298] | 0.0001 |
| M-R2 + R2 del ajuste | 0.430 [0.146, 0.784] | 0.0053 |

**Tests recíprocos** (cada candidata como variable dependiente):

| modelo | dz [IC95] | p_perm |
|---|---|---|
| REC amplitud N450 ~ 1 | 0.358 [0.054, 0.702] | 0.0178 |
| REC amplitud N450 ~ 1 + AIS_resid | 0.373 [0.070, 0.757] | 0.0155 |
| REC exponente 1/f ~ 1 | 0.405 [0.127, 0.706] | 0.0079 |
| REC exponente 1/f ~ 1 + AIS_resid | 0.488 [0.158, 0.963] | 0.0017 |
| REC AIS crudo ~ 1 | 0.365 [0.086, 0.673] | 0.0200 |
| REC AIS crudo ~ 1 + AIS_resid | 0.560 [0.121, 1.049] | 0.0006 |
| REC amplitud en el pico ~ 1 | 0.329 [0.043, 0.677] | 0.0292 |
| REC amplitud en el pico ~ 1 + AIS_resid | 0.345 [0.050, 0.763] | 0.0238 |
| REC varianza del residuo ~ 1 | 0.634 [0.436, 0.874] | 0.0001 |
| REC varianza del residuo ~ 1 + AIS_resid | 0.791 [0.538, 1.169] | 0.0001 |

## 2.4 Cierre: ¿es lineal-espectral el efecto del EEG crudo?

*Nota:* las p de esta tabla proceden de un sorteo de permutación distinto al de la familia C
(mismo generador, distinto punto de consumo), de ahí que el AIS crudo frontal aparezca como
p = 0.0155 aquí y p = 0.0196 en §2.1. La diferencia es ruido de Monte Carlo con 10 000
permutaciones y no afecta a ninguna decisión. La cifra autoritativa es la de la familia declarada.

#### Cierre — surrogados IAAFT sobre el EEG CRUDO

| ROI | cantidad | media INIT | media INHIB | dz | p_perm |
|---|---|---|---|---|---|
| frontal | AIS crudo (real) | 1.2703 | 1.2256 | 0.365 | 0.0155 |
| frontal | AIS crudo en surrogados IAAFT | 1.2535 | 1.2097 | 0.359 | 0.0172 |
| frontal | dA_NL crudo (excedente no lineal) | 0.0168 | 0.0160 | 0.024 | 0.8704 |
| cacc | AIS crudo (real) | 1.2862 | 1.3132 | -0.403 | 0.0092 |
| cacc | AIS crudo en surrogados IAAFT | 1.2863 | 1.3138 | -0.394 | 0.0114 |
| cacc | dA_NL crudo (excedente no lineal) | -0.0000 | -0.0006 | 0.026 | 0.8582 |
| lh | AIS crudo (real) | 1.0947 | 1.0899 | 0.049 | 0.7529 |
| lh | AIS crudo en surrogados IAAFT | 1.0861 | 1.0859 | 0.003 | 0.9862 |
| lh | dA_NL crudo (excedente no lineal) | 0.0086 | 0.0040 | 0.088 | 0.7737 |
| rh | AIS crudo (real) | 1.0642 | 1.0580 | 0.060 | 0.6866 |
| rh | AIS crudo en surrogados IAAFT | 1.0572 | 1.0542 | 0.027 | 0.8591 |
| rh | dA_NL crudo (excedente no lineal) | 0.0070 | 0.0038 | 0.126 | 0.3980 |

## 2.5 Barridos de sensibilidad completos

#### Barridos de sensibilidad — AIS (residuo corregido)

| barrido | valor | ROI | dz | p_perm | q (dentro del barrido, 4 ROIs) |
|---|---|---|---|---|---|
| bins | 4 | frontal | 0.428 | 0.0056 | 0.0224 |
| bins | 4 | cacc | -0.169 | 0.2578 | 0.3437 |
| bins | 6 | frontal | 0.431 | 0.0041 | 0.0164 |
| bins | 6 | cacc | -0.252 | 0.0920 | 0.1840 |
| bins | 8 | frontal | 0.391 | 0.0106 | 0.0424 |
| bins | 8 | cacc | -0.261 | 0.0815 | 0.1630 |
| bins | 10 | frontal | 0.432 | 0.0033 | 0.0132 |
| bins | 10 | cacc | -0.227 | 0.1231 | 0.1796 |
| window | +/-200ms (ajuste) | frontal | 0.391 | 0.0100 | 0.0400 |
| window | +/-200ms (ajuste) | cacc | -0.261 | 0.0774 | 0.1548 |
| window | +/-100ms | frontal | 0.266 | 0.0694 | 0.1388 |
| window | +/-100ms | cacc | -0.225 | 0.1267 | 0.1689 |
| window | +/-50ms | frontal | 0.113 | 0.4389 | 0.9155 |
| window | +/-50ms | cacc | -0.078 | 0.5978 | 0.9155 |

#### Barridos de sensibilidad — TE (residuo corregido)

| lag | bins | dirección | dz | p_perm | q (dentro del barrido, 2 direcciones) |
|---|---|---|---|---|---|
| 1 | 4 | frontal->cacc | 0.167 | 0.2757 | 0.2757 |
| 1 | 4 | cacc->frontal | -0.310 | 0.0372 | 0.0744 |
| 1 | 8 | frontal->cacc | 0.237 | 0.1146 | 0.1146 |
| 1 | 8 | cacc->frontal | -0.392 | 0.0077 | 0.0154 |
| 2 | 4 | frontal->cacc | 0.403 | 0.0072 | 0.0144 |
| 2 | 4 | cacc->frontal | -0.184 | 0.2065 | 0.2065 |
| 2 | 8 | frontal->cacc | 0.239 | 0.1058 | 0.1058 |
| 2 | 8 | cacc->frontal | -0.525 | 0.0005 | 0.0010 |
| 4 | 4 | frontal->cacc | 0.101 | 0.5008 | 0.9460 |
| 4 | 4 | cacc->frontal | -0.011 | 0.9460 | 0.9460 |
| 4 | 8 | frontal->cacc | 0.163 | 0.2634 | 0.2634 |
| 4 | 8 | cacc->frontal | -0.481 | 0.0016 | 0.0032 |

## 2.6 Extras exploratorios

#### Variabilidad trial-a-trial (EXPLORATORIO, sin familia FDR pre-declarada)

| ROI | medida | n | dz | p_perm | q (BH exploratorio, 16) |
|---|---|---|---|---|---|
| frontal | SD entre trials de ais_resid | 43 | 0.530 | 0.0011 | 0.0176 |
| frontal | SD entre trials de absA2 | 43 | 0.424 | 0.0026 | 0.0187 |
| frontal | SD entre trials de absA1 | 43 | 0.426 | 0.0035 | 0.0187 |
| frontal | SD entre trials de gamma1 | 43 | 0.233 | 0.1347 | 0.5387 |
| cacc | SD entre trials de f1 | 43 | 0.192 | 0.2182 | 0.6525 |
| cacc | SD entre trials de resid_var | 43 | 0.182 | 0.2447 | 0.6525 |
| frontal | SD entre trials de f2 | 43 | -0.133 | 0.4021 | 0.8051 |
| cacc | SD entre trials de absA2 | 43 | -0.136 | 0.4412 | 0.8051 |
| cacc | SD entre trials de absA1 | 43 | -0.133 | 0.4529 | 0.8051 |
| frontal | SD entre trials de f1 | 43 | 0.092 | 0.5501 | 0.8620 |
| frontal | SD entre trials de resid_var | 43 | 0.080 | 0.6119 | 0.8620 |
| cacc | SD entre trials de gamma1 | 43 | 0.064 | 0.6926 | 0.8620 |
| cacc | SD entre trials de ais_resid | 43 | -0.058 | 0.7015 | 0.8620 |
| cacc | SD entre trials de gamma2 | 43 | -0.031 | 0.8431 | 0.8620 |
| frontal | SD entre trials de gamma2 | 43 | 0.029 | 0.8545 | 0.8620 |
| cacc | SD entre trials de f2 | 43 | 0.029 | 0.8620 | 0.8620 |

#### Discriminación validada (EXPLORATORIO) — AIS frontal, residuo corregido

| Métrica | Valor | Criterio pre-declarado |
|---|---|---|
| AUC agrupado, AIS estandarizado within-subject (nivel de trial) | **0.5322**, p_perm = 0.0004 (nulo 0.5003) | — |
| **AUC LOSO** (logística ajustada en 46 sujetos, predicha en el excluido) | **0.5075**, **IC95 [0.474, 0.539]** | **IC95 inferior > 0.55 → NO SE CUMPLE** |
| Nivel de sujeto: fracción con INIT > INHIB | 0.660 (31/47), binomial p = 0.040 | — |

#### Modelo mixto a nivel de trial (robustez de la escalera; 4 082 trials frontales, intercepto aleatorio por sujeto)

| modelo | β(condición) | p |
|---|---|---|
| `AIS_resid ~ cond` | **+0.0332** | 1.0·10⁻⁵ |
| `+ amplitud N450` | +0.0370 | 1.1·10⁻⁶ |
| `+ pendiente 1/f (por trial)` | **+0.0113** (−66 %) | 0.043 |
| `+ AIS crudo` | **+0.0125** (−62 %) | 0.057 |
| `+ amplitud + 1/f + AIS crudo` | **+0.0053** (−84 %) | **0.31** |

#### Diagnóstico de colinealidad (correlaciones entre las diferencias INIT−INHIB, frontal, n = 47)

|  | AIS_resid | AIS crudo | amplitud | 1/f (trial) | var. residuo | 1/f (specparam) |
|---|---|---|---|---|---|---|
| **AIS_resid** | 1.000 | 0.764 | 0.313 | **0.820** | 0.610 | 0.570 |
| **AIS crudo** | 0.764 | 1.000 | 0.384 | 0.722 | 0.472 | 0.452 |
| **1/f (trial)** | 0.820 | 0.722 | 0.456 | 1.000 | 0.489 | 0.745 |

Calidad del ajuste aperiódico: specparam (`aperiodic_mode='fixed'`, 2–40 Hz — el rango se elevó de
1 a 2 Hz porque la ventana de 0.4 s da una resolución de 2.475 Hz), **GOF R² mediana 0.949**
(p05 0.663) sobre los 376 espectros sujeto × condición × ROI.

---

# 3. TABLA DE DECISIÓN

| Pregunta | Resultado | q en la familia declarada | Veredicto |
|---|---|---|---|
| ¿Sobrevive el AIS frontal a reparar el residuo? | dz = 0.391 [0.114, 0.705], p_perm = 0.010 | **q = 0.031** (familia P, 6) | **Sí** — es el único efecto AIS que pasa FDR |
| ¿Sobrevive el AIS cACC? | dz = −0.261, p_perm = 0.074 | q = 0.148 (P) | **No** |
| ¿Sobrevive el TE frontal→cACC? | dz = 0.237, p_perm = 0.109 | q = 0.163 (P) | **No** |
| ¿Hay algún TE robusto? | cACC→frontal dz = −0.392, estable en lags 1/2/4 (bins 8) | **q = 0.031** (P) | **Sí, en dirección OPUESTA a la publicada** |
| **(c)** ¿El AIS del EEG **crudo** muestra el efecto frontal? | dz = 0.365 [0.089, 0.674], p_perm = 0.020 | **q = 0.039** (familia C, 4) | **Sí — mismo signo, mismo tamaño, IC solapados** |
| **(c)** ¿Y en cACC? | crudo dz = −0.403, p = 0.009, **q = 0.037**; residual q = 0.148 | C / P | **Sí en crudo, no en residual** |
| **(c) VEREDICTO frontal** | crudo ≈ residual, ambos significativos, IC solapados | — | **(c1) EL DDS ES INERTE** |
| **(c) VEREDICTO cACC** | crudo significativo, residual no | — | **(c3) EL DDS DESTRUYE SEÑAL** |
| **(a-i)** ¿Se reduce a la amplitud N450? | subject-level dz 0.391 → 0.407 (p 0.010 → 0.007); trial-level β +0.033 → +0.037 | exploratorio | **No** |
| **(a-ii)** ¿Se reduce a la pendiente 1/f? | subject-level dz → 0.471 (p 0.003); **trial-level β −66 %, p → 0.043** | exploratorio | **Parcialmente (solo a nivel de trial)** |
| **(a-iii)** ¿Aporta sobre el AIS crudo? | subject-level dz → 0.600; **trial-level β −62 %, p → 0.057** | exploratorio | **Contradictorio entre niveles** |
| **(a)** Modelo completo | trial-level β **−84 %, p = 0.31** | exploratorio | **El efecto se absorbe por completo** |
| **(a)** Tests recíprocos | **cada** covariable retiene su efecto controlando por AIS_resid, y todas se refuerzan mutuamente (r = 0.76–0.82) | exploratorio | **Ninguna absorbe a la otra → supresión; no resoluble con n = 47** |
| **(b)** ¿Algún parámetro DDS separa condición? | \|A1\| y \|A2\| frontales q = 0.028 con la agregación pre-declarada (media) | **q = 0.028** (familia B, 24) | **Positivo por el test pre-declarado…** |
| **(b)** …¿es robusto? | mediana q = 0.38; media recortada 10 % q = 0.38; log₁₀\|A\| q = 0.072 | barrido S | **…pero NO robusto: lo llevan unos pocos ajustes desbocados (\|A\| hasta 0.246 V; `curve_fit` sin cotas)** |
| **(extra)** ¿Discrimina en la práctica? | AUC LOSO 0.5075, **IC95 [0.474, 0.539]** | criterio pre-declarado > 0.55 | **No — falla el criterio** |
| **(extra)** ¿Señal en la variabilidad trial-a-trial? | frontal: SD del AIS dz = 0.53 (q = 0.018), SD de \|A1\|/\|A2\| (q = 0.019); cACC nada | exploratorio | **Señal exploratoria, con el mismo problema de outliers en \|A\|** |
| **Cierre** ¿es lineal-espectral el efecto del EEG crudo? | frontal: real dz = 0.365 / **surrogados IAAFT dz = 0.359**; ΔA_NL contraste p = 0.87. cACC: real −0.403 / surrogados −0.394; p = 0.86. Excedente no lineal = **0.5 %** del AIS | — | **Enteramente lineal-espectral** |

---

# 4. MATRIZ RETIRAR-vs-REENCAMINAR

Regla de precedencia pre-declarada: **RETIRAR > REENCUADRAR > REENCAMINAR**.

| Ruta | Condición pre-declarada | ¿Se cumple? |
|---|---|---|
| **RETIRAR** | (c1 **o** c3) **Y** familia B nula **Y** frontal se reduce a amplitud o 1/f | **c1 ✔ y c3 ✔** (ambas). **B**: positivo por el test pre-declarado pero **no robusto** — no es «nulo» en sentido estricto. **a-ii**: se cumple a nivel de trial, no a nivel de sujeto. → **2 de 3 limpias, 1 con matiz** |
| **REENCAMINAR** (marcador lineal honesto) | **(c2)** **Y** varianza única sobre amplitud/1f/AIS-crudo **Y** AUC con IC95 inferior > 0.55 | **NO.** El veredicto es (c1), no (c2) — falla la primera condición. Y el AUC LOSO da IC95 [0.474, 0.539] — falla la tercera. **Descartada por dos vías independientes** |
| **REENCUADRAR — paper de pendiente aperiódica** | el efecto se explica por 1/f **y** χ muestra él mismo el contraste | **Parcialmente.** χ sí muestra el contraste (dz = 0.405, p_perm = 0.008) y absorbe el 66 % del efecto a nivel de trial; pero no lo disuelve a nivel de sujeto |
| **REENCUADRAR — paper dinámico (params DDS)** | familia B con q < 0.05 **robusto a los barridos** | **NO.** El único candidato colapsa al cambiar media→mediana o al recortar el 10 % |

## VEREDICTO

> ### RETIRAR el pipeline DDS-residual como contribución del paper.
>
> La razón decisiva no es la escalera de modelos (que con n = 47 no es identificable), sino el
> test más simple, mejor potenciado y pre-declarado de todos — la **familia C**:
>
> - En **frontal**, el AIS del **EEG crudo single-trial** da el **mismo efecto** que el AIS del
>   residuo DDS: dz 0.365 vs 0.391, mismo signo, IC95 solapados, ambos q < 0.05 en sus familias.
>   **La residualización DDS no añade nada.** → **(c1) DDS inerte.**
> - En **cACC**, el EEG crudo da un efecto **significativo** (dz = −0.403, q = 0.037) que
>   **desaparece** tras residualizar (q = 0.148). **La residualización destruye señal.** → **(c3).**
>
> El argumento central del manuscrito —que descomponer el N450 en sinusoides amortiguadas revela
> dinámica de la información que el ERP no muestra— queda falsado en sus dos direcciones a la vez:
> donde hay efecto, el DDS es prescindible; donde el DDS actúa, elimina un efecto real.
>
> A ello se suma que **lo que queda no es un hallazgo de información**: sobre el EEG crudo, los
> surrogados IAAFT —que preservan solo el espectro de potencia y la distribución de amplitudes—
> **reproducen el contraste completo** (frontal 0.359 vs 0.365; cACC −0.394 vs −0.403), y el
> excedente no lineal es el **0.5 %** del AIS. Y **no discrimina**: AUC LOSO 0.5075, IC95
> [0.474, 0.539], por debajo del umbral pre-declarado de 0.55.

> ### Ruta constructiva, si se quiere salvar algo
>
> Lo único defendible es un **efecto espectral de banda ancha en el EEG crudo peak-aligned**:
> el exponente aperiódico χ del segmento N450 difiere entre condiciones (dz = 0.405,
> p_perm = 0.008; specparam GOF R² = 0.95) y absorbe el 66 % del efecto AIS a nivel de trial.
> Eso es un paper **mucho** más pequeño, sin DDS, sin teoría de la información y sin
> «reorganización fronto-cingulada» — y con n = 47 y AUC ≈ 0.51 **requiere replicación previa**,
> no publicación directa.

> ### Lo que NO se sostiene, explícitamente
>
> 1. «Reorganización de la información» / «information storage & transfer»: falsado (IAAFT).
> 2. Asimetría direccional frontal→cACC: no sobrevive la reparación, y el TE robusto va en
>    **dirección contraria** (cACC→frontal, q = 0.031, estable en lags 1/2/4).
> 3. Efecto en cACC: existe en el EEG crudo, **no** en el residuo DDS.
> 4. Dinámica DDS (γ, f, A) como marcador: el único candidato no resiste cambiar la media por la
>    mediana.

---

# 5. Advertencias metodológicas que deben corregirse pase lo que pase

1. **`fit_dds_blind` corre `curve_fit` sin cotas** y produce amplitudes de hasta **0.246 V**
   (mediana 1.1·10⁻⁵ V). Cualquier estadístico sobre \|A\| basado en medias está dominado por
   esos ajustes. Debe acotarse el ajuste o pre-registrarse una regla de agregación robusta.
2. **Discrepancia código/datos**: el script de ajuste llama a `fit_dds_blind` (modelo con offset
   `C`), pero `C` no se guarda y `03b` reconstruye sin `C`. Los datos vigentes proceden de
   `fit_dds_phi` (la línea comentada). Debe sellarse.
3. **El AIS con 8 bins sobre ventanas cortas está cerca de la degeneración**: al pasar de ±200 ms
   (101 muestras) a ±50 ms (25 muestras) el AIS sube a ~1.66 bits y todos los contrastes se anulan.
   No es evidencia contra el efecto: es el estimador saturando. Cualquier análisis de ventana debe
   declararlo.
4. **Desbalance de trials** (mediana 64 INIT vs 36 INHIB) — controlado aquí por submuestreo en el
   análisis de variabilidad, pero debe declararse en el manuscrito.
5. **La familia FDR del paper deja el TE fuera.** Con la familia completa (6) los resultados no
   cambian de signo, pero la práctica debe corregirse.

---

# 6. Reproducibilidad

```
outputs/forensic/
├── hayling_fix/CHANGELOG.md          defecto, evidencia, arreglo, verificación
└── hayling_incremental/
    ├── PREREG.md                     familias y criterios, sellados antes de FASE 2
    ├── report.md                     este documento
    ├── code/   build2.py fase2.py ladder.py sweeps.py famB_robust.py closing.py core.py
    └── data/   fase2_families.csv ladder.csv ladder_diffs.csv auc.csv trialvar.csv
                sweep_ais.csv sweep_te.csv famB_robust.csv closing_raw_iaaft.csv
                fase2_chi.csv fix_verify.csv + tablas regeneradas del paper
```

Scripts del pipeline modificados (respaldos `*.pre_fix.bak`):
`DDS-Hayling/02_dds_modeling/03b_build_residual_trial_roi_timeseries.py` (arreglo),
`04_info_dynamics/04_extract_ais_n450_residual.py` y `05_extract_te_n450_residual.py`
(**solo** parametrización de E/S y bins/lag; lógica de los estimadores intacta),
`DDS-Hayling-paper/pipeline/08_stats_alt_alignment.py` y `config/analysis.yaml` (rutas `_FIXED`).

Semillas: IAAFT `seed = 42` (congelada de MDD); permutaciones/bootstrap `default_rng(2026)`.

**El manuscrito no se ha editado.**

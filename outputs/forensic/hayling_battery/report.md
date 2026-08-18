# Batería de falsación sobre AIS/TE residual del N450 (Hayling)

**Postura:** falsación-primero. El objetivo no es descubrir un efecto sino determinar si el
efecto AIS/TE residual INIT vs INHIB es una cantidad **lineal-espectral** (como el marcador de
MDD) o sobrevive a un nulo lineal. "Es lineal" es un desenlace válido y esperable.

Fecha: 2026-08-18 · Auditor: análisis forense automatizado · 47 sujetos, 20 cores, 10 000 permutaciones.

---

## 0. Inventario y parámetros congelados

### 0.1 Maquinaria congelada reutilizada (rutas exactas, sin modificar la lógica)

| Pieza | Ruta |
|---|---|
| AIS Hayling (Shannon discreto, qcut) | `DDS-HAYLING/DDS-Hayling/04_info_dynamics/04_extract_ais_n450_residual.py` |
| TE Hayling (discreta) | `DDS-HAYLING/DDS-Hayling/04_info_dynamics/05_extract_te_n450_residual.py` |
| Construcción del residuo | `DDS-HAYLING/DDS-Hayling/02_dds_modeling/03b_build_residual_trial_roi_timeseries.py` |
| Ajuste DDS peak-aligned | `DDS-HAYLING/DDS-Hayling/02_dds_modeling/03_dds_peak_aligned_trial_by_tiral.py` |
| IAAFT (N=100, 20 iter, seed 42) | `PortableSSD/DEP_DECO/reproduce/nonlinearity/iaaft_analysis.py` |
| AIS canónico MDD | `PortableSSD/DEP_DECO/scripts/utils/ais_estimator.py` |
| EPR (irreversibilidad) | `PortableSSD/DEP_DECO/scripts/utils/epr_estimator.py` |
| Pendiente aperiódica 1/f | `PortableSSD/DEP_DECO/reproduce/nonlinearity/aperiodic_slope.py` |
| Dip test (Hartigan) | `PortableSSD/DEP_DECO/reproduce/bimodality/12_test2_distribution.py` (`diptest`) |
| Estadística primaria del paper | `DDS-Hayling-paper/pipeline/05_stats_main.py`, `config/analysis.yaml` |

Copia congelada de las funciones efectivamente ejecutadas: `code/core.py`.

### 0.2 Series temporales del residuo

`DDS-Hayling/derivatives/trial_roi_timeseries_residual_r2pos/trial_roi_timeseries_residual_r2pos.csv`
(590 MB, 4 907 504 filas).

- **Sujetos:** 47 (P4, P5, P19 ya excluidos en config).
- **Unidades trial×ROI:** 16 304 (INIT 10 188 / INHIB 6 116; frontal 4 082, cACC 4 089, LH 4 060, RH 4 073).
- **Trials por condición y sujeto (ROI frontal):** mediana 64 INIT / 36 INHIB (IQR 25–78 y 26–41; mín. 2 y 1 — diseño desbalanceado por tarea, con sujetos de muy pocos trials supervivientes al filtro R²>0).
- **fs = 250 Hz** (Δt = 4 ms exacto).
- **Ventana de la serie:** −200 … +1000 ms → **301 muestras** por unidad (constante).
- **ROIs** (`config/analysis.yaml`): frontal = F3,F4,AF3,AF4,Fp1,Fp2,FC3,FC4 · **cACC = FC2, AFz, F2** · LH = F3,F5,FC3,FC5 · RH = F4,F6,FC4,FC6.

### 0.3 Parámetros congelados del pipeline **Hayling** (no importados de MDD)

Declarados y usados sin cambio en toda la batería:

| Parámetro | Valor Hayling | (valor MDD, NO usado) |
|---|---|---|
| Discretización | **cuantiles, 8 bins** (`pd.qcut`) | 4 bins percentil |
| Lag AIS/TE | **1 muestra = 4 ms** | 1 muestra @128–256 Hz |
| Ventana | **época completa −200…1000 ms, 301 pts** | 2 s de resting |
| Banda | **ninguna** (residuo del ERP filtrado 90 Hz en preproc.) | 1–40 Hz |
| MIN_SAMPLES | 20 | 10 ventanas |
| Ventana de ajuste DDS | pico N450 ± 200 ms → 101 pts | — |

Excepción declarada: el estimador **EPR** exige `n_bins²·5 + 2` muestras, de modo que con 8 bins
(322 > 301) devolvería NaN siempre. Se usa su valor congelado por defecto **n_bins = 4, lag = 1**.
Queda explícito que el EPR no comparte la discretización del AIS de Hayling.

### 0.4 Verificación de la maquinaria

El AIS recomputado con `core.py` reproduce **exactamente** el CSV publicado
(`max |Δ| = 4.4·10⁻¹⁶`, n = 16 304) y los tres tests primarios del paper:
frontal t(46) = 2.135, p = 0.0382 · cACC t(46) = −2.582, p = 0.0131 · TE frontal→cACC t(46) = 2.461, p = 0.0177.
La batería opera, por tanto, sobre la cantidad exacta que sostiene el manuscrito.

---

## 1. HALLAZGO BLOQUEANTE — el "residuo" publicado no es el residuo del ajuste DDS

Se trata de una inconsistencia de base temporal entre el ajuste y la resta, no de una cuestión
interpretativa.

- **El ajuste** (`03_dds_peak_aligned_trial_by_tiral.py`, línea `t_win = times[mask] - tmin`)
  define **t = 0 en el borde izquierdo de la ventana** pico±200 ms, y ajusta sobre 101 muestras.
- **La resta** (`03b_build_residual_trial_roi_timeseries.py`, línea `t = t_ms / 1000.0`)
  evalúa el mismo modelo en **tiempo absoluto de época**, t ∈ [−0.2, +1.0] s, sobre las 301 muestras.

Consecuencias medidas sobre las 16 304 unidades:

| Diagnóstico | Valor |
|---|---|
| R² del ajuste tal como quedó guardado, **dentro de su propia ventana** | mediana **−0.891** (IQR −1.740 … −0.371) |
| Unidades con R² < 0 dentro de su ventana | **92.2 %** |
| var(residuo publicado) / var(ERP), dentro de la ventana | mediana **1.38**; 76.5 % de unidades > 1 |
| Reconstrucción con la base temporal correcta vs. R² guardado | mediana `|Δ|` = **2.1·10⁻⁹**, r = 0.9998 |

La última fila es la prueba: al reevaluar el modelo con `t` relativo a la ventana se recupera el
R² almacenado con precisión numérica, lo que confirma que **los parámetros DDS son correctos y el
error está exclusivamente en el paso de resta**. El término `A·exp(−γ·t)` con γ hasta 150 s⁻¹ y
t = −0.2 s diverge (exp(30) ≈ 10¹³): fuera de la ventana el "residuo" está dominado por la
extrapolación del modelo, no por el EEG (valores de hasta 10³⁴).

Nota secundaria: el ajuste que se ejecuta es `fit_dds_blind` (modelo **con offset C**), pero `C`
no se guarda en `dds_n450_results.csv` y el constructor del residuo usa el modelo **sin C**. En la
práctica es inocuo — el CSV vigente reproduce exactamente el modelo sin offset — pero es una
divergencia código/datos que conviene sellar.

**Por eso la batería se corre sobre dos variantes:**

- **`pub`** — el residuo *tal como está publicado* (301 pts, base temporal absoluta). Es la cantidad que sostiene el manuscrito.
- **`cor`** — el residuo *corregido*: ERP − DDS evaluado con `t` relativo a la ventana, restringido a la ventana de ajuste (101 pts). Es el residuo que el paper dice usar.

---

## 2. Chequeo de confound de calidad de ajuste (R²), within-subject

R² del ajuste DDS por sujeto y condición (t pareado, n = 47):

| ROI | R² INIT | R² INHIB | t(46) | p | dz |
|---|---|---|---|---|---|
| frontal | 0.587 | 0.585 | 0.175 | 0.862 | 0.025 |
| cACC | 0.547 | 0.570 | −1.901 | 0.064 | −0.277 |
| LH | 0.521 | 0.511 | 0.705 | 0.484 | 0.103 |
| RH | 0.493 | 0.486 | 0.488 | 0.628 | 0.071 |

- **Frontal: no difiere** (como se esperaba en un diseño within-subject). El efecto AIS frontal no puede atribuirse a calidad de ajuste.
- **cACC: sí muestra una tendencia**, y **en la misma dirección que el efecto AIS de cACC** (INHIB > INIT en ambos). El contraste AIS de cACC está parcialmente confundido con la bondad del ajuste.

Se aplicó igualmente el control formal a todos los contrastes: ANCOVA within-subject sobre
diferencias, `ΔAIS ~ 1 + ΔR²`, con el intercepto como efecto de condición a diferencia nula de
ajuste, y permutación de signo (10 000) sobre el intercepto. Columnas `d tras control R²` de la tabla.

**Ringing / borde de ventana:** en la variante `pub` el residuo está dominado por extrapolación del
modelo (sección 1), no por ringing del filtro — el problema es más grave que el que se pedía
comprobar. En la variante `cor` el residuo se restringe a la ventana de ajuste, sin extrapolación;
la única contribución de borde es la propia truncación a 101 muestras, que afecta por igual a ambas
condiciones (diseño pareado) y por tanto no puede generar el contraste.

---

## 3. Tabla de resultados de la batería

Todos los contrastes son INIT − INHIB, t pareado sobre medias por sujeto (n = 47).
`p_perm` = permutación de signo within-subject, 10 000 permutaciones (test exacto del diseño pareado).
`q` = Benjamini-Hochberg sobre la familia de 6 contrastes de cada variante.
`ΔA_NL` = AIS real − media de 100 surrogados IAAFT (bits); IC95 bootstrap sobre sujetos (10 000).
Para TE, el nulo son **surrogados bivariados IAAFT con fase común (Prichard & Theiler 1994)**, que
preservan ambos auto-espectros **y el cross-espectro** — es decir, toda la predictibilidad cruzada
lineal / estructura de Granger — y destruyen solo el acoplamiento de orden superior.
La última columna es el test decisivo: **el mismo contraste calculado sobre los surrogados lineales**.

### variant pub
| contraste | d_bruto [IC95] | p_perm | q_FDR(6) | d tras control R² | q_FDR R² | ΔA_NL INIT [IC] | ΔA_NL INHIB [IC] | d(ΔA_NL) | q(ΔA_NL) | d en surrogados lineales (p) |
|---|---|---|---|---|---|---|---|---|---|---|
| AIS frontal | 0.311 [0.029, 0.607] | 0.0352 | 0.0704 | 0.485 | 0.0096 | 0.4377 [0.3673, 0.5115] | 0.3891 [0.3250, 0.4597] | 0.242 | 0.252 | -0.079 (0.592) |
| AIS cACC | -0.377 [-0.699, -0.100] | 0.0122 | 0.0549 | -0.428 | 0.0174 | 0.3141 [0.2768, 0.3503] | 0.3170 [0.2793, 0.3533] | -0.028 | 0.889 | -0.197 (0.182) |
| AIS LH | -0.064 [-0.448, 0.192] | 0.6898 | 0.6898 | -0.121 | 0.4209 | 0.2639 [0.2278, 0.3032] | 0.2606 [0.2195, 0.3056] | 0.022 | 0.889 | -0.107 (0.460) |
| AIS RH | -0.199 [-0.539, 0.083] | 0.1790 | 0.2148 | -0.274 | 0.0840 | 0.2042 [0.1706, 0.2382] | 0.2499 [0.2142, 0.2857] | -0.402 | 0.053 | 0.181 (0.233) |
| TE frontal->cACC | 0.359 [0.073, 0.685] | 0.0183 | 0.0549 | 0.375 | 0.0344 | -0.1249 [-0.1354, -0.1144] | -0.1203 [-0.1330, -0.1072] | -0.111 | 0.688 | 0.324 (0.031) |
| TE cACC->frontal | -0.211 [-0.523, 0.071] | 0.1576 | 0.2148 | -0.284 | 0.0838 | -0.1828 [-0.2101, -0.1579] | -0.1629 [-0.1879, -0.1409] | -0.227 | 0.252 | 0.162 (0.281) |

### variant cor
| contraste | d_bruto [IC95] | p_perm | q_FDR(6) | d tras control R² | q_FDR R² | ΔA_NL INIT [IC] | ΔA_NL INHIB [IC] | d(ΔA_NL) | q(ΔA_NL) | d en surrogados lineales (p) |
|---|---|---|---|---|---|---|---|---|---|---|
| AIS frontal | 0.397 [0.119, 0.713] | 0.0078 | 0.0297 | 0.436 | 0.0234 | 0.0218 [0.0143, 0.0296] | 0.0153 [0.0068, 0.0230] | 0.239 | 0.482 | 0.345 (0.023) |
| AIS cACC | -0.253 [-0.497, 0.023] | 0.0828 | 0.1656 | -0.261 | 0.1650 | 0.0045 [0.0002, 0.0089] | 0.0058 [0.0008, 0.0112] | -0.055 | 0.708 | -0.246 (0.101) |
| AIS LH | -0.212 [-0.449, 0.078] | 0.1618 | 0.1941 | -0.224 | 0.1826 | 0.0079 [0.0010, 0.0139] | 0.0107 [0.0050, 0.0164] | -0.109 | 0.567 | -0.199 (0.175) |
| AIS RH | -0.100 [-0.360, 0.213] | 0.5048 | 0.5048 | -0.101 | 0.5127 | 0.0082 [0.0029, 0.0136] | 0.0025 [-0.0045, 0.0089] | 0.210 | 0.482 | -0.197 (0.190) |
| TE frontal->cACC | 0.216 [-0.072, 0.470] | 0.1433 | 0.1941 | 0.214 | 0.1826 | -0.0214 [-0.0270, -0.0159] | -0.0180 [-0.0223, -0.0136] | -0.151 | 0.567 | 0.297 (0.039) |
| TE cACC->frontal | -0.387 [-0.684, -0.127] | 0.0099 | 0.0297 | -0.394 | 0.0234 | -0.0117 [-0.0170, -0.0065] | -0.0062 [-0.0166, 0.0064] | -0.123 | 0.567 | -0.283 (0.054) |

### 3.1 Tests auxiliares de la batería (EPR, 1/f, autocorrelación lag-1, dip)

Contrastes INIT − INHIB, dz pareado y p de permutación (10 000). `corr(AIS, r₁)` es la correlación
entre AIS y la autocorrelación a lag 1 sobre las 4 082–4 089 unidades trial×ROI.

| variante | ROI | dz EPR (p) | dz pendiente 1/f (p) | dz r₁ lag-1 (p) | corr(AIS, r₁) | dip p INIT | dip p INHIB |
|---|---|---|---|---|---|---|---|
| pub | frontal | −0.119 (0.421) | 0.253 (0.093) | −0.080 (0.613) | **0.096** | 0.963 | 0.740 |
| pub | cACC | −0.139 (0.357) | −0.113 (0.453) | 0.022 (0.882) | 0.096 | 0.967 | 0.825 |
| pub | LH | 0.081 (0.609) | 0.033 (0.835) | −0.215 (0.151) | 0.242 | 0.986 | 0.903 |
| pub | RH | 0.209 (0.165) | −0.114 (0.452) | 0.037 (0.807) | 0.297 | 0.804 | 0.951 |
| cor | frontal | −0.294 (0.049) | **0.342 (0.023)** | 0.277 (0.061) | **0.883** | 0.425 | 0.913 |
| cor | cACC | 0.307 (0.036) | −0.165 (0.306) | −0.282 (0.058) | 0.880 | 0.926 | 0.991 |
| cor | LH | 0.023 (0.876) | −0.067 (0.651) | −0.170 (0.261) | 0.867 | 0.972 | 0.791 |
| cor | RH | −0.104 (0.481) | 0.060 (0.689) | −0.114 (0.448) | 0.865 | 0.514 | 0.938 |

Ninguna p de EPR, de 1/f ni de r₁ sobrevive BH-FDR sobre su familia de 4 ROIs dentro de cada
variante: los mínimos son 1/f frontal `cor` q = 0.090, EPR `cor` q = 0.099 (frontal y cACC) y
r₁ `cor` q = 0.122. **Dip test: unimodalidad en los 24 tests** (todas las p ≥ 0.42) — no hay
bimodalidad ni subgrupos latentes en ninguna distribución.

---

## 4. Lectura

### 4.1 Un dato previo que hay que decir en voz alta

El brief asume que el efecto "sobrevive FDR". **La salida del propio repositorio dice lo contrario**
(`outputs/tables/stats_main_ais.csv`, generada por `05_stats_main.py`):

```
frontal  t=2.1346  p_raw=0.03815  p_fdr=0.0763   n.s.   d=0.311
cacc     t=-2.5821 p_raw=0.01307  p_fdr=0.05228  n.s.   d=-0.377
```

Los dos efectos AIS están marcados `n.s.` tras BH-FDR sobre las 4 ROIs. La batería reproduce esto
con permutación exacta (q = 0.070 y 0.049 sobre la familia de 4). No se re-litiga el efecto —
simplemente el punto de partida es más débil de lo que el brief supone, y eso condiciona el peso
que puede soportar cualquier veredicto posterior.

### 4.2 Qué sobrevive a corregir el residuo

| Efecto del paper | `pub` (publicado) | `cor` (residuo corregido) |
|---|---|---|
| AIS frontal INIT > INHIB | dz = 0.31, p = 0.035, q = 0.070 | **dz = 0.40, p = 0.008, q = 0.030 — sobrevive y se refuerza** |
| AIS cACC INHIB > INIT | dz = −0.38, p = 0.012, q = 0.055 | dz = −0.25, p = 0.083, q = 0.166 — **se disuelve** |
| TE frontal→cACC INIT > INHIB | dz = 0.36, p = 0.018, q = 0.055 | dz = 0.22, p = 0.143, q = 0.194 — **se disuelve** |
| (no reportado) TE cACC→frontal | dz = −0.21, p = 0.158 | dz = −0.39, p = 0.010, q = 0.030 — **aparece un efecto opuesto** |

El efecto **frontal** es el único que es robusto al hallazgo de la sección 1. Los otros dos
resultados primarios del manuscrito son, en grado sustancial, artefactos de la construcción
defectuosa del residuo — y la dirección dominante del TE se **invierte** al corregirlo.

### 4.3 ¿El efecto frontal es lineal-espectral?

Cuatro pruebas independientes, todas sobre la variante corregida `cor`:

1. **El nulo lineal reproduce el efecto.** Recalculado sobre los 100 surrogados IAAFT —que preservan
   solo el espectro de potencia y la distribución de amplitudes— el contraste frontal da
   **dz = 0.345, p_perm = 0.023**, frente a dz = 0.397 en los datos reales: el nulo lineal recupera
   el **87 %** del tamaño del efecto. Este es el test decisivo y su resultado es inequívoco.
2. **El excedente no lineal es trivial y no distingue condiciones.** ΔA_NL = 0.0218 [0.0143, 0.0296]
   bits en INIT y 0.0153 [0.0068, 0.0230] en INHIB — apenas el **2.3 %** del AIS frontal medio
   (0.874 bits); 1.2 % promediando las cuatro ROIs.
   Los IC no cruzan 0, pero el **contraste** entre condiciones sí: dz = 0.239, p_perm = 0.111,
   **q = 0.482**. La diferencia INIT/INHIB no reside en la parte no lineal.
3. **AIS a lag 1 ≈ autocorrelación a lag 1.** corr(AIS, r₁) = **0.883** sobre 4 082 unidades
   frontales, exactamente lo predicho por Wiener–Khinchin. La propia r₁ muestra el mismo contraste
   en la misma dirección (dz = 0.277, p = 0.061).
4. **Hay una firma espectral de banda ancha coherente.** La pendiente aperiódica 1/f del residuo
   frontal difiere entre condiciones en la misma dirección (dz = 0.342, p = 0.023; no sobrevive
   FDR sobre las 4 ROIs, q = 0.090). Un 1/f más pronunciado en INIT ⇒ más autocorrelación ⇒
   más AIS. Es el mecanismo lineal que explica el efecto.

Contraste con la variante `pub`: allí ΔA_NL era grande (0.44 y 0.39 bits, **25 %** del AIS frontal) y
corr(AIS, r₁) caía a 0.096. Eso **no** era evidencia de no linealidad neural: era la firma de la
exponencial divergente inyectada por el bug — una señal determinista, fuertemente no gaussiana,
que ningún surrogado IAAFT puede imitar. El "excedente no lineal" del residuo publicado es
excedente del artefacto.

**EPR:** ningún contraste sobrevive FDR en ninguna variante. No hay evidencia de que la
irreversibilidad temporal del residuo distinga INIT de INHIB.

**Dip test:** unimodalidad en los 24 tests. Ninguna distribución sugiere subgrupos.

---

## 5. VEREDICTO

> ### CAJA (B) — para el efecto **AIS frontal**: robusto pero lineal-espectral
>
> El contraste AIS frontal INIT > INHIB **sobrevive** la corrección del residuo (se refuerza:
> dz 0.31 → 0.40), el control por R² (dz 0.44, q = 0.023), y 10 000 permutaciones de signo con FDR
> (q = 0.030). Pero **no sobrevive el nulo lineal**: surrogados IAAFT que preservan únicamente el
> espectro de potencia reproducen el 87 % del efecto (dz = 0.345, p = 0.023); el excedente no
> lineal es el 1.2 % del AIS y su contraste entre condiciones es nulo (q = 0.482); el AIS correla
> r = 0.88 con la autocorrelación a lag 1.
>
> **Es un marcador lineal honesto.** La narrativa debe pasar de "reorganización de la información"
> a **"cambio en la predictibilidad temporal lineal del residuo del N450"**. No debe usarse
> "information storage" en sentido interpretativo. Desenlace calcado del de MDD.

> ### CAJA (A) — para **AIS cACC** y **TE frontal→cACC**: matizar la existencia
>
> Ninguno de los dos sobrevive corregir la construcción del residuo (q = 0.166 y 0.194). El de
> cACC además está confundido con la calidad del ajuste DDS, que difiere entre condiciones en la
> misma dirección (dz = −0.277, p = 0.064). Y al corregir el residuo la dirección dominante del TE
> se **invierte** (cACC→frontal, dz = −0.39, q = 0.030) — un efecto que el manuscrito no reporta y
> que la batería tampoco puede defender como no lineal (los surrogados de cross-espectro lo
> reproducen: dz = −0.283, p = 0.054).
>
> No hay base para afirmar una asimetría direccional frontal→cACC.

> ### CAJA (C) — no procede
>
> Ningún contraste sobrevive simultáneamente el control por R² **y** la batería de nulos lineales.
> No hay candidato genuinamente no lineal.

### 5.1 Acción prioritaria antes de cualquier reescritura

El bug de la sección 1 debe repararse en el pipeline (`03b_build_residual_trial_roi_timeseries.py`:
usar `t = (time_ms − window_tmin_ms)/1000` y restringir a la ventana de ajuste; y sellar la
discrepancia `fit_dds_blind`/offset `C`), y regenerarse `ais_n450_residual_r2pos` y
`te_n450_residual_r2pos`. Todos los números de las Tablas 2/3 y de las Figuras 3 y 4 cambian.
La batería ya adelanta cuál será el resultado: sobrevive frontal, desaparecen cACC y TE frontal→cACC.

---

## 6. Reproducibilidad

```
outputs/forensic/hayling_battery/
├── report.md
├── code/       core.py (estimadores congelados), build_arrays.py, step1.py,
│               check_residual.py, quant_bug.py, battery.py, infer.py
└── data/       battery_results.csv (tabla completa), unit_battery.csv (16 304 unidades ×
                AIS/EPR/r1/χ/ΔA_NL, ambas variantes), unit_te.csv (3 836 pares × TE + nulo
                bivariado), keys.csv
```

Semillas: IAAFT `seed=42` (congelada de MDD), permutaciones/bootstrap `default_rng(2026)`.
Coste: 16 304 × 100 surrogados IAAFT × 2 variantes + 3 836 × 100 surrogados bivariados × 2
variantes ≈ 3.9 M evaluaciones de AIS/TE, 20 min sobre 18 cores.

**No se ha editado el manuscrito.**

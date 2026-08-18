# PRE-REGISTRO — batería de validez incremental (Hayling / N450)

**Escrito y sellado ANTES de ejecutar cualquier análisis de FASE 2.**
Fecha de sellado: 2026-08-18. Ningún resultado de FASE 2 existía en disco en este momento.
Verificación: `outputs/forensic/hayling_incremental/data/` estaba vacío; solo existían los
artefactos de la batería de falsación previa (`hayling_battery/`).

---

## 0. Contexto congelado

- Contraste primario en TODO el documento: **INIT vs INHIB, within-subject**, medias por sujeto,
  n = 47. Estadístico: **dz pareado + IC95 bootstrap (10 000) + p por permutación de signo (10 000)**.
- Estimadores congelados, sin modificar: AIS/TE de Hayling (`04_extract_ais_n450_residual.py`,
  `05_extract_te_n450_residual.py`): **discretización por cuantiles, 8 bins, lag = 1 muestra (4 ms)**.
- Señal: residuo **corregido** (ERP − DDS en la base temporal del ajuste, ventana pico ± 200 ms,
  101 muestras, fs = 250 Hz). Filtro R² > 0 salvo donde se indique.
- ROIs (`config/analysis.yaml`): frontal = F3,F4,AF3,AF4,Fp1,Fp2,FC3,FC4 · **cACC = FC2, AFz, F2** ·
  LH = F3,F5,FC3,FC5 · RH = F4,F6,FC4,FC6.

---

## 1. FAMILIAS FDR DECLARADAS (fijadas antes de ver resultados)

Toda corrección es Benjamini–Hochberg **dentro** de la familia indicada. Frontal se reporta
**dentro de la familia P**, no en ninguna otra.

| Familia | Contenido | nº contrastes |
|---|---|---|
| **P** — primaria, información | AIS × 4 ROIs (frontal, cACC, LH, RH) **+** TE × 2 direcciones (frontal→cACC, cACC→frontal) | **6** |
| **C** — AIS sin DDS (FASE 2c) | AIS sobre EEG crudo × 4 ROIs | **4** |
| **B** — parámetros DDS (FASE 2b) | {\|A1\|, γ1, f1, \|A2\|, γ2, f2} × 4 ROIs | **24** |
| **S** — sensibilidad | cada barrido (bins, umbral R², lag, ventana) se corrige dentro de su propio barrido | variable, declarado in situ |

**Nota de transparencia sobre P:** el paper corrige AIS sobre 4 ROIs y deja el TE sin corregir.
Aquí se declara la familia más completa y más exigente (6 = todas las medidas de teoría de la
información efectivamente calculadas). Se reportará también q sobre la familia de 4 (solo AIS)
**únicamente como nota informativa**; la decisión se toma con q sobre 6.

**No pertenecen a ninguna familia FDR y se reportan como exploratorios, explícitamente marcados:**
la escalera de modelos (FASE 2a: tests anidados del *mismo* efecto, no contrastes independientes),
las fases DDS (φ1, φ2: circulares), el AUC validado, y el análisis de variabilidad trial-a-trial.

---

## 2. CRITERIOS DE DECISIÓN (fijados antes de ver resultados)

### 2c — AIS sin DDS (¿el DDS aporta algo?)
Se compara AIS sobre **EEG single-trial crudo** en la misma ventana peak-aligned, mismos ROIs,
mismo estimador, mismo lag/bins, contra AIS sobre **residuo corregido**.
Umbral de "significativo" = **q < 0.05 en su familia declarada** (C para crudo, P para residual).

- **(c1) DDS INERTE** — crudo ≈ residual y ambos significativos. Operacionalización: ambos
  q < 0.05, mismo signo, y IC95 de dz solapados.
- **(c2) LA RESIDUALIZACIÓN APORTA** — crudo no significativo (q ≥ 0.05), residual sí (q < 0.05).
- **(c3) EL DDS DESTRUYE SEÑAL** — crudo significativo (q < 0.05), residual no (q ≥ 0.05).

### 2a — Escalera de modelos (solo frontal)
ANCOVA within-subject sobre diferencias: `Δ(AIS_resid) ~ 1 + Δ(covariable)`, donde el
**intercepto** es el efecto de condición a diferencia nula de la covariable; p por permutación de
signo (10 000). Covariables, una por modelo:
  i. amplitud N450 peak-aligned (media del ERP crudo en pico ± 50 ms);
  ii. exponente aperiódico χ del residuo (specparam/FOOOF, 1–40 Hz, `aperiodic_mode='fixed'`);
  iii. AIS crudo (misma ventana, mismo estimador).
Robustez: modelo mixto a nivel de trial con intercepto aleatorio por sujeto.

- El efecto frontal **"se disuelve"** ante una covariable si su p_perm pasa de < 0.05 a **≥ 0.05**
  y |dz| cae **más de un 50 %**.
- **Explicación más parsimoniosa** = la cantidad que (a) muestra ella misma el contraste con
  p_perm < 0.05, y (b) absorbe el efecto del AIS-residual sin que el AIS-residual absorba el suyo
  (test recíproco). Si ninguna absorbe a la otra, se declara **no resoluble con n = 47**.

### 2b — Parámetros DDS como DV
Contraste pareado por parámetro y ROI, FDR dentro de la familia B (24). Positivo si **algún**
contraste da q < 0.05. Expectativa previa declarada: **nulo**.

### Extras (exploratorios)
- **AUC**: (i) AUC a nivel de trial del AIS frontal estandarizado within-subject, nulo por
  permutación de etiquetas dentro de sujeto (5 000); (ii) AUC LOSO (regresión logística ajustada
  en 46 sujetos, predicha en el excluido). Criterio de "discrimina en la práctica" declarado a
  priori: **límite inferior del IC95 > 0.55**.
- **Variabilidad trial-a-trial**: SD entre trials de los parámetros DDS y de la varianza del
  residuo, por sujeto × condición, submuestreando a n trials iguales por condición (20 sorteos)
  para neutralizar el desbalance INIT/INHIB. Exploratorio, sin FDR.

---

## 3. MATRIZ RETIRAR-vs-REENCAMINAR (fijada antes de ver resultados)

| Decisión | Condición |
|---|---|
| **RETIRAR** | (c1 **o** c3) **Y** familia B nula **Y** frontal se reduce a amplitud (a-i) o a 1/f (a-ii) |
| **REENCAMINAR** — marcador lineal honesto estilo MDD | (c2) **Y** frontal retiene varianza única sobre amplitud, 1/f **y** AIS-crudo **Y** AUC con IC95 inferior > 0.55 |
| **REENCUADRAR** — paper de pendiente aperiódica | el efecto frontal se explica por 1/f (a-ii) y χ muestra él mismo el contraste |
| **REENCUADRAR** — paper dinámico (parámetros DDS) | familia B da al menos un contraste con q < 0.05 robusto a los barridos de sensibilidad |

Si se cumple más de una condición, prevalece la más restrictiva (RETIRAR > REENCUADRAR > REENCAMINAR).
Si no se cumple ninguna limpiamente, el veredicto será **"no concluyente con n = 47"** y se
recomendará replicación, no publicación.

---

## 4. Lo que NO se hará

- No se buscará la familia FDR que favorezca a frontal: la familia P queda fijada arriba en 6.
- No se probarán ventanas, lags ni bins adicionales para rescatar un efecto: los barridos de
  sensibilidad son los que el paper ya tenía (bins, umbral R², lag, alineación), y se reportan
  completos, no selectivamente.
- No se editará el manuscrito.

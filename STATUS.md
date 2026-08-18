STATUS: ARCHIVED
RESUME: NO
DATE: 2026-08-18
REASON: DDS-residual/AIS/TE line withdrawn — no incremental validity (see CLOSURE_REPORT.md)

Esta línea de análisis (descomposición DDS del N450 → residuo → AIS / Transfer Entropy en la tarea
Hayling) está **retirada y archivada**. No es trabajo pendiente: es trabajo terminado con resultado
negativo. El pipeline se reparó, se verificó y se regeneró íntegramente antes de cerrarlo.

El motivo decisivo: el AIS del EEG crudo single-trial reproduce el efecto frontal con el mismo
tamaño que el AIS del residuo DDS (dz 0.365 vs 0.391, IC95 solapados), y en cACC el crudo es
significativo donde el residuo no lo es — la residualización DDS no aporta nada donde hay efecto y
elimina señal donde actúa. Además el efecto es lineal-espectral (surrogados IAAFT reproducen el
98 %) y no discrimina (AUC LOSO 0.5075, IC95 [0.474, 0.539]).

Antes de proponer cualquier reanálisis de este entorno, lee **CLOSURE_REPORT.md §5 (Rutas
descartadas y por qué NO reabrir ninguna)**: las cinco rutas plausibles están ya ejecutadas,
documentadas y cerradas en negativo.

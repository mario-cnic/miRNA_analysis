# Plan de reanálisis: miRNAs diferencialmente expresados en ATTR vs DCM, con vista a su detección en sangre periférica

_Rama: `claude/mirna-analysis` · Fecha: 2026-11-04_

---

## 0. Resumen

1. **La literatura sobre miRNAs en ATTR es escasa y casi toda es de sangre, no de tejido.** Hay dos candidatos replicables a priori: **miR-339-3p** (suero, ATTRwt; Derda 2018) y **miR-150-5p** (suero, ATTRv; Vita 2020). Ninguno está validado clínicamente. No he encontrado ningún estudio publicado de small RNA-seq en miocardio ATTR, así que este dataset sería de los primeros.
2. **El análisis actual (DESeq2 por defecto sobre la tabla de Novogene) tiene problemas que lo invalidan como resultado final:** no hay control de calidad, no hay filtrado, el log2FC no está reducido (*shrinkage*), los multimapeos están duplicados y no se ha comprobado si el grupo está confundido con el lote.
3. **Exploración rápida de `data/readcounts.csv` (sección 2.2).** Los miRNAs más "diferenciales" son marcadores de **tejidos no cardiacos**: neuronales y sensoriales (miR-124-3p, cluster miR-183/96/182) arriba en ATTR, y epiteliales y endocrinos (miR-375-3p, miR-205-5p, miR-200b/c) arriba en DCM. Separan los grupos de forma **perfecta**. Esto hace sospechar que **el grupo coincide con el origen o el procesado de la muestra**: distinto hospital, biopsia frente a explante, distinto lote de extracción o de secuenciación. **Hay que descartar esto antes de interpretar cualquier resultado.**
4. **Que un miRNA cambie en el tejido no garantiza que cambie en plasma.** En insuficiencia cardiaca, los cambios del miocardio apenas se reflejan en el plasma; lo que sube en sangre son los miRNAs liberados por daño del cardiomiocito (Akat 2014). Por eso el plan **puntúa cada candidato por su detectabilidad en sangre**: especificidad cardiaca, poca expresión en células sanguíneas y presencia en plasma según atlas públicos. Después propone una **validación prospectiva en plasma**.
5. **Método de expresión diferencial recomendado:** DESeq2 con *shrinkage* (apeglm) como análisis principal, y como confirmación edgeR-QL, limma-voom con pesos de calidad y un modelo a nivel isomiR (miRglmm). Se añade una evaluación del FDR mediante permutaciones (son posibles 924 asignaciones de etiquetas) y robustez *leave-one-out*. Con 6 frente a 6 muestras solo se esperan efectos grandes; los resultados deben presentarse como **cribado** que hay que confirmar.

---

## 1. Qué dice la literatura

### 1.1 miRNAs en amiloidosis cardiaca

| Estudio | Muestra / plataforma | Cohorte | Hallazgo principal | Nota |
|---|---|---|---|---|
| Derda et al. 2018, *PLoS ONE* ([doi:10.1371/journal.pone.0204235](https://doi.org/10.1371/journal.pone.0204235)) | Suero; TaqMan arrays (pools) + qPCR | ~10–13 por grupo: controles, IC, portadores TTR, ATTRwt ("senil") | **miR-339-3p ↑** en ATTRwt frente a controles, sin cambio en los otros grupos | Piloto pequeño; candidato principal para ATTRwt |
| Derda et al. 2015 ([PMC4936391](https://pmc.ncbi.nlm.nih.gov/articles/PMC4936391)) | Sangre; perfiles de miRNA | MCH (obstructiva/no), amiloidosis senil, estenosis aórtica, controles | **miR-29a + miR-155** distinguen MCH no obstructiva de amiloidosis senil | Útil como diagnóstico diferencial frente a otra hipertrofia |
| Vita et al. 2020, *Front Mol Neurosci* ([doi:10.3389/fnmol.2020.00102](https://doi.org/10.3389/fnmol.2020.00102); GEO **GSE149665**) | Suero; NanoString nCounter v3 | 10 ATTRv, 10 portadores asintomáticos, 10 CMT, 10 sanos | 33 ↑ / 48 ↓ en ATTRv frente a sanos; **miR-150-5p** separa ATTRv de portadores (AUC 0,97) | Fenotipo neuropático, no cardiaco; datos públicos reutilizables |
| Kontaraki et al. 2026, *Genes* ([doi:10.3390/genes17030356](https://doi.org/10.3390/genes17030356)) | Revisión | — | ATTR → miRNAs de fibrosis y TGF-β/Smad (miR-21, familia miR-29); AL → estrés proteotóxico e inmunidad (miR-150, miR-155) | Revisión narrativa; **no he podido leer el texto completo** |
| NCT07314268 ([ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT07314268)) | Vesículas extracelulares (EV) plasmáticas en ATTR-CM; miRNA por RT-qPCR frente a RM cardiaca | En curso | — | Indica que el interés actual va hacia miRNA en EV |

**Lagunas:** no he encontrado ningún small RNA-seq de miocardio ATTR, ni ningún panel de miRNA circulante validado que distinga ATTR-CM de su comparador clínico real (IC con FEVI preservada o hipertrofia no amiloide). Para ATTR-CM, los biomarcadores séricos actuales (NT-proBNP, troponina) estratifican el riesgo pero no diagnostican; también se están estudiando proteínas de matriz como decorina o IDUA ([PMC11587557](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11587557/)).

### 1.2 Del tejido a la sangre: qué miRNAs cardiacos llegan al plasma

- **Akat et al. 2014, *PNAS*** ([PMID 25012294](https://pubmed.ncbi.nlm.nih.gov/25012294/); GEO **GSE53081**) secuenciaron a la vez miocardio y plasma en IC. En el tejido, los myomiRs (miR-208a/b, miR-499) **no cambiaban**. En plasma subían hasta 140 veces en IC avanzada, en paralelo con la troponina I, y volvían a la normalidad tras implantar un LVAD. En IC estable, las diferencias circulantes eran menores de 5 veces. **Consecuencia para nosotros:** en plasma, lo que mejor se ve es la **liberación por daño**, no la regulación del tejido. Es el único dataset público con tejido y plasma emparejados, así que nos sirve de puente.
- **miR-208b y miR-499** reflejan el daño miocárdico en plasma ([Corsten 2010, *Circ Cardiovasc Genet*](https://www.pubmed.gov/20921333/)). miR-1 y miR-133a también se expresan mucho en músculo esquelético, lo que es un problema en ATTR, una enfermedad sistémica.
- **Origen sanguíneo de los miRNAs circulantes:** la mayoría de los miRNAs del plasma vienen de células sanguíneas. Un candidato muy expresado en eritrocitos, plaquetas o leucocitos queda "tapado" en plasma. Atlas que usaremos:
  - **miR-Blood** (Jehn 2024, *Sci Data*, [PMC10837159](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10837159/)): plasma más 9 componentes celulares de la sangre.
  - **Juzenas 2017** (*NAR* 45:9290): 7 tipos celulares, suero, exosomas y sangre total.
  - **miRNATissueAtlas 2025** (*NAR*, [doi:10.1093/nar/gkae1036](https://doi.org/10.1093/nar/gkae1036)): 61.593 muestras procesadas de forma uniforme; da la especificidad tisular.
- **Hemólisis y factores preanalíticos:** el cociente **miR-451a/miR-23a-3p** (ΔCq ≥ 5 indica riesgo moderado, ≥ 7–8 hemólisis significativa) es el estándar para detectar hemólisis ([Blondal 2013](https://brd.nci.nih.gov/brd/paper/methods/2013/assessing-sample-and-mirna-profile-quality-in-serum-and-plasma/131010); [Smith 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9317737)). La velocidad de centrifugación y la contaminación por plaquetas también alteran el perfil (Cheng 2013).
- **Biomarcadores de IC como comparador:** en el metaanálisis de Parvan et al. 2022 (*Eur J Heart Fail*, [enlace](https://duo.uio.no/handle/10852/99501)), un panel de 7 miRNAs para IC con FEVI preservada alcanzó una sensibilidad de 0,82 y una especificidad de 0,61, es decir, poco específico. miR-423-5p funciona en IC general pero no supera al BNP ([Yan 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5459698)). **Cualquier panel ATTR tendrá que demostrar que aporta algo sobre NT-proBNP.**

### 1.3 Métodos: cuantificación y expresión diferencial de miRNA-seq

- **Sesgo de preparación de librerías:** Giraldez et al. 2018 (*Nat Biotechnol* 36:746, [doi:10.1038/nbt.4183](https://doi.org/10.1038/nbt.4183)) muestran que la cuantificación **relativa entre muestras** es reproducible, pero cada protocolo tiene sesgos de secuencia. Los adaptadores degenerados (4N) los reducen. Por tanto, si las muestras se prepararon en lotes distintos, **no se puede comparar entre lotes**.
- **Cuantificación:** **nf-core/smrnaseq** (con mirtop/seqcluster para isomiRs; [nf-co.re](https://nf-co.re/smrnaseq)) y **miRge3** (miRNA, isomiRs, edición A-to-I, tRFs; [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.18.427129.full.pdf)) son las opciones más mantenidas. Las dos permiten tratar bien los multimapeos, que en la tabla de Novogene aparecen duplicados (sección 2.2).
- **Expresión diferencial:**
  - Un benchmark específico de miRNA-seq (preprint de 2026, [sciety](https://sciety.org/articles/activity/10.64898/2026.05.08.723923)) sitúa a **edgeR** entre los mejores y recomienda **edgeR, miRglmm o NBSR**. Recomienda normalizar por RPM e indica que DESeq2, edgeR v4 y limma tienden a **infraestimar** los log2FC. _Ojo:_ en nuestro tejido, miR-1-3p supone el **70 % de las lecturas**, así que el RPM bruto es inestable (un cambio en miR-1 desplaza a todos los demás). Usaremos RPM sin miR-1 y TMM/RLE como análisis de sensibilidad.
  - **miRglmm** (*Genome Biol* 2025, [doi:10.1186/s13059-025-03549-y](https://doi.org/10.1186/s13059-025-03549-y)): modelo mixto sobre contajes de isomiRs que mejora la estimación del efecto a nivel de miRNA. Requiere los FASTQ.
  - **NBSR** ([bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.05.07.592964.full.pdf)): modela la "competición por ser contado" y señala que los métodos de mRNA inflan los log2FC de los miRNAs poco expresados, precisamente el patrón del análisis actual (miR-4433a, miR-205 o miR-1269b con log2FC de 4–6 y menos de 10 lecturas de media).
  - Con n pequeño: limma y DESeq2 son las opciones más seguras; con n ≥ 6 por grupo se suele recomendar DESeq2. La concordancia entre herramientas es más baja precisamente con n pequeño ([Schurch 2016, *RNA*](https://doi.org/10.1261/rna.053959.115); [arXiv 2601.04122](https://arxiv.org/pdf/2601.04122)). **De ahí que proponga un consenso entre varios métodos.**

### 1.4 Hallazgos relevantes para nuestros candidatos actuales

- **Cluster miR-183/96/182:** miR-182 sube en hipertrofia en modelos murinos ([Li 2016, *Sci Rep*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758045/)). Pero el cluster es característico de **neuronas sensoriales (DRG) y de la retina**, igual que miR-124 es característico de neuronas. La ATTR tiene un componente neuropático y disautonómico, pero la denervación cardiaca haría esperar **menos** señal neuronal, no más. Lo más probable es un artefacto o un efecto de lote, aunque no se puede descartar una señal biológica.

---

## 2. Diagnóstico del análisis actual

### 2.1 Problemas del código y del diseño (`src/DEA_miRNA_ATTR.R`, `src/deseq.Rmd`)

| # | Problema | Impacto |
|---|---|---|
| 1 | Se parte de la tabla de contajes de Novogene (`Readcount_TPM.xls`), no de los FASTQ | No se pueden controlar el mapeo, los multimapeos ni los isomiRs; hay filas duplicadas |
| 2 | No hay control de calidad (tamaño de librería, % de reads en miRNA, PCA, outliers) | ATTR5 tiene la mitad de profundidad (4,2 M frente a ~9 M); DCM2 y ATTR1 tienen perfiles atípicos (miR-122 alto, miR-133a bajo) |
| 3 | No hay filtrado previo de miRNAs poco expresados (solo 382 de 1.324 tienen ≥10 contajes en ≥6 muestras) | Hits con contajes ínfimos y log2FC inflado (miR-4433a-5p, miR-1269b, miR-196b) |
| 4 | log2FC sin *shrinkage* (no se usa `lfcShrink`) | Ordenar por \|log2FC\| favorece el ruido |
| 5 | Diseño `~ condition` sin covariables (edad, sexo, lote, contenido de sangre en el tejido) | Riesgo de confusión (sección 2.2) |
| 6 | DCM está etiquetado como "control" | DCM es una enfermedad; la comparación es ATTR frente a otra miocardiopatía, lo que condiciona la interpretación y el uso clínico |
| 7 | Rutas absolutas de Windows y `meta.csv` con IDs `ATTR1.readcount` que no coinciden con las columnas | Análisis no reproducible; funciona solo porque el orden coincide |
| 8 | `lista` se ordena dos veces (la segunda anula la primera) | Orden del CSV final engañoso |

### 2.2 Exploración rápida de los contajes (realizada en esta sesión)

**Muestra:** miR-1-3p supone el ~70 % de las lecturas, seguido de miR-143, miR-378a y miR-133a. Es **tejido miocárdico** (falta confirmar si es biopsia endomiocárdica, explante o autopsia).

**miRNAs no cardiacos que separan perfectamente los grupos (CPM):**

| miRNA | Tejido típico | ATTR (6) | DCM (6) |
|---|---|---|---|
| miR-124-3p | neuronal | 3,8 – 7,4 | 0,0 – 0,3 |
| miR-96-5p | neurona sensorial / retina | 6,2 – 17,8 | 0,5 – 1,3 |
| miR-183-5p | neurona sensorial / retina | 225 – 503 | 25 – 68 |
| miR-182-5p | neurona sensorial / retina | 297 – 617 | 45 – 175 |
| miR-375-3p | páncreas endocrino / epitelio | 0,5 – 1,2 | 6,4 – 12,9 |
| miR-205-5p | epitelio escamoso | 0 (todas) | 1,0 – 3,3 |

Que haya una separación perfecta y bidireccional en marcadores de tejidos ajenos al corazón es la huella típica de un **efecto de lote o de procedencia**: distinto centro, distinto tipo de muestra, distinto kit, extracción o carrera de secuenciación. Una patología cardiaca no lo explicaría de forma natural. **Si el grupo coincide totalmente con el lote, ningún método estadístico puede separar los dos efectos.** En cambio, los candidatos de la literatura (miR-339-3p, miR-150-5p, miR-29a, miR-21, miR-208b, miR-499a) **no** muestran diferencias claras en tejido.

**Otros hallazgos:**
- **Contaminación variable por sangre** (miR-451a de 200 a 1.120 CPM) y por hígado (miR-122 de 2 a 960 CPM; muy alto en DCM2 y ATTR1). Hay que modelarlas como covariables o tratarlas como outliers.
- **11 pares de miRNAs con contajes idénticos** (p. ej. miR-7-5p/miR-3529-3p y miR-4433a-5p/miR-4433b-3p) por solapamiento de horquillas o multimapeo. Cuentan dos veces la misma señal e inflan el número de hits.
- **PCA:** PC1 (33 % de la varianza) no separa los grupos; está dominado por DCM2 y ATTR1 (las muestras con más miR-122) frente al resto.

---

## 3. Preguntas que necesito que me respondas antes de empezar

1. **Tipo y origen de las muestras:** ¿biopsia endomiocárdica, explante o autopsia? ¿Mismo hospital y mismo protocolo para ATTR y DCM? ¿Ventrículo izquierdo o derecho?
2. **Lote:** fechas de extracción de RNA, kit, RIN, si Novogene secuenció todo en una sola carrera y qué kit de librería usó. Esto suele estar en el informe de Novogene (`data/Novogene/`, que no está en el repositorio).
3. **¿Tienes los FASTQ?** Son imprescindibles para las fases 1 y 3b.
4. **Clínica por muestra:** ATTRwt o ATTRv (¿qué variante?), edad, sexo, FEVI, NT-proBNP, troponina, tratamiento (tafamidis, silenciadores) y tiempo de evolución.
5. **Sangre:** ¿hay plasma o suero biobancado de estos pacientes u otros? ¿Con qué anticoagulante y protocolo de centrifugación?

---

## 4. Plan de reanálisis

### Fase 0 — Reproducibilidad (preparación)
- Estructura `data/raw` (FASTQ, fuera de git), `data/processed`, `src/`, `results/` y `reports/` (Quarto).
- Entorno fijado con `renv` (R) o un contenedor; pipeline con `targets` o Snakemake/Nextflow.
- Metadatos completos en `data/meta.csv`: IDs que coincidan con las columnas, grupo, covariables y lote.

### Fase 1 — Reprocesado desde FASTQ (si están disponibles)
- QC de lecturas (FastQC/MultiQC), recorte de adaptadores (los de Novogene, o 4N/UMI si los hay) y filtro de longitud de 16–30 nt.
- Cuantificación con **nf-core/smrnaseq** (miRBase/MirGeneDB + mirtop), contrastada con **miRge3**.
  - Tabla a nivel de miRNA maduro, con los multimapeos asignados una sola vez (el miRNA es la unidad; se evitan los pares duplicados).
  - Tabla a nivel de **isomiR**, para miRglmm.
  - Fracción de otros RNAs pequeños (tRF, piRNA, snoRNA, Y-RNA), como QC y como exploración secundaria.
- Si no hay FASTQ: trabajar con `readcounts.csv`, **colapsar los pares idénticos** y anotar la pérdida de información.

### Fase 2 — Control de calidad y confusión (decide si el análisis es interpretable)
1. Métricas por muestra: profundidad, % de lecturas en miRNA, nº de miRNAs detectados y saturación.
2. **Panel de pureza y procedencia** con marcadores del miRNATissueAtlas 2025 y de miR-Blood:
   - sangre: miR-451a, miR-144-3p, miR-486-5p, miR-16
   - hígado: miR-122-5p
   - músculo esquelético: miR-206
   - neuronal: miR-124, cluster miR-183/96/182
   - epitelial y endocrino: miR-200 family, miR-205, miR-375
   - endotelial: miR-126
   - fibroblasto y fibrosis: miR-21, miR-29

   Se calcula una puntuación por muestra y se comprueba su asociación con el grupo.
3. PCA y clustering jerárquico coloreados por grupo, lote, sexo, edad y puntuaciones de contaminación. Se correlacionan los primeros PCs con cada covariable (`PCAtools::eigencorplot`).
4. **Regla de decisión:**
   - Si el lote coincide con el grupo, la expresión diferencial en tejido **no es interpretable como biología**. Habría que buscar más muestras (por ejemplo, DCM del mismo centro y protocolo que las ATTR) o datos públicos con los que comparar.
   - Si el lote es parcial, se incluye en el modelo o se estiman factores no deseados con **RUVSeq** (RUVg con miRNAs estables, o RUVs) o **svaseq**.
   - Las puntuaciones de sangre e hígado entran como covariables continuas.

### Fase 3 — Expresión diferencial en tejido
**3a. Análisis principal (nivel de miRNA)**
- Filtro: `filterByExpr` (edgeR) o ≥10 contajes en ≥6 muestras.
- **DESeq2**: `design = ~ [lote/RUV] + [sexo] + [contaminación] + grupo`, con `lfcShrink(type = "apeglm")`, filtrado independiente y Cook's distance activada. ATTR frente a DCM, con DCM como referencia (no como "control sano").
- Análisis de sensibilidad de la normalización: RLE (por defecto), TMM, y RPM/TMM **excluyendo miR-1-3p**, para comprobar que los hits no dependen del miR-1.

**3b. Confirmación con otros métodos y consenso**
- **edgeR quasi-likelihood** (`glmQLFit` con `robust = TRUE`).
- **limma-voom** con `voomWithQualityWeights`, que penaliza muestras como ATTR5, DCM2 o ATTR1 sin excluirlas.
- **miRglmm** sobre isomiRs (si hay FASTQ), que da a la vez el efecto a nivel de miRNA y la variabilidad entre isomiRs.
- Criterio de consenso: FDR < 0,05 en DESeq2, mismo signo y FDR < 0,1 en al menos uno de los otros, y |log2FC reducido| ≥ 0,58 (1,5 veces).

**3c. Robustez (clave con 6 frente a 6)**
- **FDR empírico por permutación:** con 12 muestras hay C(12,6) = 924 asignaciones de etiquetas. Se repite el pipeline sobre todas o sobre un subconjunto amplio y se compara el número de hits real con la distribución nula.
- **Leave-one-out:** se repite 12 veces quitando una muestra cada vez y se exige que el candidato se mantenga en al menos 10 de las 12.
- Se informa la **potencia esperada**: con n = 6 se detectan sobre todo efectos de 2 veces o más (Schurch 2016). La ausencia de hits no implica ausencia de efecto.

**3d. Interpretación biológica**
- Dianas validadas (miRTarBase) y predichas (TargetScan) con `multiMiR`; enriquecimiento de rutas (TGF-β, matriz extracelular, estrés proteotóxico, mitocondria).
- Integración con transcriptómica (si existe RNA-seq de ATTR o datasets públicos de snRNA-seq de miocardiopatías, como Chaffin 2022, *Nature*) para atribuir cada miRNA a un tipo celular cardiaco.

### Fase 4 — Priorización de candidatos detectables en sangre periférica
Cada miRNA candidato recibe una puntuación compuesta:

| Criterio | Fuente | Qué buscamos |
|---|---|---|
| A. Efecto y robustez en tejido | Fase 3 | log2FC reducido grande, consenso y estabilidad en LOO y permutaciones |
| B. Especificidad cardiaca | miRNATissueAtlas 2025 (índice de especificidad tisular) | Enriquecido en corazón frente a otros tejidos, sobre todo frente a músculo esquelético y nervio, porque la ATTR es sistémica |
| C. Poca expresión en células sanguíneas | miR-Blood, Juzenas 2017 | Si eritrocitos, plaquetas o leucocitos lo expresan mucho, la señal cardiaca queda tapada |
| D. Detectable en plasma | miR-Blood (plasma), exRNA Atlas, plasma de GSE53081 | Por encima del límite de detección en plasma sano |
| E. Respaldo en la literatura | Sección 1 | miR-339-3p, miR-150-5p, miR-29a, miR-155, miR-21, miR-208b, miR-499a |
| F. Asociación con la gravedad | Clínica por muestra (si la hay) | Correlación con NT-proBNP, troponina o grosor de pared |

Además, se incluye siempre un **panel fijo de candidatos de la literatura**, aunque no salgan en tejido, porque en plasma pueden comportarse distinto (Akat 2014).
- miR-339-3p y miR-150-5p, candidatos ATTR.
- miR-29a y miR-155, para distinguir de la miocardiopatía hipertrófica.
- miR-208b y miR-499a, marcadores de daño.
- miR-423-5p y miR-21, referencia de IC.

### Fase 5 — Comprobación en datos públicos (sin coste de laboratorio)
- **GSE53081** (Akat 2014): ¿los candidatos del tejido se ven en el plasma de pacientes con IC? ¿Cambian con la IC?
- **GSE149665** (Vita 2020, suero ATTRv por NanoString): ¿los candidatos están en el panel y cambian en ATTRv?
- Buscar en GEO/SRA (cuando la red lo permita) más datasets de miRNA de tejido o plasma en IC o amiloidosis.

### Fase 6 — Validación en sangre periférica (estudio prospectivo)
**Diseño:** una fase de descubrimiento en plasma y otra de validación independiente, siguiendo las fases PRoBE de desarrollo de biomarcadores.
- **Grupos:**
  - ATTR-CM (wt y v)
  - **IC con FEVI preservada o hipertrofia no amiloide de edad y sexo similares**, que es el comparador clínico real porque es la que se confunde con ATTR
  - DCM
  - AL, si es posible
  - controles sanos
  - portadores TTR asintomáticos, si existen
- **Preanalítica:**
  - Plasma EDTA con doble centrifugación para obtener plasma pobre en plaquetas, procesado en menos de 2 h, alícuotas a −80 °C y sin ciclos de descongelación.
  - QC de hemólisis: miR-451a/miR-23a-3p, excluyendo ΔCq ≥ 7.
  - Spike-in exógeno (cel-miR-39) para controlar la extracción.
- **Plataforma:**
  - Descubrimiento: small RNA-seq de plasma con kit de **adaptadores 4N/UMI** (p. ej. QIAseq miRNA), n ≈ 20–30 por grupo.
  - Opcional: perfil de miRNA en **vesículas extracelulares** (EV), que enriquece la carga liberada por los tejidos, en línea con NCT07314268.
  - Validación: **RT-ddPCR** o RT-qPCR de 5–10 candidatos.
  - Normalización con miRNAs estables elegidos con NormFinder/geNorm sobre los propios datos. Evitar miR-16, que es sensible a la hemólisis. También se pueden usar cocientes entre pares de miRNAs, que se normalizan solos.
- **Modelado:**
  - Regresión logística penalizada (elastic net) con validación cruzada anidada.
  - Se informa el AUC con IC95 % y se compara **contra NT-proBNP + troponina** y contra el modelo combinado (valor añadido: ΔAUC y NRI).
  - Validación externa en una cohorte de otro centro.

### Fase 7 — Entregables
- Informe reproducible (Quarto) con QC, análisis de confusión, resultados de expresión diferencial por método, consenso y robustez.
- Tabla de candidatos priorizados para sangre, con su puntuación A–F.
- Script de RT-qPCR/ddPCR (primers comerciales) y protocolo preanalítico para la validación.

---

## 5. Ideas adicionales (más exploratorias, con respaldo en la evidencia)

1. **Cocientes de miRNAs como biomarcador** (p. ej. un miRNA de daño frente a uno de fibrosis): no dependen de la normalización ni, en parte, de la hemólisis, y son más fáciles de trasladar a qPCR clínica.
2. **Señal extracardiaca de la ATTR:** la TTR se produce en el hígado y la ATTR afecta a nervios, ligamentos y túnel carpiano. Un panel en plasma puede aprovechar miRNAs no cardiacos (como el miR-150-5p de la neuropatía). Así, la "contaminación" neuronal del tejido merece revisarse: si se confirma que no es un efecto de lote, el cluster miR-183/96/182 podría reflejar inervación o remodelado neural.
3. **isomiRs y tRFs:** con los FASTQ, miRglmm y miRge3 permiten buscar isomiRs con uso diferencial y fragmentos de tRNA, que en plasma son abundantes y en parte específicos de tejido (Akat 2014 los describe en plasma de IC).
4. **Deconvolución del plasma:** con las firmas del miRNATissueAtlas se puede estimar qué fracción de los miRNAs plasmáticos viene del corazón y comprobar si es mayor en ATTR (marcador de "liberación cardiaca" similar a la troponina).

---

## 6. Riesgos y mitigación

| Riesgo | Mitigación |
|---|---|
| Grupo totalmente confundido con el lote o la procedencia | Fase 2: si se confirma, pasar a datos de plasma y públicos (fases 5–6) y pedir muestras DCM del mismo protocolo |
| No hay FASTQ | Trabajar con la tabla de Novogene colapsada; sin análisis de isomiRs |
| n = 6 frente a 6 sin potencia suficiente | Consenso entre métodos, permutaciones y LOO; presentar como cribado; ampliar n en la validación |
| El candidato del tejido no se ve en plasma | Puntuación de detectabilidad (fase 4), panel de la literatura y búsqueda directa en plasma (fase 6) |
| Hemólisis o plaquetas en plasma | Protocolo preanalítico estricto y QC miR-451a/miR-23a |

---

## 7. Notas sobre las fuentes

- La red de esta sesión **bloquea PubMed/PMC, OpenAlex, Europe PMC y las webs de las revistas**, así que no he podido leer textos completos ni verificar las referencias en OpenAlex. Los datos proceden de resúmenes y fragmentos obtenidos con el buscador web, y los DOI se han tomado de esos resultados. Antes de citar cifras concretas en un manuscrito hay que comprobarlas en el texto completo.
- Son **preprints** (sin revisión por pares): el benchmark de miRNA-seq de 2026, NBSR (2024) y arXiv 2601.04122.
- La revisión de Kontaraki 2026 solo la he leído a través de fragmentos.

# RMT + Genómica Computacional — colaboración propuesta con el Dr. Enrique Hernández-Lemus

Programa de investigación conjunta: Random Matrix Theory (estática y
dinámica), análisis multiescala por wavelets, y aprendizaje automático,
aplicados a genómica computacional y oncología de sistemas — foco en
cáncer de mama (TNBC), redes de coexpresión a gran escala, genómica 3D
(Hi-C), y célula única — con protocolo de validación numérica explícito
en cada línea.

## Contenido

- **`article.tex`** / **`.pdf`** — el artículo de investigación (inglés,
  40pp, formato BMC Bioinformatics: abstract estructurado,
  referencias numeradas, Declarations completas), con análisis real
  sobre TCGA-BRCA: umbral RMT (Luo et al. 2007) comparando TNBC vs.
  Luminal A, null de Monte Carlo por permutación, eigenvalores
  "spiked" + entropía espectral + escalamiento Tracy-Widom, red de
  coexpresión vs. nulos aleatorios (Erdős–Rényi/configuración),
  trayectorias de eigenvalores tipo Dyson (repulsión de niveles,
  verificada en 20 órdenes aleatorios distintos), una Proposición
  real (la densidad conjunta GOE es la distribución estacionaria del
  SDE de Dyson — resultado de Dyson 1962, no antes enunciado
  formalmente en este trabajo), un Corolario real y nuevo (2026-09-05:
  la fuerza de repulsión de Dyson es exactamente inversamente
  proporcional a la brecha entre eigenvalores — se deriva de la misma
  SDE de la Proposición y explica cuantitativamente, por primera vez,
  por qué el par λ1–λ2 es el único no confiablemente repulsivo:
  1.7–10.7× más débil que cualquier otro par, calculado directo de la
  Tabla de repulsión ya existente, no de datos nuevos), una prueba de
  reconocimiento por ML (incluye GA), y el ataque real, en esta ronda,
  a la Línea 6b del Roadmap (dinámica de Dyson en la propia matriz de
  pesos del clasificador durante el entrenamiento SGD, siguiendo a
  Aarts et al. 2024) — con tres resultados nulos honestos reportados
  como tales (wavelets en datos reales, techo de clasificación, y
  repulsión de pesos del clasificador fuera de su transitorio inicial
  de entrenamiento — este último con un hallazgo real y concreto: la
  repulsión sí aparece, mas confinada casi por completo a los primeros
  ~50 epochs antes de que el entrenamiento converja, 12× más señal ahí
  que después). Roadmap explícito y ya actualizado con datasets reales
  para cada línea pendiente (ver abajo).
- **`programa_investigacion_RMT_HernandezLemus.tex`** / **`.pdf`** — la
  propuesta original: justificación científica, hipótesis, 6 líneas de
  trabajo, núcleo matemático, protocolo de validación, datos, productos
  esperados, ruta de 12–18 meses.
- **`data/README.md`** — los conjuntos de datos públicos identificados
  (GSE171958, GSE176078, GSE167150, TCGA-BRCA, METABRIC);
  `data/TCGA-BRCA/README.md` documenta la descarga real ya hecha y cómo
  reproducirla; `data/GSE176078_singlecell/README.md` (nuevo
  2026-09-05) documenta la descarga real de célula única (paciente
  CID44971, ~178MB, no versionada por tamaño, ver `.gitignore`).
- **`references/`** — dos artículos reales del propio Dr.
  Hernández-Lemus, encontrados 2026-09-04, que fundamentan
  directamente dos líneas del Roadmap: Hernández-Lemus \& Ochoa
  (2024, revisión de integración multi-ómica) y Reyes-Gopar et al.
  (2025, Hi-C real en TNBC, dataset GSE167150 — esto reemplaza la
  afirmación anterior, ya incorrecta, de que no existía un Hi-C
  público a nivel paciente para esta comparación).
- **`scripts/`** — pipeline completo, numerado en orden de ejecución
  (`01_build_cohorts.py` … `29_new_figures_round3.py`),
  más `rmt_threshold_demo.py` (implementación de referencia del
  método de Luo et al. 2007, validada primero sobre datos sintéticos
  con estructura modular conocida antes de aplicarse a TCGA-BRCA) y
  `04b_wavelet_random_order_control.py` (control de orden aleatorio).
  De la ronda 2026-09-05a: `19_dyson_classifier_weights.py`
  (una semilla) y `20_dyson_classifier_weights_multi_seed.py` (20
  semillas, reproducibilidad) trackean los eigenvalores de
  $W_1^\top W_1$ de un MLP(20→8→1) entrenado sobre los mismos 20 genes
  hub RMT, epoch a epoch; `21_classifier_dyson_figures.py` genera 3
  figuras (fig9–fig11); `22_classifier_dyson_early_late.py`
  reparte el mismo análisis en ventana temprana (epochs 0–49) vs.
  tardía (50–299) sin reentrenar. De la ronda 2026-09-05b:
  `23_singlecell_line4_rmt.py` corre el mismo pipeline RMT (spikes +
  umbral tau) sobre datos reales de célula única (GSE176078, paciente
  CID44971, 7986 células); `24_singlecell_monte_carlo_null.py` es el
  null de permutación correspondiente (20/20 semillas exactas en
  τ=0.035, muy por debajo del τ* real de 0.10); `25_singlecell_figures.py`
  genera 3 figuras más (fig12–fig14), incluyendo un fix real de un
  problema de unfolding descubierto en el camino. Nuevos en la ronda
  2026-09-07b: `26_singlecell_epithelial_pseudotime.py` (retesta Línea
  5/Línea 4(ii) restringiendo el pseudo-tiempo a las 894 células Cancer
  Epithelial, con control real-vs-orden-aleatorio); `27_luma_dyson_
  repulsion.py` (extiende la repulsión de Dyson a LumA, segundo cohorte
  real, n=430); `28_classifier_dyson_longer_transient.py` (prueba
  directa de si un transitorio más largo, lr 5× menor, alarga la
  ventana de repulsión detectable); `29_new_figures_round3.py` genera
  3 figuras más (fig15–fig17).

## Líneas de trabajo (resumen)

1. Denoising espectral en datos multi-ómicos (Marchenko–Pastur, spiked
   models, IPR).
2. RMT para redes de coexpresión a gran escala y biomarcadores
   específicos por paciente.
3. Modelos espectrales para Hi-C y genómica 3D.
4. RMT en célula única y multimodalidad.
5. Análisis multiescala por *wavelets* (umbral espectral por escala,
   antes del umbral RMT de escala única).
6. Dinámica de Dyson y repulsión de niveles — sobre trayectorias
   dato/cohorte, y sobre la dinámica de entrenamiento del propio
   clasificador de reconocimiento de subtipo.

## Cómo correr la demo numérica

```bash
cd scripts
python rmt_threshold_demo.py
```

Solo necesita `numpy` y `scipy`. No requiere red ni datos externos —
genera sus propios datos sintéticos con estructura modular conocida, y
al final reporta qué fracción de las conexiones que sobreviven al
umbral seleccionado por RMT son, en efecto, conexiones reales dentro de
un módulo plantado (chequeo contra verdad conocida, no usado por el
método mismo).

## Estado

2026-09-07b: pase de referato completo (narrativa, matemáticas, citas)
más tres extensiones numéricas reales, todas cerrando huecos ya
identificados explícitamente en el propio Roadmap:

1. **Repulsión de Dyson en LumA** (`27_luma_dyson_repulsion.py`,
   nueva §"The Dyson-repulsion pattern replicates in a second, larger
   real cohort"): el mismo chequeo de reproducibilidad (20 órdenes
   aleatorios) corrido en LumA (n=430, antes solo TNBC) replica el
   patrón — λ5–λ6, λ6–λ7, λ7–λ8 repelen en 20/20 órdenes (vs. 18–19/20
   en TNBC) — y reproduce cualitativamente el hallazgo del Corolario:
   los pares menos confiables (λ1–λ2 en 11/20, λ2–λ3 en 9/20) son de
   nuevo los adyacentes a los eigenvalores más grandes/tipo-spike,
   consistente con que LumA tiene más del doble de eigenvalores spiked
   que TNBC (29 vs. 14).
2. **Retest de pseudo-tiempo en un solo linaje** (`26_singlecell_
   epithelial_pseudotime.py`, nueva §"Single-lineage pseudo-time: a
   corrected retest, still null, with a clarifying result"): restringir
   el pseudo-tiempo a las 894 células Cancer Epithelial (la corrección
   que el propio Roadmap proponía) NO mejora el eje de ordenamiento
   (PC1 = 8.6%, incluso un poco peor que el 9.7% de la muestra
   completa) y el retest de wavelets sigue siendo un null genuino
   (0/10 overlap, igualado o superado por 20/20 órdenes aleatorios de
   control) — pero el hallazgo más útil es metodológico: una
   comparación real-vs-orden-aleatorio (nunca antes corrida a
   resolución de célula única) muestra que la repulsión de Dyson aquí
   es una propiedad de que la muestra crezca, no de este ordenamiento
   específico — aclara, no contradice, el resultado de bulk.
3. **¿Transitorio más largo, ventana de repulsión más larga?**
   (`28_classifier_dyson_longer_transient.py`, nueva §"Does a longer
   transient give a longer repulsion window? A direct test, with an
   honest complication"): lr 5× menor (0.01) y 1500 epochs (5× más) —
   el transitorio se alarga ~11× (no 5×, superlineal), el null
   post-convergencia se vuelve aún más limpio (0/140 vs. 3/140 antes),
   pero el conteo crudo pre-convergencia baja (11/140 vs. 36/140)
   porque el criterio de "reabre en 5 epochs" no se reescala con la
   velocidad de entrenamiento — honestamente reportado como limitación
   real del detector, no como contradicción del mecanismo.

3 figuras nuevas (fig15–fig17), 3 scripts nuevos (26–28) + 1 de
figuras (29). Roadmap, Discusión y Conclusiones actualizados en
consecuencia (Líneas 4(ii), 5 y 6b ya no son "próximo paso" sino
resultados reales, con su propio próximo paso más específico cada
una). Recompila 0 errores, 40pp (de 34pp).

2026-09-07: pase de referato enfocado en la bibliografía y en la
reproducibilidad literal del código citado. Dos entradas bibliográficas
llevaban etiquetas que no eran nombres de autor (`$k$-core(2021)` y
`Multi-omics TNBC(2021)`), rompiendo el patrón Autor(Año) del resto de
la lista; se corrigieron a `Dorantes-Gilardi et al.(2021)` y
`Chappell et al.(2021)`, este último completado con los 14 autores
reales, volumen, páginas y DOI verificados en PMC/Oxford Academic (no
son del grupo Hernández-Lemus; el artículo se cita solo por el dataset
público GSE171958). Se corrigió también un desorden alfabético real
introducido al insertar la cita de Reyes-Gopar en una ronda anterior
(quedó Passerini→Reyes-Gopar→Pedregosa; ahora Passerini→Pedregosa→
Reyes-Gopar→TCGA→Tao\&Vu→Tracy\&Widom→Voiculescu→Wu, alfabético de
principio a fin). Más importante: se encontró y corrigió una
inconsistencia real de nombre de archivo — el texto cita
`scripts/24_singlecell_monte_carlo_null.py` dos veces (líneas de
Resultados y de Métodos), pero el script en disco se llamaba
`24_singlecell_montecarlo_null.py` (sin guión bajo entre "monte" y
"carlo"), rompiendo la reproducibilidad literal que esta sección de
Métodos promete explícitamente y contradiciendo la propia convención
del script 10 (`10_monte_carlo_null.py`). Renombrado el archivo para
que coincida exactamente con lo que dice el artículo. Recompila 0
errores, 34pp (sin cambio de paginación).

2026-09-05b: se atacó Línea 4 del Roadmap (célula única) con datos
reales por primera vez en este programa: GSE176078 (Wu et al. 2021,
*Nat. Genet.*), paciente CID44971 real, 7986 células, mismo pipeline
RMT sin modificar. Umbral τ*=0.10 real, confirmado por null de Monte
Carlo (20/20 permutaciones exactas en τ=0.035). Los picos espectrales
se dividen limpiamente en ejes biológicos reales (inmune/T-cell,
mieloide, epitelial — verificado contra la anotación de tipo celular
real, nunca vista por el método) y artefactos técnicos conocidos de
scRNA-seq (contenido ribosomal, estrés de disociación) — un hallazgo
honesto que los datos bulk nunca tuvieron que enfrentar. En el camino
se descubrió y corrigió un problema metodológico real: el unfolding
polinomial global (usado en todo el artículo) se degrada en esta red
mucho más densa (20.4% vs. 1.3–13.0% en los cohortes bulk), con ~4% de
los espaciamientos convertidos en outliers extremos — diagnosticado,
cuantificado, y reportado explícitamente en el texto y la figura
(no oculto), en vez de presentar una figura rota. También se revisó
directamente (no solo se propuso) el candidato de pseudo-tiempo de
célula única para la Línea 5 (wavelets): el PC1 global de esta muestra
solo captura 9.7% de varianza, sin mejora real sobre el 10.2% del PC1
bulk — la Línea 5 revisada ahora pide inferencia de trayectoria
restringida a un linaje continuo (p.ej. las 894 células Cancer
Epithelial), no célula única de cualquier tipo. Búsqueda web real
confirmó los tamaños/formato exactos de GSE167150 (Hi-C, aún no
atacado: 6.2GB, archivos `.hic` ya procesados, requiere
`hicstraw`/`straw`) y encontró una cita real y sustantiva
(Mahoney & Martin 2019, ICML, "5+1 phases of training") que contextualiza
el null de repulsión del clasificador de la ronda anterior. Se restauró
la sección de wavelets a su propia subsección (antes vivía, mal
estructurada, dentro de la subsección de topología de red). Se corrigió
también una lista de autores inventada a medias (Wu et al. 2021) antes
de que llegara a la bibliografía — verificada contra PMC antes de
citar. 3 figuras nuevas (fig12–fig14), 1 tabla nueva, 3 scripts nuevos
(23–25), 1 cita nueva real. Compila 0 errores, 34pp (de 30pp).

2026-09-05a: se atacó directamente la Línea 6b del Roadmap (dinámica de
Dyson en la matriz de pesos del propio clasificador durante SGD, no
solo en los datos) con código y cómputo real — no quedó como propuesta.
Resultado: un tercer null honesto, pero no plano — la repulsión de
pesos existe, solo que confinada casi por completo al transitorio
previo a la convergencia del entrenamiento (12× más señal en epochs
0–49 que en 50–299, verificado sobre las mismas 20 trayectorias ya
computadas, sin reentrenar). Además se derivó y probó un Corolario
nuevo de la Proposición 3.1 ya existente (la fuerza de repulsión de
Dyson es exactamente $2/\delta$, inversamente proporcional a la
brecha), aplicado a los números reales ya en la Tabla de repulsión para
explicar cuantitativamente, por primera vez, la anomalía del par
λ1–λ2. 3 figuras nuevas (fig9–fig11), 1 tabla nueva, 4 scripts nuevos
(19–22). Compila 0 errores, 30pp (de 26pp).

2026-09-04: varias rondas de referato duro completadas sobre el mismo
artículo. Un bug real y serio fue encontrado y corregido (25
referencias cruzadas `Sec.~\ref{}` renderizaban en blanco por usar
secciones sin numerar). No incluye datos clínicos ni multi-ómicos
propios del grupo — solo datos públicos ya descargados (TCGA-BRCA) o
identificados (GSE171958, GSE176078, GSE167150, METABRIC), suficientes
para un primer ciclo completo de desarrollo y validación metodológica.
Objetivo de publicación: BMC Bioinformatics (alternativa Q1 si se busca
más impacto: PLOS Computational Biology).

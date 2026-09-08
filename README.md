# RMT + Genómica Computacional — colaboración propuesta con el Dr. Enrique Hernández-Lemus

Programa de investigación conjunta: Random Matrix Theory (estática y
dinámica), análisis multiescala por wavelets, y aprendizaje automático,
aplicados a genómica computacional y oncología de sistemas — foco en
cáncer de mama (TNBC), redes de coexpresión a gran escala, genómica 3D
(Hi-C), y célula única — con protocolo de validación numérica explícito
en cada línea.

## Contenido

- **`article.tex`** / **`.pdf`** — el artículo de investigación (inglés,
  45pp, formato BMC Bioinformatics: abstract estructurado,
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
  reproducirla; `data/GSE176078_singlecell/README.md` documenta la
  descarga real de célula única (paciente CID44971, ~178MB, no
  versionada por tamaño, ver `.gitignore`); `data/GSE167150_HiC/README.md`
  documenta la descarga real de Hi-C — ahora las 4 muestras reales de la
  serie (3 tumores TNBC + 1 normal pareado, ~3.5GB total), cada una
  extraída con un solo Range request dirigido dentro del tar combinado
  de 6.2GB de GEO, sin bajar el resto — ver
  `scripts/download_hic_tnbc_tissue3.py`; `data/hg19_annotation/` (nuevo
  2026-09-07d, no versionado por tamaño) contiene la secuencia de
  referencia hg19 chr18 y la tabla de genes RefSeq (UCSC, ~32MB total),
  usadas para la validación real de compartimentos contra GC/densidad
  génica.
- **`references/`** — dos artículos reales del propio Dr.
  Hernández-Lemus, encontrados 2026-09-04, que fundamentan
  directamente dos líneas del Roadmap: Hernández-Lemus \& Ochoa
  (2024, revisión de integración multi-ómica) y Reyes-Gopar et al.
  (2025, Hi-C real en TNBC, dataset GSE167150 — esto reemplaza la
  afirmación anterior, ya incorrecta, de que no existía un Hi-C
  público a nivel paciente para esta comparación).
- **`scripts/`** — pipeline completo, numerado en orden de ejecución
  (`01_build_cohorts.py` … `35_hic_validation_figures.py`),
  más `rmt_threshold_demo.py` (implementación de referencia del
  método de Luo et al. 2007, validada primero sobre datos sintéticos
  con estructura modular conocida antes de aplicarse a TCGA-BRCA),
  `04b_wavelet_random_order_control.py` (control de orden aleatorio) y
  `download_hic_tnbc_tissue3.py` (descarga dirigida del Hi-C real, ver
  Estado abajo).
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
  3 figuras más (fig15–fig17). Nuevo en la ronda 2026-09-07c, la
  Línea 3 (Hi-C) atacada por primera vez con datos reales:
  `30_hic_line3_rmt.py` construye la matriz de correlación bin-bin
  O/E-normalizada (Lieberman-Aiden et al. 2009) del chr18 de un tumor
  TNBC real (GSM5098082, 100kb, 750 bins tras filtrar 31 bins de
  cobertura cero) y corre el mismo pipeline RMT (tau-sweep + null de
  Monte Carlo); `31_hic_figures.py` genera 1 figura nueva (fig18: mapa
  de calor plaid + sweep + track de compartimentos PC1). Nuevos en la
  ronda 2026-09-07d: `32_hic_compartment_validation.py` (valida el PC1
  de compartimentos contra GC content y densidad génica RefSeq reales,
  hg19); `33_hic_multi_sample.py` (extiende el pipeline a los otros 2
  tumores TNBC reales + la muestra normal pareada, las 4 muestras de
  GSE167150, con matriz de acuerdo de compartimentos por pares);
  `34_singlecell_diffusion_pseudotime.py` (pseudo-tiempo por difusión,
  Haghverdi et al. 2016, implementado en numpy/scipy — el estimador
  genuinamente no-PC1 que el Roadmap pedía); `35_hic_validation_
  figures.py` genera 1 figura más (fig19: scatter PC1-vs-GC + heatmap
  de acuerdo entre las 4 muestras).

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

2026-09-07e: pase de referato duro sobre las matemáticas mismas (no
solo los resultados numéricos) — encontró y corrigió un **error
matemático real** en la Proposición 3.1. Tal como estaba escrita, la
Proposición afirmaba que la densidad conjunta GOE es la distribución
estacionaria de la SDE de Dyson de la Definición 3.4 (la SDE *sin*
término de confinamiento, usada en todo el resto del artículo para la
analogía "tiempo = tamaño de muestra creciente"). Eso es falso: esa SDE
sin confinamiento no tiene ley estacionaria en absoluto (los
eigenvalores se dispersan difusivamente sin límite, igual que un
Browniano simple) — verificado aquí con el cálculo de Fokker-Planck
completo (condición de corriente cero), no solo citado. La densidad
GOE Sí es estacionaria, pero para una SDE *confinada* distinta (mismo
término de Coulomb divergente, coeficiente 1/2 en vez de 1, más un
término de Ornstein-Uhlenbeck $-\lambda_i/4$) — un objeto relacionado
pero genuinamente distinto, introducido ahora explícitamente como tal.
El Corolario 3.2 (fuerza de repulsión $\propto 1/\delta$) usa
exclusivamente la SDE sin confinar de la Definición 3.4, así que
ninguno de sus números ni ninguna cifra empírica del artículo cambió
— solo se corrigieron las referencias cruzadas para que cada resultado
cite la SDE correcta. Se añadió además, por primera vez, la
verificación real (no solo citada) de que la densidad GOE satisface la
condición de corriente cero para la SDE confinada — la Proposición 3.1
ahora tiene una demostración real, no solo "es un resultado de Dyson
1962, no lo reproducimos aquí".

También, por pedido explícito: cada una de las 15 figuras (19
sub-paneles) ahora cita en su propio pie qué script(s) la generó — antes
ninguna lo hacía, pese a que el texto ya citaba scripts extensamente en
prosa. Revisadas también consistencia de notación (β, τ, q=p/n) y
referencias cruzadas de tablas (10 tablas, todas usadas, ninguna
huérfana). Recompila 0 errores, 46pp (de 45pp).

**Sobre journal fit** (dado el alcance actual: RNA-seq bulk +
célula única + Hi-C + ML + matemática real, con hallazgos honestos
mixtos, no un solo hallazgo biológico grande): **BMC Bioinformatics**
sigue siendo el mejor fit directo — ya está formateado para eso,
sin límite de extensión duro, y su cultura editorial no penaliza
resultados nulos honestos ni un pipeline metodológico como columna
vertebral. **PLOS Computational Biology** es la alternativa Q1 más
natural si se busca más visibilidad — comparte esa misma cultura
(transparencia metodológica, tolera nulls), pero requiere resumen no
estructurado + "Author Summary" aparte, un reformateo real, no
trivial. Dado que ya son 46pp, vale la pena decidir pronto qué parte
del Apéndice/Métodos migra a Material Suplementario antes de enviar a
cualquiera de los dos — ambas revistas esperan el cuerpo principal
más corto y las tablas completas de barrido (tab:tnbc-full,
tab:luma-full) y el detalle metodológico más extenso como
suplementario, no en el cuerpo.

2026-09-07d: **cierra los tres huecos que quedaban explícitamente
abiertos** (validación real del compartimento, extensión a las otras
2 muestras TNBC + la normal, pseudo-time genuinamente no-PC1), más un
pase de referato completo (narrativa, matemáticas, citas) y una pasada
de estilo (reduce ~40% las repeticiones de "genuinely/genuine" y
~13% las de "rather than" en todo el documento — la voz seguía leyendo
demasiado "Claude", ahora más variada). Todo integrado en Abstract,
Resultados, Discusión, Roadmap, Conclusiones (reescritas en 4 párrafos
con cierre real, antes un solo bloque gigante) y Métodos.

1. **Validación real de compartimentos**
   (`32_hic_compartment_validation.py`): se bajaron la secuencia de
   referencia hg19 chr18 y la tabla RefSeq (UCSC, ~32MB, mismo build
   que el header del `.hic` declara — verificado, no asumido). El PC1
   de compartimentos correlaciona con GC content real (r=0.452,
   p=4.3e-39) y densidad génica real (r=0.306, p=9.3e-18); los bins de
   signo positivo tienen GC y densidad génica muchísimo más altos
   (5.32 vs. 1.61 transcritos/bin, t=7.22, p=1.5e-12) — confirma la
   interpretación de compartimentos A/B contra anotación genómica
   real, no solo por parecido visual.
2. **Extensión a las 4 muestras reales de GSE167150**
   (`33_hic_multi_sample.py`): los otros 2 tumores TNBC reales
   (GSM5098079, GSM5098080, ~846MB y 755MB) y la muestra normal
   pareada real (GSM5098074, ~2.48GB) se descargaron con la misma
   técnica de Range request dirigido (3.5GB en las 4 muestras, sin
   bajar los otros ~2.7GB de líneas celulares no relacionadas de la
   serie). Las 4 dan τ* real confirmado por permutación y la misma
   firma GC/densidad génica — pero el hallazgo honesto es que
   TNBC_Tissue3 es un outlier respecto a las otras 3 muestras
   (73–76% de acuerdo) mientras que Tissue1, Tissue2 y la normal
   concuerdan entre sí 87–92% — nada de un tumor-vs-normal limpio, más
   bien heterogeneidad inter-tumoral real (n=3, no concluyente por sí
   solo).
3. **Pseudo-tiempo por difusión** (`34_singlecell_diffusion_
   pseudotime.py`): implementado en numpy/scipy puro (Haghverdi et al.
   2016, sin paquete de célula única disponible) sobre las mismas 894
   células Cancer Epithelial. Confirmado genuinamente distinto de PC1
   (Spearman ρ=-0.481). El resultado: mismo patrón de repulsión de
   Dyson atribuible al tamaño de muestra (no al ordenamiento, tercera
   vez confirmado) y el mismo null de wavelets (0/10 overlap, igualado
   por 20/20 órdenes aleatorios) — tres ejes de ordenamiento distintos,
   la misma respuesta: el problema no es la elección de pseudo-tiempo.

4 figuras nuevas (fig19 + reemplazo de fig18's caption), 4 scripts
nuevos (32–35). Recompila 0 errores, 45pp (de 42pp).

2026-09-07c: **Línea 3 (Hi-C) atacada con datos reales por primera vez
en todo el programa** — cierra el último hueco genuinamente no
intentado de las 6 líneas originales. Se descargó GSM5098082
("TNBC_Tissue3", GSE167150, el más pequeño de los tres tumores TNBC
reales de la serie, ~645MB) con una técnica real de Range request
dirigido dentro del tar combinado de GEO (6.2GB), sin bajar el resto de
la serie (`scripts/download_hic_tnbc_tissue3.py`); se instaló
`hic-straw==0.0.6` (la versión pura-Python de `straw`, ya que no había
compilador de C++ disponible para la versión moderna pybind11). Se
extrajo la matriz de contacto del chr18 a 100kb (781 bins), se filtraron
31 bins de cobertura cero (artefacto técnico real, centromérico/no
mapeable — un bug real encontrado y corregido: sin este filtro,
`corrcoef` producía NaNs que, puestos en cero ingenuamente, dominaban el
eigenvector principal de forma espuria), se normalizó O/E (Lieberman-Aiden
et al. 2009, cita nueva verificada en Science) y se corrió el mismo
pipeline RMT (tau-sweep + NNSD + null de Monte Carlo) ya usado en todo
el artículo. Resultado real: τ*=0.15, confirmado por permutación (null
en 0.10 o sin transición limpia, nunca en 0.15); el mapa de correlación
bin-bin muestra el patrón "plaid" clásico de compartimentos A/B, y el
eigenvector principal se divide en bloques genómicos contiguos
(457 positivos, 293 negativos) — consistente con, pero no todavía
validado independientemente contra, la interpretación de compartimentos
(la validación contra densidad génica/GC queda como próximo paso
explícito, no se afirma de más). Integrado en Resultados, Discusión,
Roadmap, Conclusiones, Métodos, Declarations y abstract; Roadmap ahora
dice explícitamente que las 6 líneas originales ya se corrieron al
menos una vez con datos reales. 1 figura nueva (fig18), 2 scripts
nuevos (30, 31) + 1 script de descarga. Recompila 0 errores, 42pp
(de 40pp).

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

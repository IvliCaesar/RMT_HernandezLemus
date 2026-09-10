# RMT + Genómica Computacional — con Enrique Hernández-Lemus

Colaboración con el Dr. Enrique Hernández-Lemus (INMEGEN): teoría de
matrices aleatorias (RMT), dinámica de Dyson-Brownian-motion, wavelets
multiescala y aprendizaje automático aplicados a cáncer de mama. El
artículo compara redes de coexpresión génica en las cohortes TNBC y
Luminal A de TCGA-BRCA, y extiende el mismo pipeline a célula única y a
genómica 3D (Hi-C).

## Hallazgos principales

- El umbral RMT (Luo et al. 2007) separa TNBC de Luminal A en
  TCGA-BRCA; un null de Monte Carlo por permutación confirma que la
  separación no es un artefacto del umbral.
- Los eigenvalores "spiked" que exceden el borde de Marchenko–Pastur
  corresponden a programas génicos biológicamente coherentes
  (colágeno/ECM, interferón, queratina), no a ruido.
- Las trayectorias de eigenvalores muestran repulsión de niveles tipo
  Dyson al crecer la muestra, reproducible en 20 órdenes aleatorios
  distintos, en TNBC y en LumA. Un Corolario nuevo — la fuerza de
  repulsión es inversamente proporcional a la brecha entre eigenvalores
  — explica por qué el par λ1–λ2 es el único no confiablemente
  repulsivo.
- El mismo pipeline corre también sobre célula única (GSE176078) y
  sobre Hi-C (GSE167150, 4 muestras): ambos dan umbrales reales
  confirmados por permutación. La interpretación de compartimentos A/B
  en Hi-C está validada de forma independiente contra contenido de GC y
  densidad génica reales, no solo por parecido visual.
- La matriz de pesos del propio clasificador de subtipo también repele
  en el sentido de Dyson, pero solo durante el transitorio previo a la
  convergencia del entrenamiento (12× más señal ahí que después): un
  hallazgo sobre la dinámica de entrenamiento, no solo sobre los datos.
- Tres resultados nulos, reportados como tales en vez de descartados:
  la descomposición por wavelets falla en datos reales, pero un control
  sintético confirma que el método funciona — el problema es la señal
  disponible, no el método; el clasificador TNBC-vs-LumA está cerca de
  un techo que ninguna selección de variables logra romper; la
  repulsión de pesos del clasificador no aparece fuera de su
  transitorio inicial.
- Las seis líneas de trabajo propuestas originalmente para el programa
  ya corrieron sobre datos reales al menos una vez, con código y datos
  reproducibles.

## Contenido

- **`article.tex`** / **`.pdf`** — el artículo (inglés, formato PLOS
  Computational Biology), con el análisis completo descrito arriba.
- **`programa_investigacion_RMT_HernandezLemus.tex`** / **`.pdf`** — la
  propuesta original: justificación científica, hipótesis, seis líneas
  de trabajo, núcleo matemático, protocolo de validación, datos,
  productos esperados.
- **`data/`** — conjuntos de datos, cada uno con su propio README
  documentando la descarga real y cómo reproducirla:
  `TCGA-BRCA/`, `GSE176078_singlecell/` (célula única, paciente
  CID44971), `GSE167150_HiC/` (4 muestras Hi-C reales: 3 tumores TNBC y
  1 normal pareada), `hg19_annotation/` (referencia y tabla RefSeq
  usadas para validar compartimentos).
- **`references/`** — dos artículos del propio Dr. Hernández-Lemus que
  fundamentan directamente dos líneas del programa: Hernández-Lemus \&
  Ochoa (2024, integración multi-ómica) y Reyes-Gopar et al. (2025,
  Hi-C en TNBC, GSE167150).
- **`scripts/`** — pipeline completo, numerado en orden de ejecución
  (`01_build_cohorts.py` … `35_hic_validation_figures.py`), más
  `rmt_threshold_demo.py` (implementación de referencia de Luo et al.
  2007, validada primero sobre datos sintéticos con estructura modular
  conocida).

## Líneas de trabajo

1. Denoising espectral en datos multi-ómicos (Marchenko–Pastur, spiked
   models, IPR).
2. RMT para redes de coexpresión a gran escala y biomarcadores
   específicos por paciente.
3. Modelos espectrales para Hi-C y genómica 3D.
4. RMT en célula única y multimodalidad.
5. Análisis multiescala por wavelets (umbral espectral por escala,
   antes del umbral RMT de escala única).
6. Dinámica de Dyson y repulsión de niveles, en las trayectorias de
   datos y en la dinámica de entrenamiento del propio clasificador.

## Cómo correr la demo numérica

```bash
cd scripts
python rmt_threshold_demo.py
```

Solo necesita `numpy` y `scipy`. No requiere red ni datos externos:
genera sus propios datos sintéticos con estructura modular conocida y
reporta qué fracción de las conexiones que sobreviven al umbral
seleccionado por RMT son, en efecto, conexiones reales dentro de un
módulo plantado.

## Estado

Manuscrito completo, formateado para PLOS Computational Biology.
Código y datos públicos en este repositorio; los conjuntos de datos
grandes (célula única, Hi-C, anotación hg19) no están versionados por
tamaño — cada uno tiene su propio README con instrucciones de
descarga.

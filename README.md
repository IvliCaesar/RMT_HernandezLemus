# RMT + Genómica Computacional — colaboración propuesta con el Dr. Enrique Hernández-Lemus

Programa de investigación conjunta: Random Matrix Theory (estática y
dinámica), análisis multiescala por wavelets, y aprendizaje automático,
aplicados a genómica computacional y oncología de sistemas — foco en
cáncer de mama (TNBC), redes de coexpresión a gran escala, genómica 3D
(Hi-C), y célula única — con protocolo de validación numérica explícito
en cada línea.

## Contenido

- **`article.tex`** / **`.pdf`** — el artículo de investigación (inglés,
  26pp, formato BMC Bioinformatics: abstract estructurado,
  referencias numeradas, Declarations completas), con análisis real
  sobre TCGA-BRCA: umbral RMT (Luo et al. 2007) comparando TNBC vs.
  Luminal A, null de Monte Carlo por permutación, eigenvalores
  "spiked" + entropía espectral + escalamiento Tracy-Widom, red de
  coexpresión vs. nulos aleatorios (Erdős–Rényi/configuración),
  trayectorias de eigenvalores tipo Dyson (repulsión de niveles,
  verificada en 20 órdenes aleatorios distintos), una Proposición
  real (la densidad conjunta GOE es la distribución estacionaria del
  SDE de Dyson — resultado de Dyson 1962, no antes enunciado
  formalmente en este trabajo), y una prueba de reconocimiento por ML
  (incluye GA) — con dos resultados nulos honestos reportados como
  tales (wavelets en datos reales, techo de clasificación). Roadmap
  explícito y ya actualizado con datasets reales para cada línea
  pendiente (ver abajo).
- **`programa_investigacion_RMT_HernandezLemus.tex`** / **`.pdf`** — la
  propuesta original: justificación científica, hipótesis, 6 líneas de
  trabajo, núcleo matemático, protocolo de validación, datos, productos
  esperados, ruta de 12–18 meses.
- **`data/README.md`** — los conjuntos de datos públicos identificados
  (GSE171958, GSE176078, GSE167150, TCGA-BRCA, METABRIC);
  `data/TCGA-BRCA/README.md` documenta la descarga real ya hecha y cómo
  reproducirla.
- **`references/`** — dos artículos reales del propio Dr.
  Hernández-Lemus, encontrados 2026-09-04, que fundamentan
  directamente dos líneas del Roadmap: Hernández-Lemus \& Ochoa
  (2024, revisión de integración multi-ómica) y Reyes-Gopar et al.
  (2025, Hi-C real en TNBC, dataset GSE167150 — esto reemplaza la
  afirmación anterior, ya incorrecta, de que no existía un Hi-C
  público a nivel paciente para esta comparación).
- **`scripts/`** — pipeline completo, numerado en orden de ejecución
  (`01_build_cohorts.py` … `18_repulsion_reproducibility_figure.py`),
  más `rmt_threshold_demo.py` (implementación de referencia del
  método de Luo et al. 2007, validada primero sobre datos sintéticos
  con estructura modular conocida antes de aplicarse a TCGA-BRCA) y
  `04b_wavelet_random_order_control.py` (control de orden aleatorio).

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

2026-09-04: varias rondas de referato duro completadas sobre el mismo
artículo. Compila 0 errores, 26pp. Un bug real y serio fue encontrado
y corregido (25 referencias cruzadas `Sec.~\ref{}` renderizaban en
blanco por usar secciones sin numerar). No incluye datos clínicos ni
multi-ómicos propios del grupo — solo datos públicos ya descargados
(TCGA-BRCA) o identificados (GSE171958, GSE176078, GSE167150,
METABRIC), suficientes para un primer ciclo completo de desarrollo y
validación metodológica. Objetivo de publicación: BMC Bioinformatics
(alternativa Q1 si se busca más impacto: PLOS Computational Biology).

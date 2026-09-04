# RMT + Genómica Computacional — colaboración propuesta con el Dr. Enrique Hernández-Lemus

Programa de investigación conjunta: Random Matrix Theory (estática y
dinámica), análisis multiescala por wavelets, y aprendizaje automático,
aplicados a genómica computacional y oncología de sistemas — foco en
cáncer de mama (TNBC), redes de coexpresión a gran escala, genómica 3D
(Hi-C), y célula única — con protocolo de validación numérica explícito
en cada línea.

## Contenido

- **`article.tex`** / **`.pdf`** — el artículo de investigación (inglés,
  14pp, formato tipo BMC Bioinformatics), con análisis real sobre
  TCGA-BRCA: umbral RMT (Luo et al. 2007) comparando TNBC vs. Luminal A,
  chequeo de robustez a tamaño de muestra, descomposición multiescala
  por wavelets, trayectorias de eigenvalores tipo Dyson (repulsión de
  niveles), y una prueba de reconocimiento por ML — incluyendo un
  resultado nulo honesto (ver abstract). Roadmap explícito para lo que
  aún no está corrido con datos reales (Hi-C, célula única, dinámica de
  Dyson del propio clasificador).
- **`programa_investigacion_RMT_HernandezLemus.tex`** / **`.pdf`** — la
  propuesta original: justificación científica, hipótesis, 6 líneas de
  trabajo, núcleo matemático, protocolo de validación, datos, productos
  esperados, ruta de 12–18 meses.
- **`data/README.md`** — los 4 conjuntos de datos públicos identificados
  (GSE171958, GSE176078, TCGA-BRCA, METABRIC); `data/TCGA-BRCA/README.md`
  documenta la descarga real ya hecha y cómo reproducirla.
- **`scripts/`** — pipeline completo, numerado en orden de ejecución
  (`01_build_cohorts.py` … `08_spiked_eigenvalues.py`), más
  `rmt_threshold_demo.py` (implementación de referencia del método de
  Luo et al. 2007, validada primero sobre datos sintéticos con
  estructura modular conocida antes de aplicarse a TCGA-BRCA).

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

2026-09-04: primer borrador de artículo completo, corrido sobre datos
reales de TCGA-BRCA (no sintéticos), con Enrique Hernández-Lemus
confirmado como coautor. Compila 0 errores, 14pp. No incluye datos
clínicos ni multi-ómicos propios del grupo — solo datos públicos ya
descargados (TCGA-BRCA) o identificados (GSE171958, GSE176078,
METABRIC), suficientes para un primer ciclo completo de desarrollo y
validación metodológica.

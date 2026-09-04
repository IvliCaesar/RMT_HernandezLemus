# RMT + Genómica Computacional — colaboración propuesta con el Dr. Enrique Hernández-Lemus

Programa de investigación conjunta: Random Matrix Theory (estática y
dinámica), análisis multiescala por wavelets, y aprendizaje automático,
aplicados a genómica computacional y oncología de sistemas — foco en
cáncer de mama (TNBC), redes de coexpresión a gran escala, genómica 3D
(Hi-C), y célula única — con protocolo de validación numérica explícito
en cada línea.

## Contenido

- **`programa_investigacion_RMT_HernandezLemus.tex`** / **`.pdf`** — la
  propuesta completa: justificación científica, hipótesis, 6 líneas de
  trabajo, núcleo matemático, protocolo de validación, datos, productos
  esperados, ruta de 12–18 meses.
- **`data/README.md`** — los 4 conjuntos de datos públicos identificados
  (GSE171958, GSE176078, TCGA-BRCA, METABRIC), con instrucciones de
  acceso y qué línea de trabajo usa cada uno.
- **`scripts/rmt_threshold_demo.py`** — implementación de referencia,
  real y verificada, del método clásico de umbral espectral de Luo et
  al. (BMC Bioinformatics 2007): construye datos sintéticos con
  estructura modular conocida, calcula el borde de Marchenko–Pastur, el
  *inverse participation ratio*, y localiza la transición
  Poisson→Wigner-Dyson (GOE) que selecciona el umbral de correlación —
  sin tocar datos reales, como primera prueba de correctitud.

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

Primera versión del programa, lista para servir como nota conceptual de
arranque (ver "Primeros pasos sugeridos" en la propuesta). No incluye
datos clínicos ni multi-ómicos propios del grupo del Dr. Hernández-Lemus
— solo datos públicos ya identificados, suficientes para un primer ciclo
completo de desarrollo y validación metodológica.

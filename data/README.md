# Datos públicos identificados para el programa

Ninguno de estos conjuntos pertenece al grupo del Dr. Hernández-Lemus ni
requiere permiso especial: los cuatro son de acceso público y bastan para
arrancar las Líneas 1, 2, 4, 5 y 6(a) del programa sin esperar datos propios.
Ver `programa_investigacion_RMT_HernandezLemus.tex`, sección "Datos", para
cómo se usa cada uno.

## GSE171958 — TNBC multi-ómica (metilación + RNA-seq + proteómica)
- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171958>
- Proteómica pareada: ProteomeXchange **PXD025238**
- Líneas celulares: MCF10A (no tumoral), MDA-MB-231, HCC1937 (TNBC)
- Capas: metilación de DNA, RNA-seq, proteína, fosfoproteómica, marcas de
  histonas
- Fuente: *Multi-omics data integration reveals correlated regulatory
  features of triple negative breast cancer*, Molecular Omics (2021),
  DOI: 10.1039/D1MO00117E
- Uso en el programa: Línea 1 (denoising multi-ómico), Línea 5 (wavelets,
  matriz de correlación entre capas)

## GSE176078 — TNBC scRNA-seq (10 pacientes)
- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078>
- Uso en el programa: Línea 4 (célula única), Línea 6(a) (dinámica de
  Dyson a través de pacientes/pseudotiempo)

## TCGA-BRCA
- Portal: <https://portal.gdc.cancer.gov/projects/TCGA-BRCA>
- Expresión, metilación, mutación, variables clínicas; cohorte de
  referencia grande (~1100 pacientes)
- Uso en el programa: cohorte primaria de validación cruzada (§Validación
  y numérica, punto 3)

## METABRIC
- Acceso: cBioPortal, <https://www.cbioportal.org/study/summary?id=brca_metabric>
- Cohorte pública independiente de TCGA-BRCA, con subtipificación clínica
- Uso en el programa: cohorte de validación cruzada secundaria

## GSE167150 — TNBC Hi-C (paciente + tejido contralateral sano)
- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167150>
- Procesado con HiC-Pro, resolución 40kb, redes intracromosómicas
  por cromosoma (modelo hipergeométrico no central para significancia)
- Fuente: Reyes-Gopar, Pérez-Fuentes, Bendall y Hernández-Lemus,
  *Integration of chromosome conformation and gene expression networks
  reveals regulatory mechanisms in triple negative breast cancer*,
  Frontiers in Cell and Developmental Biology 13:1597245 (2025) — del
  propio grupo del Dr. Hernández-Lemus, encontrado 2026-09-04
  (`references/ReyesGopar_etal_2025_HiC_TNBC_networks.pdf`)
- Uso en el programa: Línea 3 (Hi-C/3D) — ver "Roadmap" en
  `article.tex`; esto **reemplaza** la afirmación anterior (ya
  incorrecta) de que no existía un dataset Hi-C público a nivel
  paciente para esta comparación

## Estado: TCGA-BRCA ya descargado y analizado (2026-09-04)

`article.tex` (raíz del repo) usa datos reales de TCGA-BRCA, no
sintéticos. Ver `data/TCGA-BRCA/README.md` para las URLs exactas y
cómo reproducir la descarga — los archivos crudos (69MB) están en
`.gitignore` y no viven en el repositorio, solo el código que los
reconstruye (`scripts/01_build_cohorts.py`).

GSE171958, GSE176078, GSE167150 y METABRIC siguen sin descargar — ver
la sección "Roadmap" de `article.tex` para el plan concreto de cada
una (Líneas 1-extensión, 3, 4), y no como un pendiente vago: GSE176078
es el siguiente paso directo para Línea 4, GSE167150 para Línea 3.

## Nota sobre descarga
Ninguno de estos archivos está incluido en este repositorio (son
demasiado grandes y algunos requieren aceptar términos de uso en el
portal correspondiente). Este directorio es el punto de destino
convenido para cuando se descarguen: `data/GSE171958/`,
`data/GSE176078/`, `data/TCGA-BRCA/`, `data/METABRIC/` — cada uno con
su propio `README.md` de procedencia y fecha de descarga, siguiendo la
misma disciplina de todo el proyecto: cada afirmación numérica debe
poder rastrearse hasta un archivo de datos real, nunca a un número
inventado o recordado.

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

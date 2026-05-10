# DNA Methylation Analysis Pipeline for *Heliconius erato* (v1.0)

Este repositorio contiene el flujo de trabajo bioinformático desarrollado para la identificación y caracterización de modificaciones de ADN (**5mC** y **5hmC**) en cerebros de *Heliconius erato*. El pipeline procesó ~90 millones de sitios CpG únicos, evaluando el efecto del inhibidor **RG108** sobre el paisaje epigenético.

## 1. Especificaciones Técnicas y Versiones
Para garantizar la reproducibilidad científica, se detallan las versiones validadas en este estudio:

* **Llamada de Modificaciones:** `modkit v0.5.0` (Oxford Nanopore Technologies).
* **Manejo de Alineamientos:** `samtools v1.16.1` y `minimap2 v2.24`.
* **Aritmética Genómica:** `bedtools v2.30.0`.
* **Análisis Estadístico:** `R v4.2.2` con el paquete `DSS v2.46.0` (Bioconductor).
* **Visualización:** `ggplot2`, `scales` y `RColorBrewer`.

## 2. Lógica de Filtrado y Parámetros Críticos
* **Umbral Probabilístico (0.79):** Definido a partir del percentil 10 de la distribución empírica de probabilidades de modkit. Este valor asegura un equilibrio entre sensibilidad y especificidad, alineándose con los umbrales de ~0.81 recomendados por ONT para citosinas no modificadas.
* **Cobertura Mínima:** 10x por sitio CpG en ambas condiciones.
* **Modelo Estadístico:** Test de Wald basado en una distribución Beta-binomial con estimación de dispersión compartida (`equal.disp = TRUE`), optimizado para pooling de muestras sin réplicas.

## 3. Arquitectura del Proyecto (Scripts)

### Procesamiento Core (Bash)
- `scripts/bash/modkit_pileup.sh`: Ejecución de pileup para detección dual de 5mC y 5hmC.
- `scripts/bash/filter_coverage.sh`: Filtrado técnico por profundidad y probabilidad.
- `scripts/bash/prepare_dss_input.sh`: Conversión de datos al formato `(chr, pos, N, X)`.
- `scripts/bash/run_pipeline_core.sh`: Orquestador principal del flujo de trabajo.

### Análisis y Visualización (R)
- `scripts/R/dss_analysis.R`: Detección estadística de DML y DMR.
- `scripts/R/genomic_context.R`: Anotación funcional contra GFF3 (Exones, Intrones e Intergénico).
- `scripts/R/dmr_visualization.R`: Generación de perfiles de metilación para genes candidatos (*dnmt1*, *pak*, *sin3a*).

## 4. Instrucciones de Uso
1. **Ambiente:** Recree el entorno con `mamba env create -f bioinfo.yml`.
2. **Configuración:** Actualice las rutas en `scripts/bash/config.sh`.
3. **Ejecución:** Ejecute `bash scripts/bash/run_pipeline_core.sh`.

## 5. Instituciones y Créditos
Investigación para el Máster en Bioinformática de la **Universidad Internacional de Valencia (VIU)** en colaboración con **Universidad Regional Amazónica IKIAM**.
- **Autor:** Angel Ojeda
- **Tutor:** Pablo Marín, PhD.

## 6. Licencia
Este proyecto se distribuye bajo la Licencia MIT.

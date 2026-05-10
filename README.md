# DNA Methylation Pipeline for *Heliconius erato* (v1.0)

This repository provides a complete workflow for the identification and characterization of 5mC and 5hmC modifications in *Heliconius erato* using Oxford Nanopore Technologies (ONT) long-read sequencing data.

## 1. Technical Specifications

### Bioinformatics Environment
* **Basecalling/Mapping:** High-accuracy modBAM files (Minimap2/Dorado).
* **Modification Calling:** `modkit v0.5.0`
* **Downstream Analysis:** `R v4.2+`
* **Key R Libraries:** `DSS` (Dispersion Shrinkage for Sequencing), `bsseq`, `ggplot2`.

### Processing Parameters
* **Modification Threshold:** 0.79 (Determined by the 10th percentile of the empirical probability distribution).
* **Coverage Filter:** Minimum 10x per CpG site across conditions.
* **Reference Genome:** *H. erato lativitta* v1 (Lepbase).
* **Statistical Model:** Beta-binomial distribution with shared dispersion estimation (`equal.disp = TRUE`) for non-replicated designs.

## 2. Pipeline Architecture

The pipeline is organized into modular scripts for reproducibility:

### Bash Modules (`scripts/bash/`)
1.  `config.sh`: Centralized environment variables and binary paths.
2.  `modkit_pileup.sh`: Extraction of methylation frequencies and coverage.
3.  `filter_coverage.sh`: Post-pileup filtering based on depth and probability.
4.  `prepare_dss_input.sh`: Format conversion for compatibility with the DSS R package.
5.  `run_pipeline_core.sh`: Orchestration script for the entire processing block.

### R Modules (`scripts/R/`)
1.  `dss_analysis.R`: Differential methylation testing (DML/DMR detection).
2.  `genomic_context.R`: Annotation of methylated sites against GFF3 features (Exons, Introns, Intergenic).
3.  `dmr_visualization.R`: Generation of high-resolution tracks and candidate gene figures.

## 3. Data Requirements

To execute the pipeline, the following directory structure is required:
* `datos_secuenciacion_ori/`: Input modBAM files (`.bam` and `.bai`).
* `reference_genome/`: FASTA sequences and GFF3 annotation.
* `analisis/`: (Auto-generated) Output directory for processed data and figures.

## 4. Execution

1.  Clone the repository and ensure all `.sh` files have execution permissions.
2.  Update the absolute paths in `scripts/bash/config.sh`.
3.  Run the core processing:
    ```bash
    bash scripts/bash/run_pipeline_core.sh
    ```

## 5. Acknowledgments
This research was conducted as part of the Master’s Program in Bioinformatics at **Universidad Internacional de Valencia (VIU)**, in collaboration with **Universidad Regional Amazónica IKIAM**. 

**Project Supervision:** Pablo Marín, PhD.

## 6. License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

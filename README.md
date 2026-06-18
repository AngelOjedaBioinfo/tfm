# ONT Methylation Pipeline for *Heliconius erato lativitta*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0008--7799--9525-green)](https://orcid.org/0009-0008-7799-9525)
[![DOI](https://zenodo.org/badge/1234355389.svg)](https://zenodo.org/badge/latestdoi/1234355389)

Reproducible workflow for simultaneous detection of 5-methylcytosine (5mC) and
5-hydroxymethylcytosine (5hmC) in insect brain tissue using Oxford Nanopore
Technology (ONT). The pipeline processes modBAM files from Dorado duplex basecalling
through empirical threshold calibration, coverage filtering, and beta-binomial DMR
detection via DSS.

Designed for non-model Lepidoptera with no prior methylome characterization.
Validated on *Heliconius erato lativitta* central brain tissue (20 individuals,
~44.5M CpG sites, 33–52x genome-wide coverage).

---

## Technical Specifications

| Tool | Version | Role |
|---|---|---|
| Dorado | duplex mode | Basecalling and modification detection |
| modkit | 0.5.0 | modBAM pileup and empirical threshold calibration |
| minimap2 | 2.24 | Alignment to reference genome |
| samtools | 1.16.1 | BAM file handling |
| bedtools | 2.30.0 | Genomic arithmetic and intersect operations |
| DSS | 2.58.0 (Bioconductor) | Beta-binomial DMR detection |
| R | 4.2.2 | Statistical analysis and visualization |

---

## Pipeline Overview

Full workflow: Dorado duplex basecalling → modkit pileup → coverage filtering →
bedtools intersect → DSS → CAT GFF3 annotation → candidate genes.

```mermaid
---
config:
  theme: base
  themeVariables:
    primaryColor: '#E8F4F8'
    primaryBorderColor: '#2C7DA0'
    primaryTextColor: '#1B3A4B'
    lineColor: '#5A7A8A'
    fontSize: 12px
---
flowchart LR
    subgraph P1 ["I: Raw Data & Quality"]
        direction TB
        A["<b>modBAM</b><br/>Dorado duplex<br/>+ minimap2"]:::input
        B["<b>QC & Threshold</b><br/>modkit probs<br/>p10 = 0.793"]
        A --> B
    end

    subgraph P2 ["II: Processing (Bash)"]
        direction TB
        D["<b>modkit pileup</b><br/>~44.5M CpG"]
        E["<b>Filters & Intersect</b><br/>≥10× | ~12M Common"]
        D --> E
    end

    subgraph P3 ["III: Statistics (R/DSS)"]
        direction TB
        H["<b>DSS v2.58.0</b><br/>Wald Test + DMRs"]
        I["<b>Final DMRs</b><br/>6,197 5mC | 205 5hmC"]:::result
        H --> I
    end

    subgraph P4 ["IV: Output"]
        direction TB
        J["<b>Annotation</b><br/>GFF3 + BLASTx"]
        L["<b>Figures</b><br/>ggplot2 + patchwork"]:::output
        J --> L
    end

    P1 --> P2
    P2 --> P3
    P3 --> P4

    classDef input fill:#D1E8F2,stroke:#2C7DA0,stroke-width:2px,color:#1B3A4B
    classDef result fill:#FFE0B2,stroke:#E65100,stroke-width:2px,color:#3E2723
    classDef output fill:#C8E6C9,stroke:#2E7D32,stroke-width:2px,color:#1B5E20
```

---

## Methodological Decisions

Each non-default analytical choice is documented with a decision diagram to ensure
transparency and reproducibility.

### Probability Threshold Calibration

The modkit default threshold (0.5) was not applied directly. An empirical calibration
was performed using `modkit sample-probs` (10,042 reads sampled per file). The 10th
percentile of the resulting probability distribution (p10 = 0.793) was used as a
symmetric threshold for both modifications, yielding global rates consistent with
reported lepidopteran baselines.

```mermaid
---
config:
  theme: base
  themeVariables:
    primaryColor: '#E8F4F8'
    primaryBorderColor: '#2C7DA0'
    primaryTextColor: '#1B3A4B'
    lineColor: '#5A7A8A'
    fontSize: 14px
---
flowchart TB
 subgraph Calibration["NO — Calibration"]
    direction LR
        D["<b>Empirical calibration</b><br>modkit sample-probs<br>10,042 reads sampled<br>per file"]
        E["<b>Probability distribution</b><br>5mC: peak ~1.0<br>5hmC: peak ~0.5"]
        F["<b>10th percentile</b><br>= 0.793 in both samples"]
        G["<b>Symmetric application</b><br>--mod-thresholds m:0.79<br>--mod-thresholds h:0.79"]
  end
 subgraph Direct["YES — Compatible"]
    direction LR
        H["<b>Global rates obtained</b><br>5mC: 0.86–0.89%<br>5hmC: 0.21–0.23%"]
        I["✓ Compatible with<br>lepidopteran literature"]
  end
    A["<b>modBAM ONT data</b><br>5mC + 5hmC with MM/ML tags"] --> B["<b>Default threshold</b><br>modkit threshold = 0.5"]
    B --> C{"Global rates<br>compatible with<br>literature?"}
    C -- NO --> D
    C -- YES --> H
    D --> E
    E --> F
    F --> G
    G --> H
    H --> I

     D:::process
     F:::decision
     G:::process
     H:::result
     I:::output
     A:::input
    classDef input fill:#D1E8F2,stroke:#2C7DA0,stroke-width:2px,color:#1B3A4B
    classDef process fill:#FFF3E0,stroke:#E65100,stroke-width:2px,color:#3E2723
    classDef decision fill:#FFE082,stroke:#F57F17,stroke-width:2.5px,color:#3E2723
    classDef result fill:#FFE0B2,stroke:#E65100,stroke-width:2px,color:#3E2723
    classDef output fill:#C8E6C9,stroke:#2E7D32,stroke-width:2px,color:#1B5E20
```

### Dispersion Estimation in DSS

Three options were evaluated for dispersion parameter estimation under a pooled
design (n = 1 per condition). `smoothing = FALSE` was rejected by DSS at runtime.
`smoothing = TRUE` (Park & Wu, 2016) was selected based on software recommendation
and is the established approach for designs without biological replicates. A
sensitivity analysis with `equal.disp = TRUE` is documented in the thesis (§4.9).

```mermaid
---
config:
  theme: base
  themeVariables:
    primaryColor: '#E8F4F8'
    primaryBorderColor: '#2C7DA0'
    primaryTextColor: '#1B3A4B'
    lineColor: '#5A7A8A'
    fontSize: 14px
---
flowchart TD
    A["<b>Experimental design</b><br/>n=1 per condition<br/>(10 pooled brains per group)"]:::input
    A --> B{"How to estimate<br/>dispersion parameter<br/>phi in DSS?"}

    B --> C["<b>Option 1</b><br/>smoothing=FALSE"]
    C --> C2["❌ Rejected by DSS<br/><i>'There is no biological replicates.<br/>Please set smoothing=TRUE<br/>or equal.disp=TRUE'</i>"]:::rejected

    B --> D["<b>Option 2</b><br/>smoothing=TRUE<br/>(local kernel, 500 bp default)"]:::process
    D --> D2["6,197 DMRs 5mC<br/>205 DMRs 5hmC"]:::result

    B --> E["<b>Option 3</b><br/>equal.disp=TRUE<br/>(shared variance)"]:::process
    E --> E2["62 DMRs 5mC<br/>0 DMRs 5hmC"]:::sensitivity

    D2 --> F["<b>Decision</b><br/>Adopt smoothing=TRUE<br/>(Park & Wu 2016 recommendation<br/>for designs without replicates)"]:::output
    E2 --> G["<b>Sensitivity analysis</b><br/>reported in thesis §4.9<br/>absolute DMR count<br/>depends on method"]:::warning

    classDef input fill:#D1E8F2,stroke:#2C7DA0,stroke-width:2px,color:#1B3A4B
    classDef process fill:#FFF3E0,stroke:#E65100,stroke-width:2px,color:#3E2723
    classDef result fill:#FFE0B2,stroke:#E65100,stroke-width:2px,color:#3E2723
    classDef sensitivity fill:#F3E5F5,stroke:#6A1B9A,stroke-width:2px,color:#311B92
    classDef rejected fill:#FFCDD2,stroke:#C62828,stroke-width:2px,color:#B71C1C
    classDef warning fill:#FFE082,stroke:#F57F17,stroke-width:2px,color:#3E2723
    classDef output fill:#C8E6C9,stroke:#2E7D32,stroke-width:2px,color:#1B5E20
```

### Genomic Context Classification

CpG sites were classified through hierarchical intersections with the
`Hlat.v1.1.CAT.gff3` annotation file. Introns are not explicitly annotated in GFF3
and are identified by exclusion (gene body regions minus annotated exons).

```mermaid
---
config:
  theme: base
  themeVariables:
    primaryColor: '#E8F4F8'
    primaryBorderColor: '#2C7DA0'
    primaryTextColor: '#1B3A4B'
    lineColor: '#5A7A8A'
    fontSize: 14px
---
flowchart TD
    A["<b>CpG site</b><br/>BEDMethyl filtered cov ≥10×"]:::input
    A --> B{"Overlaps an<br/>'exon' feature<br/>in GFF3?"}
    B -->|"Yes<br/>bedtools intersect -u"| C["<b>Exon</b><br/>includes CDS and UTR"]:::exon
    B -->|"No"| D{"Overlaps a<br/>'gene' feature<br/>in GFF3?"}
    D -->|"Yes<br/>bedtools intersect -u"| E["<b>Intron</b><br/>= within gene body<br/>but outside exon"]:::intron
    D -->|"No"| F["<b>Intergenic</b><br/>bedtools intersect -v<br/>vs genes.gff3"]:::inter

    G["<i>Note: introns are not explicitly annotated<br/>in GFF3; identified by exclusion</i>"]:::note

    classDef input fill:#D1E8F2,stroke:#2C7DA0,stroke-width:2px,color:#1B3A4B
    classDef exon fill:#FFB74D,stroke:#E65100,stroke-width:2px,color:#3E2723
    classDef intron fill:#BCAAA4,stroke:#5D4037,stroke-width:2px,color:#3E2723
    classDef inter fill:#E0E0E0,stroke:#616161,stroke-width:2px,color:#212121
    classDef note fill:none,stroke:none,color:#757575
```

---

## Installation & Usage

**1. Recreate the environment**

```bash
mamba env create -f bioinfo.yml
conda activate bioinfo
```

**2. Configure paths**

Edit `scripts/bash/config.sh` with your local paths:

```bash
# Required inputs
GENOME="/path/to/Hlat.v1.1.fasta"
MODBAM_DIR="/path/to/modbam/files/"
GFF3="/path/to/Hlat.v1.1.CAT.gff3"
OUTDIR="/path/to/output/"
```

**3. Run the pipeline**

```bash
bash scripts/bash/run_pipeline_core.sh
```

The pipeline runs in four sequential stages. Intermediate files are written to
`$OUTDIR` at each stage. The final DMR tables and annotated gene lists are generated
by the R scripts in `scripts/R/`.

**Computational requirements**

Tested on Ubuntu 22.04 with 32 GB RAM and 16 cores. The modkit pileup step is the
most resource-intensive (~6 h per sample at 40x coverage). DSS analysis runs in R
with ~8 GB RAM for ~12M sites.

---

## Repository Structure

```
tfm/
├── scripts/
│   ├── bash/
│   │   ├── config.sh               # Path configuration (edit before running)
│   │   ├── modkit_pileup.sh        # Dual 5mC/5hmC pileup from modBAM
│   │   ├── filter_coverage.sh      # Depth and probability filtering
│   │   ├── prepare_dss_input.sh    # Format conversion to (chr, pos, N, X)
│   │   └── run_pipeline_core.sh    # Main pipeline orchestrator
│   └── R/
│       ├── dss_analysis.R          # DML and DMR detection (DSS)
│       ├── genomic_context.R       # Functional annotation and classification
│       └── dmr_visualization.R     # Methylation profiles for candidate genes
├── assets/
│   └── figures/                    # Decision diagrams (source .mmd files)
├── bioinfo.yml                     # Conda environment specification
├── CITATION.cff                    # Machine-readable citation metadata
└── LICENSE                         # MIT
```

---

## Case Study: *H. erato lativitta* + RG108

The pipeline was validated on the following experimental system.

**Experimental design**

| Parameter | Detail |
|---|---|
| Species | *Heliconius erato lativitta* |
| Individuals | 20 total: 10 control, 10 experimental (5 males + 5 females each) |
| Source | Ikiam insectary, Tena, Ecuador |
| Treatment | RG108: 2 mM topical (2x/day) + 3 mM dietary, 7 days |
| Control | Vehicle: PBS + 1% DMSO |
| Tissue | Central brain (optical lobes excluded) |
| Sequencing | Oxford Nanopore Technology (ONT), duplex mode |
| Reference genome | *H. erato lativitta* v1.1 (CAT annotation) |

**Results summary**

| Modification | Baseline sites | Baseline rate | DMRs under RG108 | Change |
|---|---|---|---|---|
| 5mC | 8,431 | 0.228% | 6,197 | -6.48% |
| 5hmC | 32,843 | 0.888% | 205 | -3.26% |

Methylation is concentrated in exons (~17-fold enrichment over intergenic regions).
Top candidate genes: *pak*, *sin3a*, *drp1* (synaptic plasticity and epigenetic
regulatory machinery). This constitutes the first simultaneous 5mC and 5hmC
characterization in *Heliconius* brain tissue.

**Associated manuscript**

Ojeda A, Marín P, Bacquet C (in preparation). Brain methylome of *Heliconius erato
lativitta* via Oxford Nanopore Technology: first simultaneous 5mC and 5hmC atlas.

---

## Data Availability

Raw sequencing data will be deposited in the European Nucleotide Archive (ENA) upon
manuscript submission. The reference genome (*H. erato lativitta* v1.1) is available
from [Lepbase](http://lepbase.org).

---

## Citation

If you use this pipeline, please cite this repository:

```
Ojeda A (2026). ONT Methylation Pipeline for Heliconius erato lativitta.
GitHub: https://github.com/AngelOjedaBioinfo/tfm
DOI: [Zenodo DOI — assigned at v1.0.0 release]
```

A `CITATION.cff` file is included for automated citation tools (GitHub, Zotero,
Mendeley).

---

## Authors & Acknowledgments

**Ángel Andrés Ojeda Montesdeoca**
[ORCID 0009-0008-7799-9525](https://orcid.org/0009-0008-7799-9525) |
[GitHub @AngelOjedaBioinfo](https://github.com/AngelOjedaBioinfo)
Laboratorio de Biología Molecular de Docencia (LBMD),
Universidad Regional Amazónica Ikiam, Tena, Ecuador.

**Supervisors:**
Pablo Marín, PhD (Universidad Internacional de Valencia, VIU) |
Caroline Bacquet, PhD (Ikiam / Jiggins Group, Cambridge–Sanger Institute)

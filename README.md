Uncovering the temporal dynamics and regulatory networks of thermal
stress response in a hyperthermophile using transcriptomics and
proteomics
================


This repository contains bioinformatic analysis of the manuscript titled
**Uncovering the temporal dynamics and regulatory networks of thermal
stress response in a hyperthermophile using transcriptomics and
proteomics**.

You can find the preprint here:
<https://doi.org/10.1101/2023.05.02.539125>

This document provides a general workflow and overview of the tools we
have for analysis, including:  
- Differential gene expression using `DESeq2`  
- Differential protein expression using `limma`  
- Identification of transcript 3’ ends - Analysis of operons using
Nanopore PCR-cDNA sequencing  
- Downstream analysis, including arCOG enrichment & cluster analysis

------------------------------------------------------------------------

#### Table of Contents

- <a href="#repository-organization"
  id="toc-repository-organization">Repository organization</a>
- <a href="#software-versions" id="toc-software-versions">Software
  versions</a>
- <a href="#citation" id="toc-citation">Citation</a>
- <a href="#license" id="toc-license">License</a>

------------------------------------------------------------------------

### Repository organization

This repository is organized as follows:

- `raw_data/`: Contains raw data files for transcriptomics, proteomics,
  and Nanopore sequencing.
- `processed_data/`: Stores processed data and intermediate files
  generated during the analysis.
- `scripts/`: Contains the R code and other scripts used in the
  analysis.
- `documentation/`: Includes detailed documentation, such as method
  descriptions and supplementary material.
- `results/`: Stores final results, plots, and tables.

### Software versions

### Citation

If you find this work useful, please cite:

### License

This project is licensed under the \[Name of License\]. See the
`LICENSE` file for more details.

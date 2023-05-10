Uncovering the temporal dynamics and regulatory networks of thermal
stress response in a hyperthermophile using transcriptomics and
proteomics
================
<a href="https://orcid.org/0000-0001-7444-2408">Felix
Grünberger<sup>a#</sup></a>, Georg Schmid<sup>a</sup></a>, Zubeir El
Ahmad<sup>a</sup></a>, Martin Fenk<sup>a\*</sup></a>, Katharina
Vogl<sup>a</sup></a>, Robert Reichelt<sup>a</sup></a>, Winfried
Hausner<sup>a</sup></a>, Henning Urlaub<sup>b,c</sup></a>, Christof
Lenz<sup>b,c</sup></a>, Dina Grohmann<sup>a#</sup></a>  

<sup>a</sup> Institute of Biochemistry, Genetics and Microbiology,
Institute of Microbiology and Archaea Centre, Single-Molecule
Biochemistry Lab and Regensburg Center for Biochemistry, University of
Regensburg, Regensburg, Germany

<sup>b</sup> Bioanalytical Mass Spectrometry Group, Max Planck Institute
for Biophysical Chemistry, Göttingen, Germany

<sup>c</sup> Department of Clinical Chemistry, University Medical Center
Göttingen, Göttingen, Germany

<sup>\#</sup> Corresponding authors

<sup>\*</sup> Present address: Max Planck Institute for Infection
Biology, 10117 Berlin, Germany

------------------------------------------------------------------------

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
- <a href="#data-availability" id="toc-data-availability">Data
  availability</a>
- <a href="#funding" id="toc-funding">Funding</a>
- <a href="#citation" id="toc-citation">Citation</a>
- <a href="#license" id="toc-license">License</a>

------------------------------------------------------------------------

The repository is currently actively developed.

[![Active
Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)

### Repository organization

This repository is organized as follows:

- [`Rscripts/`](Rscripts): R code used for the analysis. Description in
  `documentation`
- [`data/`](data): Contains intermediate data stored in `Rdata` objects,
  and additional sets used during analysis (genome, TSS, TTS)
- [`documentation/`](documentation): Includes detailed documentation,
  such as method descriptions and description of R code  
- [`supplemental tables/`](supplemental%20tables): Contains tables S1-S4
  used for downstream differential gene and protein expression analysis,
  3’ end analysis and final summary table

### Data availability

RNA sequencing data are available at the European Nucleotide Archive
(ENA, <https://www.ebi.ac.uk/ena>) under project accession numbers
[PRJEB61174](https://www.ebi.ac.uk/ena/browser/view/PRJEB61174) (RNA-seq
data used for differential gene expression analysis) and
[PRJEB61177](https://www.ebi.ac.uk/ena/browser/view/PRJEB61177)
(Term-seq and Nanopore data).  
The mass spectrometry proteomics data have been deposited to the
ProteomeXchange Consortium via the PRIDE partner repository with the
dataset identifier PXD041262.

### Funding

Work in the Grohmann lab was supported by the Deutsche
Forschungsgemeinschaft (DFG funding schemes SPP2141 “Much more than
defence: the multiple functions and facets of CRISPR-Cas” and SFB960
TPA7 to D.G.). H.U. received funding via the DFG SPP2141 “Much more than
defence: the multiple functions and facets of CRISPR-Cas”.

### Citation

If you find this work useful, please cite:

Uncovering the temporal dynamics and regulatory networks of thermal
stress response in a hyperthermophile using transcriptomics and
proteomics

Felix Gruenberger, Georg Schmid, Zubeir El-Ahmad, Martin Fenk, Katharina
Vogl, Robert Reichelt, Winfried Hausner, Henning Urlaub, Christof Lenz,
Dina Grohmann

bioRxiv 2023.05.02.539125; doi:
<https://doi.org/10.1101/2023.05.02.539125>

### License

This project is licensed under the \[Name of License\]. See the
`LICENSE` file for more details.

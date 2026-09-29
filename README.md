# Transcriptomic analysis of male cone development in *Ginkgo biloba* and *Cycas circinalis*

This repository contains the R scripts associated with the manuscript:

**“Male cone development in *Ginkgo biloba* reveals recurrent reproductive modules shared with *Cycas circinalis*”**

The scripts reproduce the transcriptomic analyses used to investigate male cone development in *Ginkgo biloba* and to perform the comparative analysis with *Cycas circinalis*.

## Scripts

### `Supplementary_Script_S1_Ginkgo_publication_ready_minimal.R`

R workflow for the transcriptomic analysis of *Ginkgo biloba* male cone development.

The script includes:

- data import and low-count filtering;
- DESeq2 normalization and differential expression analysis;
- quality control and PCA;
- Mfuzz soft clustering of temporal expression profiles;
- Gene Ontology (GO) mapping and enrichment analyses;
- reproductive candidate-gene identification;
- generation of selected heatmaps and figures used in the manuscript;
- export of summary statistics and R session information.

Required input files:

- `count_RNAseq_MaleCones.xlsx`
- `annotazione_geni.xlsx`

### `Supplementary_Script_S2_Cycas_circinalis_publication_ready.R`

R workflow for the targeted comparative transcriptomic analysis of *Cycas circinalis* male cone development.

The script includes:

- data import and organization;
- DESeq2 normalization and differential expression analysis;
- quality control and PCA;
- UniProt-based homolog annotation and GO retrieval;
- targeted GO-category enrichment;
- reproductive candidate-gene identification;
- generation of a Ginkgo-informed representative reproductive candidate heatmap;
- export of summary statistics and R session information.

Mfuzz clustering is intentionally not applied to the *Cycas* dataset because the sampled developmental stages are not strict one-to-one equivalents of the *Ginkgo* M1–M3 developmental series.

Required input file:

- `Cycas_circinalis_Novogene_counts_annotations.csv`

When available, the *Cycas* workflow also reuses the expanded GO-category list generated during the *Ginkgo* analysis (`GO_categories_expanded_list.rds`) to maintain comparability between the two species.

## R packages

The analyses use R and packages including:

`DESeq2`, `Mfuzz`, `Biobase`, `ggplot2`, `dplyr`, `pheatmap`, `readr`, `readxl`, `stringr`, `tidyr`, `tibble`, `RColorBrewer`, `ontologyIndex`, `openxlsx`, `reshape2`, and related dependencies.

Additional package requirements are specified directly within each script.

The workflows retrieve updated annotation information from UniProt and Gene Ontology resources where required.

## Data availability

The input transcriptomic datasets are not included in this repository. Information on the availability and accession of the corresponding datasets is provided in the Data Availability Statement of the associated manuscript.

## Usage

Place the required input files in the R working directory and run the corresponding script.

Both workflows generate their analysis outputs automatically. R session information is also exported to facilitate reproducibility.

## Authors

Elisabetta Offer and co-authors  
Department of Biology, University of Padova  
Botanical Garden of Padova

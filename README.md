# AnalyteExplorer

<!-- badges: start -->
[![R-CMD-check](https://github.com/RGLab/AnalyteExplorer/workflows/R-CMD-check/badge.svg)](https://github.com/RGLab/AnalyteExplorer/actions)
<!-- badges: end -->

The goal of AnalyteExplorer is to pre-process data for the `AnalyteExplorer` module in [ImmuneSpace](https://www.immunespace.org/).

## Installation

You can install the development version of AnalyteExplorer from [GitHub](https://github.com/RGLab/AnalyteExplorer) with:

``` r
# install.packages("remotes")
remotes::install_github("RGLab/AnalyteExplorer")
```

## Workflow with ImmuneSpace

``` r
library(AnalyteExplorer)
library(UpdateAnno)
options(debug_dir = tempdir())
labkey.url.base <- "https://www.immunespace.org"
labkey.url.path <- "/AnalyteExplorer"

btm <- process_blood_transcription_modules()
check_table(btm)
res <- update_table("blood_transcript_modules", btm)

ge <- process_gene_expression()
check_table(ge)
res <- update_table("gene_expression", ge)
```

## Data Processing

### `process_blood_transcription_modules()`

- Based on [Molecular signatures of antibody responses derived from a systems biological study of 5 human vaccines](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3946932/)
- The function updates the gene symbols with Hugo.
- This table is used to display metadata about the selected module in the app.

|id   |name                                           |genes                                                                                                                                                                                                                                                                                                                      |matched_gene_ontology_terms                                                                      | number_of_genes|module_category    |
|:----|:----------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------|---------------:|:------------------|
|M0   |targets of FOSL1/2 (M0)                        |CCL2, COL1A2, DCN, IL6, CXCL8, LIF, MGP, MMP1, MMP2, MMP9, PLAU, THBD                                                                                                                                                                                                                                                      |extracellular space (11), extracellular region (11), protein binding (9)                         |              12|TF targets         |
|M1.0 |integrin cell surface interactions (I) (M1.0)  |COL1A1, COL1A2, COL5A1, DPYSL3, MYH10, NRP1, PTK2, RHOC, RRAS, SEMA6A                                                                                                                                                                                                                                                      |protein binding (26), axon guidance (20), extracellular region (16)                              |              29|molecular function |
|M1.1 |integrin cell surface interactions (II) (M1.1) |AHSP, ALAD, ALAS2, CPOX, E2F2, FECH, GATA1, HEMGN, HMBS, PLEK2, TMOD1                                                                                                                                                                                                                                                      |protein binding (12), extracellular region (12), extracellular matrix structural constituent (9) |              12|molecular function |
|M2.0 |extracellular matrix (I) (M2.0)                |CD1D, HLA-DMA, HLA-DMB, HLA-DPA1, HLA-DPB1, HLA-DQA2, METTL7A, WDFY4                                                                                                                                                                                                                                                       |protein binding (27), extracellular region (26), extracellular matrix (21)                       |              30|location           |


### `process_gene_expression()`

- This function creates a gene expression table by cohort, timepoint, and analyte type (gene, blood transcript module, or gene signature).
- Processing steps:
  1. Fetch gene expression matrices and metadata from ImmuneSpace and [combine](https://github.com/RGLab/ImmuneSpaceR/blob/main/R/ISCon-geneExpression.R#L913) them into one `ExpressionSet` object.
      - Remove genes that are not available in all expression matrices.
  1. Remove samples that have negative timepoint and select one timepoint if sample has multiple baseline timepoints.
  1. Create gene expression table summarized by analyte type.
      1. In sample level, when summarizing by blood transcript module or gene signature, compute geometric mean of the expression values of the genes in the module or signature
      2. In cohort level, compute the fold change of the expression values for all combinations of timepoints comparing to the baseline timepoint.
      3. Compute mean and standard deviation of those fold change values by analyte type
  1. Merge the three summarized tables

|cohort          |sample_type |study_accession |condition    | timepoint|analyte_id |analyte_type | mean_fold_change| sd_fold_change| id|
|:---------------|:-----------|:---------------|:------------|---------:|:----------|:------------|----------------:|--------------:|--:|
|healthy aldults |Whole blood |SDY1529         |Yellow_Fever |         0|A1CF       |gene         |                0|              0|  1|
|healthy aldults |Whole blood |SDY1529         |Yellow_Fever |         0|A2M        |gene         |                0|              0|  2|
|healthy aldults |Whole blood |SDY1529         |Yellow_Fever |         3|M0         |blood transcription module |       -0.0385090|      0.0875805| 3894315|
|healthy aldults |Whole blood |SDY1529         |Yellow_Fever |         3|M1.0       |blood transcription module |       -0.0181510|      0.1051648| 3894316|
|healthy aldults |Whole blood |SDY1529         |Yellow_Fever |         3|21357945_1_8  |gene signature |        0.7102219|      0.5402085| 4028003|
|healthy aldults |Whole blood |SDY1529         |Yellow_Fever |         3|21357945_2_9  |gene signature |        0.1105955|      0.5026520| 4028004|

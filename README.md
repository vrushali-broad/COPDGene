# Comprehensive Analysis of COPDGene Data

This repository contains a Python script designed for the processing, analysis, and visualization of clinical phenotypic data related to pulmonary diseases. It integrates various data handling and preprocessing techniques, statistical analysis, and visualization to provide insights into the impact of diseases on patient phenotypes.

## Features

- **Data Loading and Preprocessing:** Load and preprocess metadata, gene expression counts, and phenotypic information.
- **Subject Classification:** Classify subjects by disease status and filter out non-smokers.
- **Data Visualization:** Generate plots and heatmaps to visualize data distributions, correlations, and disease impacts.
- **Confounder Correction:** Correct pulmonary function and CT scan features for confounders using regression analysis.
- **Statistical Analysis:** Perform statistical analysis to understand the significance of differences across disease states.

## Dependencies

- `pandas`
- `seaborn`
- `matplotlib`
- Custom modules: `visualization_helpers`, `data_processing`, `patient_classification`, `patient_statistics`, `analyze_symptoms`, `analyze_cbc`, `correct_confounders_pft`, `correct_confounders_ct`

## Installation

To run this script, you will need Python 3.x and the above-mentioned dependencies installed on your system. You can install the dependencies using pip:

```bash
pip install pandas seaborn matplotlib
```


## Differential Expression Analysis for COPD and ACO

This repository hosts a crucial upstream script within the `COPDGene cohort analysis` package, designed for researchers and bioinformaticians focusing on Chronic Obstructive Pulmonary Disease (COPD) and Asthma-COPD Overlap (ACO). Utilizing a robust set of R packages, this script facilitates comprehensive differential expression analysis, from data preprocessing and normalization to advanced statistical testing and enrichment analysis.

### Key Features

- **Comprehensive Data Preprocessing**: Prepares gene expression data for analysis, supporting gene annotation and filtering.
- **Advanced Statistical Testing**: Leverages `limma` and `edgeR` for robust differential expression analysis.
- **Batch Effect Correction**: Implements methods to adjust for batch effects, ensuring accurate biological interpretation.
- **Enrichment Analysis**: Conducts KEGG pathway and GO term enrichment analysis to identify biologically significant pathways.
- **Detailed Logging**: Features custom logging for enhanced traceability and reproducibility of analysis results.
- **Visualization Support**: Provides tools for the visualization of analysis outcomes, aiding in the interpretation of complex datasets.

### Prerequisites

Before running this script, ensure that you have R installed (version 3.6.0 or later recommended) along with the following packages: `limma`, `edgeR`, `org.Hs.eg.db`, `data.table`, and `biomaRt`.

#### Installation Instructions

To install the necessary R packages, run the following commands in your R environment:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "edgeR", "org.Hs.eg.db", "data.table", "biomaRt"))
```
#### Usage
Ensure your workflow is set up correctly with these steps:

#### Prepare Your Data
Ensure differential expression data is in CSV format with gene IDs and log fold changes.

#### Configure the Script
Adjust script parameters, including data file path and analysis settings, as needed.

#### Execute the Analysis
Run the script in R to conduct the analysis and generate visualizations.

##### Example

```r
# Example code demonstrating script configuration and execution

# Define the data file path and analysis parameters
PATH <- "path/to/your/data"
lfc_threshold <- 0.5  # Log fold change threshold

# Commands to load data and perform analysis
load_data(PATH)
# Additional analysis commands
```

## Enrichment Analysis Script for COPD and ACO

This repository hosts the `Enrichment Analysis Script`, a key component of the `COPDGene cohort analysis` package, designed for genomic researchers focusing on Chronic Obstructive Pulmonary Disease (COPD) and Asthma-COPD Overlap (ACO). This script leverages a suite of Bioconductor packages to perform comprehensive KEGG pathway and GO term enrichment analysis, visualization, and data logging for differential expression data.

### Features

- **KEGG Pathway Analysis**: Identifies significant pathways related to the disease condition using the `clusterProfiler` package.
- **GO Term Enrichment**: Analyzes gene ontology terms associated with differentially expressed genes.
- **Visualization**: Generates chord and tree diagrams for enriched pathways and terms, facilitating intuitive interpretation of results.
- **Data Preprocessing and Mapping**: Prepares data for analysis, including mapping between gene IDs.
- **Logging**: Detailed logging of actions and package versions for reproducibility.

### Prerequisites

Ensure you have R installed (version 3.6.0 or later recommended) along with the following packages: `limma`, `edgeR`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `GOplot`, `data.table`, and `dplyr`.

#### Installation Instructions

To install the necessary R packages, run the following commands in your R environment:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "edgeR", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "GOplot", "data.table", "dplyr"))
```
### Usage

#### Prepare Your Data
Ensure your differential expression data is in CSV format with gene IDs and log fold changes.

#### Configure the Script
Adjust the script's parameters, including the data file path and KEGG pathway IDs, to suit your analysis needs.

#### Execute the Analysis
Run the script within the parent package environment to perform the analysis and generate visualizations.

##### Example

```r
# Define color settings for visualization
colors_copd <- list(high = '#A6CDC7', low = '#e89fb4', mid = 'gray88')

# Perform KEGG pathway analysis for COPD
plot_kegg_pathway(data_file = 'path/to/your/data.csv', kegg_id = 'hsa04217', colors = colors_copd)
```

> **Note:** Make sure to replace placeholder text (e.g., `path/to/your/data.csv`, `hsa04217`, `[Your Organization Name]`) with the actual information related to your project. This format helps users to easily understand how to prepare their data, configure, and execute the analysis, along with guiding contributors on how to contribute to the project.



## KEGG Pathway Analysis Subscript

### Overview
This document provides details on the KEGG Pathway Analysis Subscript, a component of the larger `COPDGene cohort analysis` package designed for comprehensive analysis of genomic data. This subscript specifically focuses on performing KEGG pathway analysis and visualization, leveraging high-throughput differential expression data to elucidate the mechanisms underlying conditions such as COPD and ACO.

### Dependencies
This subscript requires R and several Bioconductor packages. Ensure you have the following installed:
- R (Version 3.6.0 or later recommended)
- Bioconductor packages:
  - AnnotationDbi
  - org.Hs.eg.db
  - pathview
  - gage
  - gageData

To install Bioconductor and the required packages, execute the following in R:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "pathview", "gage", "gageData"))
```

### Usage

This subscript is part of the `ACO Manuscript` package, focusing on KEGG pathway analysis. To use this script effectively:

#### Prepare your data
Ensure your differential expression data is in a CSV format with gene IDs and log fold changes.

#### Configure the script
Modify the script's parameters, including the data file path and KEGG pathway IDs, to match your analysis requirements.

#### Execute the analysis
Run the script from the parent package's environment to generate pathway plots.

##### Example usage

```r
# Configuring color settings for the analysis
colors_copd <- list(high = '#A6CDC7', low = '#e89fb4', mid = 'gray88')

# Performing KEGG pathway analysis
plot_kegg_pathway(data_file = 'path/to/your/data.csv', kegg_id = 'hsa04217', colors = colors_copd)
```

### Contributing

Contributions to improve this script and the larger package are welcome. Please fork the repository, make your changes, and submit a pull request.

### License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

### Acknowledgements

Special thanks to Channing Division of Network Medicine, Harvard Medical School, Boston, MA for supporting this project, and to all contributors for their invaluable input.


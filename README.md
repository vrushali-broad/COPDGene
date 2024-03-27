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

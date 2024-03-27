"""
This script encompasses a comprehensive workflow for the analysis of COPDGene data, integrating
various aspects of data handling, preprocessing, analysis, and visualization. It is designed to process 
and analyze clinical phenotypic data, including patient metadata, gene expression counts, and detailed 
phenotypic information related to pulmonary function tests (PFT) and CT scan features. Key functionalities include:

1. Loading and preprocessing of metadata, counts data, and phenotypic information.
2. Classification of subjects by disease status and filtering of non-smokers.
3. Visualization of data distributions and statistical analysis to understand disease impact on patient phenotypes.
4. Correction of pulmonary function and CT scan features for confounders using regression analysis.
5. Generation of plots and heatmaps to visualize correlations and differences in data across disease states.

Dependencies include pandas for data manipulation, seaborn and matplotlib for visualization, and custom modules for
specific analytical and visualization tasks. The script adheres to best practices in the field, with emphasis on reproducibility, readability,
and modular design. Logging is extensively used for tracking the analysis process and troubleshooting.

Author: Vrushali D. Fangal, 
        Channing Division of Network Medicine, 
        Harvard Medical School, 
        Boston, MA.
Date: March 22, 2024
License: MIT license
"""

import os
import logging
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

# Assuming these are your modular functions from different modules
from visualization_helpers import plot  
from data_processing import load_and_preprocess_meta, load_and_preprocess_counts, load_and_preprocess_pheno
from patient_classification import classify_subjects, plot_subject_distribution
from patient_statistics import generate_statistics_table

from analyze_symptoms import analyze_and_plot_symptoms, analyze_and_plot_hospitalizations
from analyze_cbc import analyze_and_plot_cbc, create_colorbar, create_correlation_heatmap, create_correlation_heatmap_per_disease
from correct_confounders_pft import perform_regression_analysis
from correct_confounders_ct import correct_ct_scan_features, correct_contrast_ct_scan_features

## Import feature names
from pheno_features import CT_FEATURES, CT_NAMES, CT_order_names, CT_rename
from pheno_features import PFT_FEATURES, PFT_NAMES, PFT_order_names, PFT_rename

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

def create_directory_if_not_exists(directory):
    """Create a directory if it does not exist."""
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
            logging.info(f"Created directory: {directory}")
    except Exception as e:
        logging.error(f"Error creating directory {directory}: {e}")
        raise
        
def plot_features(pheno, features, order, names, plot_type, dir_name, pvalues_lst=None, figsize = (3,3), palette = None, sig=False):
    """
    Plot features from the given phenotypic data.

    Args:
        pheno (DataFrame): The DataFrame containing phenotypic data.
        features (list): A list of features to plot.
        order (list): The order of categories to be used in the plot.
        names (list): Names for each feature to be used in the plot.
        plot_type (str): Type of plot to create (e.g., 'bar', 'box').
        dir_name (str): Directory to save the plots.
        pvalues_lst (dict, optional): Dictionary of p-values for statistical annotation.
        figsize (tuple, optional): Size of the figure. Defaults to (3, 3).
        palette (list, optional): Color palette for the plots. Defaults to a predefined palette.
        sig (bool, optional): Whether to include additional statistical annotations. Defaults to False.
    
    Raises:
        ValueError: If input lists have different lengths.
        
    This function uses the `plot` function from the `visualization_helpers` module 
    to create box plots for each feature in `features`.
    """
    # Set default palette if not provided
    if palette is None:
        palette = ['#D3D3D3', '#AFE4DE', '#FAA07A', '#FCFBCB']
        
    # Plot each feature
    for feature in features:
        try:
            # Extract p-values for the current feature, if available
            feature_pvalues = pvalues_lst.get(feature, None) if pvalues_lst else None
            
            plot(pheno=pheno,
                 features=[feature],  # Wrapping feature in a list for compatibility
                 order=order,
                 names=[names[features.index(feature)]],  # Corresponding name for the feature
                 plot_type=plot_type,
                 dir_name=dir_name,
                 palette=palette,
                 pvalues=feature_pvalues,  # Use p-values if available
                 figsize=figsize,
                 sig=sig)  # Pass the sig argument
            logging.info(f"Plot generated for feature: {feature}")
        except Exception as e:
            logging.error(f"Error occurred while plotting feature {feature}: {e}")
            raise
            
def load_and_prepare_data(data_dir, processed_data_dir, figures_dir):
    """
    Load and prepare the base data for analysis.

    Args:
        data_dir (str): Directory containing the raw data files.
        processed_data_dir (str): Directory to save processed data.
        figures_dir (str): Directory to save figures.

    This function handles the initial data loading, including reading in meta, counts, 
    and phenotypic data. It also ensures necessary directories are created.
    """
    try:
        logging.info("#######################################################################")
        logging.info("################# Starting data loading and preparation... ############")
        logging.info("#######################################################################")

        # Paths for specific datasets and output files
        meta_file_path = os.path.join(data_dir, 'master.file.freeze4.txt')
        counts_file_path = os.path.join(data_dir, 'counts_raw_from_lengthScaledTPM.tsv')
        pheno_file_path = os.path.join(data_dir, 'COPDGene_P1P2P3_Flat_SM_NS_Nov21_PhenoData.txt')

        # Processing and analysis steps
        logging.info("Loading and preprocessing meta data...")
        meta = load_and_preprocess_meta(meta_file_path)

        logging.info("Loading and preprocessing counts data...")
        counts = load_and_preprocess_counts(counts_file_path)

        logging.info("Loading and preprocessing pheno data...")
        pheno = load_and_preprocess_pheno(pheno_file_path, counts.columns)

        logging.info("Classifying subjects by disease...")
        pheno = classify_subjects(pheno)
        
        # Remove smokers before processing
        logging.info("Removing non-smokers and filtering data...")
        pheno = pheno[pheno.Disease.notna() & (pheno.Disease != 'Non-smoker')]
        logging.info(f"Filtered pheno data of smokers: {pheno.shape}")

        # Updating subjects, counts, and meta based on filtered pheno data
        subjects = pheno.index
        logging.info(f"Number of subjects after filtering: {len(subjects)}")

        counts = counts[subjects]
        logging.info(f"Updated counts data shape: {counts.shape} (Number of genes: {counts.shape[0]})")

        meta = meta.loc[subjects]
        logging.info(f"Updated meta data shape: {meta.shape}")
        
        # Log total number of subjects and disease distribution
        # logging.info(f"Total subjects: {len(pheno)}")
        # disease_distribution = pheno['Disease'].value_counts().to_string()
        # logging.info(f"Disease distribution:\n{disease_distribution}")
        
        # Saving processed data
        logging.info("Saving processed data...")
        counts_file_save_path = os.path.join(processed_data_dir, 'counts.csv')
        meta_file_save_path = os.path.join(processed_data_dir, 'meta.tsv')
        pheno_file_save_path = os.path.join(processed_data_dir, 'pheno.tsv')

        counts.to_csv(counts_file_save_path, sep='\t')
        meta.to_csv(meta_file_save_path, sep='\t')
        pheno.merge(meta, how='outer', left_index=True, right_index=True).to_csv(pheno_file_save_path, sep='\t')

        logging.info("Data loading and preparation completed successfully.")
        return meta, counts, pheno
    
    except Exception as e:
        logging.error(f"Error during data loading and preparation: {e}")
        raise

def analyze_disease_distribution(pheno, figures_dir):
    """
    Analyze disease distribution, generate patient distribution plot, and create a statistics table.

    Args:
        pheno (pd.DataFrame): DataFrame containing phenotypic data with a 'Disease' column.
        patient_dist_plot_path (str): File path to save the patient distribution plot.
        processed_data_dir (str): Directory path for saving the statistics table.

    This function logs the total number of subjects and their disease distribution, generates a plot 
    for patient distribution, and creates a statistics table.
    """
    try:
        # Log total number of subjects and disease distribution
        logging.info(f"Total subjects: {len(pheno)}")
        disease_distribution = pheno['Disease'].value_counts().to_string()
        logging.info(f"Disease distribution:\n{disease_distribution}")

        # Generate and save patient distribution plot
        logging.info("Generating patient distribution plot...")
        patient_dist_plot_path = os.path.join(figures_dir, 'patient_distribution.pdf')
        plot_subject_distribution(pheno, patient_dist_plot_path)

        # Generate and save the statistics table
        logging.info("Generating statistics table...")
        stats_table_path = os.path.join(figures_dir, 'statistics_table.tsv')
        generate_statistics_table(pheno, stats_table_path)

        logging.info("Disease distribution analysis completed.")

    except Exception as e:
        logging.error(f"Error in analyzing disease distribution: {e}")
        raise

def correct_and_plot_pft(pheno, features, names, disease_order, order_names, rename, dir_name, formula):
    """
    Corrects pulmonary function test (PFT) features for confounders using regression models and plots the adjusted values.

    Args:
        pheno (pd.DataFrame): DataFrame containing the phenotypic data.
        features (list): List of lung function features to correct.
        names (list): List of names for each feature for plotting.
        order_names (list): List of confounders in the model, in order.
        rename (dict): Dictionary to rename confounders in the plot.
        dir_name (str): Directory to save the plots.
        formula_template (str): Formula template for OLS regression.

    Returns:
        pd.DataFrame, dict: Adjusted phenotypic data and p-values list.
    """
    logging.info("#################################################################################")
    logging.info("### Correcting PFT features for confounders age, gender, race, and height... ####")
    logging.info("#################################################################################")

    ## Correct PFT for confounders
    adj_pheno, pvalues_lst = perform_regression_analysis(pheno, features, names, 
                                                         order_names, rename, formula,
                                                         dir_name = dir_name)
        
    ## Plot adjusted values
    adj_features = ['Adjusted_'+col for col in features]

    # Plotting the lung function data
    plot_features(pheno = adj_pheno, features = adj_features, 
                  order = disease_order, names = names, plot_type = 'box', 
                  dir_name = dir_name, 
                  # pvalues_lst = pvalues_lst, 
                  figsize = (3,3),
                  sig = False
                 )
    
    logging.info("PFT features correction and plotting completed.")
    return adj_pheno, pvalues_lst        
    
def correct_and_plot_ct_scan_features(pheno, outcomes, names, disease_order,
                                      order_names, rename, dir_path, formula,
                                      plot_type='box', figsize=(3, 3), sig=False):
    """
    Correct CT scan features for confounders and plot the adjusted features.

    Args:
        pheno (pd.DataFrame): DataFrame containing phenotypic data including CT scan values.
        outcomes (list): List of outcome variables (CT scan features) to adjust.
        names (list): Names corresponding to each outcome for plotting.
        dir_path (str): Directory path to save the plots.
        formula (str): Template formula for the OLS regression model.
        order_names (list): List of names defining the order of factors in the model.
        rename (dict): Dictionary for renaming factors in the model.
        plot_type (str, optional): Type of plot to create. Defaults to 'box'.
        figsize (tuple, optional): Size of the figure. Defaults to (3, 3).
        sig (bool, optional): Whether to perform significance testing for the plots. Defaults to False.

    Returns:
        tuple: Adjusted phenotypic data and a dictionary of p-values for each feature.
    """
    logging.info("########################################################################################")
    logging.info("############# Starting correction of CT scan features for confounders... ###############")
    logging.info("########################################################################################")
    try:
        logging.info("Starting correction of CT scan features for confounders...")

        # Correct CT scan features for confounders
        # adj_pheno, pvalues_lst = correct_ct_scan_features(pheno, outcomes, names, dir_path, formula, order_names, rename)

        # Alternatively, use correct_contrast_ct_scan_features if needed
        adj_pheno, pvalues_lst = correct_contrast_ct_scan_features(pheno, outcomes, names, dir_path, formula, order_names, rename)
        
 
        # Plot adjusted values
        adj_features = ['Adjusted_' + col for col in outcomes]
        
        # print('xxxxxxxx')
        # print(pvalues_lst)
        # print('yyyyyyyyy')

        # Plotting the CT scan data
        plot_features(pheno=adj_pheno,
                      features=adj_features,
                      order=disease_order,  # Order of disease categories
                      names=names,
                      plot_type=plot_type,
                      dir_name=dir_path,
                      # pvalues_lst=pvalues_lst,
                      figsize=figsize,
                      sig=sig)

        logging.info("Correction and plotting of CT scan features completed successfully.")
        return adj_pheno, pvalues_lst
    
    except Exception as e:
        logging.error(f"An error occurred while correcting and plotting CT scan features: {e}")
        raise e

def main():
    """
    The main function to orchestrate the processing and analysis of data.
    It involves loading data, analyzing symptoms, hospitalizations, CBC data,
    correcting lung function features and CT scan features, and plotting results.
    """
    try:
        # Initial setup: Define base directories and other necessary variables
        # Define the base directories for raw data and processed results
        data_dir = '/udd/revfa/0.COPDGene/Clinical_Phenotypes/Data'
        processed_data_dir = 'Processed_Data' #os.path.join(data_dir, 'Processed_Data')
        figures_dir = 'Figures'
        
        # Directories for storing figures related to specific analyses
        symptoms_dir = os.path.join(figures_dir, 'Symptoms')
        cbc_dir = os.path.join(figures_dir, 'CBC')
        pft_dir = os.path.join(figures_dir, 'PFT')
        ct_dir = os.path.join(figures_dir, 'CT_scan')
        suppl_dir = os.path.join(figures_dir, 'Supplementary')

        # Ensure the creation of necessary directories
        for directory in [processed_data_dir, figures_dir, symptoms_dir, cbc_dir, ct_dir, suppl_dir]:
            create_directory_if_not_exists(directory)

        # Load and prepare the data for analysis
        meta, counts, pheno = load_and_prepare_data(data_dir, processed_data_dir, figures_dir)
        # This function handles the loading of meta, counts, and pheno data,
        # classifies subjects by disease, removes non-smokers, and updates subjects, counts, and meta.
        
        # 1. Log the total number of subjects and their disease distribution.
        # 2. Generate and save a plot illustrating the distribution of patients.
        # 3. Create and save a statistics table summarizing the phenotypic data.
        analyze_disease_distribution(pheno, figures_dir)
        
        # Define plotting order for diseases
        order = ['Control', 'Asthma', 'COPD', 'ACO']
        
        # Perform analysis and plotting for symptoms
        logging.info("Starting analysis and plotting of symptoms...")
        analyze_and_plot_symptoms(pheno, order, symptoms_dir)
        # This function focuses on analyzing and visualizing the symptom data.

        # Perform analysis and plotting for hospitalizations
        logging.info("Starting analysis and plotting of hospitalizations...")
        analyze_and_plot_hospitalizations(pheno, order, symptoms_dir)
        # This function is responsible for analyzing and plotting hospitalization data.
        
        # Analyze and plot Complete Blood Count (CBC) data
        analyze_and_plot_cbc(pheno, order, cbc_dir)
        
        # Generate correlation heatmap
        cmap = sns.diverging_palette(0, 255, sep=77, as_cmap=True)
        create_colorbar(cbc_dir, cmap = cmap)
        create_correlation_heatmap(pheno, cbc_dir, cmap = cmap)
        create_correlation_heatmap_per_disease(pheno, 'Disease', cbc_dir, cmap = cmap)

        # This function deals with the analysis and visualization of CBC data.

        # Define the formula for PFT regression model
        formula_pft = "{outcome} ~ Age_P2 +  C(gender) + C(race)  + Height_CM_P2 + C(Disease, Treatment(reference='Control'))"

        adj_pheno, pvalues_lst = correct_and_plot_pft(pheno = pheno, 
                                                      features = PFT_FEATURES,
                                                      names = PFT_NAMES,
                                                      disease_order = order,
                                                      order_names = PFT_order_names, 
                                                      rename = PFT_rename,
                                                      dir_name = pft_dir,
                                                      formula = formula_pft)
        # This function adjusts lung function features considering various confounders.

        # Plot features that are already adjusted and do not require adjustment
        pheno['TLC_pp_P2'] = pheno['TLC_Thirona_P2'] * 100/ pheno['TLC_pred_plethy_P2']
        pheno['FEV1_FVC_pp_P2'] = pheno['FEV1_FVC_post_P2']* 100/pheno['Pred_FEV1_FVC_P2']
        features = ['TLC_pp_P2', 'DLCO_GLI_pp_PbHb_adj_P2', 'FEV1pp_post_P2', 'FEV1_FVC_pp_P2']#,  'FVCpp_post_P2'] ##'FEV1pp_post_P2',
        # names = ['TLC', r'$\mathbf{DL_{CO}}$', 'FEV1pp', 'FEV1_FVC_pp_P2']#, 'FVCpp_post_P2'] #r'FEV$_1$ (%)',
        names = ['TLC', r'$DL_{CO}$', r'$FEV_{1}$', r'$FEV_{1}/FVC$']
        # names = labels = [
        #                     r'$\mathbf{TLC}$',
        #                     r'$\mathbf{DL_{CO}}$',
        #                     r'$\mathbf{FEV_{1}}$',
        #                     r'$\mathbf{FEV_{1}/FVC}$'
        #                 ]


        plot_features(pheno = pheno, features = features, order = order,
                      names = names, plot_type = 'box',
                      dir_name = pft_dir, 
                      figsize = (3,3),
                     sig = False)
        # This function handles the plotting of lung function data post-adjustment.

        # Correct CT scan features for confounders
        formula_CT = "{outcome} ~ Age_P2 + BMI_P2 + C(gender) + C(race) + C(smoking_status_P2) + C(scannerId_P2) + C(Disease, Treatment(reference='Control'))"

        adj_pheno, pvalues_lst = correct_and_plot_ct_scan_features(pheno = pheno, outcomes = CT_FEATURES, 
                                                                   names = CT_NAMES, disease_order = order,
                                                                   order_names = CT_order_names, rename = CT_rename,  
                                                                   dir_path = ct_dir, formula = formula_CT, 
                                                                   sig = True)
        # This function adjusts CT scan features for confounders and generates relevant plots.

        logging.info("Data processing and analysis completed successfully.")

    except Exception as e:
        logging.error(f"An error occurred in the main function: {e}")
        raise


if __name__ == "__main__":
    main()



    
    
    
    
    
    
    
    
 
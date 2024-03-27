import pandas as pd
import logging
import os
from visualization_helpers import plot
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def perform_chi_squared_test(df, disease_group, control_group, feature):
    # Filter data for the disease group and control group
    filtered_df = df[df['Disease'].isin([disease_group, control_group])]
    
    # Create a contingency table
    contingency_table = pd.crosstab(filtered_df['Disease'], filtered_df[feature])
    
    # Perform the Chi-squared test
    chi2, p, dof, expected = chi2_contingency(contingency_table)
    
    # # Print the results
    # logging.info(f"Chi-squared test for {feature} comparing {disease_group} vs. {control_group}:")
    # logging.info(f"Chi-squared statistic: {chi2}, p-value: {p}\n")
    return p

def analyze_and_plot_symptoms(pheno, order, dir_name = 'Figures/Symptoms/', features=None, names=None, palette = ['#D3D3D3','#AFE4DE','#FAA07A','#FCFBCB']):
    """
    Analyze symptoms and plot the data.

    This function prepares the pheno data for specific symptom features and uses the 'plot' 
    function from 'visualization_helpers' to generate visualizations.

    Parameters:
    pheno (pandas.DataFrame): The pheno data containing symptom information.
    order (list): The order in which to plot the categories.
    dir_name (str): Directory name where the plots will be saved.

    Example:
    >>> analyze_and_plot_symptoms(pheno_data, order, 'Figures/Symptoms/')
    """

    logging.info("Starting analysis and plotting of symptoms...")
    logging.info("#######################################################################")
    logging.info("############  Starting analysis and plotting of symptoms...  ##########")
    logging.info("#######################################################################")

    try:
        # Ensure the directory exists
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            logging.info(f"Created directory: {dir_name}")

        # Analyzing first set of symptoms
        features = [ "MMRCDyspneaScor_P2",  "CAT_score_P2" ]
        names = ["Dyspnea Score", "CAT Score"]

        plot(pheno=pheno,
             features=features, 
             order=order, 
             names=names, 
             plot_type='bar', 
             dir_name=dir_name, 
             palette=palette)
        
        # Analyzing first set of symptoms
        features = ["Chronic_Bronchitis_P2", "Emphysema_slv_P2",  "ShrtBrthAttk_P2", "ChestWheez12mo_P2"]
        # names = ['Chronic Bronchitis (%)', "Emphysema (%)",'Shortness of Breath Attack (%)', 'Wheeze (%)']
        names = ['Subjects with \nChronic Bronchitis (%)', 
                 "Subjects with Emphysema (%)",
                 'Subjects with \nShortness of Breath Attack (%)', 
                 'Subjects with Wheeze (%)']
        
        disease_groups = order[1:]
        control_group = 'Control'
        # Initialize a dictionary to store p-values
        p_values = {feature: [] for feature in features}
        # adjusted_p_values = {}
        # Collect p-values
        # print('===>', features)
        for feature in features:
            for disease_group in disease_groups:
                p = perform_chi_squared_test(pheno, disease_group, control_group, feature)
                p_values[feature].append(p)
        # print(p_values)
        
            # # Apply the Bonferroni correction
            # _, p_adjusted, _, _ = multipletests(p_values[feature], alpha=0.05, method='bonferroni') #fdr_bh
            # adjusted_p_values[feature] = p_adjusted
        # print(p_values, adjusted_p_values)

        # for feature in features:
        #     plot(pheno=pheno,
        #          features=[feature], 
        #          pvalues = p_values[feature],
        #          order=order, 
        #          names=names, 
        #          plot_type='bar', 
        #          dir_name=dir_name, 
        #          palette=palette, 
        #         pert = True,
        #         apply_corrections = True)
        
        # print(features, p_values)    
        plot(pheno=pheno,
             features=features, 
             pvalues = p_values,
             order=order, 
             names=names, 
             plot_type='bar', 
             dir_name=dir_name, 
             palette=palette, 
            pert = True,
            apply_corrections = True)

        logging.info("Symptoms data analysis and plotting completed successfully.")
    except Exception as e:
        logging.error(f"Error in analyzing and plotting symptoms: {e}")
        raise

def analyze_and_plot_hospitalizations(pheno, order, dir_name='Figures/Symptoms/', features=None, names=None, palette = ['#D3D3D3','#AFE4DE','#FAA07A','#FCFBCB']):
    """
    Analyze hospitalization data from the pheno dataset and create bar plots.

    Parameters:
    pheno (pandas.DataFrame): DataFrame containing the pheno data.
    order (list): List defining the order of categories for plotting.
    dir_name (str): Directory name where the plots will be saved.
    
    Example:
    >>> analyze_and_plot_hospitalization(pheno_data, order, 'Figures/Symptoms/')
    """
    logging.info("#######################################################################")
    logging.info("######### Starting analysis and plotting of hospitalizations... #######")
    logging.info("#######################################################################")
    
    try:
        # Ensure the directory exists
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            logging.info(f"Created directory: {dir_name}")
            
        # Processing hospitalization data
        # pheno['ExFr'] = (pheno["Exacerbation_Frequency_P2"] >= 2).astype(int)
        # Classifying exacerbation frequency
        pheno.loc[pheno["Exacerbation_Frequency_P2"] < 2, 'ExFr'] = 0
        pheno.loc[pheno["Exacerbation_Frequency_P2"] >= 2, 'ExFr'] = 1
        
        # Define hospitalization features and their respective names
        features = ["ExFr", "Severe_Exacerbations_P2"]
        names = ['Subjects with \nFrequent Exacerbations (%)', 'Subjects with \nSevere Exacerbations (%)']
        
        disease_groups = order[1:]
        control_group = 'Control'
        # Initialize a dictionary to store p-values
        p_values = {feature: [] for feature in features}
        adjusted_p_values = {}
        # Collect p-values
        for feature in features:
            for disease_group in disease_groups:
                p = perform_chi_squared_test(pheno, disease_group, control_group, feature)
                p_values[feature].append(p)
            # # Apply the Bonferroni correction
            # _, p_adjusted, _, _ = multipletests(p_values[feature], alpha=0.05, method='bonferroni')
            # adjusted_p_values[feature] = p_adjusted
        
        plot(pheno=pheno,
             features=features, 
             pvalues=p_values,
             order=order, 
             names=names, 
             plot_type='bar', 
             dir_name=dir_name, 
             palette=palette, 
             pert = True,
             apply_corrections = True)
        
        # for feature in features:
        #     # Plotting hospitalization data
        #     # print(p_values[feature])
        #     plot(pheno = pheno,
        #          features=[feature], 
        #          order=order, 
        #          pvalues = p_values[feature],
        #          names=names, 
        #          plot_type='bar', 
        #          dir_name=dir_name, 
        #          palette=palette,
        #         pert = True,
        #         apply_corrections = True)

        logging.info("Hospitalization data analysis and plotting completed successfully.")
    except Exception as e:
        logging.error(f"Error in analyzing hospitalizations: {e}")
        raise

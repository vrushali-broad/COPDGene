import logging
import matplotlib.pyplot as plt

def classify_subjects(pheno):
    """
    Classify subjects into different categories based on asthma and COPD criteria.

    Parameters:
    pheno (pandas.DataFrame): DataFrame containing patient data with the following key columns:
      - Asthma_P1: Have you ever had asthma (1=Yes, 0=No, 3=Don't know)
      - finalGold_P2: Final GOLD stage, post-BD Phase 2
        (0 = GOLD 0 Control (FEV1 >= 80%, FEV1/FVC >= 0.7),
         1 = GOLD 1 (FEV1 >= 80%, FEV1/FVC < 0.7),
         2 = GOLD 2 (50% <= FEV1 < 80%, FEV1/FVC < 0.7),
         3 = GOLD 3 (30% <= FEV1 < 50%, FEV1/FVC < 0.7),
         4 = GOLD 4 (FEV1 < 30%, FEV1/FVC < 0.7),
         -1 = PRISm (Preserved Ratio Impaired Spirometry),
         -2 = Never-smoked Normal, -9 = ineligible)
      - AsthmaDxByDr_P1: Diagnosed by doctor or other health professional (1=Yes, 0=No, 3=Don't know)
      - AsthmaAge_P1: Age at asthma onset
      - AsthmaAgeDK_P1: Asthma started in childhood, age not known (1=Don't Know)
      - AsthmaStillHave_P1: Current status of asthma (1=Yes, 0=No, 3=Don't know)
      - smoking_status_P2: Smoking status (0 = Never-smoked, 1 = Former smoker, 2 = Current smoker)

    Returns:
    pandas.DataFrame: DataFrame with an additional column 'Disease' indicating the patient classification.

    The function classifies patients into 'Non-smoker', 'Control', 'Asthma', 'COPD', and 'ACO' categories based on
    their asthma history, GOLD stage, diagnosis by a doctor, age of asthma onset, and smoking status.
    
    Parameters:
    pheno (pandas.DataFrame): The pheno data to be classified.
    
    Example:
    >>> classified_pheno = classify_patients(pheno_data)
    """
    try:
        # Create a copy of the DataFrame to avoid SettingWithCopyWarning
        pheno = pheno.copy()
        logging.info(f'Total subjects: {len(pheno)}')
        
        logging.info(f"PRISm: {pheno[pheno['finalGold_P2'] == -1].shape[0]}")
        logging.info(f"Never Smokers: {pheno[pheno['finalGold_P2'] == -2].shape[0]}")
        logging.info(f"GOLD stage 1: {pheno[pheno['finalGold_P2'] == 1].shape[0]}")

        # Non-smoking controls
        pheno.loc[((pheno['smoking_status_P2'] == 0) & (pheno['finalGold_P2'] <= 0)), 'Disease'] = 'Non-smoker'

        # Filter out non-smokers for further classification
        pheno = pheno[pheno['smoking_status_P2'] != 0]
        logging.info(f"Total Smokers: {pheno.shape[0]}")

        # No asthma No COPD
        pheno.loc[(pheno['Asthma_P1'] == 0) & (pheno['finalGold_P2'] == 0), 'Disease'] = 'Control'

        # Asthma and no COPD
        pheno.loc[(pheno['Asthma_P1'] == 1) & (pheno['AsthmaDxByDr_P1'] == 1) & 
                  ((pheno['AsthmaAge_P1'] <= 40) | (pheno['AsthmaAgeDK_P1'] == 1)) & 
                  (pheno['finalGold_P2'] == 0), 'Disease'] = 'Asthma'

        # No asthma and COPD
        pheno.loc[(pheno['Asthma_P1'] == 0) & (pheno['finalGold_P2'].isin([2, 3, 4])), 'Disease'] = 'COPD'

        # Asthma and COPD both
        pheno.loc[(pheno['Asthma_P1'] == 1) & (pheno['AsthmaDxByDr_P1'] == 1) & 
                  ((pheno['AsthmaAge_P1'] <= 40) | (pheno['AsthmaAgeDK_P1'] == 1)) & 
                  (pheno['finalGold_P2'].isin([2, 3, 4])), 'Disease'] = 'ACO'


        # Filter out patients with no disease classification
        pheno = pheno[pheno['Disease'].notna()]

        logging.info(f"Classified patients into diseases: {pheno['Disease'].value_counts()}")
        return pheno
    except Exception as e:
        logging.error(f"Error in classifying patients: {e}")
        raise


# def classify_subjects(pheno):
#     """
#     Classify subjects into different categories based on asthma and COPD criteria.

#     Parameters:
#     pheno (pandas.DataFrame): DataFrame containing patient data with the following key columns:
#       - Asthma_P1: Have you ever had asthma (1=Yes, 0=No, 3=Don't know)
#       - finalGold_P2: Final GOLD stage, post-BD Phase 2
#         (0 = GOLD 0 Control (FEV1 >= 80%, FEV1/FVC >= 0.7),
#          1 = GOLD 1 (FEV1 >= 80%, FEV1/FVC < 0.7),
#          2 = GOLD 2 (50% <= FEV1 < 80%, FEV1/FVC < 0.7),
#          3 = GOLD 3 (30% <= FEV1 < 50%, FEV1/FVC < 0.7),
#          4 = GOLD 4 (FEV1 < 30%, FEV1/FVC < 0.7),
#          -1 = PRISm (Preserved Ratio Impaired Spirometry),
#          -2 = Never-smoked Normal, -9 = ineligible)
#       - AsthmaDxByDr_P1: Diagnosed by doctor or other health professional (1=Yes, 0=No, 3=Don't know)
#       - AsthmaAge_P1: Age at asthma onset
#       - AsthmaAgeDK_P1: Asthma started in childhood, age not known (1=Don't Know)
#       - AsthmaStillHave_P1: Current status of asthma (1=Yes, 0=No, 3=Don't know)
#       - smoking_status_P2: Smoking status (0 = Never-smoked, 1 = Former smoker, 2 = Current smoker)

#     Returns:
#     pandas.DataFrame: DataFrame with an additional column 'Disease' indicating the patient classification.

#     The function classifies patients into 'Non-smoker', 'Control', 'Asthma', 'COPD', and 'ACO' categories based on
#     their asthma history, GOLD stage, diagnosis by a doctor, age of asthma onset, and smoking status.
    
#     Parameters:
#     pheno (pandas.DataFrame): The pheno data to be classified.
    
#     Example:
#     >>> classified_pheno = classify_patients(pheno_data)
#     """
#     try:
#         # Create a copy of the DataFrame to avoid SettingWithCopyWarning
#         pheno = pheno.copy()
#         # Convert the 'finalGold_P2' column to an integer
#         pheno['finalGold_P2'] = pheno['finalGold_P2'].astype(int)
        
#         logging.info(f'Total subjects: {len(pheno)}')
        
#         logging.info(f"PRISm: {pheno[pheno['finalGold_P2'] == -1].shape[0]}")
#         logging.info(f"Never Smokers: {pheno[pheno['finalGold_P2'] == -2].shape[0]}")
#         logging.info(f"GOLD stage 1: {pheno[pheno['finalGold_P2'] == 1].shape[0]}")

#         # Initialize the 'Disease' column to avoid any pre-existing data affecting results
#         pheno['Disease'] = None

#         # Non-smoking controls
#         pheno.loc[((pheno['smoking_status_P2'] == 0) & (pheno['finalGold_P2'] <= 0)), 'Disease'] = 'Non-smoker'

#         # Filter out non-smokers for further classification
#         pheno = pheno[pheno['smoking_status_P2'] != 0]
#         logging.info(f"Total Smokers: {pheno.shape[0]}")

#         # No asthma, GOLD stage 0 (controls)
#         pheno.loc[(pheno['Asthma_P1'] == 0) & (pheno['finalGold_P2'] == 0), 'Disease'] = 'Control'

#         # Asthma and no COPD
#         pheno.loc[(pheno['Asthma_P1'] == 1) & (pheno['AsthmaDxByDr_P1'] == 1) & 
#                   ((pheno['AsthmaAge_P1'] <= 40) | (pheno['AsthmaAgeDK_P1'] == 1)) & 
#                   (pheno['finalGold_P2'] == 0), 'Disease'] = 'Asthma'

#         # No asthma and COPD (GOLD stages 2, 3, 4 only)
#         pheno.loc[(pheno['Asthma_P1'] == 0) & (pheno['finalGold_P2'].isin([2, 3, 4])), 'Disease'] = 'COPD'

#         # Asthma and COPD both (GOLD stages 2, 3, 4 only)
#         pheno.loc[(pheno['Asthma_P1'] == 1) & (pheno['AsthmaDxByDr_P1'] == 1) & 
#                   ((pheno['AsthmaAge_P1'] <= 40) | (pheno['AsthmaAgeDK_P1'] == 1)) & 
#                   (pheno['finalGold_P2'].isin([2, 3, 4])), 'Disease'] = 'ACO'

#         # Remove patients with no disease classification
#         pheno = pheno[pheno['Disease'].notna()]

#         logging.info(f"Classified patients into diseases: {pheno['Disease'].value_counts()}")
#         return pheno
#     except Exception as e:
#         logging.error(f"Error in classifying patients: {e}")
#         raise


def classify_subjects_LLN(pheno):
    """
    Classify subjects into different categories based on asthma and COPD criteria.

    Parameters:
    pheno (pandas.DataFrame): DataFrame containing patient data with the following key columns:
      - Asthma_P1: Have you ever had asthma (1=Yes, 0=No, 3=Don't know)
      - finalGold_P2: Final GOLD stage, post-BD Phase 2
        (0 = GOLD 0 Control (FEV1 >= 80%, FEV1/FVC >= 0.7),
         1 = GOLD 1 (FEV1 >= 80%, FEV1/FVC < 0.7),
         2 = GOLD 2 (50% <= FEV1 < 80%, FEV1/FVC < 0.7),
         3 = GOLD 3 (30% <= FEV1 < 50%, FEV1/FVC < 0.7),
         4 = GOLD 4 (FEV1 < 30%, FEV1/FVC < 0.7),
         -1 = PRISm (Preserved Ratio Impaired Spirometry),
         -2 = Never-smoked Normal, -9 = ineligible)
      - AsthmaDxByDr_P1: Diagnosed by doctor or other health professional (1=Yes, 0=No, 3=Don't know)
      - AsthmaAge_P1: Age at asthma onset
      - AsthmaAgeDK_P1: Asthma started in childhood, age not known (1=Don't Know)
      - AsthmaStillHave_P1: Current status of asthma (1=Yes, 0=No, 3=Don't know)
      - smoking_status_P2: Smoking status (0 = Never-smoked, 1 = Former smoker, 2 = Current smoker)
      - FEV1_FVC_post_P2: Post-bronchodilator FEV1/FVC ratio
      - Pred_FEV1_FVC_LLN_P2: Lower limit of normal for FEV1/FVC

    Returns:
    pandas.DataFrame: DataFrame with an additional column 'Disease' indicating the patient classification.

    The function classifies patients into 'Non-smoker', 'Control', 'Asthma', 'COPD', and 'ACO' categories based on
    their asthma history, GOLD stage, diagnosis by a doctor, age of asthma onset, smoking status, and post-bronchodilator 
    spirometry values.

    - Non-smoker: Individuals who have never smoked and have normal spirometry.
    - Control: Smokers without asthma and normal spirometry.
    - Asthma: Individuals diagnosed with asthma by a doctor, with early onset (age <= 40 or unknown), and normal spirometry.
    - COPD: Individuals without asthma but with COPD defined by FEV1/FVC post-bronchodilator ratio less than the lower limit of normal (Pred_FEV1_FVC_LLN_P2).
    - ACO: Individuals with both asthma (diagnosed by a doctor, with early onset) and COPD as defined by the same spirometry criteria.

    Example:
    >>> classified_pheno = classify_patients(pheno_data)
    """
    try:
        # Create a copy of the DataFrame to avoid SettingWithCopyWarning
        pheno = pheno.copy()
        logging.info(f'Total subjects: {len(pheno)}')
        
        logging.info(f"PRISm: {pheno[pheno['finalGold_P2'] == -1].shape[0]}")
        logging.info(f"Never Smokers: {pheno[pheno['finalGold_P2'] == -2].shape[0]}")
        logging.info(f"GOLD stage 1: {pheno[pheno['finalGold_P2'] == 1].shape[0]}")

        # Non-smoking controls
        pheno.loc[((pheno['smoking_status_P2'] == 0) & (pheno['finalGold_P2'] <= 0)), 'Disease_LLN'] = 'Non-smoker'

        # Filter out non-smokers for further classification
        pheno = pheno[pheno['smoking_status_P2'] != 0]
        logging.info(f"Total Smokers: {pheno.shape[0]}")

        # No asthma No COPD
        pheno.loc[(pheno['Asthma_P1'] == 0) & (pheno['finalGold_P2'] == 0), 'Disease_LLN'] = 'Control'

        # Asthma and no COPD
        pheno.loc[(pheno['Asthma_P1'] == 1) & (pheno['AsthmaDxByDr_P1'] == 1) & 
                  ((pheno['AsthmaAge_P1'] <= 40) | (pheno['AsthmaAgeDK_P1'] == 1)) & 
                  (pheno['finalGold_P2'] == 0), 'Disease_LLN'] = 'Asthma'

        # No asthma and COPD
        pheno.loc[(pheno['Asthma_P1'] == 0) & (pheno['FEV1_FVC_post_P2'] < pheno['Pred_FEV1_FVC_LLN_P2']), 'Disease_LLN'] = 'COPD'

        # Asthma and COPD both
        pheno.loc[(pheno['Asthma_P1'] == 1) & (pheno['AsthmaDxByDr_P1'] == 1) & 
                  ((pheno['AsthmaAge_P1'] <= 40) | (pheno['AsthmaAgeDK_P1'] == 1)) & 
                  (pheno['FEV1_FVC_post_P2'] < pheno['Pred_FEV1_FVC_LLN_P2']), 'Disease_LLN'] = 'ACO'

        # Filter out patients with no disease classification
        pheno = pheno[pheno['Disease'].notna()]
        pheno.loc[pheno['Disease_LLN'].isna(), 'Disease_LLN'] = 'Excluded'

        logging.info(f"Classified patients into diseases based on LLN: {pheno['Disease_LLN'].value_counts()}")
        return pheno
    except Exception as e:
        logging.error(f"Error in classifying patients: {e}")
        raise

#####################################################################
########### Generates donut plot of patient distribution ############
#####################################################################

def plot_subject_distribution(pheno, output_path):
    """
    Plot the distribution of subject by diseases.
    
    :param data: DataFrame containing patient data
    :param column: Column name in DataFrame that contains disease categories
    :param filename: Name of the file to save the plot

    Parameters:
    pheno (pandas.DataFrame): DataFrame containing patient data with a 'Disease' column.
    output_path (str): File path where the plot will be saved.

    The function creates a pie chart showing the distribution of different diseases.
    """

    try:
        logging.info("Plotting patient disease distribution...")

        # Calculate the distribution of diseases
        disease_counts = pheno.Disease.value_counts()
        # logging.info(f"Disease counts:\n{disease_counts}")

        # Calculate the distribution of GOLD stages
        # gold_counts = pheno.finalGold_P2.value_counts()
        # logging.info(f"GOLD stage counts:\n{gold_counts}")

        # Prepare data for the pie chart
        labels = ['Control', 'Asthma', 'COPD', 'ACO']
        data = disease_counts[labels] / len(pheno[pheno.Disease.isin(labels)])
        explode = [0] * len(labels)
        colors = ['#d3d3d3', '#AFE4DE', '#ffa07a', '#fff157']

        # Create and customize the pie chart
        fig, ax = plt.subplots(figsize=(4, 4))
        wedges, texts, autotexts = ax.pie(data, colors=colors, labels=labels, 
                                          wedgeprops=dict(width=0.5), startangle=-50, 
                                          autopct='%1.1f%%', pctdistance=0.75, explode=explode,
                                          textprops={'fontsize': 16})

        plt.setp(autotexts, size=12)
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

        logging.info(f"Patient distribution plot saved to {output_path}")
    except Exception as e:
        logging.error(f"Error in plotting patient distribution: {e}")
        raise


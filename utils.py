import os
import logging
import numpy as np
import pandas as pd
from itertools import combinations
from visualization_helpers import plot
from scipy.stats import chi2_contingency, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind

## Import transformation functions
from scipy.stats import boxcox, norm
from scipy.special import logit
from sklearn.preprocessing import PowerTransformer, MinMaxScaler, StandardScaler, QuantileTransformer

## Plotting imports
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.distributions.empirical_distribution import ECDF

def format_pvalue(p):
    if p < 0.0001:
        return '****'
    elif p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'

def create_directory_if_not_exists(directory):
    """Create a directory if it does not exist."""
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
            logging.info(f"Created directory: {directory}")
    except Exception as e:
        logging.error(f"Error creating directory {directory}: {e}")
        raise

def perform_chi_squared_test(df, feature, order):
    """
    Calculate p-values for all possible pairs of groups in the 'Disease' column using the Chi-squared test.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data.
    - feature (str): The name of the feature/column for which to calculate the p-values.
    - order (list): The order of groups to consider for pairwise comparisons.

    Returns:
    - corrected_pvalues (dict): A dictionary with pairs as keys and corrected p-values as values.
    """
    # Generate all possible pairs of the groups
    # all_pairs = list(combinations(order, 2))
    
    # Define the specific pairs to consider
    all_pairs = specific_pairs = [
                                    ('Asthma', 'Control'),
                                    ('COPD', 'Control'),
                                    ('ACO', 'Control'),
                                    ('Asthma', 'COPD'),
                                    ('Asthma', 'ACO'),
                                    ('COPD', 'ACO')
                                ]
    
    # Initialize dictionary to store p-values
    pvalues = {}
    
    # Calculate p-values for each pair
    for pair in all_pairs:
        filtered_df = df[df['Disease'].isin([pair[0], pair[1]])]
        
        # Create a contingency table
        contingency_table = pd.crosstab(filtered_df['Disease'], filtered_df[feature])
        
        # Perform the Chi-squared test
        chi2, p, dof, expected = chi2_contingency(contingency_table)
        pvalues['_vs_'.join(pair)] = p
    
    corrected_pvalues = multipletests(list(pvalues.values()), alpha=0.05, method='bonferroni')[1] #method = 'fdr_bh')[1] #

    # Create a Pandas Series for corrected p-values
    corrected_pvalues_series = pd.Series(corrected_pvalues, index=pvalues.keys())
    # print('=========', feature, '===========\n', corrected_pvalues_series.apply(lambda x: x if x < 0.05 else 'n.s.'))
    print('=========', feature, '===========\n', feature, '\n', corrected_pvalues_series.apply(format_pvalue))
    
    pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
    pvalues = [corrected_pvalues_series.get(pair, None) for pair in pairs_order] 

    return pvalues

def perform_mwu_test(pheno, feature, order, ct=False):
    """
    Calculate p-values for all possible pairs of groups in the 'Disease' column.

    Parameters:
    - pheno (pd.DataFrame): DataFrame containing the data.
    - feature (str): The name of the feature/column for which to calculate the p-values.
    - order (list): The order of groups to consider for pairwise comparisons.

    Returns:
    - pvalues (dict): A dictionary with pairs as keys and p-values as values.
    """
    # Generate all possible pairs of the groups
    # all_pairs = list(combinations(order, 2))
    
    # Define the specific pairs to consider
    all_pairs = specific_pairs = [
                                    ('Asthma', 'Control'),
                                    ('COPD', 'Control'),
                                    ('ACO', 'Control'),
                                    ('Asthma', 'COPD'),
                                    ('Asthma', 'ACO'),
                                    ('COPD', 'ACO')
                                ]
    
    # Initialize dictionary to store p-values
    pvalues = {}
    
    # Calculate p-values for each pair
    for pair in all_pairs:
        group1 = pheno[pheno['Disease'] == pair[0]][feature].dropna()
        group2 = pheno[pheno['Disease'] == pair[1]][feature].dropna()
        
        _, pval = mannwhitneyu(group1, group2)
        pvalues['_vs_'.join(pair) ] = pval
        
    corrected_pvalues = multipletests( list(pvalues.values()), alpha=0.05, method='bonferroni')[1]

    # Create a Pandas Series for corrected p-values
    corrected_pvalues_series = pd.Series(corrected_pvalues, index=pvalues.keys())
    # print('=========', feature, '===========\n', feature, '\n', corrected_pvalues_series.apply(lambda x: x if x < 0.05 else 'n.s.'))

    
    
    print('=========', feature, '===========\n', feature, '\n', corrected_pvalues_series.apply(format_pvalue))

    
    if ct:
        pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control', 'COPD_vs_ACO']
    else:
        pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
    pvalues = [corrected_pvalues_series.get(pair, None) for pair in pairs_order] 
    
    return pvalues

def perform_ttest(pheno, feature, order, ct=False):
    """
    Calculate p-values for all possible pairs of groups in the 'Disease' column using t-tests.

    Parameters:
    - pheno (pd.DataFrame): DataFrame containing the data.
    - feature (str): The name of the feature/column for which to calculate the p-values.
    - order (list): The order of groups to consider for pairwise comparisons.

    Returns:
    - pvalues (dict): A dictionary with pairs as keys and p-values as values.
    """
    # Define the specific pairs to consider
    all_pairs = [
        ('Asthma', 'Control'),
        ('COPD', 'Control'),
        ('ACO', 'Control'),
        ('Asthma', 'COPD'),
        ('Asthma', 'ACO'),
        ('COPD', 'ACO')
    ]
    
    # Initialize dictionary to store p-values
    pvalues = {}
    
    # Calculate p-values for each pair
    for pair in all_pairs:
        group1 = pheno[pheno['Disease'] == pair[0]][feature].dropna()
        group2 = pheno[pheno['Disease'] == pair[1]][feature].dropna()
        
        _, pval = ttest_ind(group1, group2)
        pvalues['_vs_'.join(pair)] = pval
        
    corrected_pvalues = multipletests(list(pvalues.values()), alpha=0.05, method='bonferroni')[1]

    # Create a Pandas Series for corrected p-values
    corrected_pvalues_series = pd.Series(corrected_pvalues, index=pvalues.keys())

    # Print formatted p-values
    print('=========', feature, '===========\n', feature, '\n', corrected_pvalues_series.apply(lambda x: f'{x:.5f}' if x < 0.05 else 'n.s.'))
    
    if ct:
        pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control', 'COPD_vs_ACO']
    else:
        pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
    
    ordered_pvalues = [corrected_pvalues_series.get(pair, None) for pair in pairs_order]
    
    return ordered_pvalues


###########################################################################
################### Multivariate Regression Analysis ######################
###########################################################################

def save_kde_ecdf_plot(df, outcome, title, file_path):
    """
    Creates and saves a KDE plot and the original ECDF plot for the specified outcome as subplots.

    Parameters:
    - df (DataFrame): The data frame containing the data.
    - outcome (str): The column name of the outcome variable to plot.
    - title (str): The title of the plot.
    - file_path (str): The full file path where the plot will be saved.
    """
    # Set up the figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(6, 3))
    plt.rcParams['axes.edgecolor'] = '#D3D3D3'

    # KDE plot
    sns.histplot(df[outcome], kde=True, ax=axs[0])
    axs[0].set_title(f'KDE Plot: {title}', fontsize=10)
    axs[0].set_xlabel(title, fontsize=8)
    axs[0].set_ylabel('Density', fontsize=8)
    axs[0].tick_params(axis='both', which='major', labelsize=8)
    
    # Original ECDF plot
    data = df[outcome]
    ecdf = ECDF(data)
    x = np.linspace(min(data), max(data), num=1000)
    y = ecdf(x)
    
    axs[1].step(x, y, where="post")
    axs[1].set_title(f'ECDF Plot: {title}', fontsize=10)
    axs[1].set_xlabel(title, fontsize=8)
    axs[1].set_ylabel('ECDF', fontsize=8)
    axs[1].tick_params(axis='both', which='major', labelsize=8)

    # Adjust layout and save the plot
    plt.tight_layout()
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(file_path, bbox_inches='tight')
    plt.close()

def scale_to_unit_interval(data, column_name, epsilon=0.0001):
    """
    Scale the data to fit within the (0, 1) interval.
    
    Parameters:
    - data: pandas DataFrame containing the data
    - column_name: str, name of the column to be transformed
    - epsilon: float, small constant to avoid boundary issues (default is 0.0001)
    
    Returns:
    - transformed_data: pandas Series with transformed values
    """
    min_val = data[column_name].min()
    max_val = data[column_name].max()
    transformed_data = (data[column_name] - min_val + epsilon) / (max_val - min_val + 2 * epsilon)
    return transformed_data

# %%%%%%%%% Resting_SaO2_P2 8.472135811722177
# %%%%%%%%% FEV1_FVC_post_P2 8.472135811722177
def cubic_root_transform(data):
    return np.cbrt(data)

def fractional_polynomial_transform(data, powers):
    transformed_data = np.zeros_like(data, dtype=float)
    for power in powers:
        if power == 0:
            transformed_data += np.log(data)
        else:
            transformed_data += data ** power
    return transformed_data


def apply_transformations(df, features, transformation):
    """
    Apply the specified transformation to the given features in the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the data.
        features (list): List of feature names to transform.
        transformation (str): Type of transformation ('log', 'sqrt', 'inverse', 'boxcox', 'yeojohnson').

    Returns:
        pd.DataFrame: DataFrame with transformed features.
    """
    transformer = {}
    for feature in features:
        if transformation == 'log':
            df[f'log_{feature}'] = np.log(df[feature] + 100)  # Adding 1 to avoid log(0)
        elif transformation == 'log_shift':
            df[f'log_shift_{feature}'] = np.log(df[feature] + abs(df[feature].min()) + 1)  # Shifting to handle zero/negative values
        elif transformation == 'sqrt':
            df[f'sqrt_{feature}'] = np.sqrt(df[feature])
        elif transformation == 'inverse':
            df[f'inverse_{feature}'] = 1 / df[feature]
        elif transformation == 'quantile':
            qt = QuantileTransformer(output_distribution='normal', random_state=0)
            df[f'quantile_{feature}'] = qt.fit_transform(df[[feature]])
            transformer[feature] = qt
        elif transformation == 'boxcox':
            # Ensure all values are positive before applying Box-Cox
            df[f'boxcox_{feature}'], lambda_optimal = boxcox(df[feature] + 1)  # Adding 1 to avoid non-positive values
            # print('%%%%%%%%%', feature, lambda_optimal)
            # df[f'boxcox_{feature}'] = boxcox(df[feature] + 1, lmbda= 0.5263)  # Adding 1 to avoid non-positive values
        elif transformation == 'yeojohnson':
            pt = PowerTransformer(method='yeo-johnson')
            df[f'yeojohnson_{feature}'] = pt.fit_transform(df[[feature]])
            transformer[feature] = pt
        elif transformation == 'probit':
            # Scale to [0, 1] before applying probit
            epsilon = 1e-5  # Small constant to avoid issues with 0 and 1
            scaled_feature = (df[feature] - df[feature].min()) / (df[feature].max() - df[feature].min())
            scaled_feature = np.clip(scaled_feature, epsilon, 1 - epsilon)
            df[f'probit_{feature}'] = norm.ppf(scaled_feature)
        elif transformation == 'logit_scaled':
            # Scale to [0, 1] before applying logit
            epsilon = 1e-5  # Small constant to avoid issues with 0 and 1
            scaled_feature = (df[feature] - df[feature].min()) / (df[feature].max() - df[feature].min())
            scaled_feature = np.clip(scaled_feature, epsilon, 1 - epsilon)
            df[f'logit_scaled_{feature}'] = np.log(scaled_feature / (1 - scaled_feature))
            # Scale percentages between 0 and 1 before applying logit transformation
            # scaled_feature = df[feature] / 100.0
            # df[f'logit_{feature}'] = logit(scaled_feature.clip(1e-10, 1-1e-10))  # Clipping to avoid logit(0) and logit(1)
        elif transformation == 'minmax':
            scaler = MinMaxScaler()
            df[f'minmax_{feature}'] = scaler.fit_transform(df[[feature]])
        elif transformation == 'standard':
            scaler = StandardScaler()
            df[f'standard_{feature}'] = scaler.fit_transform(df[[feature]])
        elif transformation == 'standard_minmax':
            # Standard Scaling
            standard_scaler = StandardScaler()
            df[f'standard_{feature}'] = standard_scaler.fit_transform(df[[feature]])
            
            # Min-Max Scaling to (0, 1)
            min_max_scaler = MinMaxScaler()
            df[f'minmax_{feature}'] = min_max_scaler.fit_transform(df[[f'standard_{feature}']])
        elif transformation == 'logit':
            epsilon = 1e-5  # Small constant to avoid issues with 0 and 1
            # Ensure feature is scaled to (0, 1) before applying logit
            scaler = MinMaxScaler()
            scaled_feature = scaler.fit_transform(df[[feature]])
            df[f'logit_{feature}'] = logit(np.clip(scaled_feature, epsilon, 1 - epsilon))
        elif transformation == 'scale_to_unit_interval':
            df[f'scale_to_unit_interval_{feature}'] = scale_to_unit_interval(df, feature)
        elif transformation == 'cubic_root':
            df[f'cubic_root_{feature}'] = cubic_root_transform(df[feature])
        elif transformation == 'fractional_polynomial':
            # Example powers, this should be configured as needed
            powers = [0.5, -1, 2]
            df[f'fractional_polynomial_{feature}'] = fractional_polynomial_transform(df[feature], powers)
        else:
            raise ValueError(f"Unsupported transformation: {transformation}")
    return df, transformer

def inverse_transform(df, feature, transformed_feature, transformer):
    """
    Inverse transform the specified transformed feature in the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the data.
        feature (str): Original feature name.
        transformed_feature (str): Transformed feature name.
        transformer: Transformer object used for the transformation.

    Returns:
        pd.DataFrame: DataFrame with the inverse transformed feature.
    """
    if f'Adjusted_{transformed_feature}' in df.columns:
        df[f'Adjusted_{feature}'] = transformer.inverse_transform(df[[f'Adjusted_{transformed_feature}']])
    return df









    
def auto_detect_bounds(data, buffer_ratio=0.05):
    """
    Automatically detect lower and upper bounds with a buffer.
    :param data: Input data array.
    :param buffer_ratio: Proportion of data range to be added as buffer.
    :return: Tuple (lower_bound, upper_bound).
    """
    min_val = np.min(data)
    max_val = np.max(data)
    data_range = max_val - min_val
    
    lower_bound = min_val - buffer_ratio * data_range
    upper_bound = max_val + buffer_ratio * data_range
    
    return lower_bound, upper_bound

import numpy as np
from scipy.special import logit
def transform_data(data, lower_bound, upper_bound):
    """
    Transforms the data using logit transformation for bounded data.
    Adjusts boundary values slightly to avoid infinity issues.
    """
    epsilon = 1e-5
    # Scale data to [0, 1]
    scaled_data = (data - lower_bound) / (upper_bound - lower_bound)
    # Adjust boundary values
    scaled_data = np.clip(scaled_data, epsilon, 1 - epsilon)
    # Apply logit transformation
    transformed_data = logit(scaled_data)
    return transformed_data

def inverse_transform_data(transformed_data, lower_bound, upper_bound):
    """
    Applies the inverse logit transformation and scales the data back to original bounds.
    """
    # Inverse logit transformation
    scaled_data = expit(transformed_data)
    # Scale back to original bounds
    original_data = scaled_data * (upper_bound - lower_bound) + lower_bound
    return original_data
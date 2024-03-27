import os
import logging
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from visualization_helpers import forest_plot, get_significance_colors, get_significance_stars  

def perform_regression_analysis(pheno, outcomes, names, order_names, rename, formula, dir_name, contrast = False):
    """
    Conducts regression analysis on pulmonary function test outcomes while adjusting for confounders.

    This function prepares the input data by normalizing specific pulmonary function test outcomes and then performs 
    ordinary least squares (OLS) regression for each outcome specified. It adjusts the data for confounders as defined
    in the provided regression formula and stores the results. Additionally, it generates forest plots for each 
    outcome with the corresponding regression coefficients, confidence intervals, and p-values.

    Parameters:
    - pheno (DataFrame): The input DataFrame containing phenotypic and pulmonary function test data.
    - outcomes (list): A list of outcome variable names to be included in the regression analysis.
    - names (list): List of names corresponding to the outcomes, used in plotting.
    - order_names (list): Order of the names for plotting purposes.
    - rename (dict): Dictionary for renaming the coefficients in the forest plot.
    - formula (str): A string representing the regression formula to be used in the analysis.
    - dir_name (str): Directory name where the output forest plot images will be saved.
    - contrast (bool): Whether to calculate and include the contrast p-value between COPD and ACO in the results.

    Returns:
    - df (DataFrame): The original DataFrame augmented with residuals and adjusted outcome variables.
    - pvalues_lst (dict): Dictionary containing p-values for the regression coefficients of each outcome and, if 
      requested, the contrast p-value between COPD and ACO.

    The function adds two new columns for each outcome in the DataFrame: one for the residuals ('corrected_') 
    and one for the adjusted values ('Adjusted_'). Forest plots are saved as PDF files in the specified directory.
    """
    # Prepare the dataframe
    # pheno['TLC_pp_P2'] = pheno['TLC_Thirona_P2'] / pheno['TLC_pred_plethy_P2']
    # pheno['FEV1_FVC_pp_P2'] = pheno['FEV1_FVC_post_P2'] / pheno['Pred_FEV1_FVC_P2']
    
    df = pd.DataFrame(pheno)
    df['intercept'] = 1

    # Perform regression and store results
    results = {}
    contrast_results = {}
    for outcome in outcomes:
        current_formula = formula.format(outcome=outcome)
        model = smf.ols(formula=current_formula, data=df).fit()
        df['corrected_' + outcome] = model.resid
        df['Adjusted_' + outcome] = model.fittedvalues
        results[outcome] = model
        
        # # Calculate contrast between COPD and ACO
        # contrast_pvalue, contrast = calculate_contrast(model)
        # contrast_results[outcome] = {'p_value': contrast_pvalue, 'contrast': contrast}

    # Analysis and plotting
    pvalues_lst = {}
    for i, (key, result) in enumerate(results.items()):
        params = result.params
        
        pvalues = result.pvalues
        conf = result.conf_int()
        conf['Estimate'] = params
        conf.columns = ['Lower', 'Upper', 'Estimate']
        conf = conf.drop("Intercept")
        pvalues = pvalues.drop("Intercept")
        colors = get_significance_colors(pvalues, conf['Estimate'])
        
        fig_name = os.path.join(dir_name, 'Forest_Adj_' + outcomes[i] + '.pdf')
        forest_pvalues = forest_plot(conf.index,
                                     conf['Estimate'], conf['Lower'], conf['Upper'], 
                                     colors,
                                     pvalues=pvalues, 
                                     order=order_names, 
                                     rename=rename, 
                                     title=names[i], 
                                     file_path = fig_name)
        
        key = 'Adjusted_'+key
        pvalues_lst[key] = forest_pvalues
    
        # Calculate and store the contrast p-value if contrast analysis is requested
        if contrast:
            contrast_pvalue, _ = calculate_contrast(result)
            pvalues_lst[key].append(contrast_pvalue)
            # pvalues_lst[key]['contrast_pvalue'] = contrast_pvalue
        

    return df, pvalues_lst

def calculate_contrast(model):
    """
    Calculate the contrast between COPD and ACO from the OLS model.
    
    This function computes the contrast as the difference in the regression coefficients of the COPD and ACO groups. 
    It also calculates the p-value associated with this contrast using a t-test. The standard error of the contrast 
    is computed under the assumption of independence of the coefficients.

    Args:
        model (RegressionResults): The fitted regression model.

    Returns:
        tuple: A tuple containing the p-value and the contrast value.
        
    Raises:
        Exception: If any error occurs during the computation.
    """
    try:
        # Extract the coefficients for COPD and ACO from the model parameters
        coef_copd_control = model.params.get("C(Disease, Treatment(reference='Control'))[T.COPD]", 0)
        coef_aco_control = model.params.get("C(Disease, Treatment(reference='Control'))[T.ACO]", 0)
        
        # Calculate the contrast as the difference in coefficients
        contrast = coef_copd_control - coef_aco_control

        # Extract the standard errors for these coefficients from the model
        se_copd_control = model.bse.get("C(Disease, Treatment(reference='Control'))[T.COPD]", 0)
        se_aco_control = model.bse.get("C(Disease, Treatment(reference='Control'))[T.ACO]", 0)

        # Compute the standard error of the contrast
        # Assumes that the coefficients are independent
        se_contrast = np.sqrt(se_copd_control**2 + se_aco_control**2)

        # Calculate the t-statistic for the contrast
        t_statistic = contrast / se_contrast

        # Retrieve the degrees of freedom from the model
        doff = model.df_resid

        # Calculate the p-value associated with the t-statistic
        # This is a two-tailed test
        p_value = 2 * (1 - stats.t.cdf(np.abs(t_statistic), doff))

        # logging.info(f"Calculated contrast: {contrast}, p-value: {p_value}")
        return p_value, contrast
    except Exception as e:
        logging.error(f"Error in calculating contrast: {e}")
        raise



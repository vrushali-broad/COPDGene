import logging
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from visualization_helpers import forest_plot, get_significance_colors, get_significance_stars  

def correct_ct_scan_features(pheno, outcomes, names, dir_path, formula, order_names, rename):
    """
    Adjusts CT scan values for confounders using OLS regression.

    Args:
        pheno (pd.DataFrame): DataFrame containing phenotypic data including CT scan values.
        outcomes (list): List of outcome variables (CT scan features) to adjust.
        names (list): Names corresponding to each outcome for plotting.
        dir_path (str): Directory path to save the plots.
        formula (str): Formula for the OLS regression model.
        order_names (list): List of names defining the order of factors in the model.
        rename (dict): Dictionary for renaming factors in the model.

    Returns:
        dict: Dictionary containing regression results for each outcome.
    """
    # Create DataFrame from pheno and add intercept
    df = pd.DataFrame(pheno)
    df['intercept'] = 1

    # Initialize dictionary to store regression results
    results = {}
    pvalues_lst = {}

    # Iterate over each outcome to perform OLS regression
    for i, outcome in enumerate(outcomes):
        try:
            # Fit the model
            current_formula = formula.format(outcome=outcome)
            model = smf.ols(current_formula, data=df).fit()

            # Store residuals and fitted values
            df['corrected_' + outcome] = model.resid
            df['Adjusted_' + outcome] = model.predict(df)

            # Store model in results
            results[outcome] = model

            # Plotting and saving results
            pvalues = plot_adjusted_results(model, names[i], order_names, rename, dir_path, outcome)
            pvalues_lst[outcome] = pvalues
            logging.info(f"Adjusted and plotted outcome: {outcome}")
            
        except Exception as e:
            logging.error(f"Error in adjusting {outcome}: {e}")

    return df, pvalues_lst

def plot_adjusted_results(model, title, order_names, rename, dir_path, outcome):
    """
    Plots the adjusted results, including coefficients and confidence intervals, from the OLS regression model.

    This function generates a forest plot to visually represent the regression coefficients along with their 
    confidence intervals and significance. The plot helps in understanding the impact of each variable in the model.

    Args:
        model (RegressionResults): The fitted regression model from statsmodels.
        title (str): Title of the plot.
        order_names (list): Order of the factors as they should appear in the plot.
        rename (dict): Mapping for renaming the factors for better readability in the plot.
        dir_path (str): Directory path where the plot will be saved.
        outcome (str): Name of the outcome variable for which the model was fitted.

    Returns:
        list: List of p-values associated with each coefficient in the plot.

    Raises:
        Exception: If any error occurs during the plotting process.

    Note:
        This function assumes that the 'Intercept' is not required in the plot and removes it from the coefficients.
    """
    try:
        # Extract coefficients, p-values, and confidence intervals from the model
        params = model.params
        pvalues = model.pvalues
        conf = model.conf_int()
        conf['Estimate'] = params
        conf.columns = ['Lower', 'Upper', 'Estimate']

        # Remove the 'Intercept' if present, as it's usually not needed in the plot
        conf = conf.drop("Intercept", errors='ignore')
        pvalues = pvalues.drop("Intercept", errors='ignore')

        # Determine the color coding for significance of each coefficient
        colors = get_significance_colors(pvalues, conf['Estimate'])

        # File path for saving the plot
        fig_name = f"{dir_path}/Adj_{outcome}.pdf"

        # Create and save the forest plot
        # forest_plot function will internally handle the plotting and saving of the figure
        pvalues = forest_plot(conf.index, conf['Estimate'], conf['Lower'], conf['Upper'],
                                     colors, pvalues=pvalues, order=order_names, rename=rename,
                                     title=title, y=42, file_path=fig_name)

        logging.info(f"Forest plot created for {outcome}")
        return pvalues
    
    except Exception as e:
        logging.error(f"Error occurred while creating forest plot for {outcome}: {e}")
        raise

# def perform_ols_regression_and_contrast(df, outcomes, formula, names, order_names, rename, dir_path):
def correct_contrast_ct_scan_features(pheno, outcomes, names, dir_path, formula, order_names, rename, contrast = True):
    """
    Perform OLS regression for each outcome variable, calculate contrasts, and create forest plots.

    Args:
        pheno (pd.DataFrame): DataFrame containing the pheno data.
        outcomes (list): List of outcome variables for regression.
        formula (str): Template formula for OLS regression.
        names (list): Names corresponding to each outcome for plotting.
        order_names (list): Order of factors in the model.
        rename (dict): Dictionary for renaming factors in the model.
        dir_path (str): Directory path to save the forest plots.

    Returns:
        dict: Dictionary containing regression results and additional computed values.
    """
    # Create DataFrame from pheno and add intercept
    df = pd.DataFrame(pheno)
    df['intercept'] = 1
    
    results = {}
    pvalues_lst = {}
    
    for i, outcome in enumerate(outcomes):
        try:
            # Perform OLS regression
            current_formula = formula.format(outcome=outcome)
            model = smf.ols(current_formula, data=df).fit()

            # Calculate contrast between COPD and ACO
            # contrast_pvalue, contrast = calculate_contrast(model)

            # Update DataFrame with corrected and adjusted values
            df['corrected_' + outcome] = model.resid
            df['Adjusted_' + outcome] = model.fittedvalues

            # Extract coefficients, confidence intervals, and p-values
            params = model.params
            conf = model.conf_int()
            conf['Estimate'] = params
            conf.columns = ['Lower', 'Upper', 'Estimate']
            conf = conf.drop("Intercept", errors='ignore')
            pvalues = model.pvalues.drop("Intercept", errors='ignore')

            # Prepare data for forest plot
            labels = conf.index
            estimates = conf['Estimate']
            lower = conf['Lower']
            upper = conf['Upper']
            colors = get_significance_colors(pvalues, estimates)

            # Generate forest plot
            file_path = f"{dir_path}/Forest_Adj_{outcome}.pdf"
            forest_pvalues = forest_plot(labels, estimates, lower, upper, colors, pvalues, order=order_names, rename=rename, title=names[i], file_path=file_path)

        
            outcome = 'Adjusted_'+outcome
            pvalues_lst[outcome] = forest_pvalues
            
            # Calculate contrast between COPD and ACO
            # Calculate and store the contrast p-value if contrast analysis is requested
            if contrast:
                contrast_pvalue, _ = calculate_contrast(model)
                pvalues_lst[outcome].append(contrast_pvalue)
            
            logging.info(f"Processed and plotted outcome: {outcome}")

        except Exception as e:
            logging.error(f"Error in processing outcome {outcome}: {e}")
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


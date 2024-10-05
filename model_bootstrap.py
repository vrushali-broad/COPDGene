import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
import os
import matplotlib.pyplot as plt
from visualization_helpers import forest_plot, get_significance_colors, get_significance_stars  
from assumption_checks_regression import check_assumptions
from utils import save_kde_ecdf_plot
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import t, studentized_range

def pvalue_to_stars(pvalue):
    if pvalue < 0.0001:
        return '****'
    if pvalue < 0.001:
        return '***'
    elif pvalue < 0.01:
        return '**'
    elif pvalue < 0.05:
        return '*'
    else:
        return 'n.s.'
    
def bootstrap_regression_estimates_hc3(data, formula, n_bootstraps=1000, seed=42):
    """
    Performs bootstrap resampling with replacement to estimate regression coefficients,
    standard errors, confidence intervals, and residuals using HC3 robust standard errors.

    Parameters:
    - data (DataFrame): The input data frame.
    - formula (str): The regression formula.
    - n_bootstraps (int): Number of bootstrap iterations.
    - seed (int): Random seed for reproducibility.

    Returns:
    - results (DataFrame): DataFrame containing bootstrap estimates, standard errors, confidence intervals, and p-values.
    - mean_residuals (ndarray): Array of mean residuals from the bootstrap samples.
    """
    # Set the random seed for reproducibility
    np.random.seed(seed)

    # Fit the initial model with HC3 robust standard errors
    model = smf.ols(formula=formula, data=data).fit(cov_type='HC3')
    
    # Extract variable names
    variable_names = model.params.index

    # Initialize arrays for bootstrap estimates and residuals
    bootstrap_estimates = np.zeros((n_bootstraps, len(model.params)))

    # Perform bootstrap resampling
    for i in range(n_bootstraps):
        # Sample with replacement
        bootstrap_sample = data.sample(n=len(data), replace=True, random_state=seed + i)
        
        # Fit the model to the bootstrap sample with HC3 robust standard errors
        bootstrap_model = smf.ols(formula=formula, data=bootstrap_sample).fit(cov_type='HC3')
        
        # Align the parameters to the original model's design matrix
        bootstrap_params = pd.Series(bootstrap_model.params, index=variable_names)#.fillna(0)
        
        # Store the parameter estimates
        bootstrap_estimates[i, :] = bootstrap_params
        
    ################################################################################     
    # Convert to DataFrame to handle NaNs
    bootstrap_estimates = pd.DataFrame(bootstrap_estimates, columns=variable_names)
    # print('>>>>>', bootstrap_estimates)
    # Impute NaNs with the mean of the non-NaN values for each coefficient
    bootstrap_estimates = bootstrap_estimates.apply(lambda col: col.fillna(col.mean()), axis=0)
    # Impute NaNs with the median of the non-NaN values for each coefficient
    # bootstrap_estimates_df = bootstrap_estimates_df.apply(lambda col: col.fillna(col.median()), axis=0)
    ################################################################################ 

    # Calculate bootstrap standard errors, means, and confidence intervals
    bootstrap_means = np.mean(bootstrap_estimates, axis=0)
    bootstrap_se = np.std(bootstrap_estimates, axis=0)
    bootstrap_ci_lower = np.percentile(bootstrap_estimates, 2.5, axis=0)
    bootstrap_ci_upper = np.percentile(bootstrap_estimates, 97.5, axis=0)
    
    # Calculate p-values based on bootstrap standard errors
    z_scores = bootstrap_means / bootstrap_se
    p_values = 2 * (1 - norm.cdf(np.abs(z_scores)))

    # Create a results DataFrame
    results = pd.DataFrame({
        'Estimate': bootstrap_means,
        'Bootstrap SE': bootstrap_se,
        'Bootstrap 95% CI Lower': bootstrap_ci_lower,
        'Bootstrap 95% CI Upper': bootstrap_ci_upper,
        'p-value': p_values
    }, index=variable_names)

    # Calculate residuals for original data using mean estimates
    X = model.model.exog
    y = model.model.endog
    mean_fitted_values = np.dot(X, bootstrap_means)
    mean_residuals = y - mean_fitted_values
    
    used_indices = model.model.data.row_labels
    
    # Create a DataFrame for raw bootstrap estimates
    bootstrap_estimates_df = pd.DataFrame(bootstrap_estimates, columns=variable_names)
    
    return results, mean_residuals, mean_fitted_values, bootstrap_estimates_df, used_indices

def calculate_all_contrasts(model, results_df, bootstrap_samples_df, requested_contrasts):
    """
    Calculate all possible contrasts between specified groups in a regression model and return the requested ones.

    Parameters:
    - results_df (DataFrame): DataFrame containing regression estimates and statistics.
    - requested_contrasts (list of str): List of requested contrasts in the format 'Group1_vs_Group2'.

    Returns:
    - contrast_pvalues (Series): Series containing p-values for each requested contrast.
    - contrast_estimates (Series): Series containing contrast estimates.
    """
    # Define all possible contrasts
    groups = ['Asthma', 'COPD', 'ACO', 'Control']
    all_contrasts = requested_contrasts # [(g1, g2) for i, g1 in enumerate(groups) for g2 in groups[i+1:]]

    # Initialize dictionaries to store all possible contrast results
    all_contrast_estimates = {}
    all_contrast_se = {}
    all_contrast_t_statistic = {}
    all_contrast_p_values = {}
    
    # Calculate covariance matrix from bootstrap samples
    cov_matrix = bootstrap_samples_df.cov()
    
    # # Calculate covariance matrix from bootstrap samples
    # cov_matrix = results_df['Estimate'].cov()

    # Calculate all possible contrasts
    # for group1, group2 in all_contrasts:
    for contrast in requested_contrasts:
        group1, group2 = contrast.split('_vs_')
        if group2 == 'Control':
            var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
            contrast_estimate = results_df.loc[var1, 'Estimate']
            se = results_df.loc[var1, 'Bootstrap SE']
        elif group1 == 'Control':
            var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
            contrast_estimate = -results_df.loc[var2, 'Estimate']
            se = results_df.loc[var2, 'Bootstrap SE']
        else:
            var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
            var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
            contrast_estimate = results_df.loc[var1, 'Estimate'] - results_df.loc[var2, 'Estimate']
            # Assuming independence between coefficients
            se = np.sqrt(results_df.loc[var1, 'Bootstrap SE']**2 + results_df.loc[var2, 'Bootstrap SE']**2)
            # Assuming coefficients are not independent
            # se = np.sqrt(results_df.loc[var1, 'Bootstrap SE']**2 + results_df.loc[var2, 'Bootstrap SE']**2 - 2 * cov_matrix.loc[var1, var2])
            
        contrast_se_value = se
        t_statistic = contrast_estimate / contrast_se_value
        # Two-sided
        # p_value = 2 * (1 - norm.cdf(np.abs(t_statistic)))
        # Two-sided
        # p_value = 2 * norm.cdf(-np.abs(t_statistic))        
        # # Two-sided
        # Extract degrees of freedom from the model
        doff = model.df_resid ## ddof of the bootstrapped model is the same as the original model (the #observations and #parameters are the same)
        p_value = 2 * (1 - stats.t.cdf(np.abs(t_statistic), doff))
        
        # # Calculate p-value using the studentized range distribution
        # q = t_statistic * np.sqrt(2)
        # p_value = studentized_range.cdf(q, num_groups, doff, loc=0, scale=1)
        # p_value = 2 * (1 - p_value)  # Convert to two-sided p-value

        # Store the results
        contrast_key = f"{group1}_vs_{group2}"
        all_contrast_estimates[contrast_key] = contrast_estimate
        all_contrast_se[contrast_key] = contrast_se_value
        all_contrast_t_statistic[contrast_key] = t_statistic
        all_contrast_p_values[contrast_key] = p_value
        
    # Filter the requested contrasts
    contrast_estimates = pd.Series({k: all_contrast_estimates[k] for k in requested_contrasts if k in all_contrast_estimates})
    contrast_pvalues = pd.Series({k: all_contrast_p_values[k] for k in requested_contrasts if k in all_contrast_p_values})

    return contrast_pvalues, contrast_estimates

def perform_regression_analysis_bootstrap(pheno, outcomes, names, order_names, rename, formula, dir_name, contrast=False, ct=None):
    df = pd.DataFrame(pheno)
    df['Intercept'] = 1

    results = {}
    contrast_results = {}
    bootstrap_samples_dict = {}
    models = {}
    tukeys = {}
    
    for i, outcome in enumerate(outcomes):
        
        save_kde_ecdf_plot(df, outcome, names[i], os.path.join(dir_name, f'KDE_{outcome}.pdf'))

        current_formula = formula.format(outcome=outcome)
        model = smf.ols(formula=current_formula, data=df).fit(cov_type='HC3')
        
        # infl = model.get_influence()
        # print(infl.summary_table())
        
        # cov_matrix = model.cov_params()
        
        # df['corrected_' + outcome] = model.resid
        # df['Adjusted_' + outcome] = model.fittedvalues
        # results[outcome] = model
        models[outcome] = model

        output_dir = os.path.join('Figures', 'Assumption_Checks')
        residuals = model.resid
        fitted_values = model.fittedvalues
        check_assumptions(model, residuals, fitted_values, outcome, output_dir, figure_title = names[i])
        
        # Perform bootstrapping with HC3
        bootstrap_results, bootstrapped_residuals, bootstrapped_fittedvalues, bootstrap_estimates, used_indices = bootstrap_regression_estimates_hc3(df, current_formula, n_bootstraps=1000, seed=42)
        # bootstrap_results, bootstrapped_residuals, bootstrapped_fittedvalues, bootstrap_estimates = bca_bootstrap(df, current_formula, n_bootstraps=1000, seed=42)
        results[outcome] = bootstrap_results  # Update results with bootstrap results
        # Store bootstrap samples for this outcome
        bootstrap_samples_dict[outcome] = bootstrap_estimates
        
        df.loc[used_indices,'corrected_' + outcome] = bootstrapped_residuals
        df.loc[used_indices,'Adjusted_' + outcome] = bootstrapped_fittedvalues
        
#         #########>>>>>>>>>>>>>>>>>>################
        # Set pandas options to display full DataFrame
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.expand_frame_repr', False)
        pd.set_option('display.max_colwidth', None)
        
        df.loc[used_indices, 'bootstrapped_fitted_values'] = bootstrapped_fittedvalues
        # Remove rows with NaNs in the relevant columns
        df = df.dropna(subset=['bootstrapped_fitted_values', 'Disease'])

        # Perform ANOVA
        formula_anova = f"bootstrapped_fitted_values ~ C(Disease, Treatment(reference='Control'))"
        anova_model = smf.ols(formula_anova, data=df.loc[used_indices]).fit() 
        anova_results  = sm.stats.anova_lm(anova_model, typ=2)
        print(f'ANOVA results for {outcome}:')
        print('*****', anova_results , '******')
        
        # If ANOVA results indicate significant differences, perform Tukey's HSD test
        if anova_results['PR(>F)'][0] < 0.05:
            tukey = pairwise_tukeyhsd(endog=df['bootstrapped_fitted_values'], groups=df['Disease'], alpha=0.05)

            # Convert Tukey results to a DataFrame for better readability
            tukey_results_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
            
            # Set a tolerance level for considering values as zero
            tolerance = 1e-10
            # Apply the tolerance: any value within [-tolerance, tolerance] is set to zero
            tukey_results_df['p-adj'] = tukey_results_df['p-adj'].apply(lambda x: 0.0 if abs(x) < tolerance else x)

            # # Convert p-adj to float for accurate inspection and format to 30 decimal places
            # tukey_results_df['p-adj'] = tukey_results_df['p-adj'].astype(float).map(lambda x: f"{x:.5f}")

            # Print Tukey's HSD results
            print('Tukey HSD results:')
            print(tukey_results_df)
            
            # # Check for negative p-values
            # negative_pvalues = tukey_results_df[tukey_results_df['p-adj'] < 0]
            # print(negative_pvalues, '>>>>>>>>>')

            # Create a custom index for the p-values
            index = [f"{row['group1']}_vs_{row['group2']}" for _, row in tukey_results_df.iterrows()]
            tukeys[outcome] = pd.Series(tukey_results_df['p-adj'].values, index=index) 
            # print(tukeys[outcome])
            stars = tukeys[outcome].apply(pvalue_to_stars)
            stars_series = pd.Series(stars, index=tukeys[outcome].index, name='P-value')
            print('######################', outcome,'\n',stars_series,'#########################')
        else:
            print('No significant differences found in ANOVA. Tukey HSD test will not be performed.')

#         # Calculate pairwise p-values for differences using Tukey's HSD test
#         mc_results = pairwise_tukeyhsd(endog=df.loc[used_indices, 'bootstrapped_fitted_values'], groups=df.loc[used_indices,'Disease'], alpha=0.05)
#         # Convert results to DataFrame and display p-values up to 5 decimal points
#         tukey_results_df = pd.DataFrame(data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
#         tukey_results_df['p-adj'] = tukey_results_df['p-adj'].astype(float).round(5)
#         print(f'Tukey HSD results for {outcome}:')
#         print(tukey_results_df)
        
#         # Create a custom index for the p-values
#         index = [f"{row['group1']}_vs_{row['group2']}" for _, row in tukey_results_df.iterrows()]
#         tukeys[outcome] = pd.Series(tukey_results_df['p-adj'].values, index=index) 
#         print(tukeys[outcome])
#         print('===========================================================')
        ###########################################

        # Check assumptions on bootstrapped estimates
        check_assumptions(model, bootstrapped_residuals, bootstrapped_fittedvalues, outcome, output_dir, names[i], prefix='bootstrapped_')

    pvalues_lst = {}
    for i, (key, result) in enumerate(results.items()):
        params = result['Estimate']
        pvalues = result['p-value']
        conf = result[['Bootstrap 95% CI Lower', 'Bootstrap 95% CI Upper', 'Estimate']]
        conf.columns = ['Lower', 'Upper', 'Estimate']
        
        if "Intercept" in conf.index:
            conf = conf.drop("Intercept")
        if "Intercept" in pvalues.index:
            pvalues = pvalues.drop("Intercept")
            
        corrected_pvalues = multipletests(pvalues, method='fdr_bh')[1] #method='bonferroni')[1] #

        colors = get_significance_colors(corrected_pvalues, conf['Estimate'])
        
        fig_name = os.path.join(dir_name, 'Forest_Adj_' + outcomes[i] + '.pdf')
        forest_pvalues = forest_plot(conf.index,
                                     conf['Estimate'], conf['Lower'], conf['Upper'], 
                                     colors,
                                     pvalues=corrected_pvalues, 
                                     order=order_names, 
                                     rename=rename, 
                                     title=names[i], 
                                     file_path=fig_name)
        
        adj_key = 'Adjusted_' + key
        pvalues_lst[adj_key] = corrected_pvalues 
    
        if contrast:
#             # Define the requested contrasts
#             requested_contrasts = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control', 'Asthma_vs_COPD', 'Asthma_vs_ACO', 'COPD_vs_ACO']
            
#             # Calculate contrasts based on bootstrap estimates
#             bootstrap_samples_df = bootstrap_samples_dict[key]
#             contrast_pvalues, contrasts = calculate_all_contrasts(models[key], result, bootstrap_samples_df, requested_contrasts)
#             corrected_contrast_pvalues = multipletests(contrast_pvalues, method='fdr_bh')[1] #method='bonferroni')[1] #
#             corrected_contrast_pvalues = pd.Series(corrected_contrast_pvalues, index=contrast_pvalues.index, name="corrected p-value")
#             # print('######################', adj_key, contrast_pvalues, corrected_contrast_pvalues,'#########################')

#             # Convert p-values to stars and save as a series
#             stars = corrected_contrast_pvalues.apply(pvalue_to_stars)
#             stars_series = pd.Series(stars, index=corrected_contrast_pvalues.index, name='Stars')
#             print('######################', adj_key,'\n',stars_series,'#########################')

#             pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
#             if ct:
#                 pairs_order.append('COPD_vs_ACO')
#             pvalues = [corrected_contrast_pvalues.get(pair, None) for pair in pairs_order]
            
            ## Tukey Pvalues
            pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
            if ct:
                pairs_order.append('ACO_vs_COPD')
            pvalues = [tukeys[key].get(pair, None) for pair in pairs_order]
            pvalues_lst[adj_key] = pvalues 

    return df, pvalues_lst

def calculate_copd_aco_contrast(results_df, bootstrap_samples_df):
    var_copd = "C(Disease, Treatment(reference='Control'))[T.COPD]"
    var_aco = "C(Disease, Treatment(reference='Control'))[T.ACO]"

    contrast_estimate = results_df.loc[var_aco, 'Estimate'] - results_df.loc[var_copd, 'Estimate']
    
    cov_matrix = bootstrap_samples_df.cov()
    se = np.sqrt(
        results_df.loc[var_aco, 'Bootstrap SE']**2 +
        results_df.loc[var_copd, 'Bootstrap SE']**2 -
        2 * cov_matrix.loc[var_aco, var_copd]
    )

    t_statistic = contrast_estimate / se
    p_value = 2 * (1 - norm.cdf(np.abs(t_statistic)))

    return p_value, contrast_estimate#, se, t_statistic, p_value



















# def perform_regression_analysis_bootstrap(pheno, outcomes, names, order_names, rename, formula, dir_name, contrast=False, ct=None):
#     df = pd.DataFrame(pheno)
#     df['Intercept'] = 1

#     results = {}
#     contrast_results = {}
#     bootstrap_samples_dict = {}
#     models = {}
#     tukeys = {}
    
#     for i, outcome in enumerate(outcomes):
        
#         save_kde_ecdf_plot(df, outcome, names[i], os.path.join(dir_name, f'KDE_{outcome}.pdf'))

#         current_formula = formula.format(outcome=outcome)
#         model = smf.ols(formula=current_formula, data=df).fit(cov_type='HC3')
       
#         models[outcome] = model

#         output_dir = os.path.join('Figures', 'Assumption_Checks')
#         residuals = model.resid
#         fitted_values = model.fittedvalues
#         check_assumptions(model, residuals, fitted_values, outcome, output_dir, figure_title = names[i])
        
        
        
#         # Add fitted values to the DataFrame
#         df['fitted_values'] = fitted_values
#         # Perform ANOVA
#         anova_model = smf.ols("fitted_values ~ C(Disease, Treatment(reference='Control'))", data=df).fit()
#         anova_table = sm.stats.anova_lm(anova_model, typ=2)
#         print('ANOVA results:')
#         print(anova_table)
#         # Perform Tukey's HSD test
#         df = df[['fitted_values', 'Disease']].dropna()
#         mc_results = pairwise_tukeyhsd(endog=df['fitted_values'], groups=df['Disease'], alpha=0.05)

#         # Convert to DataFrame
#         tukey_results_df = pd.DataFrame(data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
#         tukey_results_df['p-adj'] = tukey_results_df['p-adj'].astype(float).round(5)

#         # Check for any negative p-values
#         if (tukey_results_df['p-adj'] < 0).any():
#             raise ValueError("Negative p-values detected in Tukey HSD test results.")

#         print('Tukey HSD results:')
#         print(tukey_results_df)

        
        
        
        

#         # Perform bootstrapping with HC3 on clean data
#         bootstrap_results, bootstrapped_residuals, bootstrapped_fittedvalues, bootstrap_estimates, used_indices = bootstrap_regression_estimates_hc3(df, current_formula, n_bootstraps=10, seed=42)
#         results[outcome] = bootstrap_results  # Update results with bootstrap results
#         bootstrap_samples_dict[outcome] = bootstrap_estimates
        
#         df.loc[used_indices, 'corrected_' + outcome] = bootstrapped_residuals
#         df.loc[used_indices, 'Adjusted_' + outcome] = bootstrapped_fittedvalues
#         df.loc[used_indices, 'bootstrapped_fitted_values'] = bootstrapped_fittedvalues
        
# #         # Clean data before ANOVA
# #         clean_df = df.dropna(subset=['bootstrapped_fitted_values', 'Disease'])
  
# #         # Perform ANOVA on clean data
# #         anova_model = smf.ols("bootstrapped_fitted_values ~ C(Disease, Treatment(reference='Control'))", data=clean_df.loc[used_indices]).fit() 
# #         anova_table = sm.stats.anova_lm(anova_model, typ=2)
# #         print(f'ANOVA results for {outcome}:')
# #         print('*****', anova_table, '******')
        
# #         # Check for NaNs and off values before Tukey's HSD test
# #         print("Checking data for Tukey's HSD test...")
# #         print("Number of NaNs in 'bootstrapped_fitted_values':", clean_df.loc[used_indices,'bootstrapped_fitted_values'].isna().sum())
# #         print("Number of NaNs in 'Disease':", clean_df.loc[used_indices,'Disease'].isna().sum())
        
        
# #         # Ensure the data used for Tukey's test is clean
# #         clean_data = clean_df.loc[used_indices, ['bootstrapped_fitted_values', 'Disease']].dropna()
# #         # Perform Tukey's HSD test
# #         mc_results = pairwise_tukeyhsd(endog=clean_data['bootstrapped_fitted_values'], groups=clean_data['Disease'], alpha=0.05)
# #         # Check Tukey HSD results directly
# #         print(mc_results)

# #         # # Calculate pairwise p-values for differences using Tukey's HSD test
# #         # mc_results = pairwise_tukeyhsd(endog=clean_df.loc[used_indices, 'bootstrapped_fitted_values'], groups=clean_df.loc[used_indices,'Disease'], alpha=0.05)
# #         tukey_results_df = pd.DataFrame(data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
        
# #         # Ensure p-values are positive and round them
# #         tukey_results_df['p-adj'] = tukey_results_df['p-adj'].astype(float).round(5)
# #         print(f'Tukey HSD results for {outcome}:')
# #         print(tukey_results_df)
        
# #         index = [f"{row['group1']}_vs_{row['group2']}" for _, row in tukey_results_df.iterrows()]
# #         tukeys[outcome] = pd.Series(tukey_results_df['p-adj'].values, index=index) 
# #         print(tukeys[outcome])
# #         print('===========================================================')
        
#         check_assumptions(model, bootstrapped_residuals, bootstrapped_fittedvalues, outcome, output_dir, names[i], prefix='bootstrapped_')

#     pvalues_lst = {}
#     for i, (key, result) in enumerate(results.items()):
#         params = result['Estimate']
#         pvalues = result['p-value']
#         conf = result[['Bootstrap 95% CI Lower', 'Bootstrap 95% CI Upper', 'Estimate']]
#         conf.columns = ['Lower', 'Upper', 'Estimate']
        
#         if "Intercept" in conf.index:
#             conf = conf.drop("Intercept")
#         if "Intercept" in pvalues.index:
#             pvalues = pvalues.drop("Intercept")
            
#         corrected_pvalues = multipletests(pvalues, method='fdr_bh')[1] 

#         colors = get_significance_colors(corrected_pvalues, conf['Estimate'])
        
#         fig_name = os.path.join(dir_name, 'Forest_Adj_' + outcomes[i] + '.pdf')
#         forest_pvalues = forest_plot(conf.index,
#                                      conf['Estimate'], conf['Lower'], conf['Upper'], 
#                                      colors,
#                                      pvalues=corrected_pvalues, 
#                                      order=order_names, 
#                                      rename=rename, 
#                                      title=names[i], 
#                                      file_path=fig_name)
        
#         adj_key = 'Adjusted_' + key
#         pvalues_lst[adj_key] = corrected_pvalues 
    
#         if contrast:
#             pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
#             if ct:
#                 pairs_order.append('ACO_vs_COPD')

#             pvalues = [tukeys[key].get(pair, None) for pair in pairs_order]
#             pvalues_lst[adj_key] = pvalues 

#     return df, pvalues_lst




# def perform_regression_analysis_bootstrap_test(pheno, outcomes, names, order_names, rename, formula, dir_name, contrast=False, ct=None):
#     df = pd.DataFrame(pheno)
#     df['Intercept'] = 1

#     results = {}
#     contrast_results = {}
#     bootstrap_samples_dict = {}
#     models = {}
    
#     for i, outcome in enumerate(outcomes):
        
#         save_kde_ecdf_plot(df, outcome, names[i], os.path.join(dir_name, f'KDE_{outcome}.pdf'))

#         current_formula = formula.format(outcome=outcome)
#         model = smf.ols(formula=current_formula, data=df).fit(cov_type='HC3')
        
#         # cov_matrix = model.cov_params()

#         models[outcome] = model

#         output_dir = os.path.join('Figures', 'Assumption_Checks')
#         residuals = model.resid
#         fitted_values = model.fittedvalues
#         check_assumptions(model, residuals, fitted_values, outcome, output_dir, figure_title = names[i])
        
#         # Perform bootstrapping with HC3
#         bootstrap_results, bootstrapped_residuals, bootstrapped_fittedvalues, bootstrap_estimates, used_indices = bootstrap_regression_estimates_hc3(df, current_formula, n_bootstraps=1000, seed=42)
#         # bootstrap_results, bootstrapped_residuals, bootstrapped_fittedvalues, bootstrap_estimates = bca_bootstrap(df, current_formula, n_bootstraps=1000, seed=42)
#         results[outcome] = bootstrap_results  # Update results with bootstrap results
#         # Store bootstrap samples for this outcome
#         bootstrap_samples_dict[outcome] = bootstrap_estimates
        
#         #########>>>>>>>>>>>>>>>>>>################
#         df.loc[used_indices,'corrected_' + outcome] = bootstrapped_residuals
#         df.loc[used_indices,'Adjusted_' + outcome] = bootstrapped_fittedvalues
#         # Calculate pairwise p-values for differences using Tukey's HSD test
#         df.loc[used_indices, 'bootstrapped_fitted_values'] = bootstrapped_fittedvalues
#         mc_results = pairwise_tukeyhsd(endog=df.loc[used_indices, 'bootstrapped_fitted_values'], groups=df.loc[used_indices,'Disease'], alpha=0.05)
#         # Convert results to DataFrame and display p-values up to 5 decimal points
#         tukey_results_df = pd.DataFrame(data=mc_results._results_table.data[1:], columns=mc_results._results_table.data[0])
#         tukey_results_df['p-adj'] = tukey_results_df['p-adj'].astype(float).round(5)
#         print(f'Tukey HSD results for {outcome}:')
#         print(tukey_results_df)
        
        
#         # Perform ANOVA
#         anova_model = smf.ols(f'bootstrapped_fitted_values ~ C(Disease)', data=df.loc[used_indices]).fit()
#         anova_table = sm.stats.anova_lm(anova_model, typ=2)
#         print(f'ANOVA results for {outcome}:')
#         print(anova_table)
#         print('===========================================================')
#         ###########################################

#         # Check assumptions on bootstrapped estimates
#         check_assumptions(model, bootstrapped_residuals, bootstrapped_fittedvalues, outcome, output_dir, names[i], prefix='bootstrapped_')

#     pvalues_lst = {}
#     for i, (key, result) in enumerate(results.items()):
#         params = result['Estimate']
#         pvalues = result['p-value']
#         conf = result[['Bootstrap 95% CI Lower', 'Bootstrap 95% CI Upper', 'Estimate']]
#         # conf = result[['BCa 95% CI Lower', 'BCa 95% CI Upper', 'Estimate']]
#         conf.columns = ['Lower', 'Upper', 'Estimate']
        
#         if "Intercept" in conf.index:
#             conf = conf.drop("Intercept")
#         if "Intercept" in pvalues.index:
#             pvalues = pvalues.drop("Intercept")
            
#         corrected_pvalues = multipletests(pvalues, method='fdr_bh')[1] #method='bonferroni')[1] #

#         colors = get_significance_colors(corrected_pvalues, conf['Estimate'])
        
#         fig_name = os.path.join(dir_name, 'Forest_Adj_' + outcomes[i] + '.pdf')
#         forest_pvalues = forest_plot(conf.index,
#                                      conf['Estimate'], conf['Lower'], conf['Upper'], 
#                                      colors,
#                                      pvalues=corrected_pvalues, 
#                                      order=order_names, 
#                                      rename=rename, 
#                                      title=names[i], 
#                                      file_path=fig_name)
        
#         adj_key = 'Adjusted_' + key
#         pvalues_lst[adj_key] = corrected_pvalues 
    
#         if contrast:
#             # Define the requested contrasts
#             # requested_contrasts = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control', 'Asthma_vs_COPD', 'Asthma_vs_ACO', 'COPD_vs_ACO']
#             # requested_contrasts = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control', 'COPD_vs_ACO']
#             requested_contrasts = [ 'COPD_vs_ACO'] #'Asthma_vs_COPD', 'Asthma_vs_ACO',
            
#             # Calculate contrasts based on bootstrap estimates
#             bootstrap_samples_df = bootstrap_samples_dict[key]
#             contrast_pvalues, contrasts = calculate_all_contrasts(models[key], result, bootstrap_samples_df, requested_contrasts)
            
            
            
            
            
#             # =================================================== #
#             print('+++++', contrast_pvalues)
# #             # Mapping dictionary for the desired format
# #             mapping = {
# #                 "C(Disease, Treatment(reference='Control'))[T.ACO]": "ACO_vs_Control",
# #                 "C(Disease, Treatment(reference='Control'))[T.Asthma]": "Asthma_vs_Control",
# #                 "C(Disease, Treatment(reference='Control'))[T.COPD]": "COPD_vs_Control"
# #             }
# #             # Combine pvalues from model and contrasts
# #             contrast_pvalues = pd.concat([contrast_pvalues, pvalues])
# #             contrast_pvalues = contrast_pvalues.rename(index=mapping)
# #             # Adding Asthma_vs_COPD and Asthma_vs_ACO with p-value 0
# #             contrast_pvalues['Asthma_vs_COPD'] = 0
# #             contrast_pvalues['Asthma_vs_ACO'] = 0
# #             # =================================================== #
            
            
            
#             corrected_contrast_pvalues = multipletests(contrast_pvalues, method='fdr_bh')[1] #method='bonferroni')[1] #
#             corrected_contrast_pvalues = pd.Series(corrected_contrast_pvalues, index=contrast_pvalues.index, name="corrected p-value")
#             # print('######################', adj_key, contrast_pvalues, corrected_contrast_pvalues,'#########################')
            
            
#             # Convert p-values to stars and save as a series
#             stars = corrected_contrast_pvalues.apply(pvalue_to_stars)
#             stars_series = pd.Series(stars, index=corrected_contrast_pvalues.index, name='Stars')
#             print('######################', adj_key,'\n',stars_series,'#########################')

            
#             pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
#             if ct:
#                 pairs_order.append('COPD_vs_ACO')

#             pvalues = [corrected_contrast_pvalues.get(pair, None) for pair in pairs_order]
#             pvalues_lst[adj_key] = pvalues 

#     return df, pvalues_lst



# def bca_bootstrap(data, formula, n_bootstraps=1000, seed=42):
#     """
#     Performs BCa bootstrap resampling to estimate regression coefficients,
#     standard errors, confidence intervals, and p-values.

#     Parameters:
#     - data (DataFrame): The input data frame.
#     - formula (str): The regression formula.
#     - n_bootstraps (int): Number of bootstrap iterations.
#     - seed (int): Random seed for reproducibility.

#     Returns:
#     - results (DataFrame): DataFrame containing bootstrap estimates, standard errors, confidence intervals, and p-values.
#     - mean_residuals (ndarray): Array of mean residuals from the bootstrap samples.
#     - mean_fitted_values (ndarray): Array of mean fitted values from the bootstrap samples.
#     - bootstrap_estimates_df (DataFrame): DataFrame containing raw bootstrap estimates for each parameter.
#     """
#     np.random.seed(seed)
#     model = smf.ols(formula=formula, data=data).fit(cov_type='HC3')
#     variable_names = model.params.index
#     bootstrap_estimates = np.zeros((n_bootstraps, len(model.params)))

#     for i in range(n_bootstraps):
#         bootstrap_sample = data.sample(n=len(data), replace=True, random_state=seed + i)
#         bootstrap_model = smf.ols(formula=formula, data=bootstrap_sample).fit(cov_type='HC3')
#         bootstrap_params = pd.Series(bootstrap_model.params, index=variable_names).fillna(0)
#         bootstrap_estimates[i, :] = bootstrap_params

#     bootstrap_means = np.mean(bootstrap_estimates, axis=0)
#     bootstrap_se = np.std(bootstrap_estimates, axis=0)
#     bootstrap_ci_lower = np.percentile(bootstrap_estimates, 2.5, axis=0)
#     bootstrap_ci_upper = np.percentile(bootstrap_estimates, 97.5, axis=0)

#     # BCa adjusted intervals
#     alpha = 0.05
#     z0 = norm.ppf((bootstrap_estimates < bootstrap_means).mean(axis=0))
#     jacks = np.mean(bootstrap_estimates, axis=0)
#     accel = np.mean((bootstrap_estimates - jacks) ** 3, axis=0) / (6 * np.mean((bootstrap_estimates - jacks) ** 2, axis=0) ** 1.5)

#     lower_percentile = norm.cdf(z0 + (z0 + norm.ppf(alpha / 2)) / (1 - accel * (z0 + norm.ppf(alpha / 2))))
#     upper_percentile = norm.cdf(z0 + (z0 + norm.ppf(1 - alpha / 2)) / (1 - accel * (z0 + norm.ppf(1 - alpha / 2))))

#     bca_ci_lower = np.percentile(bootstrap_estimates, (lower_percentile * 100).flatten(), axis=0)
#     bca_ci_upper = np.percentile(bootstrap_estimates, (upper_percentile * 100).flatten(), axis=0)

#     results = pd.DataFrame({
#         'Estimate': bootstrap_means,
#         'Bootstrap SE': bootstrap_se,
#         'Bootstrap 95% CI Lower': bootstrap_ci_lower,
#         'Bootstrap 95% CI Upper': bootstrap_ci_upper,
#         'BCa 95% CI Lower': bca_ci_lower,
#         'BCa 95% CI Upper': bca_ci_upper
#     }, index=variable_names)

#     # Calculate residuals for original data using mean estimates
#     X = model.model.exog
#     y = model.model.endog
#     mean_fitted_values = np.dot(X, bootstrap_means)
#     mean_residuals = y - mean_fitted_values
    
#     # Create a DataFrame for raw bootstrap estimates
#     bootstrap_estimates_df = pd.DataFrame(bootstrap_estimates, columns=variable_names)

#     return results, mean_residuals, mean_fitted_values, bootstrap_estimates_df


        
    
# from bootstrapped import bootstrap as bs
# from bootstrapped import stats_functions as bs_stats
# from scipy.stats import norm
# def bca_bootstrap(data, formula, n_bootstraps=1000, seed=42):
#     """
#     Performs BCa bootstrap resampling to estimate regression coefficients,
#     standard errors, confidence intervals, and p-values.

#     Parameters:
#     - data (DataFrame): The input data frame.
#     - formula (str): The regression formula.
#     - n_bootstraps (int): Number of bootstrap iterations.
#     - seed (int): Random seed for reproducibility.

#     Returns:
#     - results (DataFrame): DataFrame containing bootstrap estimates, standard errors, confidence intervals, and p-values.
#     """
#     np.random.seed(seed)
#     model = smf.ols(formula=formula, data=data).fit(cov_type='HC3')
#     variable_names = model.params.index

#     bootstrap_estimates = np.zeros((n_bootstraps, len(model.params)))

#     for i in range(n_bootstraps):
#         bootstrap_sample = data.sample(n=len(data), replace=True, random_state=seed + i)
#         bootstrap_model = smf.ols(formula=formula, data=bootstrap_sample).fit(cov_type='HC3')
#         bootstrap_params = pd.Series(bootstrap_model.params, index=variable_names).fillna(0)
#         bootstrap_estimates[i, :] = bootstrap_params

#     bootstrap_means = np.mean(bootstrap_estimates, axis=0)
#     bootstrap_se = np.std(bootstrap_estimates, axis=0)
#     bootstrap_ci_lower = np.percentile(bootstrap_estimates, 2.5, axis=0)
#     bootstrap_ci_upper = np.percentile(bootstrap_estimates, 97.5, axis=0)

#     results = pd.DataFrame({
#         'Estimate': bootstrap_means,
#         'Bootstrap SE': bootstrap_se,
#         'Bootstrap 95% CI Lower': bootstrap_ci_lower,
#         'Bootstrap 95% CI Upper': bootstrap_ci_upper
#     }, index=variable_names)

#     # BCa adjusted intervals
#     for i, var in enumerate(variable_names):
#         results.loc[var, 'BCa 95% CI Lower'], results.loc[var, 'BCa 95% CI Upper'] = bs.bootstrap_ab(data[var], stat_func=bs_stats.mean, is_pivotal=False, return_distribution=False, alpha=0.05)

#     return results   
    
    
    
    
    
    
    
    
    
    
    
    




# ################################################################################################
# def perform_regression_analysis_bootstrap(pheno, outcomes, names, order_names, rename, formula, dir_name, contrast=False, ct=None):
#     df = pd.DataFrame(pheno)
#     df['Intercept'] = 1

#     results = {}
#     contrast_results = {}
    
#     for i, outcome in enumerate(outcomes):
        
#         save_kde_ecdf_plot(df, outcome, names[i], os.path.join(dir_name, f'KDE_{outcome}.pdf'))

#         current_formula = formula.format(outcome=outcome)
#         model = smf.ols(formula=current_formula, data=df).fit(cov_type='HC3')
        
#         df['corrected_' + outcome] = model.resid
#         df['Adjusted_' + outcome] = model.fittedvalues
#         results[outcome] = model

#         output_dir = os.path.join('Figures', 'Assumption_Checks')
#         residuals = model.resid
#         fitted_values = model.fittedvalues
#         check_assumptions(model, residuals, fitted_values, outcome, output_dir)
        
#         bootstrap_results, bootstrapped_residuals, bootstrapped_fittedvalues = bootstrap_regression_estimates_hc3(df, current_formula, n_bootstraps=100, seed=42)
#         results[outcome] = bootstrap_results

#         check_assumptions(model, bootstrapped_residuals, bootstrapped_fittedvalues, outcome, output_dir, prefix='bootstrapped_')

#     pvalues_lst = {}
#     for i, (key, result) in enumerate(results.items()):
#         params = result['Estimate']
#         pvalues = result['p-value']
#         conf = result[['Bootstrap 95% CI Lower', 'Bootstrap 95% CI Upper', 'Estimate']]
#         conf.columns = ['Lower', 'Upper', 'Estimate']
        
#         if "Intercept" in conf.index:
#             conf = conf.drop("Intercept")
#         if "Intercept" in pvalues.index:
#             pvalues = pvalues.drop("Intercept")
            
#         corrected_pvalues = multipletests(pvalues, method='fdr_bh')[1]

#         colors = get_significance_colors(corrected_pvalues, conf['Estimate'])
        
#         fig_name = os.path.join(dir_name, 'Forest_Adj_' + outcomes[i] + '.pdf')
#         forest_pvalues = forest_plot(conf.index,
#                                      conf['Estimate'], conf['Lower'], conf['Upper'], 
#                                      colors,
#                                      pvalues=corrected_pvalues, 
#                                      order=order_names, 
#                                      rename=rename, 
#                                      title=names[i], 
#                                      file_path=fig_name)
        
#         key = 'Adjusted_' + key
#         pvalues_lst[key] = corrected_pvalues 
    
#         if contrast:
#             contrast_pairs = [('Asthma', 'Control'), ('COPD', 'Control'), ('ACO', 'Control')]
#             if ct:
#                 contrast_pairs.append(('COPD', 'ACO'))

#             contrast_pvalues, contrasts = calculate_all_contrasts(model, contrast_pairs)
#             corrected_contrast_pvalues = multipletests(contrast_pvalues, method='fdr_bh')[1]
#             corrected_contrast_pvalues = pd.Series(corrected_contrast_pvalues, index=contrast_pvalues.index, name="corrected p-value")
            
#             pairs_order = ['Asthma_vs_Control', 'COPD_vs_Control', 'ACO_vs_Control']
#             if ct:
#                 pairs_order.append('COPD_vs_ACO')

#             pvalues = [corrected_contrast_pvalues.get(pair, None) for pair in pairs_order]
#             pvalues_lst[key] = pvalues 

#     return df, pvalues_lst

# def calculate_all_contrasts(model, contrasts):
#     """
#     Calculate contrasts between specified groups in a regression model.

#     Parameters:
#     - model: A fitted regression model (e.g., statsmodels OLS result).
#     - contrasts (list of tuples): Each tuple contains the names of the groups to be contrasted.

#     Returns:
#     - contrast_pvalues (Series): Series containing p-values for each contrast.
#     - contrast_estimates (Series): Series containing contrast estimates.
#     """
#     # Extract coefficients and covariance matrix from the model
#     params = model.params
#     cov_matrix = model.cov_params()

#     # Initialize lists to store contrast results
#     contrast_estimates = []
#     contrast_se = []
#     contrast_z_scores = []
#     contrast_p_values = []

#     # Calculate contrasts
#     for group1, group2 in contrasts:
#         if group2 == 'Control':
#             var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
#             contrast_estimate = params[var1]
#             se = cov_matrix.loc[var1, var1]
#         elif group1 == 'Control':
#             var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
#             contrast_estimate = -params[var2]
#             se = cov_matrix.loc[var2, var2]
#         else:
#             var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
#             var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
#             contrast_estimate = params[var1] - params[var2]
#             se = cov_matrix.loc[var1, var1] + cov_matrix.loc[var2, var2] - 2 * cov_matrix.loc[var1, var2]

#         contrast_estimates.append(contrast_estimate)
#         contrast_se_value = np.sqrt(se)
#         contrast_se.append(contrast_se_value)

#         # Calculate z-score and p-value for the contrast
#         z_score = contrast_estimate / contrast_se_value
#         contrast_z_scores.append(z_score)
#         p_value = 2 * (1 - norm.cdf(np.abs(z_score)))
#         contrast_p_values.append(p_value)

#     # Create Series to store contrast results
#     contrast_pvalues = pd.Series(contrast_p_values, index=[f"{group1}_vs_{group2}" for group1, group2 in contrasts])
#     contrast_estimates = pd.Series(contrast_estimates, index=[f"{group1}_vs_{group2}" for group1, group2 in contrasts])

#     return contrast_pvalues, contrast_estimates

# def calculate_all_contrasts_works(result):
#     """
#     Calculate the contrasts between COPD and ACO, Asthma and ACO, and Asthma and COPD from the bootstrap results.
    
#     This function computes the contrasts as the differences in the regression coefficients of the groups. 
#     It also calculates the p-values associated with these contrasts using t-tests. The standard errors of the contrasts 
#     are computed under the assumption of independence of the coefficients.

#     Args:
#         result (DataFrame): The bootstrap results DataFrame.

#     Returns:
#         pd.Series: A series containing the p-values for each pair (COPD vs ACO, Asthma vs ACO, Asthma vs COPD, etc.).
#         dict: A dictionary containing the contrasts for each pair.
        
#     Raises:
#         Exception: If any error occurs during the computation.
#     """
#     try:
#         contrasts = {}
#         p_values_dict = {}
        
#         # Extract the coefficients for each group from the result DataFrame
#         coef_copd = result['Estimate'].get("C(Disease, Treatment(reference='Control'))[T.COPD]", 0)
#         coef_aco = result['Estimate'].get("C(Disease, Treatment(reference='Control'))[T.ACO]", 0)
#         coef_asthma = result['Estimate'].get("C(Disease, Treatment(reference='Control'))[T.Asthma]", 0)
#         coef_control = 0  # Reference group is control, so coefficient is implicitly 0
        
#         # Extract the standard errors for these coefficients from the result DataFrame
#         se_copd = result['Bootstrap SE'].get("C(Disease, Treatment(reference='Control'))[T.COPD]", 0)
#         se_aco = result['Bootstrap SE'].get("C(Disease, Treatment(reference='Control'))[T.ACO]", 0)
#         se_asthma = result['Bootstrap SE'].get("C(Disease, Treatment(reference='Control'))[T.Asthma]", 0)
#         se_control = 0  # Reference group is control, so standard error is implicitly 0
        
#         # Calculate contrasts and standard errors
#         contrasts['COPD_vs_Control'] = (coef_copd - coef_control, np.sqrt(se_copd**2 + se_control**2))
#         contrasts['ACO_vs_Control'] = (coef_aco - coef_control, np.sqrt(se_aco**2 + se_control**2))
#         contrasts['Asthma_vs_Control'] = (coef_asthma - coef_control, np.sqrt(se_asthma**2 + se_control**2))
#         contrasts['COPD_vs_ACO'] = (coef_copd - coef_aco, np.sqrt(se_copd**2 + se_aco**2))
#         contrasts['Asthma_vs_ACO'] = (coef_asthma - coef_aco, np.sqrt(se_asthma**2 + se_aco**2))
#         contrasts['Asthma_vs_COPD'] = (coef_asthma - coef_copd, np.sqrt(se_asthma**2 + se_copd**2))
        
#         # Calculate t-statistics and p-values for each contrast
#         for key, (contrast, se_contrast) in contrasts.items():
#             t_statistic = contrast / se_contrast
#             p_value = 2 * (1 - norm.cdf(np.abs(t_statistic)))
#             p_values_dict[key] = p_value
        
#         return pd.Series(p_values_dict), contrasts
#     except Exception as e:
#         logging.error(f"Error in calculating contrasts: {e}")
#         raise
# ################################################################################################
        
        
        
        
        
        
        
        
        
        
       
        
        
        
def calculate_all_contrasts_no_bootstrap(model):
    """
    Calculate the contrasts between specified pairs from the regression model.
    
    This function computes the contrasts as the differences in the regression coefficients of the groups. 
    It also calculates the p-values associated with these contrasts using t-tests. The standard errors of the contrasts 
    are computed under the assumption of independence of the coefficients.

    Args:
        model (RegressionResults): The fitted regression model.

    Returns:
        pd.Series: A series containing the p-values for each contrast pair.
        dict: A dictionary containing the contrasts for each pair.
        
    Raises:
        Exception: If any error occurs during the computation.
    """
    try:
        # Initialize dictionaries to store contrasts and p-values
        contrasts = {}
        p_values_dict = {}

        # Extract the coefficients and standard errors for each group from the model parameters
        coef_asthma = model.params.get("C(Disease, Treatment(reference='Control'))[T.Asthma]", 0)
        coef_copd = model.params.get("C(Disease, Treatment(reference='Control'))[T.COPD]", 0)
        coef_aco = model.params.get("C(Disease, Treatment(reference='Control'))[T.ACO]", 0)
        coef_control = 0  # Reference group is control, so coefficient is implicitly 0
        
        se_asthma = model.bse.get("C(Disease, Treatment(reference='Control'))[T.Asthma]", 0)
        se_copd = model.bse.get("C(Disease, Treatment(reference='Control'))[T.COPD]", 0)
        se_aco = model.bse.get("C(Disease, Treatment(reference='Control'))[T.ACO]", 0)
        se_control = 0  # Reference group is control, so standard error is implicitly 0

        # Define pairs and compute contrasts and standard errors
        pairs = {
            'Asthma_vs_Control': (coef_asthma, coef_control, se_asthma, se_control),
            'COPD_vs_Control': (coef_copd, coef_control, se_copd, se_control),
            'ACO_vs_Control': (coef_aco, coef_control, se_aco, se_control),
            'Asthma_vs_COPD': (coef_asthma, coef_copd, se_asthma, se_copd),
            'Asthma_vs_ACO': (coef_asthma, coef_aco, se_asthma, se_aco),
            'COPD_vs_ACO': (coef_copd, coef_aco, se_copd, se_aco)
        }

        # Calculate contrasts, standard errors, t-statistics, and p-values
        for key, (coef1, coef2, se1, se2) in pairs.items():
            contrast = coef1 - coef2
            se_contrast = np.sqrt(se1**2 + se2**2)
            t_statistic = contrast / se_contrast
            doff = model.df_resid
            p_value = 2 * (1 - stats.t.cdf(np.abs(t_statistic), doff))
            contrasts[key] = contrast
            p_values_dict[key] = p_value

        # Convert p-values dictionary to a pandas Series
        p_values_series = pd.Series(p_values_dict, name="p-value")
        
        return p_values_series, contrasts
    except Exception as e:
        logging.error(f"Error in calculating contrasts: {e}")
        raise
 




















# def calculate_all_contrasts(model, requested_contrasts):
#     """
#     Calculate all possible contrasts between specified groups in a regression model and return only the requested contrasts.

#     Parameters:
#     - model: A fitted regression model (e.g., statsmodels OLS result).
#     - requested_contrasts (list of str): List of requested contrast names in the format "Group1_vs_Group2".

#     Returns:
#     - contrast_pvalues (Series): Series containing p-values for each requested contrast.
#     - contrast_estimates (Series): Series containing contrast estimates.
#     """
#     # Extract coefficients and covariance matrix from the model
#     params = model.params
#     cov_matrix = model.cov_params()

#     # Define all possible groups based on the model parameters
#     groups = ['Control', 'Asthma', 'COPD', 'ACO']
    
#     # Initialize lists to store contrast results
#     all_contrast_estimates = {}
#     all_contrast_se = {}
#     all_contrast_p_values = {}

#     # Calculate all possible pairwise contrasts
#     for i, group1 in enumerate(groups):
#         for group2 in groups[i+1:]:
#             if group1 == 'Control':
#                 var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
#                 contrast_estimate = -params[var2]
#                 se = cov_matrix.loc[var2, var2]
#             elif group2 == 'Control':
#                 var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
#                 contrast_estimate = params[var1]
#                 se = cov_matrix.loc[var1, var1]
#             else:
#                 var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
#                 var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
#                 contrast_estimate = params[var1] - params[var2]
#                 se = cov_matrix.loc[var1, var1] + cov_matrix.loc[var2, var2] - 2 * cov_matrix.loc[var1, var2]

#             se_value = np.sqrt(se)
#             z_score = contrast_estimate / se_value
#             p_value = 2 * (1 - norm.cdf(np.abs(z_score)))

#             contrast_name = f"{group1}_vs_{group2}"
#             all_contrast_estimates[contrast_name] = contrast_estimate
#             all_contrast_se[contrast_name] = se_value
#             all_contrast_p_values[contrast_name] = p_value

#     # Filter to return only the requested contrasts
#     contrast_estimates = pd.Series({k: v for k, v in all_contrast_estimates.items() if k in requested_contrasts})
#     contrast_pvalues = pd.Series({k: v for k, v in all_contrast_p_values.items() if k in requested_contrasts})

#     return contrast_pvalues, contrast_estimates




# def calculate_all_contrasts(model, contrasts):
#     """
#     Calculate contrasts between specified groups in a regression model.

#     Parameters:
#     - model: A fitted regression model (e.g., statsmodels OLS result).
#     - contrasts (list of tuples): Each tuple contains the names of the groups to be contrasted.

#     Returns:
#     - contrast_pvalues (Series): Series containing p-values for each contrast.
#     - contrast_estimates (Series): Series containing contrast estimates.
#     """
#     # Extract coefficients and covariance matrix from the model
#     params = model.params
#     cov_matrix = model.cov_params()

#     # Initialize lists to store contrast results
#     contrast_estimates = []
#     contrast_se = []
#     contrast_z_scores = []
#     contrast_p_values = []
   
#     # Calculate contrasts
#     for group1, group2 in contrasts:
#         if group2 == 'Control':
#             var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
#             contrast_estimate = params[var1]
#             se = cov_matrix.loc[var1, var1]
#         elif group1 == 'Control':
#             var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
#             contrast_estimate = -params[var2]
#             se = cov_matrix.loc[var2, var2]
#         else:
#             var1 = f"C(Disease, Treatment(reference='Control'))[T.{group1}]"
#             var2 = f"C(Disease, Treatment(reference='Control'))[T.{group2}]"
#             contrast_estimate = params[var1] - params[var2]
#             se = cov_matrix.loc[var1, var1] + cov_matrix.loc[var2, var2] - 2 * cov_matrix.loc[var1, var2]

#         # Compute the contrast estimate (difference in coefficients)
#         contrast_estimate = params[var1] - params[var2]
#         contrast_estimates.append(contrast_estimate)

#         # Compute the standard error of the contrast
#         se1 = cov_matrix.loc[var1, var1]
#         se2 = cov_matrix.loc[var2, var2]
#         cov12 = cov_matrix.loc[var1, var2]
#         contrast_se_value = np.sqrt(se1 + se2 - 2 * cov12)
#         contrast_se.append(contrast_se_value)

#         # Calculate z-score and p-value for the contrast
#         z_score = contrast_estimate / contrast_se_value
#         contrast_z_scores.append(z_score)
#         p_value = 2 * (1 - norm.cdf(np.abs(z_score)))
#         contrast_p_values.append(p_value)

#     # Create Series to store contrast results
#     contrast_pvalues = pd.Series(contrast_p_values, index=[f"{group1}_vs_{group2}" for group1, group2 in contrasts])
#     contrast_estimates = pd.Series(contrast_estimates, index=[f"{group1}_vs_{group2}" for group1, group2 in contrasts])

#     return contrast_pvalues, contrast_estimates


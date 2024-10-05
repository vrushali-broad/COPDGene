import os
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.diagnostic import linear_reset, het_breuschpagan, het_white
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy.stats import shapiro
from statsmodels.stats.stattools import durbin_watson
import numpy as np
from scipy.stats import ks_2samp, norm
from scipy.stats import jarque_bera
from scipy.stats import normaltest

def check_assumptions(model, residuals, fitted_values, feature_name, output_dir, figure_title, prefix=""):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created directory: {output_dir}")

    logging.info('=====================================================')
    
    # Create a (1, 3) grid of subplots
    fig, axs = plt.subplots(1, 3, figsize=(9, 3))
    # fig, axs = plt.subplots(1, 3, figsize=(7, 2))
    plt.rcParams['axes.edgecolor'] = '#D3D3D3'
    
    # Histogram of residuals
    sns.histplot(residuals, kde=True, color='#6fa8dc', ax=axs[0])
    axs[0].set_xlabel('Residuals', fontsize=10)
    axs[0].set_ylabel('Count', fontsize=10)
    axs[0].set_title('Histogram of Residuals', fontsize=10)
    
    # Normality - Kolmogorov-Smirnov (KS) test
    ks_statistic, ks_pvalue = ks_2samp(residuals, np.random.normal(loc=np.mean(residuals), scale=np.std(residuals), size=len(residuals)))
    axs[0].text(0.05, 0.95, f'KS test p-value: {ks_pvalue:.2f}', transform=axs[0].transAxes,
                verticalalignment='top', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    # Normality - Shapiro-Wilk test
    # shapiro_test = shapiro(residuals)
    # axs[0].text(0.05, 0.95, f'Shapiro-Wilk p-value: {shapiro_test[1]:.2f}', transform=axs[0].transAxes,
    #             verticalalignment='top', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    # if shapiro_test[1] < 0.05:
    #     logging.warning(f'Shapiro-Wilk test indicates non-normality (p-value < 0.05) for {feature_name}.')
    # else:
    #     logging.info(f'Shapiro-Wilk test indicates normality (p-value >= 0.05) for {feature_name}.')

    # # Normality - D'Agostino and Pearson's test
    # dp_statistic, dp_pvalue = normaltest(residuals)
    # axs[0].text(0.05, 0.95, f'D\'Agostino and Pearson\'s test p-value: {dp_pvalue:.2f}', transform=axs[0].transAxes,
    #             verticalalignment='top', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Q-Q Plot
    qqplot = sm.qqplot(residuals, line='45', fit=True, ax=axs[1], color='#6fa8dc', marker='o', alpha=0.5)
    # Adjust marker size and color
    for line in qqplot.gca().get_lines():
        line.set_marker('o')
        line.set_markerfacecolor('#6fa8dc')
        line.set_markeredgecolor('#6fa8dc')
        line.set_markersize(5)  # Adjust marker size
    axs[1].set_title('Q-Q Plot', fontsize=10)
    
    # Linearity - Residuals vs Fitted
    sns.residplot(x=fitted_values, y=residuals, lowess=True,
                  # line_kws={'color': 'red'}, scatter_kws={'s': 5, 'color': '#6fa8dc'}, 
                  line_kws={'color': 'red', 'zorder': 3}, 
                  scatter_kws={'s': 5, 'color': '#6fa8dc', 'zorder': 2}, 
                  ax=axs[2])
    axs[2].set_xlabel('Fitted values',fontsize=10)
    axs[2].set_ylabel('Residuals', fontsize=10)
    axs[2].set_title('Residuals vs Fitted', fontsize=10)

    # Homoscedasticity - Breusch-Pagan test
    _, bp_pval, __, f_pval = het_breuschpagan(residuals, model.model.exog)
    axs[2].text(0.05, 0.95, f'Breusch-Pagan p-value: {bp_pval:.2f}', transform=axs[2].transAxes,
                verticalalignment='top', fontsize=8 , bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    if bp_pval < 0.05:
        logging.warning(f'Breusch-Pagan test indicates heteroscedasticity (p-value < 0.05) for {feature_name}.')
    else:
        logging.info(f'Breusch-Pagan test indicates no heteroscedasticity (p-value >= 0.05) for {feature_name}.')


    # Add the main title for the figure
    fig.suptitle(figure_title, fontsize=12, fontweight='bold', y=0.95)
    # fig.subplots_adjust(top=0.90)  # Adjust 'top' as needed
    
    # Save the figure with the prefix
    plot_filename = f"{prefix}{feature_name}_diagnostic_plots.pdf"
    plot_path = os.path.join(output_dir, plot_filename)
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()
    logging.info(f'Saved diagnostic plots for {feature_name} at {plot_path}')
    
    logging.info(f'Check the Residuals vs Fitted plot for {feature_name}: The residuals should be randomly scattered around zero.')
    logging.info(f'Check the Histogram of Residuals for {feature_name}: The residuals should approximately follow a normal distribution.')
    logging.info(f'Check the Q-Q plot for {feature_name}: The residuals should lie approximately along the 45-degree line.')
    
    # Independence - Durbin-Watson Test
    dw = durbin_watson(residuals)
    logging.info(f'Durbin-Watson statistic for {feature_name}: {dw}')
    if 1.5 < dw < 2.5:
        logging.info(f'Durbin-Watson statistic is within the acceptable range (1.5-2.5), indicating no autocorrelation for {feature_name}.')
    else:
        logging.warning(f'Durbin-Watson statistic is outside the acceptable range (1.5-2.5), indicating potential autocorrelation for {feature_name}.')
        
    # Multicollinearity (excluding the intercept)
    logging.info(f'Variance Inflation Factor (VIF) for {feature_name}:')
    exog = model.model.exog
    exog_names = model.model.exog_names

    for i in range(1, exog.shape[1]):  # Skip the intercept (index 0)
        vif = variance_inflation_factor(exog, i)
        logging.info(f'{exog_names[i]}: {vif:.2f}')
        if vif > 10:
            logging.warning(f'High multicollinearity detected for {exog_names[i]} (VIF > 10).')

    logging.info('=====================================================')

    
    
    
    
# import os
# import logging
# import matplotlib.pyplot as plt
# import seaborn as sns
# import statsmodels.api as sm
# from statsmodels.stats.diagnostic import linear_reset, het_breuschpagan, het_white
# from statsmodels.stats.outliers_influence import variance_inflation_factor
# from scipy.stats import shapiro
# from statsmodels.stats.stattools import durbin_watson
# import numpy as np
# from scipy.stats import ks_2samp, norm
# from scipy.stats import jarque_bera
# from scipy.stats import normaltest

# def check_assumptions(model, residuals, fitted_values, feature_name, output_dir, figure_title, prefix=""):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#         logging.info(f"Created directory: {output_dir}")

#     logging.info('=====================================================')
    
#     # Create a (1, 3) grid of subplots
#     fig, axs = plt.subplots(1, 3, figsize=(12, 4))
#     plt.rcParams['axes.edgecolor'] = '#D3D3D3'
    
#     # Histogram of residuals
#     sns.histplot(residuals, kde=True, color='#6fa8dc', ax=axs[0])
#     axs[0].set_xlabel('Residuals', fontsize=10)
#     axs[0].set_ylabel('Count', fontsize=10)
#     axs[0].set_title('Histogram of Residuals', fontsize=10)
    
#     # Normality - Kolmogorov-Smirnov (KS) test
#     ks_statistic, ks_pvalue = ks_2samp(residuals, np.random.normal(loc=np.mean(residuals), scale=np.std(residuals), size=len(residuals)))
#     axs[0].text(0.05, 0.95, f'KS test p-value: {ks_pvalue:.2f}', transform=axs[0].transAxes,
#                 verticalalignment='top', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

#     # Q-Q Plot
#     qqplot = sm.qqplot(residuals, line='45', fit=True, ax=axs[1], color='#6fa8dc', marker='o', alpha=0.5)
#     # Adjust marker size and color
#     for line in qqplot.gca().get_lines():
#         line.set_marker('o')
#         line.set_markerfacecolor('#6fa8dc')
#         line.set_markeredgecolor('#6fa8dc')
#         line.set_markersize(5)  # Adjust marker size
#     axs[1].set_title('Q-Q Plot', fontsize=10)
    
#     # Linearity - Residuals vs Fitted
#     sns.residplot(x=fitted_values, y=residuals, lowess=True, line_kws={'color': 'red'}, scatter_kws={'s': 5, 'color': '#6fa8dc'}, ax=axs[2])
#     axs[2].set_xlabel('Fitted values', fontsize=10)
#     axs[2].set_ylabel('Residuals', fontsize=10)
#     axs[2].set_title('Residuals vs Fitted', fontsize=10)

#     # Homoscedasticity - Breusch-Pagan test
#     _, bp_pval, __, f_pval = het_breuschpagan(residuals, model.model.exog)
#     axs[2].text(0.05, 0.95, f'Breusch-Pagan p-value: {bp_pval:.2f}', transform=axs[2].transAxes,
#                 verticalalignment='top', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
#     if bp_pval < 0.05:
#         logging.warning(f'Breusch-Pagan test indicates heteroscedasticity (p-value < 0.05) for {feature_name}.')
#     else:
#         logging.info(f'Breusch-Pagan test indicates no heteroscedasticity (p-value >= 0.05) for {feature_name}.')

#     # Add the main title for the figure
#     fig.suptitle(figure_title, fontsize=12, fontweight='bold', y=0.95)
    
#     # Save the figure with the prefix
#     plot_filename = f"{prefix}{feature_name}_diagnostic_plots.pdf"
#     plot_path = os.path.join(output_dir, plot_filename)
#     plt.tight_layout()
#     plt.savefig(plot_path)
#     plt.close()
#     logging.info(f'Saved diagnostic plots for {feature_name} at {plot_path}')
    
#     logging.info(f'Check the Residuals vs Fitted plot for {feature_name}: The residuals should be randomly scattered around zero.')
#     logging.info(f'Check the Histogram of Residuals for {feature_name}: The residuals should approximately follow a normal distribution.')
#     logging.info(f'Check the Q-Q plot for {feature_name}: The residuals should lie approximately along the 45-degree line.')
    
#     # Independence - Durbin-Watson Test
#     dw = durbin_watson(residuals)
#     logging.info(f'Durbin-Watson statistic for {feature_name}: {dw}')
#     if 1.5 < dw < 2.5:
#         logging.info(f'Durbin-Watson statistic is within the acceptable range (1.5-2.5), indicating no autocorrelation for {feature_name}.')
#     else:
#         logging.warning(f'Durbin-Watson statistic is outside the acceptable range (1.5-2.5), indicating potential autocorrelation for {feature_name}.')
        
#     # Multicollinearity (excluding the intercept)
#     logging.info(f'Variance Inflation Factor (VIF) for {feature_name}:')
#     exog = model.model.exog
#     exog_names = model.model.exog_names

#     for i in range(1, exog.shape[1]):  # Skip the intercept (index 0)
#         vif = variance_inflation_factor(exog, i)
#         logging.info(f'{exog_names[i]}: {vif:.2f}')
#         if vif > 10:
#             logging.warning(f'High multicollinearity detected for {exog_names[i]} (VIF > 10).')

#     logging.info('=====================================================')

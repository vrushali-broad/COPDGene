import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from visualization_helpers import plot  

from utils import perform_mwu_test
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from assumption_checks_regression import check_assumptions
from visualization_helpers import forest_plot, get_significance_colors, get_significance_stars  
from sklearn.preprocessing import QuantileTransformer, MinMaxScaler

from utils import save_kde_ecdf_plot

###############################################################
###################### CBC Boxplots ###########################
###############################################################
def analyze_and_plot_cbc(pheno, order, dir_name='Figures/CBC/', palette = ['#D3D3D3','#AFE4DE','#FAA07A','#FCFBCB']):
    
    """
    Analyze complete blood count (CBC) data from the pheno dataset and create plots.

    Parameters:
    pheno (pandas.DataFrame): DataFrame containing the pheno data.
    order (list): List defining the order of categories for plotting.
    dir_name (str): Directory name where the plots will be saved.
    palette (list): List of colors for the plots.
    """
    logging.info("#######################################################################")
    logging.info("##############  Starting analysis and plotting of CBC data...   #######")
    logging.info("#######################################################################")

    try:
        # Ensure the directory exists
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            logging.info(f"Created directory: {dir_name}")
            
        # Define CBC features and their respective names
        features = ['neutrophl_P2', 'eosinphl_P2', 'lymphcyt_P2', 'wbc_P2', 'monocyt_P2']
        names = ['Neutrophils', 'Eosinophils', 'Lymphocytes', 'WBC', 'Monocytes']

        # Initialize a dictionary to store p-values
        p_values = {feature: [] for feature in features}
        for feature in features:
            p_values[feature] = perform_mwu_test(pheno, feature, order)
            
        # Plotting the CBC data
        plot(pheno=pheno,
             features=features, 
             order=order, 
             names=names, 
             pvalues = p_values,
             plot_type='boxen', 
             dir_name=dir_name, 
             palette=palette,
            # sig = True
            )
        
        # Perform regression and store results
        results = {}
        contrast_results = {}
        
        i = 0
        
        # for feature in features: 
        for i, feature in enumerate(['eosinphl_P2']):#: features#
            print('xxxxxxx', feature)
            names = ['Eosinophils']
            df = pd.DataFrame(pheno)
            df['intercept'] = 1
            
            save_kde_ecdf_plot(df, feature, names[i], os.path.join(dir_name, f'KDE_{feature}.pdf'))
            
            qt = QuantileTransformer(output_distribution='normal', random_state=0)
            df[f'quantile_{feature}'] = qt.fit_transform(df[[feature]])
            # print('===>', feature, df[f'quantile_{feature}'])
            
#             # Initialize the transformer for Box-Cox
#             from sklearn.preprocessing import PowerTransformer
#             pt_boxcox = PowerTransformer(method='box-cox', standardize=True)
#             # Apply Box-Cox transformation
#             df[feature] = df[feature] + 1
#             df[f'boxcox_{feature}'] = pt_boxcox.fit_transform(df[[feature]])

            ############# YeoJohnson #############
            from scipy.stats import yeojohnson
            df[f'yeojohnson_{feature}'], lambda_optimal = yeojohnson(df[feature])
    
            ############# Log Transform #############
            df[f'log_{feature}'] = np.log(df[feature]+1)
            
            ############## Boxcox Transform ##########
            from scipy.stats import boxcox, norm
            df[f'boxcox_{feature}'], lambda_optimal = boxcox(df[feature] + 1) 
            
#             ############# Logit #############
#             # Ensure the data is in the range (0, 1)
#             scaler = MinMaxScaler()
#             scaled_feature = scaler.fit_transform(df[[feature]])

#             # Apply logit transformation
#             epsilon = 1e-5  # Small constant to avoid issues with 0 and 1
#             clipped_feature = np.clip(scaled_feature, epsilon, 1 - epsilon)
#             logit_transformed_column = np.log(clipped_feature / (1 - clipped_feature))

#             # Add the transformed column to the DataFrame
#             df[f'logit_{feature}'] = logit_transformed_column
#             ###################################

            df = df[df.Disease.isin(['Asthma', 'Control'])]
            # Define the formula for CBC regression model
            # it helps determine if the effect of being a current smoker on eosinophil counts is different for individuals with COPD compared to those with asthma or ACO.
            # formula_cbc = "{feature} ~ C(smoking_status_P2) * C(Disease, Treatment(reference='Control')) + Age_P2 + BMI_P2 + C(gender)"
            formula_cbc = "yeojohnson_{feature} ~ C(smoking_status_P2) * C(Disease, Treatment(reference='Control')) + Age_P2 + BMI_P2 + C(gender)"
            current_formula = formula_cbc.format(feature=feature)

            # Fit model
            model = smf.ols(formula=current_formula, data=df).fit(cov_type='HC3')
            df['corrected_' + feature] = model.resid
            df['Adjusted_' + feature] = model.fittedvalues
            results[feature] = model
            
            # Check assumptions
            output_dir = os.path.join('Figures', 'Assumption_Checks')
            residuals = model.resid
            fitted_values = model.fittedvalues
            check_assumptions(model, residuals, fitted_values, feature, output_dir, figure_title = names[i])

            ## Generate Forest plot
            params = model.params
            pvalues = model.pvalues
            conf = model.conf_int()
            conf.columns = ['Lower', 'Upper']
            conf['Estimate'] = params

            if "Intercept" in conf.index:
                conf = conf.drop("Intercept")
            if "Intercept" in pvalues.index:
                pvalues = pvalues.drop("Intercept")

            # Apply multiple testing correction (Bonferroni/Benjamini-Hochberg)
            corrected_pvalues = multipletests(pvalues, alpha=0.05, method='bonferroni')[1] #method='fdr_bh')[1]

            colors = get_significance_colors(corrected_pvalues, conf['Estimate'])
          
            order_names = ["C(smoking_status_P2)[T.2.0]",
                            "C(Disease, Treatment(reference='Control'))[T.Asthma]",
                            "C(gender)[T.2]",
                            "C(smoking_status_P2)[T.2.0]:C(Disease, Treatment(reference='Control'))[T.Asthma]",
                            "Age_P2",
                            "BMI_P2"]
            # order_names = order_names[::-1]
            rename = {"Age_P2":"Age",
                       "BMI_P2":"BMI",
                      "C(gender)[T.2]":"Gender",
                      "C(smoking_status_P2)[T.2.0]":"Smoking Status",
                      "C(Disease, Treatment(reference='Control'))[T.Asthma]":"Asthma",
                      "C(smoking_status_P2)[T.2.0]:C(Disease, Treatment(reference='Control'))[T.Asthma]":"Smoking status in Asthma"}

            fig_name = os.path.join(dir_name, 'Forest_Adj_' + feature + '.pdf')
            forest_pvalues = forest_plot(conf.index,
                                         conf['Estimate'], conf['Lower'], conf['Upper'], 
                                         colors,
                                         pvalues=corrected_pvalues, 
                                         order=order_names, 
                                         rename=rename, 
                                         title=names[i], 
                                         file_path=fig_name)
            i += 1
        
        logging.info("CBC analysis and plotting completed successfully.")
    except Exception as e:
        logging.error(f"Error in analyzing CBC data: {e}")
        raise

###############################################################
#################### Correlation Plots ########################
###############################################################
## Plot heatmap 
def create_colorbar(output_dir, cmap='PuOr'):
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    fig, ax = plt.subplots(1, figsize=(6, 2),
                           subplot_kw=dict(xticks=[], yticks=[]))
    ax.set(frame_on=False)
    ax.imshow([colors], extent=[0, 10, 2, 0.5])
    xmarks = [0, 5, 10]
    ax.set_xticks(xmarks)
    ax.set_xticklabels(['-1', '0', '1'], fontsize=30)
    output_path = os.path.join(output_dir, 'cbar_corr.pdf')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
def create_correlation_heatmap(df, output_dir, cmap='PuOr'):
    columns = ['wbc_P2', 'eosinphl_P2', 'neutrophl_P2', 'monocyt_P2', 'lymphcyt_P2']
    df = df[columns]
    df.columns = ['WBC', 'Eosinophil', 'Neutrophil', 'Monocyte', 'Lymphocyte']
    corr = df.corr()
    corr = np.round(corr, 2)
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    output_path = os.path.join(output_dir, 'corr.pdf')
    sns.heatmap(corr, cmap=cmap, vmin=-1, vmax=1,
                annot=True, annot_kws={'size': 10}, cbar=False)
    plt.xticks(rotation=90, fontsize=12, fontweight='bold')
    plt.yticks(rotation=0, fontsize=12, fontweight='bold')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

# This function will create and save a heatmap of correlations for each disease category
def create_correlation_heatmap_per_disease(df, disease_column, output_dir, cmap='PuOr'):
    # List of leukocyte types to be considered for correlation
    columns = ['wbc_P2', 'eosinphl_P2', 'neutrophl_P2', 'monocyt_P2', 'lymphcyt_P2']
    
    # Renaming columns for better readability in the heatmap
    rename_dict = {'wbc_P2': 'WBC', 'eosinphl_P2': 'Eosinophil',
                   'neutrophl_P2': 'Neutrophil', 'monocyt_P2': 'Monocyte', 'lymphcyt_P2': 'Lymphocyte'}
    
    # Check if the output directory exists; if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Iterate over each disease category to create separate heatmaps
    for disease in df[disease_column].unique():
        disease_df = df[df[disease_column] == disease][columns].rename(columns=rename_dict)
        corr = disease_df.corr()
        corr = np.round(corr, 2)
        fig = plt.figure(figsize=(4, 4))
        sns.heatmap(corr, cmap=cmap, vmin=-1, vmax=1,
                    annot=True, annot_kws={'size': 10}, cbar=False)
        plt.xticks(rotation=90, fontsize=12, fontweight='bold')
        plt.yticks(rotation=0, fontsize=12, fontweight='bold')
        plt.title(disease, fontsize=12, fontweight='bold')
        # Define the output path with the disease name
        output_path = os.path.join(output_dir, f'{disease}_corr.pdf')
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()


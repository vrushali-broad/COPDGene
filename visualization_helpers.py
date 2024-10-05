## Import libraries
import re
import logging
import scipy.stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pymannkendall as mk
from statannot import add_stat_annotation
import matplotlib.ticker as ticker
from statsmodels.stats.multitest import multipletests

import warnings
warnings.filterwarnings('ignore')


# # Create a dictionary to hold the library versions
# import statannot
# versions = {
#     "re": re.__version__,
#     "logging": logging.__version__,
#     "scipy": scipy.__version__,
#     "numpy": np.__version__,
#     "pandas": pd.__version__,
#     "seaborn": sns.__version__,
#     "pymannkendall": mk.__version__,
#     "statannot": statannot.__version__,
# }

# print(versions)

import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
plt.rcParams['font.family'] = 'DejaVu Sans'  

# Global properties for plots
PROPS = {
    'boxprops': {'edgecolor': '#929292'},
    'medianprops': {'color': '#929292'},
    'whiskerprops': {'color': '#929292'}
}
flierprops = dict(marker='o', markerfacecolor='None', markersize=1, markeredgecolor='grey')

def generate_valid_linux_filename(name):
    # Remove invalid characters and replace spaces with underscores
    # Also avoid leading or trailing dots or spaces
    # valid_name = re.sub(r'[\/:*?"<>|]', '_', name).strip('. ')
    # Replace spaces with underscores and remove invalid characters
    valid_name = re.sub(r'[\/:*?"<>|]', '_', name.replace(' ', '_')).strip('. ')

    # Check if the name is too long
    MAX_LENGTH = 255
    if len(valid_name) > MAX_LENGTH:
        valid_name = valid_name[:MAX_LENGTH]

    return valid_name

#####################################################################
################ Generates different types of plots #################
#####################################################################

def plot(pheno, features, order, names, dir_name, plot_type, figsize=(3,3), palette=['#D3D3D3','#AFE4DE','#FAA07A','#FCFBCB'], pvalues=None, sig=False, pert=False, apply_corrections = False, hue = None):
    """
    Create various types of plots for the given features.

    Args:
    pheno: DataFrame containing the data to plot.
    features: List of features (columns in pheno) to plot.
    order: Order in which to plot the categories.
    names: Names to use for the plots.
    dir_name: Directory to save the plots.
    plot_type: Type of plot to create (e.g., 'bar', 'box', etc.).
    figsize: Size of the figure.
    palette: Color palette for the plot.
    pvalues: If provided, p-values for statistical annotation.
    sig: Whether to perform significance testing.
    pert: If True, y-ticks will be multiplied by 100.
    """
    for i, feature in enumerate(features):
        plt.figure(figsize=figsize)
        plt.rcParams['axes.edgecolor'] = '#D3D3D3'
        
        ## Check if name is valid, if not generate a valid filename
        name = 'disease_' + names[i] + '.pdf'
        name = generate_valid_linux_filename(name)
        
        
        ax = None
        if plot_type == 'bar':
            ax = sns.barplot(data=pheno, 
                             x="Disease", 
                             y=feature, 
                             order=order, 
                             palette=palette, 
                             errorbar = 'se',
                             capsize=.2, 
                             errcolor=".5", 
                             errwidth=2, 
                             # errorbar = None,
                             linewidth=1, 
                             edgecolor=".5", 
                             width=0.70)

            # Calculate means and 95% confidence intervals
            grouped_data = pheno.groupby('Disease')[feature]
            means = grouped_data.mean().reindex(order)
            cis = grouped_data.apply(lambda grp: scipy.stats.t.interval(
                0.95, len(grp.dropna())-1, loc=np.mean(grp.dropna()), 
                scale=scipy.stats.sem(grp.dropna())) if len(grp.dropna()) > 1 else (np.nan, np.nan))

            # Print means and 95% confidence intervals
            # Calculate counts for each group
            counts = pheno['Disease'].value_counts().reindex(order)
            # Calculate missing data for each group
            missing_data = pheno.groupby('Disease')[feature].apply(lambda x: x.isnull().sum()).reindex(order)

            for disease in order:
                mean = means.get(disease, None)
                ci = cis.get(disease, None)
                if pert and mean is not None and ci is not None:
                    # Multiply mean and CI by 100
                    mean *= 100
                    ci = [x * 100 for x in ci]
                nonna = pheno.groupby('Disease')[feature].apply(lambda x: x.notnull().sum()).reindex(order)
                count = nonna.get(disease, None)  # Get the count for the disease
                if mean is not None and not any(np.isnan(ci)):
                    print(f"{disease} --> {feature}: Count = {count},  Mean = {mean:.2f}, 95% CI = [{ci[0]:.2f}, {ci[1]:.2f}]")
                else:
                    print(f"{disease} --> {feature}: Count = {count},  Mean = {mean:.2f}, 95% CI = Not Available")

            if pert:
                # Convert y-tick labels to percentages
                locs, labels = plt.yticks()  # Get current y-ticks
                plt.yticks(locs, [f"{loc * 100:.0f}" for loc in locs])  # Set new y-tick labels to percentages
                
            # else:
            #     # locs, labels = plt.yticks()  # Get current y-ticks
            ylim = ax.get_ylim()
            # Calculate the new y-limits
            new_ylim = (ylim[0], ylim[1] * 1.1)
            if feature == 'ShrtBrthAttk_P2':
                    new_ylim = (ylim[0], ylim[1] * 1.15)
            ax.set_ylim(new_ylim)
                
            ## Write sample size on the plot
            # Determine the top y-limit of the plot
            y_max = ax.get_ylim()[1]

            # Annotate bars with sample sizes just below the top line of the bounding box
            # counts = pheno['Disease'].value_counts().reindex(order)
            nonna_counts = pheno.groupby('Disease')[feature].apply(lambda x: x.notnull().sum()).reindex(order)
            for j, count in enumerate(nonna_counts):
                ax.annotate(f'{count}',#f'N={count}', f'N={count}', 
                            xy=(j, y_max), 
                            xycoords='data', 
                            textcoords='offset points', 
                            va='top', 
                            ha='center', 
                            fontsize=12, 
                            color= '#494949',
                            xytext=(0, -5))  # Adjust the offset as needed

        elif plot_type == 'box':
            ax = sns.boxplot(data=pheno, 
                             x="Disease",
                             y=feature, 
                             order=order,
                             palette=palette,
                             width=0.5, 
                             flierprops=flierprops,
                             showcaps=False, 
                             # showfliers = False,
                             **PROPS,
                            )
            
            locs, labels = plt.yticks()  # Get current y-ticks
            ylim = ax.get_ylim()
            # Calculate the new y-limits
            new_ylim = (ylim[0], ylim[1] * 1.02)
            # print('---->', feature)
            if feature == 'Adjusted_yeojohnson_Resting_SaO2_P2':
                new_ylim = (ylim[0], ylim[1] + 0.5)
                
            ax.set_ylim(new_ylim)
            
            ## Write sample size on the plot
            # Determine the top y-limit of the plot
            y_max = ax.get_ylim()[1]

            # Annotate bars with sample sizes just below the top line of the bounding box
            # counts = pheno['Disease'].value_counts().reindex(order)
            nonna_counts = pheno.groupby('Disease')[feature].apply(lambda x: x.notnull().sum()).reindex(order)
            for j, count in enumerate(nonna_counts):
                ax.annotate(f'{count}',#f'N={count}', 
                            xy=(j, y_max), 
                            xycoords='data', 
                            textcoords='offset points', 
                            va='top', 
                            ha='center', 
                            fontsize=12, 
                            color= '#494949',
                            xytext=(0, -5))  # Adjust the offset as needed


            # Calculate means and 95% confidence intervals
            grouped_data = pheno.groupby('Disease')[feature]
            means = grouped_data.mean().reindex(order)
            cis = grouped_data.apply(lambda grp: scipy.stats.t.interval(
                0.95, len(grp.dropna())-1, loc=np.mean(grp.dropna()), 
                scale=scipy.stats.sem(grp.dropna())) if len(grp.dropna()) > 1 else (np.nan, np.nan))

            # Print means and 95% confidence intervals
            for disease in order:
                mean = means.get(disease, None)
                ci = cis.get(disease, None)
                # Calculate the missing data for each group
                nonna = pheno.groupby('Disease')[feature].apply(lambda x: x.notnull().sum()).reindex(order)
                count = nonna.get(disease, None)  # Get the count for the disease
                # total_count = grouped_data.count().get(disease, 0)  # Count of non-missing values

                if mean is not None and not any(np.isnan(ci)):
                    if f"{ci[0]:.2f}" == f"{ci[1]:.2f}":
                        # Print the confidence interval with more precision
                        print(f"{disease} --> {feature}: Total = {count}, Mean = {mean:.2f}, 95% CI = [{ci[0]:.4f}, {ci[1]:.4f}]")
                    else:
                        # Print the confidence interval with two decimal places
                        print(f"{disease} --> {feature}: Total = {count}, Mean = {mean:.2f}, 95% CI = [{ci[0]:.2f}, {ci[1]:.2f}]")
                    # print(f"{disease} --> {feature}: Mean = {mean:.2f}, 95% CI = [{ci[0]:.2f}, {ci[1]:.2f}]")
                else:
                    print(f"{disease} --> {feature}: Total = {count}, Mean = {mean:.2f}, 95% CI = Not Available")

        elif plot_type == 'violin':
            ax = sns.violinplot(data=pheno, 
                                x="Disease", 
                                y=feature, 
                                order=order,
                                palette=palette,
                                width=0.6, **PROPS)
            
        elif plot_type == 'boxen':
            ax = sns.boxenplot(data=pheno, 
                               x="Disease", 
                               y=feature, 
                               order=order,
                               palette=palette, 
                               width=0.6,
                               flier_kws = dict(marker='o', edgecolors='grey',linewidths = 0, c = None, s = 5),
                               line_kws = {'color':'#929292', 'linewidth': 1},
                               #showfliers = False,
                              )
            ## Write sample size on the plot
            # Determine the top y-limit of the plot
            y_max = ax.get_ylim()[1]

            # Annotate bars with sample sizes just below the top line of the bounding box
            # counts = pheno['Disease'].value_counts().reindex(order)
            nonna_counts = pheno.groupby('Disease')[feature].apply(lambda x: x.notnull().sum()).reindex(order)
            for j, count in enumerate(nonna_counts):
                ax.annotate(f'{count}',#f'N={count}', f'N={count}', 
                            xy=(j, y_max), 
                            xycoords='data', 
                            textcoords='offset points', 
                            va='top', 
                            ha='center', 
                            fontsize=12, 
                            color= '#494949',
                            xytext=(0, -5))  # Adjust the offset as needed


            # Calculate means and 95% confidence intervals
            grouped_data = pheno.groupby('Disease')[feature]
            means = grouped_data.mean().reindex(order)
            cis = grouped_data.apply(lambda grp: scipy.stats.t.interval(
                0.95, len(grp.dropna())-1, loc=np.mean(grp.dropna()), 
                scale=scipy.stats.sem(grp.dropna())) if len(grp.dropna()) > 1 else (np.nan, np.nan))

            # Print means and 95% confidence intervals
            for disease in order:
                mean = means.get(disease, None)
                ci = cis.get(disease, None)
                # total_count = grouped_data.count().get(disease, 0)  # Count of non-missing values
                nonna = pheno.groupby('Disease')[feature].apply(lambda x: x.notnull().sum()).reindex(order)
                count = nonna.get(disease, None)  # Get the count for the disease
                if mean is not None and not any(np.isnan(ci)):
                    print(f"{disease} --> {feature}: Total = {count},  Mean = {mean:.2f}, 95% CI = [{ci[0]:.2f}, {ci[1]:.2f}]")
                else:
                    print(f"{disease} --> {feature}: Total = {count}, Mean = {mean:.2f}, 95% CI = Not Available")

        elif plot_type == 'swarm':
            ax = sns.swarmplot(data=pheno, 
                               x="Disease", 
                               y=feature, 
                               order=order,
                               palette=palette, 
                               size=2)
            
        elif plot_type == 'strip':
            ax = sns.stripplot(data=pheno,
                               x="Disease",
                               y=feature, 
                               order=order,
                               palette=palette,
                               size=1.5)

        # Statistical annotations
        # Define box_pairs depending on sig
        box_pairs = [("Control", "Asthma"), 
                     ("Control", "COPD"), 
                     ("Control", "ACO"), 
                     ("COPD", "ACO")] if sig else [("Control", "Asthma"), 
                                                   ("Control", "COPD"),
                                                   ("Control", "ACO")]

        # Statistical annotations
        if pvalues:
            # print('=======>', pvalues[feature])
            ## pvalues dictionary should be corrected for bonferroni
            add_stat_annotation(ax, data=pheno, 
                                    x="Disease", y=feature, 
                                    order=order, 
                                    box_pairs=box_pairs, 
                                    perform_stat_test=False, 
                                    pvalues = pvalues[feature], 
                                    test=None, 
                                    text_format='star',
                                    loc='outside',
                                    verbose=0, 
                                    line_height=0, 
                                    text_offset=2,
                                    color='grey',
                                    fontsize=14)
        else:
            ## Defaults to Bonferroni correction
            add_stat_annotation(ax, data=pheno, 
                                x="Disease", y=feature, 
                                order=order,
                                box_pairs=box_pairs, 
                                test='Mann-Whitney', 
                                text_format='star', 
                                loc='outside',
                                verbose=0,
                                # verbose=1, 
                                line_height=0, 
                                text_offset=2, 
                                color='grey', 
                                fontsize=14)


        # Calculate and print trend information
        data = [grp[feature].sum() / len(grp) for dis, grp in pheno.groupby('Disease') if dis in ['Control', 'Asthma', 'COPD', 'ACO']]

        ylabel_fontsize = 14 if pert else 16
        
        # Enable LaTeX in Matplotlib
        # plt.rcParams['text.usetex'] = True
        plt.rcParams['axes.labelweight'] = 'bold'
        # Define font properties for bold text
        from matplotlib.font_manager import FontProperties
        bold_font = FontProperties(weight='bold')


        ax.set_ylabel(names[i], fontsize=ylabel_fontsize, fontweight='bold', fontproperties=bold_font) 
        # plt.ylabel(names[i], fontsize=ylabel_fontsize, fontweight='bold')#, color='blue', backgroundcolor='lightgrey')
        plt.xlabel('Disease', fontsize=16, fontweight='bold', fontproperties=bold_font) 
        plt.xticks(fontsize=14, rotation=30)
        plt.yticks(fontsize=14)
        
        plt.savefig(f'{dir_name}/{name}', dpi=300, bbox_inches='tight')
        # i += 1
        # plt.show()
        plt.close()
        
#####################################################################
##################### Generates forest plots ########################
#####################################################################       
def get_significance_colors(pvalues, estimates, alpha=0.05):
    """
    Generate a list of colors representing the significance of estimates based on p-values.

    Args:
    pvalues (list of float): List of p-values for each estimate.
    estimates (list of float): List of estimate values corresponding to the p-values.
    alpha (float, optional): Significance level. Defaults to 0.05.

    Returns:
    list of str: List of color codes representing the significance of each estimate.

    Raises:
    ValueError: If the length of pvalues and estimates do not match.

    Example:
    >>> get_significance_colors([0.01, 0.04, 0.2], [1.5, -2.0, 0.3])
    ['#a5d732', '#ffb353', '#929292']

    The function assigns green (#a5d732) for significant positive estimates, 
    red (#ffb353) for significant negative estimates, and grey (#929292) for non-significant estimates.
    """

    # Check if the length of pvalues and estimates match
    if len(pvalues) != len(estimates):
        logging.error("Length of pvalues and estimates do not match.")
        raise ValueError("pvalues and estimates lists must have the same length.")

    colors = []
    for p, e in zip(pvalues, estimates):
        if p < alpha:
            if e > 0:
                colors.append('#a5d732')  # Green for significant positive estimate
            else:
                colors.append('#ffb353')  # Red for significant negative estimate
        else:
            colors.append('#929292')  # Grey for non-significant
    return colors

def get_significance_stars(pvalues):
    """
    Convert a list of p-values into significance stars.

    Each p-value is converted into a string of stars ('*', '**', '***', '****') 
    based on its significance level. The function handles a list of p-values 
    and returns a corresponding list of strings representing the significance.

    Args:
    pvalues (list of float): A list of p-values.

    Returns:
    list of str: A list of strings, each representing the significance level of the corresponding p-value.

    Example:
    >>> get_significance_stars([0.00005, 0.0005, 0.05, 0.1])
    ['****', '***', '*', '']

    Raises:
    ValueError: If any element in pvalues is not a float.
    """

    # Initialize an empty list to store significance stars
    stars = []

    # Iterate over each p-value in the input list
    for p in pvalues:
        # Check if the input is a valid float
        if not isinstance(p, float):
            logging.error("Non-float value found in p-values: %s", p)
            raise ValueError(f"Invalid p-value: {p}. P-values must be floats.")

        # Append the appropriate number of stars based on the p-value
        if p < 0.0001:
            stars.append('****')
        elif p < 0.001:
            stars.append('***')
        elif p < 0.01:
            stars.append('**')
        elif p < 0.05:
            stars.append('*')
        else:
            stars.append('')

    return stars

def forest_plot(labels, estimates, lower, upper, colors, pvalues, order=None, rename=None, y = None, title="Forest Plot", xlabel="Estimate", file_path = 'Forest_plot.pdf'):

    """
    Create a forest plot and save it to a file.

    Args:
        labels (list): Labels for each row in the plot.
        estimates (list): Point estimates for each row.
        lower (list): Lower bounds of the confidence intervals.
        upper (list): Upper bounds of the confidence intervals.
        colors (list): Colors for each point estimate.
        pvalues (list): P-values for each comparison.
        file_path (str): File path to save the plot.
        order (list, optional): Order of the rows in the plot.
        rename (dict, optional): Mapping of original labels to new labels.
        y (float, optional): Y-position for the plot title.
        title (str): Title of the plot.
        xlabel (str): Label for the x-axis.

    Returns:
        list: List of p-values.
    """
    # Convert data to a DataFrame for easy manipulation
    data = pd.DataFrame({
        'labels': labels,
        'estimates': estimates,
        'lower': lower,
        'upper': upper,
        'colors': colors,
        'pvalues': pvalues
    })

    # Reorder rows if 'order' is provided
    if order:
        data['labels'] = pd.Categorical(data['labels'], categories=order, ordered=True)
        data = data.sort_values('labels')

    # Rename labels if 'rename' is provided
    if rename:
        data['labels'] = data['labels'].replace(rename)
           
    # Create plot
    fig, ax = plt.subplots(figsize=(2, len(data) * 0.3))

    stars = get_significance_stars(data['pvalues'])
    for i, (index, row) in enumerate(data.iterrows()):
        ax.errorbar(row['estimates'], i, xerr=[[row['estimates'] - row['lower']],
                                       [row['upper'] - row['estimates']]], 
                    # fmt = 'o', 
                    color = row['colors'],
                    ecolor = '#929292',
                    elinewidth = 1.5,
                    marker = "o",# "s", 
                    capsize = 0, 
                    ms = 4)
        
        ax.text(row['estimates'], i + 0.1, stars[i], color='k', ha='center', fontsize = 9)  # Use `i` here
        
        max_value = max([max(map(abs, lower)), max(map(abs, upper))])
            
        # Calculate the total range of the x-axis
        x_range = 2*max_value + 2  # This assumes ax.set_xlim([-max_value-1, max_value+1])

        # Define fixed fractions of the x_range for column spacing
        spacing_pvalue = 0.15 * x_range
        spacing_estimate = 0.50 * x_range
        spacing_ci = 0.90 * x_range

        # Calculate offsets
        offset_pvalue = max_value + 1 + spacing_pvalue
        offset_estimate = max_value + 1 + spacing_estimate
        offset_ci = max_value + 1 + spacing_ci
        
        ## P value text
        if row['pvalues'] < 0.0001:
            p_value_text = "$<10^{-4}$"
        elif row['pvalues'] < 0.001:
            p_value_text = "$<10^{-3}$"
        else:
            p_value_text = f"{row['pvalues']:.3f}"
        
        # Apply offsets in the text positioning
        ax.text(offset_pvalue, i, p_value_text, ha='left', va='center', fontsize=9)
        ax.text(offset_estimate, i, f"{row['estimates']:.2f}", ha='left', va='center', fontsize=9)
        ax.text(offset_ci, i, f"({row['lower']:.2f}, {row['upper']:.2f})", ha='left', va='center', fontsize=9)
        
    # Adding column headers
    ax.text(offset_pvalue, len(labels), 'P-value', ha='left', va='center', weight='bold', fontsize = 8)
    ax.text(offset_estimate, len(labels), 'Estimate', ha='left', va='center', weight='bold', fontsize = 8)
    ax.text(offset_ci, len(labels), '95% CI', ha='left', va='center', weight='bold', fontsize = 8)

    ax.axvline(x=0, color="#929292", linestyle='--', linewidth = 1) #6fa8dc
    ax.set_xlabel(xlabel, fontsize = 10)#, fontweight="bold")
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data['labels'], fontsize = 10)#, fontweight = 'bold')
    # Change the fontsize of x-axis tick labels
    ax.tick_params(axis='x', labelsize=10)

    
    ax.set_xlim([-max_value-1, max_value+1])
    
    # print('#################################', title)
    # if title == 'Eosinophils':
    #     ax.set_xlim([-max_value, max_value])
    
    ## Title
    title_x_position = (max_value) / 2  #(2*max_value + 2) / 2  # This centers the title by placing it in the middle of the x-axis
    title_y_position = len(labels) + 0.5  # This aligns the title vertically with the column headers
    if y:
        title_y_position = y
    ax.text(0, title_y_position, title, fontsize=10, fontweight="bold", ha='center')
    
    # # Add grid lines without losing ticks
    plt.grid(True, linestyle="--", alpha=0.7)   # Customize grid appearance as needed
    sns.despine(left=True, right=True, top=True)
    
    plt.savefig(file_path, bbox_inches='tight')
    plt.close(fig)
    
    ordered_data = data[data.labels.isin(['Asthma', 'COPD', 'ACO'])][['labels', 'pvalues']].set_index('labels').reindex(['Asthma', 'COPD', 'ACO']).reset_index()
    return list(ordered_data.pvalues.values)


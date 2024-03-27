import os
import logging
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from visualization_helpers import plot  # Ensure this is correctly imported

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

        # Plotting the CBC data
        plot(pheno=pheno,
             features=features, 
             order=order, 
             names=names, 
             plot_type='boxen', 
             dir_name=dir_name, 
             palette=palette,
            # sig = True
            )
        
        logging.info("CBC analysis and plotting completed successfully.")
    except Exception as e:
        logging.error(f"Error in analyzing CBC data: {e}")
        raise

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


import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore
from statannot import add_stat_annotation
sns.set_style('ticks')

# Setting up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_and_visualize_gene_data(file_paths, gene_lists=None, disease_val='Disease'):
    """
    Processes and visualizes genetic data for specific gene lists.

    Parameters:
    - gene_lists (dict): A dictionary of gene lists, with each key being a list name and its value being a list of genes.
    - file_paths (dict): A dictionary containing file paths for 'DE_ACO', 'DE_COPD', 'Annotations', 'Gene Counts', and 'Pheno'.
    - disease_val (str, optional): The column name in 'Pheno' file that indicates the disease. Defaults to 'Disease'.

    Returns:
    None - Generates and saves boxplot visualizations for each gene in the provided lists.
    """

    try:
        # Load Data
        de_aco = pd.read_csv(file_paths['DE_ACO']).set_index(['Unnamed: 0'])
        de_copd = pd.read_csv(file_paths['DE_COPD']).set_index(['Unnamed: 0'])
        annot = pd.read_csv(file_paths['Annotations'], sep='\t', low_memory=False).set_index('ensembl_gene_id')
        df = pd.read_csv(file_paths['Gene Counts'], low_memory=False)
        meta = pd.read_csv(file_paths['Pheno'], sep='\t', low_memory=False).sort_values(disease_val).set_index(['sid'])

        # Data Processing
        de_aco = de_aco[de_aco.logFC > 0]
        de_copd = de_copd[de_copd.logFC > 0]
        de_aco.columns = [col + '_aco' for col in de_aco.columns]
        de_copd.columns = [col + '_copd' for col in de_copd.columns]

        merged_df = pd.merge(de_aco, de_copd, left_index=True, right_index=True, how='outer')
        merged_df = merged_df[['adj.P.Val_aco', 'adj.P.Val_copd']]
        merged_df['adj.P.Val_asthma'] = 1
        merged_df = merged_df[['adj.P.Val_asthma', 'adj.P.Val_copd', 'adj.P.Val_aco']]
        merged_df = merged_df.fillna(1)

        condition = (merged_df['adj.P.Val_copd'] <= 0.05) & (merged_df['adj.P.Val_aco'] <= 0.05)
        filtered_rows = merged_df[condition]
        de_all = list(filtered_rows.index.values)
        
#         condition = (merged_df['adj.P.Val_copd'] <= 0.05) & (merged_df['adj.P.Val_aco'] > 0.05)
#         filtered_rows_copd = merged_df[condition]
#         de_copd = list(filtered_rows_copd.index.values)
#         condition = (merged_df['adj.P.Val_aco'] <= 0.05) & (merged_df['adj.P.Val_copd'] > 0.05)
#         filtered_rows_aco = merged_df[condition]
#         de_aco = list(filtered_rows_aco.index.values)
        
        # Convert DataFrame to dictionary for p-values
        pvalues = {index: row.tolist() for index, row in merged_df.iterrows()}

        # Further Data Preparation
        df.columns = ['Genes'] + list(df.columns[1:].values)
        df_all = df.set_index('Genes').T
        df_all = df_all.apply(zscore, axis=0)
            
        # Conditional execution based on whether gene_lists is provided
        if gene_lists:
            ## Plot genes provided through gene_lists
            # for gene_list in gene_lists.values():
            for gene_list_name, gene_list in gene_lists.items():
                logging.info(f"Processing gene list: {gene_list_name}")
                for gene_name in gene_list:
                    if not annot[annot.hgnc_symbol == gene_name].empty:
                        gene_ids = annot[annot.hgnc_symbol == gene_name].index.values
                        for gene_id in gene_ids:
                            if gene_id in df_all.columns:
                                logging.info(f"Plotting data for gene: {gene_name} (ID: {gene_id})")
                                tmp = pd.DataFrame(df_all[gene_id])
                                tmp[disease_val] = meta.loc[tmp.index, disease_val]

                                visualize_gene_expression(tmp, gene_id, annot, pvalues, disease_val)
                    else:
                        logging.warning(f"Gene {gene_name} not found in annotations.")
        else:
            for gene_id in de_all:
                tmp = pd.DataFrame(df_all[gene_id])
                tmp[disease_val] = meta.loc[tmp.index, disease_val]
                visualize_gene_expression(tmp, gene_id, annot, pvalues, disease_val)
            
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
        
def visualize_gene_expression(data, gene_id, annotations, pvalues, disease_val):
    """
    Plots gene expression data as a boxplot for different disease conditions.

    Parameters:
    - data (DataFrame): The gene data to be visualized.
    - gene_id (str): The ID of the gene.
    - annotations (DataFrame): The annotations DataFrame.
    - pvalues (dict): Dictionary of p-values for statistical annotation.
    - disease_val (str): The disease column name.

    Returns:
    None - Generates and saves a boxplot visualization for the gene.
    """
    # Visualization
    sns.set_style('ticks')
    plt.figure(figsize=(2.5, 2.5))
    plt.rcParams['axes.edgecolor'] = '#D3D3D3'
    
    # Flier properties for the boxplot.
    flierprops = dict(marker='o', markerfacecolor='None', markersize=0.5, markeredgecolor='grey')

    # Boxplot properties.
    props_genes = {'boxprops': {'edgecolor': '#929292', 'linewidth': 0.1},
                    'medianprops': {'color': '#929292'},
                    'whiskerprops': {'color': '#929292'}
                    }
    
    order = ['Control', 'Asthma', 'COPD', 'ACO']
#     ax = sns.boxplot(data=data, x=disease_val, y=gene_id, order=order, palette=['#d3d3d3', '#AFE4DE', '#ffa07a', '#fff157'],
#                      notch=True, width=0.5, flierprops=flierprops, showcaps=False, **props_genes)
    
    ax = sns.boxplot(data=data, x=disease_val, y=gene_id, order=order,
                     palette=['#d3d3d3', '#AFE4DE', '#ffa07a', '#fff157'],
                     notch=True, width=0.5,
                     showcaps=False,
                     flierprops=flierprops, **props_genes)
                     # flierprops=dict(marker='o', 
                     #                 markerfacecolor='None',
                     #                 markersize=1,
                     #                 markeredgecolor='grey'),
                     # boxprops={'edgecolor': '#929292', 'linewidth': 0.1},
                     # medianprops={'color': '#929292'},
                     # whiskerprops={'color': '#929292'})

    
    # Check if gene_id is in pvalues
    if gene_id in pvalues:
        # Add statistical annotation
        add_stat_annotation(ax, data=data, x=disease_val, y=gene_id,
                            order=order,
                            box_pairs=[("Control", "Asthma"), 
                                       ("Control", "COPD"), 
                                       ("Control", "ACO")],
                            perform_stat_test=False, 
                            pvalues=pvalues[gene_id],
                            test=None,
                            text_format='star', 
                            loc='outside',
                            verbose=0,
                            line_height=0,
                            text_offset=2,
                            color='grey',
                            fontsize=14)
    else:
        logging.warning(f"No p-values found for gene ID: {gene_id}")
        
    # Set border width
    for spine in ax.spines.values():
        spine.set_linewidth(1)

    title = annotations.loc[gene_id, 'hgnc_symbol']
    # title = title if isinstance(title, str) else title.values[0]
    # logging.info(f'Gene name: {title}')
    
    # Check if title is a Series, NaN, or a single value
    if isinstance(title, pd.Series):
        if title.empty:
            logging.warning(f"No gene symbol found for gene ID: {gene_id}")
            title = gene_id
            # return  # Skip this gene_id
        title = title.iloc[0]  # Take the first value if there are multiple
    elif pd.isna(title):
        logging.warning(f"No gene symbol found for gene ID: {gene_id}")
        return  # Skip this gene_id

    logging.info(f'Gene name: {title}')
    
    plt.ylabel(title, fontsize=16, fontweight='bold')
    plt.xlabel('')
    plt.xticks(np.arange(len(order)), order, fontsize=14, rotation=30)
    plt.yticks(fontsize=14)

    # Saving the figure
    fig_name = f'Figures/Genes/{title}.pdf'
    plt.savefig(fig_name, bbox_inches='tight')
    plt.close()


# Example usage
gene_sets = {
    # 'th2':['IL13RA1', 'CSF1', 'IL10RB'],
        'eos': ['PTGS2','P2RY14','HES1','DACH1', 'CAMK1', 'OLIG1', 'OLIG2', 'GFOD1', 
                'IDO1', 'SORD', 'SORD2P', 'CEBPE','CDK15',
                'PLAAT5', 'ACSM3', 'ADGRE1', 'ACACB', 'PIK3R6','CACNG6', 'SMPD3', 'AJAP1', 'ADGRE4P',
                'ADORA3','RNASE3', 'RNASE2'],
        'th1': ['IL1R2', 'IL18', 'IFNGR1', 'IFNGR2', 'IFNAR1', 'IFNAR2','HIF1A', 'CXCR2',
                'CXCR1', 'IL1RAP', 'IL1B', 'VAV1', 'CXCL8', 'TNFAIP2',
                'JAK1', 'TYK2', 'IL1R2', 'IL1RN', 'TNF', 'CCR2',
                'IL18RAP'] + ['CD8B', 'CD8A', 'CD84', 'CD83', 'CXCR3', 'STAT4'],
        'neut': ['CSF2RA', 'CSF2RB', 'CSF3R', 'SERPINA1', 'DUSP1', 'FUT7', 'VMP1','EPHB1'],
        'inflammation': ['AZU1', 'NFKBIA', 'DUSP6', 'EGLN1'],
        'defense': ['KLF4','BPI','DEFA4','PLA2G4A'],
    }

# gene_sets = { 'th1': ['IL1R2', 'IL18', 'IFNGR1']}

# gene_lists = {'eos': eos, 'th1': th1, 'neut': neut}
file_paths = {
    'DE_ACO': 'Data/Processed_Data/DE_ACO.csv',
    'DE_COPD': 'Data/Processed_Data/DE_COPD.csv',
    'Annotations': 'Data/Processed_Data/Annotations.tsv',
    'Gene Counts': 'Data/Processed_Data/normalized_batch_corrected_gene_counts.csv',
    'Pheno': 'Processed_Data/pheno.tsv'
}


process_and_visualize_gene_data(file_paths = file_paths)#, gene_lists = gene_sets)


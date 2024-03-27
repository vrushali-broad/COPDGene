import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_and_preprocess_meta(file_path):
    """
    Load and preprocess meta data from a specified file.

    This function reads a CSV file containing meta data, sets 'actual_id' as the index,
    and filters out rows where 'final.analyzed' is not 1.

    Parameters:
    file_path (str): The path to the meta data file.

    Returns:
    pandas.DataFrame: Preprocessed meta data.

    Example:
    >>> meta_data = load_and_preprocess_meta('Data/master.file.freeze4.txt')
    """

    try:
        meta = pd.read_csv(file_path, sep='\t').set_index('actual_id')
        # Select only those patients that have final.analyzed = 1
        meta = meta[meta['final.analyzed'] == 1]
        logging.info(f"Loaded and preprocessed meta data: {meta.shape}")
        return meta
    except Exception as e:
        logging.error(f"Error in loading and preprocessing meta data: {e}")
        raise

        
def load_and_preprocess_counts(file_path):
    """
    Load and preprocess counts data from a specified file.

    This function reads a TSV file containing counts data. It sets 'GENEID' as the index
    and preprocesses the index to remove any additional characters after a period in gene IDs.

    Parameters:
    file_path (str): The path to the counts data file.

    Returns:
    pandas.DataFrame: Preprocessed counts data.

    Example:
    >>> counts_data = load_and_preprocess_counts('Data/counts_raw_from_lengthScaledTPM.tsv')
    """

    try:
        counts = pd.read_csv(file_path, sep='\t', low_memory=False).set_index('GENEID')
        counts.index = [gene.split('.')[0] for gene in counts.index.values]
        logging.info(f"Loaded and preprocessed counts data: {counts.shape}")
        return counts
    except Exception as e:
        logging.error(f"Error in loading and preprocessing counts data: {e}")
        raise

def load_and_preprocess_pheno(file_path, counts_columns):
    """
    Load and preprocess pheno data from a specified file.

    This function reads a TSV file containing pheno data, sets 'sid' as the index,
    and filters the data to include only those records present in counts data columns.
    It also filters out records missing certain key measurements.

    Parameters:
    file_path (str): The path to the pheno data file.
    counts_columns (list): List of columns to filter the pheno data.

    Returns:
    pandas.DataFrame: Preprocessed pheno data.

    Example:
    >>> pheno_data = load_and_preprocess_pheno('Data/COPDGene_P1P2P3_Flat_SM_NS_Nov21_PhenoData.txt', counts_data.columns)
    """

    try:
        pheno = pd.read_csv(file_path, sep='\t', low_memory=False).set_index('sid')
        pheno = pheno[pheno.index.isin(counts_columns)]
        # print(pheno.shape)
        ## Select patients that have CBC data
        pheno = pheno[(pheno['neutrophl_P2'].notna()) & (pheno['neutrophl_pct_P2'].notna()) &
                      (pheno['eosinphl_P2'].notna()) & (pheno['eosinphl_pct_P2'].notna())]
        logging.info(f"Loaded and preprocessed pheno data: {pheno.shape}")
        return pheno
    except Exception as e:
        logging.error(f"Error in loading and preprocessing pheno data: {e}")
        raise

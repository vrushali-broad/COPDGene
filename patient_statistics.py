# File: patient_statistics.py
import pandas as pd
import logging


def calculate_mean(data, variable, group):
    """
    Calculate mean of a variable for each group.

    Parameters:
    data (pandas.DataFrame): DataFrame containing the data.
    variable (str): The name of the column for which the mean is to be calculated.
    group (str): The name of the column by which to group the data.

    Returns:
    pandas.Series: A series with the mean of the variable for each group.
    """
    try:
        mean_values = data.groupby(group)[variable].mean().round(2)
        logging.info(f"Calculated mean for '{variable}' grouped by '{group}'.")
        return mean_values
    except Exception as e:
        logging.error(f"Error in calculating mean for '{variable}' grouped by '{group}': {e}")
        raise

def calculate_percentage(data, variable, group):
    """
    Calculate percentage of a variable for each group.

    Parameters:
    data (pandas.DataFrame): DataFrame containing the data.
    variable (str): The name of the column for which the percentage is to be calculated.
    group (str): The name of the column by which to group the data.

    Returns:
    pandas.DataFrame: A DataFrame with percentages of the variable for each group, 
                      rounded to two decimal places.
    """
    try:
        percentage_values = data.groupby(group)[variable].value_counts(normalize=True).mul(100)
        # Round the percentages to two decimal places
        percentage_values_rounded = percentage_values.round(2)
        logging.info(f"Calculated and rounded percentage for '{variable}' grouped by '{group}'.")
        return percentage_values_rounded.unstack()
    except Exception as e:
        logging.error(f"Error in calculating percentage for '{variable}' grouped by '{group}': {e}")
        raise

def generate_statistics_table(pheno, output_path):
    """
    Generate a summary statistics table for patient data and save it to a file.

    This function calculates the mean age and BMI, as well as the percentage distribution 
    of gender, race, smoking status, and GOLD stage for each disease category in the pheno data. 
    The results are saved to a specified output path as a tab-separated CSV file.

    Parameters:
    pheno (pandas.DataFrame): DataFrame containing patient data. Expected columns:
                              - Age_P2: Age of patients.
                              - BMI_P2: Body Mass Index of patients.
                              - gender: Gender of patients (1=Male, 2=Female).
                              - race: Race of patients (1=White, 2=Black or African American).
                              - smoking_status_P2: Smoking status (0=Never-smoked, 1=Former smoker, 2=Current smoker).
                              - finalGold_P2: GOLD stage of patients (ranging from 0 to 4).
    output_path (str): File path where the statistics table will be saved as a CSV file.

    The function prints the generated table to the console and saves it to the given output path.

    Returns:
    pandas.DataFrame: The generated statistics table.
    
    Raises:
    Exception: If any error occurs during the generation of the table.

    The 'Disease' column in the pheno DataFrame is used to group the data for analysis.
    """
    try:
        table = pd.DataFrame()

        # Age
        table['Age (years)'] = calculate_mean(pheno, 'Age_P2', 'Disease')

        # BMI
        table['BMI'] = calculate_mean(pheno, 'BMI_P2', 'Disease')

        # Gender percentage
        gender = calculate_percentage(pheno, 'gender', 'Disease')
        for gender_type in gender.columns:
            table[f'Gender {gender_type} (%)'] = gender[gender_type]

        # Race percentage
        race = calculate_percentage(pheno, 'race', 'Disease')
        for race_type in race.columns:
            table[f'Race {race_type} (%)'] = race[race_type]

        # Smoking status percentage
        smoking = calculate_percentage(pheno, 'smoking_status_P2', 'Disease')
        for smoking_type in smoking.columns:
            table[f'Smoking {smoking_type} (%)'] = smoking[smoking_type]

        # Gold stage percentage
        gold_stage = calculate_percentage(pheno, 'finalGold_P2', 'Disease')
        for stage in gold_stage.columns:
            table[f'Gold Stage {stage} (%)'] = gold_stage[stage]
            
        # Ensure that 'Disease' is set as the index (row labels) if not already
        # if 'Disease' not in table.index.names:
        #     table = table.set_index('Disease')

        # Check if the DataFrame has a meaningful index; if not, set one
        if 'Disease' in table.columns:
            table.set_index('Disease', inplace=True)
            
        # Convert the DataFrame to a string and print it
        table_string = table.to_string()
        print(table_string)

        # Save the table as a CSV file with tab separator
        table.to_csv(output_path, sep='\t')#, index=False)
        logging.info(f"Statistics table saved to {output_path}")

        return table
    except Exception as e:
        logging.error(f"Error in generating statistics table: {e}")
        raise

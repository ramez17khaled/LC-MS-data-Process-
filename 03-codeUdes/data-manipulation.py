inputPOSFile = "D:/data/MSDial/01.1-TermoData/1.2-NICOresultsPOS-thermo/AlignPOSData.csv"
inputNEGFile = "D:/data/MSDial/01.1-TermoData/1.2-NICOresultsNEG-thermo/AlignNEGData.csv"
metadataPOS = "D:/data/MSDial/02-Nico_metadata/TermoMetaData/thermoPOS.csv"
metadataNEG = "D:/data/MSDial/02-Nico_metadata/TermoMetaData/thermoNEG.csv"
outputPath = "D:/data/MSDial/05-codeOutput/Thermo_results/"
outputPOSFile = 'POS-manipulated'
outputNEGFile = 'NEG-manipulated'

import pandas as pd
import numpy as np
import os

def read_file(file_path, sep=';', encoding='utf-8'):
    """
    Read a file from the specified path.

    Parameters:
        file_path (str): Path to the file.
        sep (str, optional): Separator to use when reading CSV files. Defaults to ','.
        encoding (str, optional): Encoding to use when reading the file. Defaults to 'utf-8'.

    Returns:
        DataFrame: DataFrame containing the data from the file, or None if the file doesn't exist.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found at the specified path: {file_path}")

    # Determine file type based on extension
    file_ext = os.path.splitext(file_path)[1].lower()

    if file_ext == ".csv":
        # Read CSV file
        data = pd.read_csv(file_path, sep=sep, encoding=encoding)
    elif file_ext in [".xls", ".xlsx"]:
        # Read Excel file
        data = pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file format. Only CSV and Excel files are supported.")
    
    return data

def add_metabolite_column(data):
    """
    Add a column called 'metabolite' to the DataFrame based on the values in the 'RT(min)' and 'MZ' columns.

    Parameters:
        data (DataFrame): Input DataFrame containing 'RT(min)' and 'MZ' columns.

    Returns:
        DataFrame: DataFrame with the 'metabolite' column added.
    """
    data['metabolite'] = 'M' + data['Mz'].astype(str).str.split('.').str[0] + 'T' + data['Rt(min)'].astype(str)
    return data


def filter_column(data, include=None, exclude=None):
    """
    Filter columns from a DataFrame based on included and excluded columns.

    Parameters:
        data (DataFrame): Input DataFrame.
        include(list, optional): List of column names to include. Default is None.
        exclude (list, optional): List of column names to exclude. Default is None.

    Returns:
        DataFrame: DataFrame containing filtered columns.
    """
    if include is None and exclude is None:
        raise ValueError("Either 'include' or 'exclude' must be provided.")

    if include is not None and exclude is not None:
        raise ValueError("Only one of 'include' or 'exclude' should be provided.")

    if include is not None:
        selected_columns = include
    else:
        selected_columns = [col for col in data.columns if col not in exclude]

    selected_data = data[selected_columns]
    return selected_data

def select_columns_with_metabolites_columns(data, list_columns=None):
    """
    Selects columns from a DataFrame that start with "M" and includes additional columns specified in the list.

    Parameters:
    - data: DataFrame
        The input DataFrame from which columns will be selected.
    - list_columns: list or None, optional
        A list of additional columns to include along with the columns starting with "M". If None, only columns starting with "M" will be selected. Default is None.

    Returns:
    - DataFrame
        A subset of the input DataFrame containing selected columns.
    """
    # Select columns that start with "M"
    m_columns = [col for col in data.columns if col.startswith("M")]
    
    # If list_columns is not provided, return only columns starting with "M"
    if list_columns is None:
        return data[m_columns]
    
    # Include additional columns specified in the list_columns list
    selected_columns = m_columns + list_columns
    
    # Return the subset of the DataFrame with selected columns
    return data[selected_columns]

def filter_rows(data, include=None, exclude=None):
    """
    Select and/or exclude rows based on specific conditions.

    Parameters:
        data (DataFrame): Input DataFrame.
        include (str or list): Condition(s) to include rows (optional) ex: list =["column1 == 'value'", "column2 <10"].
        exclude (str or list): Condition(s) to exclude rows (optional) ex: list =["column1 == 'value'", "column2 <10"].

    Returns:
        DataFrame: Filtered DataFrame.
    """
    if include is not None:
        if isinstance(include, str):
            include_condition = include
        elif isinstance(include, list):
            include_condition = ' or '.join(include)
        else:
            raise ValueError("Include parameter must be a string or a list of strings.")
    else:
        include_condition = None

    if exclude is not None:
        if isinstance(exclude, str):
            exclude_condition = exclude
        elif isinstance(exclude, list):
            exclude_condition = ' and '.join(['not (' + cond + ')' for cond in exclude])
        else:
            raise ValueError("Exclude parameter must be a string or a list of strings.")
    else:
        exclude_condition = None

    if include_condition is not None:
        selected_data = data.query(include_condition)
    else:
        selected_data = data

    if exclude_condition is not None:
        excluded_data = selected_data.query(exclude_condition)
    else:
        excluded_data = selected_data

    return excluded_data


def transpose_data(data, metabolite_column):
    """
    Transpose the data with specified column values as column names.

    Parameters:
        data (DataFrame): Input DataFrame.
        metabolite_column (str): Name of the column to use as column names.

    Returns:
        DataFrame: Transposed DataFrame with specified column values as column names.
    """
    transposed_data = data.set_index(metabolite_column).T
    transposed_data.index.name ='sample_name'
    return transposed_data

def merge_data(df1, df2):
    """
    Merge two datasets based on their indices.
    Make sue that the indices are the same in both datasets

    Parameters:
        df1 (DataFrame): First DataFrame.
        df2 (DataFrame): Second DataFrame.
    
    Returns:
        DataFrame: Merged DataFrame.
    """
    merged_data = pd.merge(df1, df2, left_index=True, right_index=True, how='outer')
    return merged_data


def blank_filter(data, condition_column='SampleType', condition_values=['blank'], operation='!=', threshold=0):
    """
    Filter columns from a DataFrame based on a condition applied to specific rows.

    Parameters:
    - data: DataFrame
        The input DataFrame from which columns will be removed.
    - condition_column: str
        The name of the column to use as the condition (default parameter = 'SampleType').
    - condition_values: list
        A list of values in the condition_column to consider for the condition (default parameter = ['blank']).
    - operation: str
        The comparison operation to apply. Supported operations are '!=', '<=', '>=', and '=' (default parameter = '!=').
    - threshold: int or float
        The threshold value for the condition (default parameter = 0).

    Returns:
    - DataFrame
        A modified DataFrame with columns removed based on the specified condition.
    """
    # Find rows where the condition is met
    rows_to_check = data[data[condition_column].isin(condition_values)].index
    special_rows = data.loc[rows_to_check]

    # Get columns based on the specified operation and threshold
    if operation == '!=':
        filtered_columns = special_rows.loc[:, (special_rows != threshold).any()]
    elif operation == '<=':
        filtered_columns = special_rows.loc[:, (special_rows <= threshold).any()]
    elif operation == '>=':
        filtered_columns = special_rows.loc[:, (special_rows >= threshold).any()]
    elif operation == '=':
        filtered_columns = special_rows.loc[:, (special_rows == threshold).any()]
    else:
        raise ValueError("Unsupported operation. Supported operations are '!=', '<=', '>=', and '='")

    # Get column names
    column_names = filtered_columns.columns.tolist()

    # Remove columns from the DataFrame
    data_filtered = data.drop(column_names, axis=1)

    return data_filtered

def cv_filter(data, raw_list, threshold =20, operation = "<=" ):
    """
    Calculate the coefficient of variation (CV) for each column between given rows and filter based on threshold.

    Requirement:
    -  filter_rows function

    Parameters:
    - data: DataFrame
        The input DataFrame from which columns will be selected.
    - raw_list: list
        liste of rows to compare (ex: raw_list = ["sample_name == '240326NCE_Globale_neg_QC3-DIL8'", "sample_name == '240326NCE_Globale_neg_QC3'" , "sample_name == '240326NCE_Globale_neg_QC3-DIL2'"]).
    - threshold: int 
        cv threasholf of filtering (defolt = 20)
    - operation: "<=", ">="
        the fitration index: drop those how huer or less than threshold 

    Returns:
    - DataFrame
        A DataFrame data filtered based on the threshold given.
    """
    #filter QCDil2, QCDil8 and QC3 raws
    CV_data = filter_rows(data, include = raw_list)
    
    # Calculate the mean and standard deviation for each column
    mean_values = CV_data.mean()
    std_values = CV_data.std()
    
    # Calculate the coefficient of variation (CV) for each column
    cvs = (std_values / mean_values) * 100

    # Identify columns with CV greater than the threshold
    if operation == '<=':
        columns_to_drop = cvs[cvs <= threshold].index.tolist()
    elif operation == '>=':
        columns_to_drop = cvs[cvs >= threshold].index.tolist()

    # Drop the identified columns from the DataFrame
    filtered_data = data.drop(columns_to_drop, axis=1)
    
    return filtered_data

def QC_filter_with_zeros(data, threshold=0.75):
    """
    Drop columns from a DataFrame where the value is 0 in at least 3/4 of the QCs.

    Parameters:
    - data: DataFrame
        The input DataFrame from which columns will be removed.
    - threshold: float
        The threshold proportion of zeros in the column to trigger removal (default parameter = 0.75).

    Returns:
    - DataFrame
        A modified DataFrame with columns removed based on the specified condition.
    """
    # Calculate the proportion of zeros in each column
    zero_proportion = (data == 0).sum() / len(data)

    # Identify columns where the proportion of zeros exceeds the threshold
    columns_to_drop = zero_proportion[zero_proportion >= threshold].index.tolist()

    # Drop identified columns from the DataFrame
    data_filtered = data.drop(columns_to_drop, axis=1)

    return data_filtered


def save_as_csv(data, output_dir, output_file, file_conflict="skip"):
    """
    Save DataFrame as a CSV file in the given directory.

    Parameters:
        data (DataFrame): DataFrame to be saved.
        output_dir (str): Directory where the CSV file will be saved.
        output_file (str): Name of the CSV file (without the extension).
        file_conflict (str): Behavior in case of a file conflict.
            - "skip": Skip saving the file (default).
            - "replace": Replace the existing file.
            - "append": Append to the existing file.

    Returns:
        str: Path to the saved CSV file, or None if saving is skipped.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct the full path to the output file
    output_path = os.path.join(output_dir, output_file + ".csv")
    
    # Check if a file with the same name already exists
    if os.path.exists(output_path):
        if file_conflict == "skip":
            print("Saving skipped.")
            return None
        elif file_conflict == "replace":
            # Replace the existing file
            data.to_csv(output_path, index=True, sep = ';')
            print("Existing file replaced.")
        elif file_conflict == "append":
            # Append to the existing file
            existing_data = pd.read_csv(output_path)
            combined_data = pd.concat([existing_data, data], ignore_index=True)
            combined_data.to_csv(output_path, index=True, sep = ';')
            print("Data appended to the existing file.")
        else:
            print("Invalid value for 'file_conflict'. Skipping saving.")
            return None
    else:
        # Save the DataFrame as a new CSV file
        data.to_csv(output_path, index=True)
        print("CSV file saved.")
    
    return output_path


## POS

#data inport
dataPOS = read_file(inputPOSFile, sep = '\t')
metadataPOSData = read_file(metadataPOS, sep = ';', encoding='latin-1' )
#creating metabolite column start with M
dataPOS = add_metabolite_column(dataPOS)
excludDataPOS = ['Rt(min)', 'Mz']
fdataPOS = filter_column(dataPOS, exclude = excludDataPOS)
#transposing => metabolites as columns name
tdataPOS = transpose_data(fdataPOS, 'metabolite')
#exclud machine blc
exludRaw = ["sample_name == 'blc'","sample_name == 'blc_20240403164953'","sample_name == 'blc_20240404121923'", "sample_name == 'blc_20240404124839'", "sample_name == 'istd_ode'"]
ftdataPOS = filter_rows (tdataPOS,exclude=exludRaw)
#metadata merging 
excludDataNEG = ['id natif', 'class', 'injectionOrder']
metadataPOSData = filter_column(metadataPOSData,exclude=excludDataNEG)
metadataPOSData.set_index('sample_name', inplace = True)
start_dataPOS = merge_data(ftdataPOS, metadataPOSData)
# blank filtering 0 besed
blank_filtered_dataPOS= blank_filter(start_dataPOS)
# QC filtering: if 0 in 3/4 QC elimined
QC_blank_filterPOS = QC_filter_with_zeros (blank_filtered_dataPOS, 0.75)
# QCDil filtering CV <=10 elimined
raw_list = ["sample_name == '240326NCE_Globale_neg_QC3-DIL8'", "sample_name == '240326NCE_Globale_neg_QC3'" , "sample_name == '240326NCE_Globale_neg_QC3-DIL2'"]
QCDil_QC_blank_filterPOS = cv_filter (QC_blank_filterPOS, raw_list, 10, "<=")
# caving in a csv
output_pathPOS = save_as_csv(QCDil_QC_blank_filterPOS, output_dir=outputPath, output_file=outputPOSFile, file_conflict="replace")
print("CSV file saved at:", output_pathPOS)




## NEG

#data inport
dataNEG = read_file(inputNEGFile, sep = ';')
metadataNEGData = read_file(metadataNEG, sep = ';', encoding='latin-1' )
#creating metabolite column start with M
dataNEG = add_metabolite_column(dataNEG)
excludData = ['Rt(min)', 'Mz']
fdataNEG = filter_column(dataNEG, exclude=excludDataNEG)
#transposing => metabolites as columns name
tdataNEG = transpose_data(fdataNEG, 'metabolite')
#exclud machine blc
exludRawNEG = ["sample_name == 'blc'","sample_name == 'blc_20240403164953'","sample_name == 'blc_20240404121923'", "sample_name == 'blc_20240404124839'", "sample_name == 'istd_ode'"]
ftdataNEG = filter_rows (tdataNEG,exclude=exludRawNEG)
#metadata merging 
excludDataNEG = ['id natif', 'class', 'injectionOrder']
metadataNEGData = filter_column(metadataNEGData,exclude=excludDataNEG)
metadataNEGData.set_index('sample_name', inplace = True)
start_dataNEG = merge_data(ftdataNEG, metadataNEGData)
# blank filtering 0 besed
blank_filtered_dataNEG= blank_filter(start_dataNEG)
# QC filtering: if 0 in 3/4 QC elimined
QC_blank_filterNEG = QC_filter_with_zeros (blank_filtered_dataNEG, 0.75)
# QCDil filtering CV <=10 elimined
raw_list = ["sample_name == '240326NCE_Globale_neg_QC3-DIL8'", "sample_name == '240326NCE_Globale_neg_QC3'" , "sample_name == '240326NCE_Globale_neg_QC3-DIL2'"]
QCDil_QC_blank_filterNEG = cv_filter (QC_blank_filterNEG, raw_list, 10, "<=")
# caving in a csv
output_pathNEG = save_as_csv(QCDil_QC_blank_filterNEG, output_dir=outputPath, output_file=outputNEGFile, file_conflict="replace")
print("CSV file saved at:", output_pathNEG)

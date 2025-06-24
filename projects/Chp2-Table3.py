### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pandas as pd
import numpy as np
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
from scipy import stats


small_df_size = 37
starting_epidemic_year = 1994
tspan = np.arange(0, small_df_size, 1)
folder_path = os.path.join(base_path, 'output', 'Chp2_scenarios')
result_path = os.path.join(base_path, 'output', 'Chp2_tables')


# Function to extract F_W1_median from an Excel file
def extract_excel_files(df, small_df_size):
    number_of_small_dfs = len(df) // small_df_size
    writing_sols = []
    
    for i in range(number_of_small_dfs):
        start_index = starting_epidemic_year
        end_index = start_index + small_df_size
        small_df = df.iloc[i*small_df_size:(i+1)*small_df_size].copy()
        small_df.set_index(pd.RangeIndex(start=start_index, stop=end_index), inplace=True)


        writing_sols.append(small_df)
    
    return writing_sols


starting_epidemic_year = 1994
year_tostart_evaluation= 2025
ending_evaluation_year = 2030

def extract_data_from_dfs(small_dfs, column_or_formula, start_idx, end_idx):
    extracted_data = []
    
    for df in small_dfs:
        # Ensure that df has enough rows to perform the operation
        if df.shape[0] >= end_idx:
            if isinstance(column_or_formula, str):  # If a column name is provided
                # Compute a single sum value from start_idx to end_idx
                # value = df[column_or_formula].iloc[start_idx:end_idx].sum()
                value = df[column_or_formula].iloc[end_idx] - df[column_or_formula].iloc[start_idx]

            else:  # If a formula (function) is provided
                # Apply the formula, then compute a single sum value from start_idx to end_idx
                # value = df.apply(column_or_formula, axis=1).iloc[start_idx:end_idx].sum()

                # do a summation first, then apply the formula
                # value = df.iloc[start_idx:end_idx].sum().pipe(column_or_formula)
                df_copy = df.copy()
                df_copy = df_copy.drop('Simulation', axis=1)
                value = (df_copy.iloc[end_idx] - df_copy.iloc[start_idx]).pipe(column_or_formula)

        else:
            # Handle cases where df has fewer rows than end_idx (you might adjust this part as per your needs)
            value = np.nan
        
        extracted_data.append(value)
    
    return np.array(extracted_data)



def custom_round(val, threshold=100, decimals=2):
    """Round a value based on a threshold."""
    return round(val) if val > threshold else round(val, decimals)


def calculate_modes(data, axis=0):
    return stats.mode(np.round(data,decimals=0), axis=axis).mode[0]

def compute_summary_statistics(data, decimals=2):

    absolute_data = np.fabs(data)
    median_abs = custom_round(np.median(absolute_data, axis=0),threshold=100, decimals=1)
    mode_abs = custom_round(calculate_modes(absolute_data, axis=0),threshold=100, decimals=1)

    P2_5_abs = custom_round(np.percentile(absolute_data, 2.5, axis=0),threshold=100, decimals=1)
    P97_5_abs = custom_round(np.percentile(absolute_data, 97.5, axis=0),threshold=100, decimals=1)


    # Format the summary
    summary =   f"{median_abs:,} ({P2_5_abs:,}-{P97_5_abs:,})", median_abs, f"({P2_5_abs:,}-{P97_5_abs:,})"

    
    return summary

def subtract_dataframes(df1, df2):
    """
    Perform cell-wise subtraction of two dataframes, excluding the "Simulation" column.
    
    Parameters:
        df1 (pd.DataFrame): The first dataframe (minuend).
        df2 (pd.DataFrame): The second dataframe (subtrahend).
    
    Returns:
        pd.DataFrame: A dataframe with the subtracted values.
    """
    # Ensure the dataframes have the same shape and columns
    assert df1.shape == df2.shape, "DataFrames must have the same shape."
    assert all(df1.columns == df2.columns), "DataFrames must have the same columns."
    
    # Create a copy of the dataframes to keep original data intact
    result_df = df1.copy()
    
    # Perform cell-wise subtraction for numeric columns, excluding "Simulation"
    for col in df1.columns:
        if col != "Simulation" and pd.api.types.is_numeric_dtype(df1[col]):
            result_df[col] = df1[col] - df2[col]
    
    return result_df


excel_files = ["Baseline_Routine_2_WithDolutegravir.xlsx",
                ##### POCVL only, without DR testing
               "Baseline_POC_2_WithDolutegravir.xlsx",
               "Baseline_POC_3_WithDolutegravir.xlsx",
               "Baseline_POC_4_WithDolutegravir.xlsx",
               ##### POCVL and pre-treatment DR Testing
               "Timevarying_POC_2_WithDolutegravir.xlsx",
               "Timevarying_POC_3_WithDolutegravir.xlsx",
               "Timevarying_POC_4_WithDolutegravir.xlsx",
               #### POCVL with only DR testing at failure
               "Baseline_POC_AcquiredDR_2_WithDolutegravir.xlsx",
               "Baseline_POC_AcquiredDR_3_WithDolutegravir.xlsx",
               "Baseline_POC_AcquiredDR_4_WithDolutegravir.xlsx",
               ##### POCVL, pre-treatment DR Testing and DR testing at failure
               "Timevarying_POC_AcquiredDR_2_WithDolutegravir.xlsx",
               "Timevarying_POC_AcquiredDR_3_WithDolutegravir.xlsx",
               "Timevarying_POC_AcquiredDR_4_WithDolutegravir.xlsx"
            ]  

excel_files_labels = ["Baseline", 
                      ##### POCVL only, without DR testing
                      "POCVL ACT-UP and No Drug Resistance Testing", 
                      "POCVL 2 times ACT-UP level and No Drug Resistance Testing", 
                      "POCVL 3 times ACT-UP level and No Drug Resistance Testing", 
                      ##### POCVL and pre-treatment DR Testing
                      "POCVL and pre-treatment DR testing at ACT-UP level but no DR testing at failure",
                      "POCVL and pre-treatment DR testing at 2 times ACT-UP level but no DR testing at failure",
                      "POCVL and pre-treatment DR testing at 3 times ACT-UP level but no DR testing at failure",
                      #### POCVL with only DR testing at failure
                      "POCVL  with DR testing at failure at ACT-UP level but No pre-treatment DR testing", 
                      "POCVL  with DR testing at failure at 2 times ACT-UP level but No pre-treatment DR testing",
                      "POCVL  with DR testing at failure at 3 times ACT-UP level but No pre-treatment DR testing",
                      ##### POCVL, pre-treatment DR Testing and DR testing at failure
                      "POCVL, pre-treatment DR testing and DR testing at failure at ACT-UP level",
                      "POCVL, pre-treatment DR testing and DR testing at failure at 2 times ACT-UP level",
                      "POCVL, pre-treatment DR testing and DR testing at failure at 3 times ACT-UP level"
                      ]


data_sources = [




###### Start of dataframes that were used in the calculation sheet 

     lambda x: np.divide(x['VLtest_POC'],x['incidence'],out=np.zeros_like(x['incidence']),where= x['incidence'] != 0),
     lambda x: np.divide(x['VLtest_POC'],   x['newAcquiredDR'] + x['newTransmittedDR'],out=np.zeros_like(x['newAcquiredDR']),where= (x['newAcquiredDR'] + x['newTransmittedDR'])!= 0),

     lambda x: np.divide((x['preDRtest'] + x['postDRtest_routine'] + x['postDRtest_POC']),x['incidence_Resist'],out=np.zeros_like(x['incidence_Resist']),where= x['incidence_Resist'] != 0),
     lambda x: np.divide((x['preDRtest'] + x['postDRtest_routine'] + x['postDRtest_POC']),   x['newMisTreatedDR'],out=np.zeros_like(x['newMisTreatedDR']),where= (x['newMisTreatedDR'])!= 0), 





]

data_sources_labels = [

    'Number of POCVL test per infection averted',
    'Number of POCVL test per new DR Failure averted',

    'Number of DR test per new DR infections averted',
    'Number of DR test per wrongly prescribed DR PLHIV averted',

]




additional_sources = [

     'incidence',
     'incidence_Resist',

     lambda x: x['newAcquiredDR'] + x['newTransmittedDR'],

     'newMisTreatedDR'

]

additional_sources_labels = [

    'Number of new infection averted',
    'Number of new DR infections averted',

    'Number of new DR Failure averted',

    'Number of wrongly prescribed DR PLHIV averted',

]





cumulative_data_sources = [
    'incidence',
    'incidence_Resist',
    # 'newMisTreatedDR',
    lambda x: x['newAcquiredDR'] + x['newTransmittedDR'],
    'newMisTreatedDR',

    'VLtest_POC',
     lambda x: (x['preDRtest'] + x['postDRtest_routine'] + x['postDRtest_POC'])
]



cumulative_data_sources_labels = [
    'Newly infected',
    'Newly infected with Drug Resistance',
    # 'Number of Mistreated DR using first-line',
    'Total number of Failure due to Acquired and Transmitted DR',
    'Number of HIV DR wrongly prescribed first-line ART',

    'Number of Viral load levels test in POC',
    'Total number of DR test'
]






# Read the Excel file into a DataFrame
baseline_file = "Baseline_Routine_2_WithDolutegravir.xlsx"
baseline_file_path = os.path.join(folder_path, baseline_file)
baseline_df = pd.read_excel(baseline_file_path)



absolute_df = pd.DataFrame(columns=['Excel File'])
interval_df = pd.DataFrame(columns=['Excel File'])



for excel_file, excel_label in zip(excel_files, excel_files_labels):
    calculation_df = pd.read_excel(os.path.join(folder_path, excel_file))
    subtracted_df = subtract_dataframes(calculation_df, baseline_df)

    small_dfs = extract_excel_files(subtracted_df, small_df_size)
    absolute_row = {'Excel File': excel_label}  
    interval_row = {'Excel File': excel_label}  


    
    for data_source, label in zip(data_sources, data_sources_labels):
        data = extract_data_from_dfs(small_dfs, data_source,
                                     year_tostart_evaluation - starting_epidemic_year,
                                     ending_evaluation_year - starting_epidemic_year)
        
        absolute,median,interval = compute_summary_statistics(data)
        column_name = f"{label}"
        absolute_row[column_name] = absolute
        interval_row[column_name] = interval

        
    
    # Use concat to add the row to the summary_df
    absolute_df = pd.concat([absolute_df, pd.DataFrame([absolute_row])], ignore_index=True)
    interval_df = pd.concat([interval_df, pd.DataFrame([interval_row])], ignore_index=True)



##### Calculate the "subtracted" dataframes 
##### Results are used to represent the relative number of VL and DR tests done per new infections, new DR.... averted 

# output_file_path = os.path.join(result_path, "Final Table 3_Absolute.xlsx")
# absolute_df.to_excel(output_file_path, index=False)
# absolute_df



absolute_df_cumulative = pd.DataFrame(columns=['Excel File'])
median_df = pd.DataFrame(columns=['Excel File'])

for excel_file, excel_label in zip(excel_files, excel_files_labels):
    calculation_df = pd.read_excel(os.path.join(folder_path, excel_file))

    small_dfs = extract_excel_files(calculation_df, small_df_size)
    absolute_row = {'Excel File': excel_label}  
    median_row = {'Excel File': excel_label} 

    
    for data_source, label in zip(cumulative_data_sources, cumulative_data_sources_labels):
        data = extract_data_from_dfs(small_dfs, data_source,
                                     year_tostart_evaluation - starting_epidemic_year,
                                     ending_evaluation_year - starting_epidemic_year)
        
        absolute,median,interval = compute_summary_statistics(data)
        column_name = f"{label}"
        absolute_row[column_name] = absolute
        median_row[column_name]  = median

        
    
    # Use concat to add the row to the summary_df
    absolute_df_cumulative = pd.concat([absolute_df_cumulative, pd.DataFrame([absolute_row])], ignore_index=True)
    median_df = pd.concat([median_df, pd.DataFrame([median_row])], ignore_index=True)

##### Calculate the cumulative number of New infections, new VL tests or new DR tests done during a period 
##### Essentially, the total cumulative number of tests and VL done 
##### These are then included in the new table. at the end 

# output_file_path = os.path.join(result_path, "Final Table 3_Absolute(cumulative).xlsx")
# absolute_df_cumulative.to_excel(output_file_path, index=False)
# absolute_df_cumulative


##### Processing of the two above dataframes 
##### The dataframes would then be used in the manuscript, 
##### after modifying, changing column names and other reviews 


#### Final display to in Manuscript 

concatenated_df = pd.concat(
    [absolute_df_cumulative,
    absolute_df.drop(absolute_df.columns[0],axis=1)],axis=1
)


### Actual table that was used in the Manuscript Processed --> other information can be updated and changed in the previous codes 
### in the above section 

# output_file_path = os.path.join(result_path, "Final Table 3_ActualManuscript.xlsx")
# concatenated_df.to_excel(output_file_path, index=False)
# concatenated_df



### revising for the new dataframe 
##### Subtracting from the baseline to print out the differences in median used in the Table 3 
##### (the calculated median is not correctly equal to the calculation in the Table)
median_df.iloc[0:, 1:] = median_df.iloc[0:, 1:].sub(median_df.iloc[0, 1:])
# converting to absolute values
median_df.iloc[0:, 1:] = median_df.iloc[0:, 1:].abs()


###### Overwrite the median estimates (as calculated directly from median of differences) to the differences in medians of reporting statistics --> then calculating the ratios in the following codes. 
result_df = pd.DataFrame()
result_df['Excel File'] = median_df.iloc[:,0]

result_df['Number of POCVL test per infection averted'] = (median_df.iloc[:, 5] / median_df.iloc[:, 1]).astype(float).round(1).astype(str)
result_df['Number of POCVL test per new DR Failure averted'] = (median_df.iloc[:, 5] / median_df.iloc[:, 3]).astype(float).round(1).astype(str)
result_df['Number of DR test per new DR infections averted'] = (median_df.iloc[:, 6] / median_df.iloc[:, 2]).astype(float).round(1).astype(str)
result_df['Number of DR test per wrongly prescribed DR PLHIV averted'] = (median_df.iloc[:, 6] / median_df.iloc[:, 4]).astype(float).round(1).astype(str)


print_df=  result_df + ' ' + interval_df


final_df = pd.concat(
    [absolute_df_cumulative,
    print_df.drop(print_df.columns[0],axis=1)],axis=1
)

output_file_path = os.path.join(result_path, "Chapter 2 - Table 3.xlsx")
final_df.to_excel(output_file_path, index=False)
final_df








####### Additional analyses 

# Read the Excel file into a DataFrame
baseline_file = "Baseline_Routine_2_WithDolutegravir.xlsx"
baseline_file_path = os.path.join(folder_path, baseline_file)
baseline_df = pd.read_excel(baseline_file_path)



absolute_df = pd.DataFrame(columns=['Excel File'])
interval_df = pd.DataFrame(columns=['Excel File'])



for excel_file, excel_label in zip(excel_files, excel_files_labels):
    calculation_df = pd.read_excel(os.path.join(folder_path, excel_file))
    subtracted_df = subtract_dataframes(calculation_df, baseline_df)

    small_dfs = extract_excel_files(subtracted_df, small_df_size)
    absolute_row = {'Excel File': excel_label}  
    interval_row = {'Excel File': excel_label}  


    
    for data_source, label in zip(additional_sources, additional_sources_labels):
        data = extract_data_from_dfs(small_dfs, data_source,
                                     year_tostart_evaluation - starting_epidemic_year,
                                     ending_evaluation_year - starting_epidemic_year)
        
        absolute,median,interval = compute_summary_statistics(data)
        column_name = f"{label}"
        absolute_row[column_name] = absolute
        interval_row[column_name] = interval

        
    
    # Use concat to add the row to the summary_df
    absolute_df = pd.concat([absolute_df, pd.DataFrame([absolute_row])], ignore_index=True)
    interval_df = pd.concat([interval_df, pd.DataFrame([interval_row])], ignore_index=True)


###### Overwrite the median estimates (as calculated directly from median of differences) to the differences in medians of reporting statistics --> then calculating the ratios in the following codes. 
result_df = pd.DataFrame()
result_df['Excel File'] = interval_df.iloc[:,0]


result_df['Number of new infection averted'] = median_df.iloc[:, 1].astype(float).round(0).apply(lambda x: f"{int(x):,}")
result_df['Number of new DR infections averted'] = median_df.iloc[:, 2].astype(float).round(0).apply(lambda x: f"{int(x):,}")
result_df['Number of new DR Failure averted'] = median_df.iloc[:, 3].astype(float).round(0).apply(lambda x: f"{int(x):,}")
result_df['Number of wrongly prescribed DR PLHIV averted'] =median_df.iloc[:, 4].astype(float).round(0).apply(lambda x: f"{int(x):,}")

print_df=  result_df + ' ' + interval_df



##### Calculate the "subtracted" dataframes 
##### Results are used to represent the relative number of VL and DR tests done per new infections, new DR.... averted 

output_file_path = os.path.join(result_path, "Chapter 2 - Additional comparative analyses.xlsx")
print_df.to_excel(output_file_path, index=False)


absolute_df
print_df
print('Remember to enter the last two dataframes into console for checking!')
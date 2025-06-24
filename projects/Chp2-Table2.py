### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pandas as pd
import numpy as np
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


small_df_size = 37
starting_epidemic_year = 1994
tspan = np.arange(0, small_df_size, 1)
# folder_path = "output/Posterior_results_B/All_Posteriors"
folder_path = os.path.join(base_path, 'output', 'Chp2_scenarios')
# result_path = "output/Posterior_results_B/Result_excel"
result_path = os.path.join(base_path, 'output', 'Chp2_tables')

year_starting_tablereference = 2020

# Function to extract F_W1_median from an Excel file
def extract_excel_files(excel_file, small_df_size):
    full_file_path = os.path.join(folder_path, excel_file)
    df = pd.read_excel(full_file_path)
    number_of_small_dfs = len(df) // small_df_size
    writing_sols = []
    
    for i in range(number_of_small_dfs):
        start_index = starting_epidemic_year
        end_index = start_index + small_df_size
        small_df = df.iloc[i*small_df_size:(i+1)*small_df_size].copy()
        small_df.set_index(pd.RangeIndex(start=start_index, stop=end_index), inplace=True)

        # Calculating the gradient of 'incidence' with respect to 'tspan' --> add new 
        small_df['newly_infected'] = np.gradient(small_df['incidence'], tspan)
        small_df['yearly_deaths'] = np.gradient(small_df['deaths'],tspan)
        small_df['newly_VFduetoDR'] = np.gradient(small_df['newAcquiredDR'],tspan) + np.gradient(small_df['newTransmittedDR'],tspan)

        small_df['newly_infectedDR'] =  np.gradient(small_df['incidence_Resist'], tspan)


        small_df['newly_DRtests'] = np.gradient(small_df['preDRtest'],tspan) + np.gradient(small_df['postDRtest_routine'],tspan) +np.gradient(small_df['postDRtest_POC'],tspan)
        small_df['newly_VLtest_routine'] = np.gradient(small_df['VLtest_routine'],tspan) 
        small_df['newly_VLtest_POC'] = np.gradient(small_df['VLtest_POC'],tspan)
        small_df['newly_VLtest'] = np.gradient(small_df['VLtest_routine'],tspan) + np.gradient(small_df['VLtest_POC'],tspan)

        writing_sols.append(small_df)
    
    return writing_sols


def extract_data_from_dfs(small_dfs, column_or_formula):
    extracted_data = []
    
    for df in small_dfs:
        if isinstance(column_or_formula, str):  # If a column name is provided
            values = df[column_or_formula].values
        else:  # If a formula (function) is provided
            values = df.apply(column_or_formula, axis=1).values
        extracted_data.append(values)
    
    return np.array(extracted_data)


def custom_round(val, threshold=100, decimals=1):
    """Round a value based on a threshold."""
    return round(val) if val > threshold else round(val, decimals)


def compute_summary_statistics(data, index, decimals=1, is_percentage=False):
    # Compute the median
    median_val = custom_round(val = np.median(data, axis=0)[index], decimals = decimals)
    
    # Compute the 2.5th and 97.5th percentiles
    P2_5 = custom_round(val = np.percentile(data, 2.5, axis=0)[index], decimals = decimals)
    P97_5 = custom_round(val = np.percentile(data, 97.5, axis=0)[index], decimals = decimals)

    # Format the summary based on the is_percentage flag
    if is_percentage:
        summary = f"{median_val}% ({P2_5}%-{P97_5}%)"
    else:
        summary = f"{median_val:,} ({P2_5:,}-{P97_5:,})"
    
    return summary



def compute_median(data,index, decimals=1):
    # Compute the median
    median_val = custom_round(val = np.median(data, axis=0)[index], decimals = decimals)

    return median_val



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
               "Timevarying_POC_AcquiredDR_4_WithDolutegravir.xlsx",
               "Baseline_Routine_2_NoDolutegravir.xlsx"
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
                      "POCVL, pre-treatment DR testing and DR testing at failure at 3 times ACT-UP level",

                      "No Dolutegravir"
                      ]


indices = [year_starting_tablereference-starting_epidemic_year, 
           2025-starting_epidemic_year, 
           2030 - starting_epidemic_year]

data_sources = [
    'Total',
        'newly_infected',
                lambda x: (x["F_W1"]+ x["F_W2"] + x["F_TDR1"] + x["F_ADR1"]),
                    lambda x: (x["D_DR"] + x["I_DR"]),
    'DResistance',

     'newly_DRtests' ,
     'newly_VLtest_routine',
     'newly_VLtest_POC',
     'newly_VLtest',


    'yearly_deaths',
    lambda x: (x["I_DR"]+x["D_DR"])/(x["I"]+x["D"])*100,

    'newly_VFduetoDR',
    lambda x: np.divide(x["T"],x["T"] + x["F"],out=np.zeros_like(x["T"] + x["F"]),where= x["T"] + x["F"] != 0)*100,

    lambda x: x["F_TDR1"] + x["F_ADR1"],

    'newly_infectedDR',
    lambda x: (x["T"]+x["F"])

]
data_sources_labels = [
    'Total Number of PLHIV in PNG',
    'Newly infected',
                'VF in Total (excluded second-line nonadherence)',
    'Number of pre-treatment (I&D) Drug Resistance',
    'Total number of PLHIV with Drug Resistant HIV',

    'Total number of DR test',
    'Number of Viral load levels test in Routine',
    'Number of Viral load levels test in POC',
    'Total number of Viral load levels test',


    'Total deaths due to HIV/AIDS', 
    'DR prevalence in Treatment-naive HIV',

    'Newly VF due to Drug Resistance',
    'Viral suppression level %',

    'VF due to drug resistance (accumulated)',
    'Newly infected with Drug Resistance',
        'Total number on Treatment in PNG'

]



percentage_labels = [
    'DR prevalence in Treatment-naive HIV',
    'Viral suppression level %'
]



summary_df = pd.DataFrame(columns=['Excel File'] + [f"{label} {starting_epidemic_year + idx}" for label in data_sources_labels for idx in indices])
median_df = pd.DataFrame(columns=['Excel File'] + [f"{label} {starting_epidemic_year + idx}" for label in data_sources_labels for idx in indices])

for excel_file, excel_label in zip(excel_files, excel_files_labels):
    small_dfs = extract_excel_files(excel_file, small_df_size)
    summary_row = {'Excel File': excel_label}  
    median_row = {'Excel File': excel_label}
    
    for data_source, label in zip(data_sources, data_sources_labels):
        data = extract_data_from_dfs(small_dfs, data_source)
        is_percentage = label in percentage_labels  # Check if the label requires percentage formatting

        for index in indices:
            summary = compute_summary_statistics(data, index, is_percentage=is_percentage)
            median = compute_median(data, index)
            column_name = f"{label} {starting_epidemic_year + index}"
            summary_row[column_name] = summary
            median_row[column_name] = median
        
    # Use concat to add the row to the summary_df
    summary_df = pd.concat([summary_df, pd.DataFrame([summary_row])], ignore_index=True)
    median_df = pd.concat([median_df, pd.DataFrame([median_row])], ignore_index=True)


# output_file_path = os.path.join(result_path, "Final Table 2_Absolute.xlsx")
# summary_df.to_excel(output_file_path, index=False)
# summary_df


# Extract the "Reference" row
reference_row = median_df.iloc[0]

string_column = median_df.iloc[:, 0] 
numeric_df = median_df.iloc[:, 1:]

percentage_decrease = (reference_row.iloc[1:] - numeric_df) / reference_row.iloc[1:] * 100
percentage_decrease = percentage_decrease.apply(pd.to_numeric, errors='coerce')

formatted_numeric_df = numeric_df.astype(str) + " (" + percentage_decrease.round(2).astype(str) + "%)"

# Concatenating the string column back
final_median_df = pd.concat([string_column, formatted_numeric_df], axis=1)

# final_median_df
# output_file_path = os.path.join(result_path, "Final Table 2_PercentageReduction.xlsx")
# final_median_df.to_excel(output_file_path, index=False)
# final_median_df




selected_column_names = [    'Total Number of PLHIV in PNG',
    'Newly infected',
    'Newly infected with Drug Resistance',
    'Viral suppression level %',
    'Newly VF due to Drug Resistance'

]

# The suffixes you are interested in
years_selected = ['2025', '2030']

# Generate the new list with the full column names
full_column_names = [f'{base} {suffix}' for base in selected_column_names for suffix in years_selected]

final_dataset = summary_df[full_column_names]
output_file_path = os.path.join(result_path, "Chapter 2 - Table 2.xlsx")
final_dataset.to_excel(output_file_path, index=False)
final_dataset



###### Line of code to calculate the difference between two years 2020 and 2030

def calculate_median_and_CI(data,digits=1):
    median_val = np.round(np.median(data),decimals=digits)
    P2_5 = np.round(np.percentile(data, 2.5),decimals=digits)
    P97_5 = np.round(np.percentile(data, 97.5),decimals=digits)
    return median_val, P2_5, P97_5



# Load and process the data
excel_file = "Baseline_Routine_2_WithDolutegravir.xlsx"
small_dfs = extract_excel_files(excel_file, small_df_size)


# Sum the columns "T" and "F" for the years 2020 and 2030
total_2020 = np.array([df.loc[2020, 'T'] + df.loc[2020, 'F'] for df in small_dfs])
total_2030 = np.array([df.loc[2030, 'T'] + df.loc[2030, 'F'] for df in small_dfs])

# Compute the differences
percentage_diff = (total_2030 - total_2020)/total_2020*100

displayed_median_diff = np.round((np.median(total_2030) - np.median(total_2020))/np.median(total_2020)*100, decimals=1)

# Calculate the median and 95% CI for the differences
median_diff, P2_5_diff, P97_5_diff = calculate_median_and_CI(percentage_diff)

# Format the output
output = f"{displayed_median_diff}% ({P2_5_diff}%-{P97_5_diff}%)"
print("Difference between year 2030 and 2020 for total of columns T and F:", output)




# Load and process the data
excel_file = "Baseline_Routine_2_NoDolutegravir.xlsx"
small_dfs = extract_excel_files(excel_file, small_df_size)


# Sum the columns "T" and "F" for the years 2020 and 2030
number_2020 = np.array([df.loc[2020, 'newly_VFduetoDR'] for df in small_dfs])
number_2030 = np.array([df.loc[2030, 'newly_VFduetoDR']  for df in small_dfs])

# Compute the differences
percentage_diff = (number_2030 - number_2020)/number_2020*100

displayed_median_diff = np.round((np.median(number_2030) - np.median(number_2020))/np.median(number_2020)*100, decimals=1)

# Calculate the median and 95% CI for the differences
median_diff, P2_5_diff, P97_5_diff = calculate_median_and_CI(percentage_diff)

# Format the output
output = f"{displayed_median_diff}% ({P2_5_diff}%-{P97_5_diff}%)"
print("Difference between year 2030 and 2020 for difference in newly_VFduetoDR:", output)

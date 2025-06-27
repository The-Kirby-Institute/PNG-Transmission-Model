### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from model.meta_parameters import *
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
# import matplotlib.image as mpimg
# from matplotlib.ticker import ScalarFormatter
# from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import os


from model.Model_init import *



import random

from matplotlib.ticker import MaxNLocator

from EconEval.cea_functions import *
from EconEval.prb_data_matrix import*
import copy
from datetime import datetime


start_time = datetime.now()
print(f"Starting  code execution at {start_time}...")


random.seed(334)


# folder_path = "CEA_results/Scenarios_results"
folder_path = os.path.join(base_path, 'output', 'Chp3_scenarios')



def extract_excel_files(excel_file, small_df_size):
    df = pd.read_excel(excel_file)
    number_of_small_dfs = len(df) // small_df_size
    writing_sols = [df.iloc[i*small_df_size:(i+1)*small_df_size] for i in range(number_of_small_dfs)]
    return writing_sols


# def setCorrectxAxis(frame, frequency_ticks =5,starting_position=0):
#     plt.figure()
#     default_x_ticks = list(range(len(frame.I_W)))
#     new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+len(frame.I_W) ))
#     plt.xticks(default_x_ticks,new_x_ticks)
#     plt.xticks(np.append(np.arange(starting_position, len(frame.I_W)+1, frequency_ticks),default_x_ticks[-1]))
#     plt.xticks(fontsize=8, rotation=45)



def report_statistics(input_list, scale='actual', fixed_value=np.nan, deci =0):
    # Compute the median
    if np.isnan(fixed_value):
        median_value = np.median(input_list)
    else:
        median_value = fixed_value
    
    # Compute the 95% CI using the 2.5th and 97.5th percentiles
    lower_bound = np.percentile(input_list, 2.5)
    upper_bound = np.percentile(input_list, 97.5)
    
    # Scaling the values
    if scale == 'millions':
        median_value /= 1_000_000
        lower_bound /= 1_000_000
        upper_bound /= 1_000_000
    elif scale == 'thousands':
        median_value /= 1_000
        lower_bound /= 1_000
        upper_bound /= 1_000
    elif scale == 'percentage':
        median_value *= 100
        lower_bound *= 100
        upper_bound *= 100


    # Rounding to 1 decimal place
    median_value = round(median_value, deci)
    lower_bound = round(lower_bound, deci)
    upper_bound = round(upper_bound, deci)
    
    # Return the formatted string with commas as thousand separators
    if deci == 0:
        return_str = f'{median_value:,.0f} ({lower_bound:,.0f} - {upper_bound:,.0f})'
    else: 
        return_str = f'{median_value:,.1f} ({lower_bound:,.1f} - {upper_bound:,.1f})'
    return return_str



import os
import numpy as np
import pandas as pd

def calculate_pairwise_icer(file_names, folder_path, p_costs, p_daly, discount_rate, year_start_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies):

    incremental_costs = []
    incremental_effects = []
    icers = []
    comparisons = []
    median_costs = []
    median_effects = []

    # Loop through the file list, comparing each file with the previous one
    for i in range(1, len(file_names)):
        # Extract data for the current and previous files
        prev_file = extract_excel_files(os.path.join(folder_path, file_names[i - 1]), year_end_evaluation - starting_epidemic_year + 1)
        curr_file = extract_excel_files(os.path.join(folder_path, file_names[i]), year_end_evaluation - starting_epidemic_year + 1)

        # Calculate cost and effect for the previous file
        prev_cost = np.array(probabilistic_all_populations(prev_file, p_costs, discount_rate, year_start_policy - starting_epidemic_year, year_end_evaluation - starting_epidemic_year, tspan=tspan))
        prev_effect = np.array(probabilistic_all_populations(prev_file, p_daly, discount_rate, year_start_policy - starting_epidemic_year, year_end_evaluation - starting_epidemic_year, tspan=tspan))

        # Calculate cost and effect for the current file
        curr_cost = np.array(probabilistic_all_populations(curr_file, p_costs, discount_rate, year_start_policy - starting_epidemic_year, year_end_evaluation - starting_epidemic_year, tspan=tspan))
        curr_effect = np.array(probabilistic_all_populations(curr_file, p_daly, discount_rate, year_start_policy - starting_epidemic_year, year_end_evaluation - starting_epidemic_year, tspan=tspan))

        # Compute incremental cost and effect
        incremental_cost = curr_cost - prev_cost
        incremental_effect = -(curr_effect - prev_effect)

        # Compute ICER
        icer = incremental_cost / incremental_effect

        # Calculate medians
        median_cost_prev = np.median(prev_cost)
        median_cost_curr = np.median(curr_cost)
        median_effect_prev = np.median(prev_effect)
        median_effect_curr = np.median(curr_effect)

        ### Store results for each pair 
        comparisons.extend([f"{naming_strategies[i - 1]}",  f"{naming_strategies[i]}", np.nan])
        median_costs.extend([report_statistics(prev_cost,scale = "millions",deci=1), 
                             report_statistics(curr_cost,scale = "millions",deci=1),np.nan])
        median_effects.extend([report_statistics(prev_effect,scale = "thousands"), 
                               report_statistics(curr_effect,scale = "thousands"),np.nan])
        incremental_costs.extend(["-", 
                                  report_statistics(incremental_cost,scale = "millions",fixed_value=median_cost_curr - median_cost_prev,deci=1),np.nan])
        incremental_effects.extend(["-", 
                                    report_statistics(incremental_effect,fixed_value= -(median_effect_curr - median_effect_prev)),np.nan])
        icers.extend(["-", 
                      report_statistics(icer,scale = "actual", fixed_value= (median_cost_curr - median_cost_prev)/-(median_effect_curr - median_effect_prev)),np.nan])

    # Create a DataFrame to store the results
    return pd.DataFrame({
        'Comparison': comparisons,
        # 'Strategy': strategies,
        'Median Cost': median_costs,
        'Incremental Cost': incremental_costs,
        'Median Effect': median_effects,
        'Incremental Effect': incremental_effects,
        'ICER': icers
    })



discount_rate = 5/100
year_start_policy = 2020
year_end_evaluation = 2100

year_end_results = 2100
years = year_end_results - starting_epidemic_year +1
tspan = np.arange(0, years, 1)


naming_strategies = [
    "Status-quo scenario"
    ,"20% access to POC VL testing and DR testing to confirm treatment failure" 
    ,"40% access to POC VL testing and DR testing to confirm treatment failure"
    ,"60% access to POC VL testing and DR testing to confirm treatment failure"
    ,"Expansion scenario"
]


###### Main strategies in the modelling 
#### 2020 - 2050 period
strategies_on_cea_frontier =  ["Routine_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_56_56.xlsx"]


# dframe_result30 = calculate_pairwise_icer(strategies_on_cea_frontier, folder_path, p_costs, p_daly, discount_rate, year_start_policy, year_end_evaluation, starting_epidemic_year, tspan,[naming_strategies[0]] + naming_strategies[-3:])

### Strategies on the CEA frontier is different dependent on the costs associated 
dframe_result30 = calculate_pairwise_icer(strategies_on_cea_frontier, folder_path, p_costs, p_daly, discount_rate, year_start_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies)



#### 2020 - 2030 period
strategies_on_cea_frontier =  ["Routine_2_Baseline_WithDolutegravir_36_36.xlsx"
                            #    ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_36_36.xlsx"
                            #    ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_36_36.xlsx"]


dframe_result10 = calculate_pairwise_icer(strategies_on_cea_frontier, folder_path, p_costs, p_daly, discount_rate, year_start_policy, year_end_evaluation, starting_epidemic_year, tspan,[naming_strategies[0]] + naming_strategies[-2:])

### Strategies on the CEA frontier is different dependent on the costs associated 
# dframe_result10 = calculate_pairwise_icer(strategies_on_cea_frontier, folder_path, p_costs, p_daly, discount_rate, year_start_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies)



#### 2020 - 2035 period
strategies_on_cea_frontier =  ["Routine_2_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_41_41.xlsx"]


# dframe_result15 = calculate_pairwise_icer(strategies_on_cea_frontier, folder_path, p_costs, p_daly, discount_rate, year_start_policy, year_end_evaluation, starting_epidemic_year, tspan,[naming_strategies[0]] + naming_strategies[-3:])

### Strategies on the CEA frontier is different dependent on the costs associated 
dframe_result15 = calculate_pairwise_icer(strategies_on_cea_frontier, folder_path, p_costs, p_daly, discount_rate, year_start_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies)


output_path = os.path.join(base_path, 'output', 'Chp3_tables', 'Chp3 - Table3.xlsx')
with pd.ExcelWriter(output_path) as writer:
    dframe_result30.to_excel(writer, sheet_name='Horizon 2020-2050', index=False)
    dframe_result10.to_excel(writer, sheet_name='Horizon 2020-2030', index=False)
    dframe_result15.to_excel(writer, sheet_name='Horizon 2020-2035', index=False)





end_time = datetime.now()
print(f"Code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")

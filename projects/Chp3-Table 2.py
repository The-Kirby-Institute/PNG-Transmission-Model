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


starting_epidemic_year  = 1994
# plt.style.use('ggplot')

reference_year = 2000
adj =  reference_year - starting_epidemic_year

year_end_results = 2100
years = year_end_results - starting_epidemic_year +1
tspan = np.arange(0, years, 1)
discount_rate = 5/100

#### Costing done based on the discounted into the future numbers 
folder_path = os.path.join(base_path, 'output', 'Chp3_scenarios')
corresponding_Strategies = ["Status-quo scenario"
,"Current implementation scenario - 20% access to POC VL testing, no DR testing"
,"40% access to POC VL testing , no DR testing"
,"60% access to POC VL testing , no DR testing"
,"20% access to POC VL testing and DR testing before treatment initiation"
,"40% access to POC VL testing and DR testing before treatment initiation"
,"60% access to POC VL testing and DR testing before treatment initiation"
,"20% access to POC VL testing and DR testing to confirm treatment failure" 
,"40% access to POC VL testing and DR testing to confirm treatment failure"
,"60% access to POC VL testing and DR testing to confirm treatment failure"
,"20% access to POC VL testing, with DR testing both before treatment initiation and confirm treatment failure"
,"40% access to POC VL testing, with DR testing both before treatment initiation and confirm treatment failure"
,"Expansion scenario - 60% access to POC VL testing, with DR testing both before treatment initiation and confirm treatment failure"
]
matching_groups = ['G1'] * 3 + ['G2'] * 3 + ['G3'] * 3 + ['G4'] * 3



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



def report_statistics(input_list, scale='actual',fixed_value = np.nan,deci =0):
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





def return_all_results_modelling(file_names, p_costs, p_daly,p_incidences,p_deaths,p_ADR,p_TDR,p_secondline,p_mistreated,p_ADR_Total,p_TDR_Total,p_incidences_DR,p_Total_TFailure,p_Viral_unsuppressed,p_new_TFailure, discount_rate, year_start_policy,year_end_evaluation,year_end_policy, starting_epidemic_year, tspan):
    median_total_costs = []
    median_total_effects = []
    report_incremental_costs = [0]
    report_incremental_effects = [0]
    percentage_cost = [0]
    ICERs = [0]

    percentage_effect = [0]
    #### Other outputs to be added in the frame

    
    comparator = extract_excel_files(os.path.join(folder_path, file_names[0]), year_end_evaluation - starting_epidemic_year + 1)
    comparator_cost = np.array(probabilistic_all_populations(comparator, p_costs, discount_rate,year_start_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
    comparator_effect = np.array(probabilistic_all_populations(comparator, p_daly, discount_rate,year_start_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
    
    #### additional outcomes to calculate in the modelling 
    #### all additional outcomes should not be discounted , unlike costs and DALYs
    
    ### Add in the values of median total cost and effects
    median_total_costs.append(report_statistics(comparator_cost,scale = "millions",deci=1))
    median_total_effects.append(report_statistics(comparator_effect,scale = 'thousands'))

    median_comparator_cost = np.median(comparator_cost)
    median_comparator_effect = np.median(comparator_effect)

    #### Add in additional outcomes in the modell


    for file_name in file_names[1:]:
        evaluation = extract_excel_files(os.path.join(folder_path, file_name), year_end_evaluation - starting_epidemic_year + 1)

        evaluation_cost = np.array(probabilistic_all_populations(evaluation, p_costs, discount_rate,year_start_policy-starting_epidemic_year,year_end_evaluation-starting_epidemic_year,tspan=tspan))
        evaluation_effect = np.array(probabilistic_all_populations(evaluation, p_daly, discount_rate,year_start_policy-starting_epidemic_year,year_end_evaluation-starting_epidemic_year,tspan=tspan))

        median_evaluation_cost = np.median(evaluation_cost)
        median_evaluation_effect = np.median(evaluation_effect)

        ### These are reporting dataframes only
        report_incremental_costs.append(report_statistics(evaluation_cost- comparator_cost,scale ='millions',fixed_value= median_evaluation_cost - median_comparator_cost,deci=1))
        report_incremental_effects.append(report_statistics(-(evaluation_effect-comparator_effect),scale = "actual", fixed_value= -(median_evaluation_effect - median_comparator_effect)))

        ### Add in the total effects
        median_total_costs.append(report_statistics(evaluation_cost,scale ="millions",deci=1))
        median_total_effects.append(report_statistics(evaluation_effect,scale ='thousands'))

        ### Percentage of costs 
        # percentage_cost.append(report_statistics((evaluation_cost- comparator_cost)/comparator_cost,scale ="percentage", fixed_value= (median_evaluation_cost - median_comparator_cost)/median_comparator_cost))
        ### Percentage of DALYs
        # percentage_effect.append(report_statistics(-(evaluation_effect-comparator_effect)/comparator_effect,scale ="percentage", fixed_value= -(median_evaluation_effect-median_comparator_effect)/median_comparator_effect))
        ### Calculate the ICER column
        ICERs.append(report_statistics((evaluation_cost- comparator_cost)/-(evaluation_effect-comparator_effect),scale ="actual", fixed_value= (median_evaluation_cost- median_comparator_cost)/-(median_evaluation_effect-median_comparator_effect)))








    return pd.DataFrame({
        'Strategies': corresponding_Strategies,
        'Total Cost (in millions $)':median_total_costs,
        'Incremental Cost (in thousands $)': report_incremental_costs,
        # '% of Total Costs':percentage_cost,
        'Total Effect (in thousands DALYs)': median_total_effects,
        'Incremental Effect': report_incremental_effects,
        # '% of Total Effect':percentage_effect,
        'ICERs': ICERs
    })






year_start_policy = 2020
year_end_evaluation = 2100
year_end_policy = 2050
### GDP per capita, PPP (constant 2021 international $) 
WTP= 4_119
file_names = ["Routine_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_3_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_4_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_2_Timevarying_WithDolutegravir_56_56.xlsx"
                               ,"POC_3_Timevarying_WithDolutegravir_56_56.xlsx"
                               ,"POC_4_Timevarying_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_2_Timevarying_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_3_Timevarying_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_56_56.xlsx"]


# file_names = ["Routine_2_Baseline_WithDolutegravir_36_36.xlsx"
#                                ,"POC_2_Baseline_WithDolutegravir_36_36.xlsx"
#                                ,"POC_3_Baseline_WithDolutegravir_36_36.xlsx"
#                                ,"POC_4_Baseline_WithDolutegravir_36_36.xlsx"
#                                ,"POC_2_Timevarying_WithDolutegravir_36_36.xlsx"
#                                ,"POC_3_Timevarying_WithDolutegravir_36_36.xlsx"
#                                ,"POC_4_Timevarying_WithDolutegravir_36_36.xlsx"
#                                ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_36_36.xlsx"
#                                ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_36_36.xlsx"
#                                ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_36_36.xlsx"
#                                ,"POC_AcquiredDR_2_Timevarying_WithDolutegravir_36_36.xlsx"
#                                ,"POC_AcquiredDR_3_Timevarying_WithDolutegravir_36_36.xlsx"
#                                ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_36_36.xlsx"]


display_frame = return_all_results_modelling(file_names,p_costs, p_daly,p_incidences,p_deaths,p_ADR,p_TDR,p_secondline,p_mistreated,p_ADR_Total,p_TDR_Total,p_incidences_DR,p_Total_TFailure,p_Viral_unsuppressed,p_new_TFailure, discount_rate, year_start_policy,year_end_evaluation, year_end_policy,
                                              starting_epidemic_year, tspan)


#### add percentage signs into two columns of data 
#### remove the percentage columns 
# display_frame['% of Total Costs'] = display_frame['% of Total Costs'].str.replace(r'(\d+\.\d+)', r'\1%', regex=True)
# display_frame['% of Total Effect'] = display_frame['% of Total Effect'].str.replace(r'(\d+\.\d+)', r'\1%', regex=True)
output = display_frame
### Output to Excel for Table 2 results 5 years
output.to_excel(os.path.join(base_path, 'output', 'Chp3_tables', 'Chp3 - Table2.xlsx'))


end_time = datetime.now()
print(f"Code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")

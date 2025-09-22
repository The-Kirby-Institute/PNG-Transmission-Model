### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

#### Efficency frontier diagrams - Figure 2 (Python)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.image as mpimg
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import os
import re


import random

from matplotlib.ticker import MaxNLocator

from EconEval.cea_functions import *
from EconEval.prb_data_matrix import*
import copy
from datetime import datetime


start_time = datetime.now()
print(f"Starting  code execution at {start_time}...")


random.seed(334)

small_df_size = 37

starting_epidemic_year  = 1994
# plt.style.use('ggplot')

reference_year = 2000
adj =  reference_year - starting_epidemic_year

year_end_results = 2100
years = year_end_results - starting_epidemic_year +1
tspan = np.arange(0, years, 1)
discount_rate = 5/100

#### Costing done based on the discounted into the future numbers 
folder_path = "output/Chp3_scenarios"


abbreviated_Strategies = [
    "",  "(1)",  "(2)",  "(3)",  
    "(4)",  "(5)",  "(6)",  "(7)",  
    "(8)",  "(9)",  "(10)",  "(11)",  
    "(12)"
]

### new names based on manuscript
abbreviated_Strategies = [
    "","S1", "S2", "S3", "S4", "S5", "S6",
    "S7", "S8", "S9", "S10", "S11", "S12"
]

shortened_Strategies = [""
,"(1) Current implementation scenario"
,"(2) 40% access to POC VL testing , no DR testing"
,"(3) 60% access to POC VL testing , no DR testing"
,"(4) 20% access to POC VL testing and DR testing \nbefore treatment initiation"
,"(5) 40% access to POC VL testing and DR testing \nbefore treatment initiation"
,"(6) 60% access to POC VL testing and DR testing \nbefore treatment initiation"
,"(7) 20% access to POC VL testing and DR testing \nto confirm treatment failure" 
,"(8) 40% access to POC VL testing and DR testing \nto confirm treatment failure"
,"(9) 60% access to POC VL testing and DR testing \nto confirm treatment failure"
,"(10) 20% access to POC VL testing, with DR testing \nboth before treatment initiation \nand to confirm treatment failure"
,"(11) 40% access to POC VL testing, with DR testing \nboth before treatment initiation \nand to confirm treatment failure"
,"(12) Expansion scenario"
]


shortened_Strategies_prest = [""
,"(1) 20% access to POC VL testing , no DR testing"
,"(2) 40% access to POC VL testing , no DR testing"
,"(3) 60% access to POC VL testing , no DR testing"
,"(4) 20% access to POC VL testing and DR testing \nbefore treatment initiation"
,"(5) 40% access to POC VL testing and DR testing \nbefore treatment initiation"
,"(6) 60% access to POC VL testing and DR testing \nbefore treatment initiation"
,"(7) 20% access to POC VL testing and DR testing \nto confirm treatment failure" 
,"(8) 40% access to POC VL testing and DR testing \nto confirm treatment failure"
,"(9) 60% access to POC VL testing and DR testing \nto confirm treatment failure"
,"(10) 20% access to POC VL testing, with DR testing \nboth before treatment initiation \nand to confirm treatment failure"
,"(11) 40% access to POC VL testing, with DR testing \nboth before treatment initiation \nand to confirm treatment failure"
,"(12) 60% access to POC VL testing, with DR testing \nboth before treatment initiation \nand to confirm treatment failure"
]



markers = ['x', '*', 'H', 'h', 'p', '+', 'd', 'D', 'v', '^', 's', 'o', 'x', '*', 'H', 'h', 'p']
colours = ['white', 'red', 'blue', 'cyan', 'orange', 'green', 'green', 'orange', 'black', 'blue', 'red', 'green', 'purple','lime',  'pink', 'grey', 'brown']



def extract_excel_files(excel_file, small_df_size):
    df = pd.read_excel(excel_file)
    number_of_small_dfs = len(df) // small_df_size
    writing_sols = [df.iloc[i*small_df_size:(i+1)*small_df_size] for i in range(number_of_small_dfs)]
    return writing_sols


def setCorrectxAxis(frame, frequency_ticks =5,starting_position=0):
    plt.figure()
    default_x_ticks = list(range(len(frame.I_W)))
    new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+len(frame.I_W) ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, len(frame.I_W)+1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)




def plot_cost_effectiveness_plane(file_names, p_costs, p_daly, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan):
    median_incremental_costs = [0]
    median_incremental_effects = [0]
    ci_cost_higher =[0]
    ci_cost_lower =[0]
    ci_effect_higher = [0]
    ci_effect_lower = [0]
    comparator = extract_excel_files(os.path.join(base_path, folder_path, file_names[0]), year_end_evaluation - starting_epidemic_year + 1)
    comparator_cost = np.array(probabilistic_all_populations(comparator, p_costs, discount_rate,year_end_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
    comparator_effect = np.array(probabilistic_all_populations(comparator, p_daly, discount_rate,year_end_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
    
    md_comparator_cost = np.median(comparator_cost)
    md_comparator_effect = np.median(comparator_effect)
    
    ### Add in the values of median total cost and effects

    for file_name in file_names[1:]:
        evaluation = extract_excel_files(os.path.join(base_path, folder_path, file_name), year_end_evaluation - starting_epidemic_year + 1)

        evaluation_cost = np.array(probabilistic_all_populations(evaluation, p_costs, discount_rate,year_end_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
        evaluation_effect = np.array(probabilistic_all_populations(evaluation, p_daly, discount_rate,year_end_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
        
        md_evaluation_cost = np.median(evaluation_cost)
        md_evaluation_effect = np.median(evaluation_effect)
        
        median_incremental_costs.append(md_evaluation_cost- md_comparator_cost)
        median_incremental_effects.append(-(md_evaluation_effect - md_comparator_effect))

        # append the upper and lower of ICERs
        ci_cost_lower.append(np.percentile(median_incremental_costs, 2.5))
        ci_cost_higher.append(np.percentile(median_incremental_costs, 97.5))
        ci_effect_lower.append(np.percentile(median_incremental_effects, 2.5))
        ci_effect_higher.append(np.percentile(median_incremental_effects, 97.5))


    return median_incremental_costs,median_incremental_effects, [ci_cost_lower,ci_cost_higher],[ci_effect_lower,ci_effect_higher]


year_end_policy = 2020
year_end_evaluation = 2100
### GDP per capita, PPP (constant 2021 international $) 
WTP= 4_119


##### Strategies of 25 years 
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

median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs, p_daly, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)



#### revising efficiency frontier
undominated_list = ["Routine_2_Baseline_WithDolutegravir_56_56.xlsx"
                    ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_56_56.xlsx"
                    ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_56_56.xlsx"
                    ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_56_56.xlsx"
                    ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_56_56.xlsx"]

efficient_mask = np.array([item in undominated_list for item in file_names])
efficient_costs = np.array(median_incremental_costs)[efficient_mask]
efficient_effects = np.array(median_incremental_effects)[efficient_mask]
# Sort the efficient points again by effectiveness for plotting
sorted_indices = np.argsort(efficient_effects)
efficient_costs = efficient_costs[sorted_indices]
efficient_effects = efficient_effects[sorted_indices]



plt.figure(figsize=(10, 8))
### enter the CI
# plt.errorbar(median_incremental_effects, median_incremental_costs, 
#              xerr=ci_effect, yerr=ci_cost, fmt='none', ecolor='gray',  elinewidth=0.5,     
#              linestyle='--',capsize=2)          

for i, file_name in enumerate(file_names):
    plt.scatter( median_incremental_effects[i],median_incremental_costs[i], label=abbreviated_Strategies[i], s=70, marker=markers[i], color=colours[i])
plt.ticklabel_format(style='plain', axis='y')
x = np.linspace(0, 1/4*max(median_incremental_effects), 10)  
y =  x*WTP

plt.plot(x, y, label='', color="red", linestyle=":")

# plt.scatter(median_incremental_costs, median_incremental_effects, marker='o', color='blue')

# # Adding labels for each evaluation
# for i, file_name in enumerate(file_names):
#     plt.annotate(file_name, (median_incremental_costs[i], median_incremental_effects[i]), fontsize=9)
# Plot the efficient frontier
plt.plot(efficient_effects, efficient_costs, color='green', linestyle='--', marker='o', markersize=1)

plt.gca().get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x / 1_000_000), ',')))
plt.gca().get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x/1_000), ',')))

plt.axhline(0, color='grey', lw=1, ls='--')
plt.axvline(0, color='grey', lw=1, ls='--')
# plt.title('Cost-Effectiveness Plane')
plt.ylabel('Incremental Cost (millions)',fontsize=15)
plt.xlabel('DALYs averted (thousands)',fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=12) 

plt.grid()
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', borderaxespad=0., ncol=6,fontsize=14, title = "Scenario",title_fontsize=12)
# plt.savefig('CEA_results/Figure 3 CEA planes/Figure 2 CEA plane_25horizon.png', dpi=500,bbox_inches='tight')
# plt.savefig(os.path.join(base_path, 'output','Chp3-figures', 'Figure 2 CEA plane_25horizon.png'), dpi=500,bbox_inches='tight')
plt.savefig(os.path.join(base_path, 'output', 'Figure 2 CEA plane_25horizon.png'), dpi=500,bbox_inches='tight')
plt.show()




##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
########## Grouped bar charts diagrams in the model 

#### Costing done based on the discounted into the future numbers 
folder_path = "output/Chp3_scenarios"
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
,"Expanded scenario - 60% access to POC VL testing, with DR testing both before treatment initiation and confirm treatment failure"
]
matching_groups = ['G1'] * 3 + ['G2'] * 3 + ['G3'] * 3 + ['G4'] * 3



def extract_excel_files(excel_file, small_df_size):
    df = pd.read_excel(excel_file)
    number_of_small_dfs = len(df) // small_df_size
    writing_sols = []
    
    for i in range(number_of_small_dfs):
        start_index = starting_epidemic_year
        end_index = start_index + small_df_size
        small_df = df.iloc[i*small_df_size:(i+1)*small_df_size].copy()
        small_df.set_index(pd.RangeIndex(start=start_index, stop=end_index), inplace=True)

        ## only use this when reading data, not when with the main functions, p_incidences is correct 
        small_df['newly_infected'] = np.gradient(small_df['incidence'], tspan)
        small_df['yearly_deaths'] = np.gradient(small_df['deaths'],tspan)


        writing_sols.append(small_df)
    
    return writing_sols


def setCorrectxAxis(frame, frequency_ticks =5,starting_position=0):
    plt.figure()
    default_x_ticks = list(range(len(frame.I_W)))
    new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+len(frame.I_W) ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, len(frame.I_W)+1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)



def report_statistics(input_list, scale='actual'):
    # Compute the median
    median_value = np.median(input_list)
    
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
    median_value = round(median_value, 1)
    lower_bound = round(lower_bound, 1)
    upper_bound = round(upper_bound, 1)
    
    # Return the formatted string with commas as thousand separators
    return f'{median_value:,.1f} ({lower_bound:,.1f} - {upper_bound:,.1f})'


def extract_data_from_dfs(small_dfs, column_or_formula):
    extracted_data = []
    
    for df in small_dfs:
        if isinstance(column_or_formula, str):  # If a column name is provided
            values = df[column_or_formula].values
        else:  # If a formula (function) is provided
            values = df.apply(column_or_formula, axis=1).values
        extracted_data.append(values)
    
    return np.array(extracted_data)



def compute_summary_statistics(data, index, decimals=0, is_percentage=False):
    # Compute the median
    median_val = np.round( np.median(data, axis=0)[index],  decimals)
    
    # Compute the 2.5th and 97.5th percentiles
    P2_5 = np.round( np.percentile(data, 2.5, axis=0)[index],  decimals)
    P97_5 = np.round( np.percentile(data, 97.5, axis=0)[index],  decimals)

    # Format the summary based on the is_percentage flag
    if is_percentage:
        summary = f"{median_val}% ({P2_5}%-{P97_5}%)"
    else:
        summary = f"{int(median_val):,} ({int(P2_5):,}-{int(P97_5):,})"
    
    return summary




def return_all_results_modelling(file_names, p_costs, p_daly,p_incidences,p_deaths,p_ADR,p_TDR,p_secondline,p_mistreated,p_ADR_Total,p_TDR_Total,p_incidences_DR,p_Total_TFailure,p_Viral_unsuppressed,p_new_TFailure, discount_rate, year_start_policy,year_end_evaluation,year_end_policy, starting_epidemic_year, tspan):
    median_total_costs = []
    median_total_effects = []
    report_incremental_costs = [0]
    report_incremental_effects = [0]
    percentage_cost = [0]
    ICERs = [0]

    percentage_effect = [0]
    #### Other outputs to be added in the frame
    incidences = []
    deaths = []
    ADR = []
    TDR =[]
    secondline = []
    mistreated =[]
    ADR_Total =[]
    TDR_Total =[]
    incidences_DR =[]
    Total_TF =[]
    Viral_unsuppressed = []
    new_TFailure = []

    percentage_incidences = [0]
    percentage_deaths = [0]
    percentage_ADR = [0]
    percentage_TDR =[0]
    percentage_secondline = [0]
    percentage_mistreated =[0]
    percentage_ADR_Total =[0]
    percentage_TDR_Total =[0]
    percentage_incidences_DR =[0]
    percentage_Total_TF = [0]
    percentage_Viral_unsuppressed = [0]
    percentage_new_TFailure = [0]

    
    comparator = extract_excel_files(os.path.join(base_path,folder_path, file_names[0]), year_end_evaluation - starting_epidemic_year + 1)
    comparator_cost = np.array(probabilistic_all_populations(comparator, p_costs, discount_rate,year_start_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
    comparator_effect = np.array(probabilistic_all_populations(comparator, p_daly, discount_rate,year_start_policy-starting_epidemic_year,
                                                             year_end_evaluation-starting_epidemic_year,tspan=tspan))
    
    #### additional outcomes to calculate in the modelling 
    #### all additional outcomes should not be discounted , unlike costs and DALYs
    comparator_incidences = np.array(probabilistic_all_populations(comparator, p_incidences, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_deaths = np.array(probabilistic_all_populations(comparator, p_deaths, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_ADR = np.array(probabilistic_all_populations(comparator, p_ADR, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_TDR = np.array(probabilistic_all_populations(comparator, p_TDR, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_secondline = np.array(probabilistic_all_populations(comparator, p_secondline, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_mistreated = np.array(probabilistic_all_populations(comparator, p_mistreated, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_ADR_Total = np.array(probabilistic_all_populations(comparator, p_ADR_Total, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_TDR_Total = np.array(probabilistic_all_populations(comparator, p_TDR_Total, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_incidences_DR = np.array(probabilistic_all_populations(comparator, p_incidences_DR, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_Total_TF = np.array(probabilistic_all_populations(comparator, p_Total_TFailure, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_Viral_unsuppressed= np.array(probabilistic_all_populations(comparator, p_Viral_unsuppressed, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
    comparator_new_TFailure= np.array(probabilistic_all_populations(comparator, p_new_TFailure, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))

    
    ### Add in the values of median total cost and effects
    median_total_costs.append(report_statistics(comparator_cost,"millions"))
    median_total_effects.append(report_statistics(comparator_effect,'thousands'))


    #### Add in additional outcomes in the modell
    incidences.append(report_statistics(comparator_incidences))
    deaths.append(report_statistics(comparator_deaths))
    ADR.append(report_statistics(comparator_ADR))
    TDR.append(report_statistics(comparator_TDR))
    secondline.append(report_statistics(comparator_secondline))
    mistreated.append(report_statistics(comparator_mistreated))
    ADR_Total.append(report_statistics(comparator_ADR_Total))
    TDR_Total.append(report_statistics(comparator_TDR_Total))
    incidences_DR.append(report_statistics(comparator_incidences_DR))
    Total_TF.append(report_statistics(comparator_Total_TF))
    Viral_unsuppressed.append(report_statistics(comparator_Viral_unsuppressed))
    new_TFailure.append(report_statistics(comparator_new_TFailure))


    for file_name in file_names[1:]:
        evaluation = extract_excel_files(os.path.join(base_path,folder_path, file_name), year_end_evaluation - starting_epidemic_year + 1)

        evaluation_cost = np.array(probabilistic_all_populations(evaluation, p_costs, discount_rate,year_start_policy-starting_epidemic_year,year_end_evaluation-starting_epidemic_year,tspan=tspan))
        evaluation_effect = np.array(probabilistic_all_populations(evaluation, p_daly, discount_rate,year_start_policy-starting_epidemic_year,year_end_evaluation-starting_epidemic_year,tspan=tspan))


        #### Additional comparisions outcomes
        #### additional outcomes like these should not be discounted 0%
        evaluation_incidences = np.array(probabilistic_all_populations(evaluation, p_incidences, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_deaths = np.array(probabilistic_all_populations(evaluation, p_deaths, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_ADR = np.array(probabilistic_all_populations(evaluation, p_ADR, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_TDR = np.array(probabilistic_all_populations(evaluation, p_TDR, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_secondline = np.array(probabilistic_all_populations(evaluation, p_secondline, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_mistreated = np.array(probabilistic_all_populations(evaluation, p_mistreated, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_ADR_Total = np.array(probabilistic_all_populations(evaluation, p_ADR_Total, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_TDR_Total = np.array(probabilistic_all_populations(evaluation, p_TDR_Total, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_incidences_DR = np.array(probabilistic_all_populations(evaluation, p_incidences_DR, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_Total_TF = np.array(probabilistic_all_populations(evaluation, p_Total_TFailure, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_Viral_unsuppressed= np.array(probabilistic_all_populations(evaluation, p_Viral_unsuppressed, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))
        evaluation_new_TFailure= np.array(probabilistic_all_populations(evaluation, p_new_TFailure, 0/100,year_start_policy-starting_epidemic_year,year_end_policy-starting_epidemic_year,tspan=tspan))


        ### These are reporting dataframes only
        report_incremental_costs.append(report_statistics(evaluation_cost- comparator_cost,'thousands'))
        report_incremental_effects.append(report_statistics(-(evaluation_effect-comparator_effect)))

        ### Add in the total effects
        median_total_costs.append(report_statistics(evaluation_cost,"millions"))
        median_total_effects.append(report_statistics(evaluation_effect,'thousands'))

        ### Percentage of costs 
        percentage_cost.append(report_statistics((evaluation_cost- comparator_cost)/comparator_cost,"percentage"))
        ### Percentage of DALYs
        percentage_effect.append(report_statistics(-(evaluation_effect-comparator_effect)/comparator_effect,"percentage"))
        ### Calculate the ICER column
        ICERs.append(report_statistics((evaluation_cost- comparator_cost)/-(evaluation_effect-comparator_effect)))

        #### additional outcomes in the modelling 
        incidences.append(report_statistics(evaluation_incidences))
        deaths.append(report_statistics(evaluation_deaths))
        ADR.append(report_statistics(evaluation_ADR))
        TDR.append(report_statistics(evaluation_TDR))
        secondline.append(report_statistics(evaluation_secondline))
        mistreated.append(report_statistics(evaluation_mistreated))
        ADR_Total.append(report_statistics(evaluation_ADR_Total))
        TDR_Total.append(report_statistics(evaluation_TDR_Total))
        incidences_DR.append(report_statistics(evaluation_incidences_DR))
        Total_TF.append(report_statistics(evaluation_Total_TF))
        Viral_unsuppressed.append(report_statistics(evaluation_Viral_unsuppressed))
        new_TFailure.append(report_statistics(evaluation_new_TFailure))


        #### Percentage averted of outcomes in the modelling 
        percentage_incidences.append(report_statistics(-(evaluation_incidences-comparator_incidences)/comparator_incidences,"percentage"))
        percentage_deaths.append(report_statistics(-(evaluation_deaths-comparator_deaths)/comparator_deaths,"percentage"))
        percentage_ADR.append(report_statistics(-(evaluation_ADR-comparator_ADR)/comparator_ADR,"percentage"))
        percentage_TDR.append(report_statistics(-(evaluation_TDR-comparator_TDR)/comparator_TDR,"percentage"))
        ### second-line and mistreated would be reduced
        percentage_secondline.append(report_statistics((evaluation_secondline-comparator_secondline)/comparator_secondline,"percentage"))
        percentage_mistreated.append(report_statistics((evaluation_mistreated-comparator_mistreated)/comparator_mistreated,"percentage"))
        percentage_ADR_Total.append(report_statistics(-(evaluation_ADR_Total-comparator_ADR_Total)/comparator_ADR_Total,"percentage"))
        percentage_TDR_Total.append(report_statistics(-(evaluation_TDR_Total-comparator_TDR_Total)/comparator_TDR_Total,"percentage"))
        percentage_incidences_DR.append(report_statistics(-(evaluation_incidences_DR-comparator_incidences_DR)/comparator_incidences_DR,"percentage"))
        percentage_Total_TF.append(report_statistics(-(evaluation_Total_TF-comparator_Total_TF)/comparator_Total_TF,"percentage"))
        percentage_Viral_unsuppressed.append(report_statistics(-(evaluation_Viral_unsuppressed-comparator_Viral_unsuppressed)/comparator_Viral_unsuppressed,"percentage"))
        percentage_new_TFailure.append(report_statistics(-(evaluation_new_TFailure-comparator_new_TFailure)/comparator_new_TFailure,"percentage"))








    return pd.DataFrame({
        'Strategies': corresponding_Strategies,
        'Total Cost (in millions $)':median_total_costs,
        'Incremental Cost (in thousands $)': report_incremental_costs,
        '% of Total Costs':percentage_cost,
        'Total Effect (in thousands DALYs)': median_total_effects,
        'Incremental Effect': report_incremental_effects,
        '% of Total Effect':percentage_effect,
        'ICERs': ICERs
    }), pd.DataFrame({
        'Strategies': corresponding_Strategies,
        'Total Cost (in millions $)':median_total_costs,
        'Incremental Cost (in thousands $)': report_incremental_costs,
        '% of Total Costs':percentage_cost,
        'Total Effect (in thousands DALYs)': median_total_effects,
        'Incremental Effect': report_incremental_effects,
        '% of Total Effect':percentage_effect,
        'ICERs': ICERs,
        'incidences':incidences,
        '% of Total incidences':percentage_incidences,
        'deaths':deaths,
        '% of Total deaths':percentage_deaths,
        'ADR':ADR,
        '% of New Failures of ADR':percentage_ADR,
        'TDR':TDR,
        '% of New Failures of TDR':percentage_TDR,
        'Second-line ART':secondline,
        '% of Total secondline ART':percentage_secondline,
        'Mistreated DR on first-line':mistreated,
        '% of Total of mistreated DR': percentage_mistreated,
        'Total ADR numbers':ADR_Total,
        '% of Total of ADR': percentage_ADR_Total,
        'Total TDR numbers':TDR_Total,
        '% of Total of TDR': percentage_TDR_Total,
        'Total incidences with DR':incidences_DR,
        '% of incidences with DR averted':percentage_incidences_DR,
        'Total Treatment Failures due to DR':Total_TF,
        '% total TF due to DR averted':percentage_Total_TF,
        'Total viral unsuppressed':Viral_unsuppressed,
        '% total viral unsuppressed averted':percentage_Viral_unsuppressed,
        'New treatment failures in total':new_TFailure,
        '% New treatment failures averted':percentage_new_TFailure,

    })



def parse_value(value):
    match = re.match(r"(\d+(?:\.\d+)?) \((\d+(?:\.\d+)?) - (\d+(?:\.\d+)?)\)", value)
    match = re.match(r"(-?\d+(?:\.\d+)?) \((-?\d+(?:\.\d+)?) - (-?\d+(?:\.\d+)?)\)", value)
    if match:
        return float(match.group(1)), float(match.group(2)), float(match.group(3))
    else:
        raise ValueError(f"Invalid format: {value}")


def plot_grouped_bar_chart(file_path, value,yaxistitle, bar_width=0.2, figsize=(12, 6)):
    df = pd.read_excel(file_path)
    df = df.iloc[1:]
    df[['Median', 'Lower_CI', 'Upper_CI']] = df[value].apply(parse_value).apply(pd.Series)

    groups = ['VL testing only, no DR testing'] * 3 + ['VL testing with \npre-treatment DR testing'] * 3 + ['VL testing with \npost-treatment DR testing'] * 3 + ['VL testing with \nboth pre- and\n post-treament DR testing'] * 3  
    df['Group'] = groups
    opacities = [0.5, 0.7, 0.9]  # Opacity values for 3 bars within each group
    group_colors = {
        'VL testing only, no DR testing': '#1f77b4',  # Blue
        'VL testing with \npre-treatment DR testing': '#ff7f0e',  # Orange
        'VL testing with \npost-treatment DR testing': '#2ca02c',  # Green
        'VL testing with \nboth pre- and\n post-treament DR testing': '#d62728'   # Red
    }

    group_labels = df['Group'].unique()
    items_per_group = len(df) // len(group_labels)

    fig, ax = plt.subplots(figsize=figsize)
    # Add percentage formatter to y-axis
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))
    ax.tick_params(axis='y', which='major', labelsize=14, length=8, width=2)  # Y-axis
    ax.tick_params(axis='x', which='major', labelsize=12, length=8, width=2)  # X-axis

    x = np.repeat(np.arange(len(group_labels)), items_per_group) + np.tile(np.arange(items_per_group), len(group_labels)) * bar_width
    group_starts = np.arange(0, len(df), items_per_group)  

    for i, group in enumerate(group_labels):
        group_data = df[df['Group'] == group]
        x_positions = x[group_starts[i]:group_starts[i] + items_per_group]

        for j in range(items_per_group):
            ax.bar(
                x_positions[j],  # Individual x position
                group_data.iloc[j]['Median'],
                width=bar_width,
                yerr=np.array([
                    [group_data.iloc[j]['Median'] - group_data.iloc[j]['Lower_CI']], 
                    [group_data.iloc[j]['Upper_CI'] - group_data.iloc[j]['Median']]
                ]),
                capsize=5,
                alpha=opacities[j],  # Unique alpha for each bar
                color=group_colors[group],  # Assign color based on group
                label=group if j == 0 else "",
                error_kw={'elinewidth': 0.8}
            )
    ax.set_ylabel(yaxistitle, fontsize=14)
    ax.set_xticks(np.arange(len(group_labels)) + (items_per_group - 1) * bar_width / 2)
    ax.set_xticklabels(group_labels)
    # ax.legend(title='Group')

    plt.tight_layout()
    # plt.savefig('output/Chp3-figures/' +value+ '.png', dpi=500)
    plt.savefig(os.path.join(base_path,'output/Chp3-figures/' +value+ '.png'), dpi=500)
    plt.show()




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
display_frame,barchart_frame = return_all_results_modelling(file_names,p_costs, p_daly,p_incidences,p_deaths,p_ADR,p_TDR,p_secondline,p_mistreated,p_ADR_Total,p_TDR_Total,p_incidences_DR,p_Total_TFailure,p_Viral_unsuppressed,p_new_TFailure, discount_rate, year_start_policy,year_end_evaluation, year_end_policy,
                                              starting_epidemic_year, tspan)


#### add percentage signs into two columns of data 
display_frame['% of Total Costs'] = display_frame['% of Total Costs'].str.replace(r'(\d+\.\d+)', r'\1%', regex=True)
display_frame['% of Total Effect'] = display_frame['% of Total Effect'].str.replace(r'(\d+\.\d+)', r'\1%', regex=True)
output = display_frame
### Output to Excel for Table 2 results 5 years
# output.to_excel("CEA_results/Output_Results/Final_CEA_Table2.xlsx")
output.to_excel(os.path.join(base_path, 'output','Chp3-figures', 'Final_CEA_Table2.xlsx'))
# barchart_frame.to_excel("CEA_results/Output_Results/Final_BarchartCEA results.xlsx")
barchart_frame.to_excel(os.path.join(base_path, 'output','Chp3-figures', 'Final_BarchartCEA results.xlsx'))


### graph_title instead
file_outout = os.path.join(base_path, 'output','Chp3-figures', 'Final_BarchartCEA results.xlsx')

plot_grouped_bar_chart(file_outout,"% New treatment failures averted","% Reduction in new people \nexperiencing treatment failures 2020-2050")
plot_grouped_bar_chart(file_outout,"% of Total incidences","% HIV incidences averted 2020-2050 \n")
plot_grouped_bar_chart(file_outout,"% of Total Costs","% Incremental Life-time Costs\n")


#### Produce results of new infections in 2035 and 2050 across the selected scenarios for modelling
extra_files = ["Routine_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_3_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_4_Baseline_WithDolutegravir_56_56.xlsx"
                               
                               ,"POC_2_Timevarying_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_56_56.xlsx"
                               ,"POC_AcquiredDR_2_Timevarying_WithDolutegravir_56_56.xlsx"]

labels = ["Baseline",
          "POC Current -No DR 20%",
          "POC VL 40%",
          "POC VL 60%",
          "POC with 20% pre-treatment",
          "POC with 20% post-treatment",
          "POC with 20%  both pre- and post-treatment"]

for file, label in zip(extra_files,labels):
    print(f"Current evaluation results at {label}.")
    comparator = extract_excel_files(os.path.join(base_path, folder_path, file), year_end_evaluation - starting_epidemic_year + 1)
    new_infections = extract_data_from_dfs(comparator,"newly_infected")
    print(f"Infections in 2035 is {compute_summary_statistics(new_infections, 2035-starting_epidemic_year)}.")
    print(f"Infections in 2050 is {compute_summary_statistics(new_infections, 2050-starting_epidemic_year)}.")





#### Combine the figure for CEA 

listB = [
    '% New treatment failures averted.png',
    '% of Total incidences.png',
    '% of Total Costs.png'
]

img_list = [mpimg.imread(os.path.join(base_path,'output/Chp3-figures', filename)) for filename in listB]

# Create a 1x3 grid for displaying the images:
fig, axes = plt.subplots(3, 1, figsize=(8.5, 12)) 

for ax, img, title in zip(axes, img_list, listB):
    ax.imshow(img)
    ax.axis('off')

labels = ['(a)', '(b)', '(c)']
for ax, img, label in zip(axes, img_list, labels):
    ax.imshow(img)
    ax.axis('off')
    ax.text(0, 1, label, transform=ax.transAxes, fontsize=12, ha='center', va='top')
# Add horizontal lines to separate the 2 rows
for ax in axes:
    ax.axhline(y=ax.get_ylim()[0], color='black', linewidth=2, linestyle='--')  


# for ax in axes[:-1]:  
    # ax.axvline(x=ax.get_xlim()[1], color='black', linewidth=1, linestyle='--')
plt.tight_layout(pad=0.1, h_pad=-0.1)

# plt.savefig('CEA_results/Barcharts Figure 3/Combined Figure 3_Final_Manuscript.png', dpi=500)  # Save the figure with high dpi
plt.savefig(os.path.join(base_path, 'output','Chp3-figures', 'Combined Figure 3_Final_Manuscript.png'), dpi=500,bbox_inches='tight')
plt.show()











##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################


def extract_excel_files(excel_file, small_df_size):
    df = pd.read_excel(excel_file)
    number_of_small_dfs = len(df) // small_df_size
    writing_sols = [df.iloc[i*small_df_size:(i+1)*small_df_size] for i in range(number_of_small_dfs)]
    return writing_sols


def setCorrectxAxis_otherposition_maindiagram(frame, frequency_ticks =5,starting_position=0,starting_year = starting_epidemic_year):
    plt.figure()
    default_x_ticks = list(range(len(frame.I_W)) )
    new_x_ticks = list(range(starting_year,starting_year+len(frame.I_W) ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, len(frame.I_W) - starting_position +1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)



#### Costing based on extending the 2100

year_end_results = 2100

years = year_end_results - starting_epidemic_year +1
tspan = np.arange(0, years, 1)
folder_path = "output/Chp3_scenarios"




file_names = [
    "Routine_2_Baseline_WithDolutegravir_56_56.xlsx",
    "POC_2_Baseline_WithDolutegravir_56_56.xlsx",
    "POC_3_Baseline_WithDolutegravir_56_56.xlsx",
    "POC_4_Baseline_WithDolutegravir_56_56.xlsx",
    "POC_2_Timevarying_WithDolutegravir_56_56.xlsx",
    "POC_3_Timevarying_WithDolutegravir_56_56.xlsx",
    "POC_4_Timevarying_WithDolutegravir_56_56.xlsx",
    "POC_AcquiredDR_2_Baseline_WithDolutegravir_56_56.xlsx",
    "POC_AcquiredDR_3_Baseline_WithDolutegravir_56_56.xlsx",
    "POC_AcquiredDR_4_Baseline_WithDolutegravir_56_56.xlsx",
    "POC_AcquiredDR_2_Timevarying_WithDolutegravir_56_56.xlsx",
    "POC_AcquiredDR_3_Timevarying_WithDolutegravir_56_56.xlsx",
    "POC_AcquiredDR_4_Timevarying_WithDolutegravir_56_56.xlsx"
]
results = [
    extract_excel_files(os.path.join(base_path, folder_path, file), year_end_results - starting_epidemic_year + 1)
    for file in file_names
]



# Calculate treatment values for each solution "On treatment, Suppressed"
treatment_results = [
    np.array([sol["T"] for sol in result])
    for result in results
]

(comparator, intervention1, intervention2, intervention3, intervention4, 
 intervention5, intervention6, intervention7, intervention8, 
 intervention9, intervention10, intervention11, intervention12) = results

# Unpack the treatment results into individual variables
(comparator_treatment, intervention1_treatment, intervention2_treatment, intervention3_treatment, 
 intervention4_treatment, intervention5_treatment, intervention6_treatment, 
 intervention7_treatment, intervention8_treatment, intervention9_treatment, 
 intervention10_treatment, intervention11_treatment, intervention12_treatment) = treatment_results


# Define groups and assign interventions to each group
groups = [
    'VL testing only, no DR testing',               # Interventions 1, 2, 3
    'VL testing with pre-treatment DR testing',        # Interventions 4, 5, 6
    'VL testing with post-treatment DR testing',  # Interventions 7, 8, 9
    'VL testing with both pre- and post-treatment DR testing'    # Interventions 10, 11, 12
]

# Line styles to cycle through for each group
linestyles = ['--', ':', '-']



setCorrectxAxis_otherposition_maindiagram(comparator[0],frequency_ticks=5, 
        starting_position=1, starting_year=starting_epidemic_year)

plt.plot(np.median(intervention1_treatment - comparator_treatment, axis =0),color='#1f77b4', linestyle=linestyles[0], linewidth=1)
plt.plot(np.median(intervention2_treatment - comparator_treatment, axis =0),color='#1f77b4', linestyle=linestyles[1], linewidth=1)
plt.plot(np.median(intervention3_treatment - comparator_treatment, axis =0),color='#1f77b4', linestyle=linestyles[2], linewidth=1)

plt.title(groups[0])
plt.axhline(y=0, color='black', linestyle='-.', linewidth=0.5)

plt.ylabel('Additional people with suppressed VL')
plt.xlim(2020 - starting_epidemic_year, 2050 - starting_epidemic_year)
plt.ylim(top = 4_500, bottom = -12_000)
plt.xticks(fontsize = 12)
# plt.savefig('CEA_results/Actual Linecharts Figure 3/First_group', dpi=500, bbox_inches='tight')
plt.savefig(os.path.join(base_path, 'output','Chp3-figures', 'First_group.png'), dpi=500,bbox_inches='tight')
plt.show()




setCorrectxAxis_otherposition_maindiagram(comparator[0],frequency_ticks=5, 
        starting_position=1, starting_year=starting_epidemic_year)

plt.plot(np.median(intervention4_treatment - comparator_treatment, axis =0),color='#ff7f0e', linestyle=linestyles[0], linewidth=1)
plt.plot(np.median(intervention5_treatment - comparator_treatment, axis =0),color='#ff7f0e', linestyle=linestyles[1], linewidth=1)
plt.plot(np.median(intervention6_treatment - comparator_treatment, axis =0),color='#ff7f0e', linestyle=linestyles[2], linewidth=1)

plt.title(groups[1])
plt.axhline(y=0, color='black', linestyle='-.', linewidth=0.5)

# plt.ylabel('Additional people with suppressed VL')
plt.xlim(2020 - starting_epidemic_year, 2050 - starting_epidemic_year)
plt.ylim(top = 4_500, bottom = -12_000)
plt.xticks(fontsize = 12)
# plt.savefig('CEA_results/Actual Linecharts Figure 3/Second_group', dpi=500, bbox_inches='tight')
plt.savefig(os.path.join(base_path, 'output','Chp3-figures', 'Second_group.png'), dpi=500,bbox_inches='tight')
plt.show()






setCorrectxAxis_otherposition_maindiagram(comparator[0],frequency_ticks=5, 
        starting_position=1, starting_year=starting_epidemic_year)

plt.plot(np.median(intervention7_treatment - comparator_treatment, axis =0),color='#2ca02c', linestyle=linestyles[0], linewidth=1)
plt.plot(np.median(intervention8_treatment - comparator_treatment, axis =0),color='#2ca02c', linestyle=linestyles[1], linewidth=1)
plt.plot(np.median(intervention9_treatment - comparator_treatment, axis =0),color='#2ca02c', linestyle=linestyles[2], linewidth=1)

plt.title(groups[2])
plt.axhline(y=0, color='black', linestyle='-.', linewidth=0.5)

plt.ylabel('Additional people with suppressed VL')
plt.xlim(2020 - starting_epidemic_year, 2050 - starting_epidemic_year)
plt.ylim(top = 4_500, bottom = -12_000)
plt.xticks(fontsize = 12)
# plt.savefig('CEA_results/Actual Linecharts Figure 3/Third_group', dpi=500, bbox_inches='tight')
plt.savefig(os.path.join(base_path, 'output','Chp3-figures', 'Third_group.png'), dpi=500,bbox_inches='tight')
plt.show()






setCorrectxAxis_otherposition_maindiagram(comparator[0],frequency_ticks=5, 
        starting_position=1, starting_year=starting_epidemic_year)

plt.plot(np.median(intervention10_treatment - comparator_treatment, axis =0),color='#d62728', linestyle=linestyles[0], linewidth=1, label = "20% access")
plt.plot(np.median(intervention11_treatment - comparator_treatment, axis =0),color='#d62728', linestyle=linestyles[1], linewidth=1, label = "40% access")
plt.plot(np.median(intervention12_treatment - comparator_treatment, axis =0),color='#d62728', linestyle=linestyles[2], linewidth=1, label = "60% access")

plt.title(groups[3])
plt.axhline(y=0, color='black', linestyle='-.', linewidth=0.5)

# plt.ylabel('Additional people with suppressed VL')
plt.xlim(2020 - starting_epidemic_year, 2050 - starting_epidemic_year)
plt.ylim(top = 4_500, bottom = -12_000)
plt.xticks(fontsize = 12)
plt.legend(loc = "lower left",fontsize = 14)
# plt.savefig('CEA_results/Actual Linecharts Figure 3/Fourth_group', dpi=500, bbox_inches='tight')
plt.savefig(os.path.join(base_path, 'output','Chp3-figures', 'Fourth_group.png'), dpi=500,bbox_inches='tight')
plt.show()






listB = [

        "First_group.png",
        "Second_group.png",

        "Third_group.png",
        "Fourth_group.png"
        ]

img_list = [mpimg.imread(os.path.join(base_path,'output/Chp3-figures',filename)) for filename in listB]

# Create a 4x3 grid for displaying the images:
fig, axes = plt.subplots(2, 2, figsize=(15,13))

for ax, img, title in zip(axes.ravel(), img_list, listB):
    ax.imshow(img)
    ax.axis('off')

plt.tight_layout(h_pad = -8)  # Added vertical padding between subplots
# plt.savefig('CEA_results/Actual Linecharts Figure 3/Figure 2_Impact of time.png', dpi=500,bbox_inches='tight')  # Save the figure with high dpi
plt.savefig(os.path.join(base_path, 'output','Chp3-figures', 'Figure 2_Impact of time.png'), dpi=500,bbox_inches='tight')
plt.show()














end_time = datetime.now()
print(f"Code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")

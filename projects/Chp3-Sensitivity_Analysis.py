### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# from Bayesian_Analysis.meta_parameters import *

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.image as mpimg
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import os


# from Bayesian_Analysis.Posteriors_model import *



import random

from matplotlib.ticker import MaxNLocator

from EconEval.cea_functions import *
from EconEval.prb_data_matrix import*
import copy
from datetime import datetime


start_time = datetime.now()
print(f"Starting  code execution at {start_time}...")


random.seed(334)
random_seed = 334

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

corresponding_Strategies = ["Baseline scenario"
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


### new names based on manuscript
abbreviated_Strategies = [
    "","S1", "S2", "S3", "S4", "S5", "S6",
    "S7", "S8", "S9", "S10", "S11", "S12"
]

shortened_Strategies = [""
,"Current implementation scenario"
,"40% access to POC VL testing , no DR testing"
,"60% access to POC VL testing , no DR testing"
,"20% access to POC VL testing and DR testing \nbefore treatment initiation"
,"40% access to POC VL testing and DR testing \nbefore treatment initiation"
,"60% access to POC VL testing and DR testing \nbefore treatment initiation"
,"20% access to POC VL testing and DR testing \nto confirm treatment failure" 
,"40% access to POC VL testing and DR testing \nto confirm treatment failure"
,"60% access to POC VL testing and DR testing \nto confirm treatment failure"
,"20% access to POC VL testing, with DR testing \nboth before treatment initiation \nand to confirm treatment failure"
,"40% access to POC VL testing, with DR testing \nboth before treatment initiation \nand to confirm treatment failure"
,"Expansion scenario"
]

naming_strategies = [
    "Status-quo scenario"
    ,"20% access to POC VL testing and DR testing to confirm treatment failure" 
    ,"40% access to POC VL testing and DR testing to confirm treatment failure"
    ,"60% access to POC VL testing and DR testing to confirm treatment failure"
    ,"Expanded scenario"
]




markers = ['x', '*', 'H', 'h', 'p', '+', 'd', 'D', 'v', '^', 's', 'o', 'x', '*', 'H', 'h', 'p']
colours = ['white', 'red', 'blue', 'cyan', 'orange', 'green', 'green', 'orange', 'black', 'blue', 'red', 'green', 'purple','lime',  'pink', 'grey', 'brown']

# from Bayesian_Analysis.read_posterior_samples import*


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

def report_statistics(input_list, scale='actual', fixed_value=np.nan, deci =1):
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


#### calculate pairwise ICERs
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
        prev_file = extract_excel_files(os.path.join(base_path, folder_path, file_names[i - 1]), year_end_evaluation - starting_epidemic_year + 1)
        curr_file = extract_excel_files(os.path.join(base_path, folder_path, file_names[i]), year_end_evaluation - starting_epidemic_year + 1)

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
        median_costs.extend([report_statistics(prev_cost,scale = "millions"), 
                             report_statistics(curr_cost,scale = "millions"),np.nan])
        median_effects.extend([report_statistics(prev_effect,scale = "thousands"), 
                               report_statistics(curr_effect,scale = "thousands"),np.nan])
        incremental_costs.extend(["-", 
                                  report_statistics(incremental_cost,scale = "millions",fixed_value=median_cost_curr - median_cost_prev),np.nan])
        incremental_effects.extend(["-", 
                                    report_statistics(incremental_effect,fixed_value= -(median_effect_curr - median_effect_prev)),np.nan])
        icers.extend(["-", 
                      report_statistics(icer,scale = "actual", fixed_value= (median_cost_curr - median_cost_prev)/-(median_effect_curr - median_effect_prev),deci=0),np.nan])

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






year_end_policy = 2020
year_end_evaluation = 2100
### GDP per capita, PPP (constant 2021 international $) 
WTP= 4_119




###########################################################################################
###########################################################################################
###########################################################################################


file_names = ["Routine_2_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_2_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_3_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_4_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_2_Timevarying_WithDolutegravir_36_36.xlsx"
                               ,"POC_3_Timevarying_WithDolutegravir_36_36.xlsx"
                               ,"POC_4_Timevarying_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_2_Timevarying_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_3_Timevarying_WithDolutegravir_36_36.xlsx"
                               ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_36_36.xlsx"]



median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs, p_daly, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)


#### revised frontier
undominated_list = ["Routine_2_Baseline_WithDolutegravir_36_36.xlsx"
                    ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_36_36.xlsx"
                    # ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_36_36.xlsx"
                    ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_36_36.xlsx"
                    ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_36_36.xlsx"]

#### CEA plane in this time, I need to remove the last strategy as the all combined strategies were not on the cost-effectiveness plane
efficient_mask = np.array([item in undominated_list  for item in file_names])
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
x = np.linspace(0, 1/2*max(median_incremental_effects), 10)
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
plt.savefig(os.path.join(base_path,'output/Chp3-figures/Chp3-SensitivityAnalysis/Figure 2 CEA plane.png'), dpi=500,bbox_inches='tight')
plt.show()


dframe_result2030 = calculate_pairwise_icer(undominated_list, folder_path, p_costs, p_daly, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan
# ,naming_strategies
,[naming_strategies[i] for i in [0, 1, 3, 4]])









file_names = ["Routine_2_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_2_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_3_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_4_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_2_Timevarying_WithDolutegravir_41_41.xlsx"
                               ,"POC_3_Timevarying_WithDolutegravir_41_41.xlsx"
                               ,"POC_4_Timevarying_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_2_Timevarying_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_3_Timevarying_WithDolutegravir_41_41.xlsx"
                               ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_41_41.xlsx"]

median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs, p_daly, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)


### revising efficiency frontier in the model
undominated_list = ["Routine_2_Baseline_WithDolutegravir_41_41.xlsx"
                    ,"POC_AcquiredDR_2_Baseline_WithDolutegravir_41_41.xlsx"
                    # ,"POC_AcquiredDR_3_Baseline_WithDolutegravir_41_41.xlsx"
                    ,"POC_AcquiredDR_4_Baseline_WithDolutegravir_41_41.xlsx"
                    ,"POC_AcquiredDR_4_Timevarying_WithDolutegravir_41_41.xlsx"]

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
x = np.linspace(0, 1/2*max(median_incremental_effects), 10)  
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
plt.savefig(os.path.join(base_path,'output/Chp3-figures/Chp3-SensitivityAnalysis/Figure 2 CEA plane_10horizon.png'), dpi=500,bbox_inches='tight')
plt.show()




dframe_result2035 = calculate_pairwise_icer(undominated_list, folder_path, p_costs, p_daly, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan
# ,naming_strategies
,[naming_strategies[i] for i in [0, 1, 3, 4]])






###########################################################################################
###########################################################################################
###########################################################################################


##### First sensitivity analysis scenario 
### revision of disability weights: disability weight of unsuppressed on treatment is equal to untreated 
p_daly_specified1 = {
    'I_W': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early, False),        
    'I_DR': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early, False),      
    'D_W': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early, False),
    'D_DR': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,False),
    #### main analyses all suppressed
    'T_W1': (p_Suppressed,False),
    'T_W2': (p_Suppressed,False),
    'T_DR1': (p_Suppressed,False),
    'T_DR2': (p_Suppressed,False), 
    'F_W1': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,False),
    'F_W2': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,False), 
    'F_TDR1': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,False),
    'F_ADR1': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,False),
    'F_DR2': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,False),
    #### deaths will be have to be accumulated, in the DALYs calculations
    'all_mortality': (u_mortality,False)
    # 'deaths': (u_mortality,True)

}
p_daly1 = generate_p_values(p_daly_specified1,n)

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

### revise the disability weights for sensitivity analysis scenario 1
median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs, p_daly1, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)



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








##### CEA plane where the legends is on the left side 
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

plt.plot(x, y, label='', color="red", linestyle="--")

# plt.scatter(median_incremental_costs, median_incremental_effects, marker='o', color='blue')

# # Adding labels for each evaluation
# for i, file_name in enumerate(file_names):
#     plt.annotate(file_name, (median_incremental_costs[i], median_incremental_effects[i]), fontsize=9)
# Plot the efficient frontier
plt.plot(efficient_effects, efficient_costs, color='green', linestyle='--', marker='o', markersize=1, label='Cost-effectiveness frontier')

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

plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', borderaxespad=0., ncol=6,fontsize=14, title = "Strategy",title_fontsize=12)
plt.savefig(os.path.join(base_path,'output/Chp3-figures/Chp3-SensitivityAnalysis/Figure 2 CEA plane_25horizon_Sens1st.png'), dpi=500,bbox_inches='tight')
plt.show()







dframe_resultsens1 = calculate_pairwise_icer(undominated_list, folder_path, p_costs, p_daly1, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies)














##### Second sensitivity analysis scenario 
### revision of proportion of Probability of prob_AIDS at 10% instead


prob_AIDS2 = beta_dist.rvs(2.0, 10.0, size=n, random_state=random_seed)
p_daly_specified2 = {
    'I_W': (prob_AIDS2 * p_Untreated_AIDS + (1-prob_AIDS2) * p_Untreated_Early, False),        
    'I_DR': (prob_AIDS2 * p_Untreated_AIDS + (1-prob_AIDS2) * p_Untreated_Early, False),      
    'D_W': (prob_AIDS2 * p_Untreated_AIDS + (1-prob_AIDS2) * p_Untreated_Early, False),
    'D_DR': (prob_AIDS2 * p_Untreated_AIDS + (1-prob_AIDS2) * p_Untreated_Early,False),
    #### main analyses all suppressed
    'T_W1': (p_Suppressed,False),
    'T_W2': (p_Suppressed,False),
    'T_DR1': (p_Suppressed,False),
    'T_DR2': (p_Suppressed,False), 
    'F_W1': (p_Suppressed,False),
    'F_W2': (p_Suppressed,False), 
    'F_TDR1': (p_Suppressed,False),
    'F_ADR1': (p_Suppressed,False),
    'F_DR2': (p_Suppressed,False),
    #### deaths will be have to be accumulated, in the DALYs calculations
    'all_mortality': (u_mortality,False)
    # 'deaths': (u_mortality,True)

}
p_daly2 = generate_p_values(p_daly_specified2,n)


### revise the disability weights for sensitivity analysis scenario 1
median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs, p_daly2, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)



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








##### CEA plane where the legends is on the left side 
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

plt.plot(x, y, label='', color="red", linestyle="--")

# plt.scatter(median_incremental_costs, median_incremental_effects, marker='o', color='blue')

# # Adding labels for each evaluation
# for i, file_name in enumerate(file_names):
#     plt.annotate(file_name, (median_incremental_costs[i], median_incremental_effects[i]), fontsize=9)
# Plot the efficient frontier
plt.plot(efficient_effects, efficient_costs, color='green', linestyle='--', marker='o', markersize=1, label='Cost-effectiveness frontier')

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

plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', borderaxespad=0., ncol=6,fontsize=14, title = "Strategy",title_fontsize=12)
plt.savefig(os.path.join(base_path,'output/Chp3-figures/Chp3-SensitivityAnalysis/Figure 2 CEA plane_25horizon_Sens2nd.png'), dpi=500,bbox_inches='tight')
plt.show()






dframe_resultsens2 = calculate_pairwise_icer(undominated_list, folder_path, p_costs, p_daly2, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan
# ,[naming_strategies[0]] + naming_strategies[-2:]
,naming_strategies)













##### Third sensitivity analysis scenario 
### revision of proportion of Probability of prob_AIDS at 72% instead


prob_AID3 = beta_dist.rvs(5.0, 2.556, size=n, random_state=random_seed)
p_daly_specified3 = {
    'I_W': (prob_AID3 * p_Untreated_AIDS + (1-prob_AID3) * p_Untreated_Early, False),        
    'I_DR': (prob_AID3 * p_Untreated_AIDS + (1-prob_AID3) * p_Untreated_Early, False),      
    'D_W': (prob_AID3 * p_Untreated_AIDS + (1-prob_AID3) * p_Untreated_Early, False),
    'D_DR': (prob_AID3 * p_Untreated_AIDS + (1-prob_AID3) * p_Untreated_Early,False),
    #### main analyses all suppressed
    'T_W1': (p_Suppressed,False),
    'T_W2': (p_Suppressed,False),
    'T_DR1': (p_Suppressed,False),
    'T_DR2': (p_Suppressed,False), 
    'F_W1': (p_Suppressed,False),
    'F_W2': (p_Suppressed,False), 
    'F_TDR1': (p_Suppressed,False),
    'F_ADR1': (p_Suppressed,False),
    'F_DR2': (p_Suppressed,False),
    #### deaths will be have to be accumulated, in the DALYs calculations
    'all_mortality': (u_mortality,False)
    # 'deaths': (u_mortality,True)

}
p_daly3 = generate_p_values(p_daly_specified3,n)


### revise the disability weights for sensitivity analysis scenario 1
median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs, p_daly3, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)



#### revising efficiency frontier
undominated_list = ["Routine_2_Baseline_WithDolutegravir_56_56.xlsx"
                    # ,"POC_2_Baseline_WithDolutegravir_56_56.xlsx"
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








##### CEA plane where the legends is on the left side 
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

plt.plot(x, y, label='', color="red", linestyle="--")

# plt.scatter(median_incremental_costs, median_incremental_effects, marker='o', color='blue')

# # Adding labels for each evaluation
# for i, file_name in enumerate(file_names):
#     plt.annotate(file_name, (median_incremental_costs[i], median_incremental_effects[i]), fontsize=9)
# Plot the efficient frontier
plt.plot(efficient_effects, efficient_costs, color='green', linestyle='--', marker='o', markersize=1, label='Cost-effectiveness frontier')

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

plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', borderaxespad=0., ncol=6,fontsize=14, title = "Strategy",title_fontsize=12)
plt.savefig(os.path.join(base_path,'output/Chp3-figures/Chp3-SensitivityAnalysis/Figure 2 CEA plane_25horizon_Sens3rd.png'), dpi=500,bbox_inches='tight')
plt.show()






dframe_resultsens3 = calculate_pairwise_icer(undominated_list, folder_path, p_costs, p_daly3, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies)

















##### Fourth sensitivity analysis scenario 
### revision of proportion of Probability on Old ART across all treatment failure compartments - regardless of new or old 


p_costs_specified4 = {
    'I_W': (0, False),        
    'I_DR': (0, False),      
    'D_W': (0, False),
    'D_DR': (0,False),
    'T_W1': (p_ART + p_ClinicVisit,False),
    'T_W2': (p_ART + p_ClinicVisit,False),
    'T_DR1': (p_ART + p_ClinicVisit,False),
    'T_DR2': (p_ART_2ndline + p_ClinicVisit,False), 
    'F_W1': ((p_ART + p_ClinicVisit)*p_probOldART,False),
    'F_W2': ((p_ART + p_ClinicVisit)*p_probOldART,False), 
    'F_TDR1': ((p_ART + p_ClinicVisit)*p_probOldART,False),
    'F_ADR1': ((p_ART + p_ClinicVisit)*p_probOldART,False),
    'F_DR2': ((p_ART_2ndline + p_ClinicVisit)*p_probOldART,False),
    'VLtest_POC': (p_VLtest + p_VLtest_labour + p_AdherenceCounsel + p_ClinicVisit,True),
    'preDRtest':(p_DRTesting + p_DRTesting_labour + p_ARTinitiation + 2*(p_VLtest + p_VLtest_labour + p_AdherenceCounsel  + p_ClinicVisit),True),
    'postDRtest_POC':(p_DRTesting + p_DRTesting_labour + p_ARTinitiation + 2*(p_VLtest + p_VLtest_labour + p_AdherenceCounsel  + p_ClinicVisit),True)
    
    ,'deaths': (0,True)
    ,'diagnoses':(p_cDiagnosis,True)
    ,'treat_inits':(p_ARTinitiation,True)
}
p_costs4 = generate_p_values(p_costs_specified4,n)


### revise the disability weights for sensitivity analysis scenario 1
median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs4, p_daly, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)



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








##### CEA plane where the legends is on the left side 
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

plt.plot(x, y, label='', color="red", linestyle="--")

# plt.scatter(median_incremental_costs, median_incremental_effects, marker='o', color='blue')

# # Adding labels for each evaluation
# for i, file_name in enumerate(file_names):
#     plt.annotate(file_name, (median_incremental_costs[i], median_incremental_effects[i]), fontsize=9)
# Plot the efficient frontier
plt.plot(efficient_effects, efficient_costs, color='green', linestyle='--', marker='o', markersize=1, label='Cost-effectiveness frontier')

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

plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', borderaxespad=0., ncol=6,fontsize=14, title = "Strategy",title_fontsize=12)
plt.savefig(os.path.join(base_path,'output/Chp3-figures/Chp3-SensitivityAnalysis/Figure 2 CEA plane_25horizon_Sens4th.png'), dpi=500,bbox_inches='tight')
plt.show()




# dframe_resultsens4 = calculate_pairwise_icer(undominated_list, folder_path, p_costs4, p_daly, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies)


dframe_resultsens4 = calculate_pairwise_icer(undominated_list, folder_path, p_costs4, p_daly, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan
# ,[naming_strategies[0]] + naming_strategies[-3:]
,naming_strategies)




















##### Fifth sensitivity analysis scenario 
### revision of proportion of Probability on Old ART across all treatment failure compartments - regardless of new or old 


p_costs_specified5 = {
    'I_W': (0, False),        
    'I_DR': (0, False),      
    'D_W': (0, False),
    'D_DR': (0,False),
    'T_W1': (p_ART + p_ClinicVisit,False),
    'T_W2': (p_ART + p_ClinicVisit,False),
    'T_DR1': (p_ART + p_ClinicVisit,False),
    'T_DR2': (p_ART_2ndline + p_ClinicVisit,False), 
    'F_W1': ((p_ART + p_ClinicVisit)*p_probNewART,False),
    'F_W2': ((p_ART + p_ClinicVisit)*p_probNewART,False), 
    'F_TDR1': ((p_ART + p_ClinicVisit)*p_probNewART,False),
    'F_ADR1': ((p_ART + p_ClinicVisit)*p_probNewART,False),
    'F_DR2': ((p_ART_2ndline + p_ClinicVisit)*p_probNewART,False),
    'VLtest_POC': (p_VLtest + p_VLtest_labour + p_AdherenceCounsel + p_ClinicVisit,True),
    'preDRtest':(p_DRTesting + p_DRTesting_labour + p_ARTinitiation + 2*(p_VLtest + p_VLtest_labour + p_AdherenceCounsel  + p_ClinicVisit),True),
    'postDRtest_POC':(p_DRTesting + p_DRTesting_labour + p_ARTinitiation + 2*(p_VLtest + p_VLtest_labour + p_AdherenceCounsel  + p_ClinicVisit),True)
    
    ,'deaths': (0,True)
    ,'diagnoses':(p_cDiagnosis,True)
    ,'treat_inits':(p_ARTinitiation,True)
}
p_costs5 = generate_p_values(p_costs_specified5,n)


### revise the disability weights for sensitivity analysis scenario 1
median_incremental_costs,median_incremental_effects,ci_cost,ci_effect = plot_cost_effectiveness_plane(file_names,
                              p_costs5, p_daly, discount_rate, year_end_policy,year_end_evaluation, starting_epidemic_year, tspan)



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








##### CEA plane where the legends is on the left side 
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

plt.plot(x, y, label='', color="red", linestyle="--")

# plt.scatter(median_incremental_costs, median_incremental_effects, marker='o', color='blue')

# # Adding labels for each evaluation
# for i, file_name in enumerate(file_names):
#     plt.annotate(file_name, (median_incremental_costs[i], median_incremental_effects[i]), fontsize=9)
# Plot the efficient frontier
plt.plot(efficient_effects, efficient_costs, color='green', linestyle='--', marker='o', markersize=1, label='Cost-effectiveness frontier')

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

plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', borderaxespad=0., ncol=6,fontsize=14, title = "Strategy",title_fontsize=12)
plt.savefig(os.path.join(base_path,'output/Chp3-figures/Chp3-SensitivityAnalysis/Figure 2 CEA plane_25horizon_Sens5th.png'), dpi=500,bbox_inches='tight')
plt.show()




dframe_resultsens5 = calculate_pairwise_icer(undominated_list, folder_path, p_costs5, p_daly, discount_rate, year_end_policy, year_end_evaluation, starting_epidemic_year, tspan,naming_strategies)























output_path = os.path.join(base_path,"output/Chp3-figures/Chp3-SensitivityAnalysis/CEA frontier_sens.xlsx")
# dframe_result25.to_excel(output_path, index=False)
with pd.ExcelWriter(output_path) as writer:
    dframe_resultsens1.to_excel(writer, sheet_name='SensSc1', index=False)
    dframe_resultsens2.to_excel(writer, sheet_name='SensSc2', index=False)
    dframe_resultsens3.to_excel(writer, sheet_name='SensSc3', index=False)
    dframe_resultsens4.to_excel(writer, sheet_name='SensSc4', index=False)
    dframe_resultsens5.to_excel(writer, sheet_name='SensSc5', index=False)
    dframe_result2030.to_excel(writer, sheet_name='2030results', index=False)
    dframe_result2035.to_excel(writer, sheet_name='2035results', index=False)










end_time = datetime.now()
print(f"Code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")

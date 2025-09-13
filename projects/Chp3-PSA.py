### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
from scipy import stats 
import os



import random

from matplotlib.ticker import MaxNLocator

from EconEval.cea_functions import *
from EconEval.prb_data_matrix import*
import matplotlib.patches as patches
from datetime import datetime


start_time = datetime.now()
print(f"Starting  code execution at {start_time}...")


random.seed(334)

small_df_size = 37
folder_path = "output/Chp3_scenarios"

starting_epidemic_year  = 1994
# plt.style.use('ggplot')



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


def generate_ceac_plot(file_names, titles, folder_path, comparator, p_costs, p_effect, discount_rate, year_end_evaluation, year_end_results, tspan):
    wtp_values = np.linspace(0, 10000, 40)
    all_nmb = []

    for file in file_names:
        evaluate = extract_excel_files(os.path.join(base_path, folder_path, file), year_end_results - starting_epidemic_year + 1)
        cost_eval = np.array(probabilistic_all_populations(evaluate, p_costs, discount_rate, year_end_evaluation - starting_epidemic_year, year_end_results - starting_epidemic_year, tspan=tspan))
        cost_comp = np.array(probabilistic_all_populations(comparator, p_costs, discount_rate, year_end_evaluation - starting_epidemic_year, year_end_results - starting_epidemic_year, tspan=tspan))
        effect_eval = np.array(probabilistic_all_populations(evaluate, p_effect, discount_rate, year_end_evaluation - starting_epidemic_year, year_end_results - starting_epidemic_year, tspan=tspan))
        effect_comp = np.array(probabilistic_all_populations(comparator, p_effect, discount_rate, year_end_evaluation - starting_epidemic_year, year_end_results - starting_epidemic_year, tspan=tspan))

        incr_cost = cost_eval - cost_comp
        incr_effect = effect_comp - effect_eval
        nmb_matrix = np.array([wtp * incr_effect - incr_cost for wtp in wtp_values])
        all_nmb.append(nmb_matrix)

    all_nmb = np.array(all_nmb)  # shape: (n_scenarios, n_wtp, n_sim)
    ceac_probs = (all_nmb == np.max(all_nmb, axis=0)).mean(axis=2)  # shape: (n_scenarios, n_wtp)

    # abbreviated_Strategies = [
    #      "(S1)", "(S2)", "(S3)", "(S4)", "(S5)", "(S6)",
    #     "(S7)", "(S8)", "(S9)", "(S10)", "(S11)", "(S12)","Status Quo",
    # ]
    ### changed the name of abbreviated strategies consistent with text 
    abbreviated_Strategies = [
    "S1", "S2", "S3", "S4", "S5", "S6",
    "S7", "S8", "S9", "S10", "S11", "S12", "Status Quo"
    ]
    # markers = ['x', '*', 'H', 'h', 'p', '+', 'd', 'D', 'v', '^', 's', 'o', 'x', '*']
    # colours = ['black', 'red', 'blue', 'cyan', 'orange', 'green', 'green', 'orange',
    #            'black', 'blue', 'red', 'green', 'purple']

    ### matching colours and shapes
    markers = [ '*', 'H', 'h', 'p', '+', 'd', 'D', 'v', '^', 's', 'o', 'x','d']
    colours = ['red', 'blue', 'cyan', 'orange', 'green', 'green', 'orange', 'black', 'blue', 'red', 'green', 'purple','brown']

    # Sort to show Status Quo first
    strategy_order = np.argsort([0 if "Status quo" in t or "Status-quo" in t else 1 for t in titles])

    plt.figure(figsize=(10, 7))
    for idx in strategy_order:
        probs = ceac_probs[idx]
        plt.plot(
            wtp_values, probs,
            label=abbreviated_Strategies[idx] if idx < len(abbreviated_Strategies) else f"({idx})",
            marker=markers[idx % len(markers)],
            color=colours[idx % len(colours)],
            linewidth=1
        )

    plt.xlabel("Willingness-to-pay threshold", fontsize=14)
    plt.ylabel("% Iterations Cost-Effective", fontsize=14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize =12)
    plt.grid(True)

    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'${int(x):,}'))
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))

    plt.legend(
        title="",
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),
        ncol=7,
        frameon=False,
        columnspacing=0.8,
        handletextpad=0.5,
        borderpad=0.5,
        fontsize=10
    )

    plt.tight_layout()
    plt.savefig(os.path.join(base_path, "output/Chp3-figures/Acceptability curve_sens.png"), dpi=800)
    plt.show()
#### Costing based on extending the 2100


year_end_results = 2100
years = year_end_results - starting_epidemic_year +1
tspan = np.arange(0, years, 1)
discount_rate = 5/100

# #### Costing done based on the discounted into the future numbers 
# folder_path = "CEA_results/Scenarios_results"





def create_psa_plot(folder_path, comparator, evaluate_file, year_end_results, 
                    starting_epidemic_year,year_end_evaluation, p_costs, p_effect, discount_rate, tspan, xlabel,
                    wtp=4119, save_path='CEA_results/PSA_Diagrams/Comparisons_PSA_baseline.png',
                    ylimax = 175_000_000,xlimax = 30_000,maximising_axes = False, title=None):
    # Extract data from Excel files
    evaluate = extract_excel_files(os.path.join(base_path, folder_path, evaluate_file), 
                                   year_end_results - starting_epidemic_year + 1)
    
    # Calculate incremental costs and effects
    incremental_costs = np.array(probabilistic_all_populations(
        evaluate, p_costs, discount_rate, year_end_evaluation-starting_epidemic_year, year_end_results-starting_epidemic_year, tspan=tspan)) - np.array(probabilistic_all_populations(
        comparator, p_costs, discount_rate, year_end_evaluation-starting_epidemic_year, year_end_results-starting_epidemic_year, tspan=tspan))
    
    incremental_effects = np.array(probabilistic_all_populations(
        comparator, p_effect, discount_rate, year_end_evaluation-starting_epidemic_year, year_end_results-starting_epidemic_year, tspan=tspan)) - np.array(probabilistic_all_populations(
        evaluate, p_effect, discount_rate, year_end_evaluation-starting_epidemic_year, year_end_results-starting_epidemic_year, tspan=tspan))
    
    # Create plot
    plt.figure(figsize=(8, 5))
    
    # Add horizontal and vertical lines at 0
    plt.axhline(0, color='black', linewidth=1)
    plt.axvline(0, color='black', linewidth=1)
    
    # Scatter plot of incremental costs vs effects
    plt.scatter(incremental_effects, incremental_costs, s=10)
    plt.gca().get_yaxis().set_major_formatter(
        ticker.FuncFormatter(lambda x, p: format(int(x / 1_000_000), ',')))
    plt.gca().get_xaxis().set_major_formatter(
        ticker.FuncFormatter(lambda x, p: format(int(x / 1_000), ',')))
    
    # Add labels and title
    plt.xlabel(xlabel,fontsize= 14)
    plt.ylabel('Discounted Cost (millions USD)',fontsize=14)
    # plt.title('Incremental Cost vs Effect between 2025 and 2100')
    
    # Plot WTP line
    # x = np.linspace(min(min(incremental_effects),0), min(max(incremental_effects), max(incremental_costs) / wtp), 10)
    if maximising_axes == False:
        x = np.linspace(0, 1/2 * max(max(incremental_effects), xlimax), 10)
    else:
        x = np.linspace(0, min(1/2 * max(incremental_effects), max(incremental_costs) / wtp), 10)

    y = x * wtp
    plt.plot(x, y, label=f'WTP {wtp:,}', color="red", linestyle="--")
    
    if title is not None:
        plt.title(title, fontsize=16)

    

    ### set the limits in the y x axes
    if maximising_axes == False:
        plt.xlim(xmax=xlimax)
        plt.ylim(ymax=ylimax)
    # Adjust ticks and save the plot
    plt.xticks(rotation=45, fontsize = 12)
    plt.yticks(fontsize = 12)
    # plt.legend(loc="lower right")
    plt.savefig(save_path, dpi=500, bbox_inches='tight')
    
    # Show the plot
    plt.show()


comparator = extract_excel_files(os.path.join(base_path, folder_path, "Routine_2_Baseline_WithDolutegravir_56_56.xlsx"), year_end_results - starting_epidemic_year +1) 



file_names = ["POC_2_Baseline_WithDolutegravir_56_56.xlsx"
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

file_names = file_names +["Routine_2_Baseline_WithDolutegravir_56_56.xlsx"] 


titles = [
    "POC VL 20% (No DR testing)",
    "POC VL 40% (No DR testing)",
    "POC VL 60% (No DR testing)",
    "POC VL 20% + pre-treatment DR testing",
    "POC VL 40% + pre-treatment DR testing",
    "POC VL 60% + pre-treatment DR testing",
    "POC VL 20% + post-treatment DR testing",
    "POC VL 40% + post-treatment DR testing",
    "POC VL 60% + post-treatment DR testing",
    "POC VL 20% + pre- and post-treatment DR testing",
    "POC VL 40% + pre- and post-treatment DR testing",
    "POC VL 60% + pre- and post-treatment DR testing"
]

titles = titles +["Status quo: Routine VL 20% (no DR testing)"] 


#### All the PSA diagrams for DALYs modelling 
for evaluate_file, title in zip(file_names, titles):
    save_path = os.path.join(base_path,f"output/Chp3-figures/Chp3-PSA/DALYs_{evaluate_file.replace('.xlsx', '')}_PSA.png")
    
    create_psa_plot(
        folder_path="output/Chp3_scenarios",
        comparator=comparator,
        evaluate_file=evaluate_file,
        year_end_results=2100, 
        starting_epidemic_year=1994, 
        year_end_evaluation=2020,
        p_costs=p_costs, 
        p_effect=p_daly,  ## p_daly to calculate DALYs averted
        discount_rate=discount_rate, 
        tspan=tspan,
        xlabel="DALYs averted (thousands)",
        save_path=save_path,
        ylimax = 200_000_000,xlimax = 350_000,maximising_axes=True,
        title=title
    )


image_path = "output/Chp3-figures"
# Generate CEAC plot for all scenarios using shared comparator
generate_ceac_plot(
    file_names=file_names,
    titles=titles,
    folder_path=folder_path,
    comparator=comparator,
    p_costs=p_costs,
    p_effect=p_daly,
    discount_rate=discount_rate,
    year_end_evaluation=2020,
    year_end_results=2100,
    tspan=tspan
)



end_time = datetime.now()
print(f"Code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")

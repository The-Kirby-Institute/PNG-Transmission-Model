# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 11:55:07 2023

@author: nguye
"""
"""
This block of code is used to generate Excel sheets for the scenarios in Chapter 2 of the thesis.
It sets up the configurations for different scenarios, runs the ODE model for each scenario,
and saves the results to Excel files.

Note that the entire folder Chp2-Scenarios are ignored by git, due to the large size of the output files.
"""

### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

from model.meta_parameters import *

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import seaborn as sns

import importlib
import json


from scipy.interpolate import UnivariateSpline, CubicSpline
import random
from datetime import datetime

from matplotlib.ticker import MaxNLocator

#### Calculation mode = Posteriors -  read Posterior Samples 
from model.read_posterior_samples import*




from model import partial_functions
from model import meta_parameters
from model import core_scenarios
from model import Model_ode
from model import Model_init

from model.Model_ode import *
from model.Model_init import *
from model.core_scenarios import *


start_time = datetime.now()
print(f"Starting  code execution at {start_time}...")

random.seed(334)



starting_epidemic_year  = 1994

reference_year = 2000
adj =  reference_year - starting_epidemic_year


def setCorrectxAxis(frame, frequency_ticks =5,starting_position=0):
    plt.figure()
    default_x_ticks = list(range(len(frame.I_W)))
    new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+len(frame.I_W) ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, len(frame.I_W)+1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)


def update_config(Monitoring_Type, scenario, DR_testing_scenario,treatment,running_mode,analysis_mode,year_OffInfs,year_OffInterv):
    config = {"Monitoring_Type": Monitoring_Type, "scenario": scenario, "DR_testing_scenario": DR_testing_scenario,"treatment":treatment,'running_mode':running_mode,'analysis_mode':analysis_mode,"year_OffInfs":year_OffInfs,"year_OffInterv":year_OffInterv}
    with open(os.path.join(base_path, 'model', 'config.json'), 'w') as file:
        json.dump(config, file)


#### Start some sanity check of the models before running the full model
update_config("Routine", "2", "Baseline", "WithDolutegravir", "Continuous", "Inference",60,60)




#### Initialise all the appropriate parameters and conditions for the UNCERTAINTY analyses

S,I_W, I_DR, D_W, D_DR, T_W1,T_W2, T_DR1,T_DR2, F_W1,F_W2, F_TDR1,F_ADR1,F_DR2,incidence,deaths,total_PLHIV,incidence_Resist,preDRtest,postDRtest_routine,postDRtest_POC,VLtest_routine,VLtest_POC,newAcquiredDR,newTransmittedDR,newMisTreatedDR,diagnoses,treat_inits,all_mortality = initial_conditions    



sols = []
for i in range(len(theta_samples)):

    params_series = theta_samples[i], beta_u_samples[i], beta_t_samples[i], beta_f_samples[i], delta_U_samples[i], delta_T_samples[i], delta_F_samples[i],delta_B_samples[i], h1_samples[i], h2_samples[i], eta_1_samples[i], eta_2_samples[i], eta_3_samples[i],  eta_4_samples[i], eta_5_samples[i], eta_6_samples[i], b_asterisk_samples[i], b_k_samples[i],  c_asterisk_samples[i], c_k_samples[i],  rho_asterisk_samples[i], counsel_samples[i], transfer_2ndline_samples[i], f1_rate_samples[i], f2_rate_samples[i], f3_rate_samples[i], f4_rate_samples[i], mu1_samples[i], qt_samples[i],year_OffInfs


    result = ode_solver(params_series, initial_conditions, tspan,treatment)


    df = pd. DataFrame(result, columns=['S','I_W', 'I_DR', 'D_W', 'D_DR', 'T_W1','T_W2', 'T_DR1','T_DR2', 'F_W1','F_W2', 'F_TDR1','F_ADR1','F_DR2','incidence','deaths','total_PLHIV','incidence_Resist','preDRtest','postDRtest_routine','postDRtest_POC','VLtest_routine','VLtest_POC','newAcquiredDR','newTransmittedDR','newMisTreatedDR','diagnoses','treat_inits','all_mortality'])
    actual_population = compute_results(df)

    actual_population[actual_population.select_dtypes(include=['number']).columns] *= corrected_population['Total_Pop'][starting_epidemic_year]


    # Store the results
    sols.append(actual_population)





HIV_incident = np.array([np.gradient(sol.incidence ,tspan) for sol in sols])
median_Total = np.median(HIV_incident, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(HIV_incident, 25, axis=0)
Q3 = np.percentile(HIV_incident, 75, axis=0)

# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(HIV_incident, 2.5, axis=0)
P97_5 = np.percentile(HIV_incident, 97.5, axis=0)


setCorrectxAxis(actual_population,5,1)

# plt.plot(actual_population.incidence,c="blue", label="Modelled incident adults")

plt.plot(median_Total,c="black", label="Modelled newly infected adults",linewidth = 0.9)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)

### Updated UNAIDS/ SPECTRUM estimates in 2022
scatter_points = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Adults (15+) newly infected with HIV']]
lower_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing31']]
higher_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing32']]

scatter_points[scatter_points =="<500"] = np.nan
lower_CI[lower_CI =="<500"] = np.nan
higher_CI[higher_CI =="<500"] = np.nan
lower_CI[lower_CI =="<200"] = np.nan
higher_CI[higher_CI =="<200"] = np.nan

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))



x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Adults (15+) newly infected with HIV',c="red")
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='o',  capsize=2, c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates")


plt.grid(True)
plt.legend(loc="upper left",prop={'size': 10})
plt.ylabel('Number of newly infected HIV cases in PNG')
plt.subplots_adjust(bottom=0.15,left=0.12)
# plt.savefig('Bayesian Predictive Pictures/Incidence in PNG (adults)_with IQR range_all(production).png', dpi=500)
plt.savefig(os.path.join(base_path, 'output','diagnostics','Chp2_diag', 'Baseline_Newly infected HIV cases.png'), dpi=500)
plt.show()











Total = np.array([sol.Total for sol in sols])
median_Total = np.median(Total, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(Total, 25, axis=0)
Q3 = np.percentile(Total, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(Total, 2.5, axis=0)
P97_5 = np.percentile(Total, 97.5, axis=0)



setCorrectxAxis(actual_population,starting_position=1)
# plt.plot(mean_Total,c="blue", label="Mean Total PLHIV",linewidth = 0.5)
plt.plot(median_Total,c="black", label="Modelled Total Adults living with HIV",linewidth = 0.8)
# plt.fill_between(tspan, mean_Total - 1.96*std_Total, mean_Total + 1.96*std_Total, color='cornflowerblue', alpha=0.5, label='Uncertainty range around Modelled Total PLHIV (~95%)')
# max_array = max([sol['Total'] for sol in sols], key=lambda x: x.iloc[-1])
# min_array = min([sol['Total'] for sol in sols], key=lambda x: x.iloc[-1])
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)


    
#### Updated UNAIDS/ SPECTRUM estimates in 2022
scatter_points = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Estimated adults (15+) living with HIV']]
lower_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing19']]
higher_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing20']]

lower_CI[lower_CI =="<500"] = 0
higher_CI[higher_CI =="<500"] = 0


scatter_points = scatter_points["Estimated adults (15+) living with HIV"].tolist()
# scatter_points = [item for item in scatter_points if isinstance(item,int)] + [int(item.replace(" ",""))  for item in scatter_points if isinstance(item,str)]
lower_CI = lower_CI["Missing19"].tolist()
# lower_CI = [item for item in lower_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in lower_CI if isinstance(item,str)]
higher_CI = higher_CI["Missing20"].tolist()
# higher_CI = [item for item in higher_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in higher_CI if isinstance(item,str)]

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))


x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Estimated adults (15+) with HIV',c="yellow")
plt.errorbar(x_values_scatters,scatter_points,yerr=np.c_[lower_CI, higher_CI].T, fmt='--o',  capsize=2, c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates")



scatter_points = dat_estimates.loc[dat_estimates.index[dat_estimates.index >= starting_epidemic_year],['ACTUAL prevalence of HIV in reports']]
scatter_points = scatter_points["ACTUAL prevalence of HIV in reports"].tolist()
scatter_points += [np.nan] * (len(actual_population.Total) - len(scatter_points))
scatter_points = scatter_points*(actual_population.Total+actual_population.S)/100
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'PNG National Estimates',c="#DC143C")

plt.grid(True)
plt.ylabel('Number of Adults (15+) living with HIV in PNG')
plt.legend(loc="upper left",prop={'size': 10})
# plt.savefig('Bayesian Predictive Pictures/Uncertainty_graph_PLHIV_(production).png', dpi=500)
plt.savefig(os.path.join(base_path, 'output','diagnostics','Chp2_diag', 'Baseline_Adults living with HIV.png'), dpi=500)
plt.show()











##### Prevalence of HIV among adults older than 15 years old 

prevalenceHIV = np.array([(sol.Total)/(sol["Total"] + sol["S"] )*100 for sol in sols])
median_Total = np.median(prevalenceHIV, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(prevalenceHIV, 25, axis=0)
Q3 = np.percentile(prevalenceHIV, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(prevalenceHIV, 2.5, axis=0)
P97_5 = np.percentile(prevalenceHIV, 97.5, axis=0)


setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled Prevalence of HIV",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)




#### Updated UNAIDS/ SPECTRUM estimates in 2022
scatter_points = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Estimated adults (15+) living with HIV']]
lower_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing19']]
higher_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing20']]

lower_CI[lower_CI =="<500"] = 0
higher_CI[higher_CI =="<500"] = 0


scatter_points = scatter_points["Estimated adults (15+) living with HIV"].tolist()
# scatter_points = [item for item in scatter_points if isinstance(item,int)] + [int(item.replace(" ",""))  for item in scatter_points if isinstance(item,str)]
lower_CI = lower_CI["Missing19"].tolist()
# lower_CI = [item for item in lower_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in lower_CI if isinstance(item,str)]
higher_CI = higher_CI["Missing20"].tolist()
# higher_CI = [item for item in higher_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in higher_CI if isinstance(item,str)]

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI)).tolist()
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points)).tolist()

### divide by the number of adults greater than 15 years old 
df2 = corrected_population.loc[starting_epidemic_year:]['Total_Pop'].tolist()
scatter_points = [a/b*100 for a,b in zip(scatter_points,df2)]
lower_CI = [a/b*100 for a,b in zip(lower_CI,df2)]
higher_CI = [a/b*100 for a,b in zip(higher_CI,df2)]


x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Adults (15-49) prevalence',c="yellow")
plt.errorbar(x_values_scatters,scatter_points,yerr=[lower_CI, higher_CI], fmt='o',  capsize=2,c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates")
########

scatter_points = dat_estimates.loc[dat_estimates.index[dat_estimates.index >= starting_epidemic_year],['ACTUAL prevalence of HIV in reports']]
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'PNG National  estimate',c="#DC143C")

plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1f}%'))
plt.xticks(fontsize=10)
plt.ylabel("Prevalence of HIV in PNG \namong adults")
plt.grid(True)
# plt.ylim(0,1.5)
plt.yticks(np.arange(0, 2.2, 0.2))
# plt.legend(loc="lower right",prop={'size': 10})
plt.savefig(os.path.join(base_path, 'output','diagnostics','Chp2_diag', 'HIV prevalence (with CIs)_adults15+_IQR_all(production).png'), dpi=500)
plt.show()








##### Set up the Baseline reference condition for Bayesian model - Continuous running mode for Bayesian time-varying functions
##### These time-varying functions implement the various interventions/ trends that happened in HIV epidemics in PNG
runnning_mode = "Continuous"
analysis_mode = "Inference"

DR_testing_scenariosss = ["Timevarying","Baseline"]
Monitoring_Typesss = ["POC_AcquiredDR","Routine","POC"]
scenariosss = ["2","3","4"]
treatments = ["WithDolutegravir","NoDolutegravir"]



counter = 0


all_combine = []

for x in DR_testing_scenariosss:
    for y in Monitoring_Typesss:
        for z in scenariosss:
            if y == "Routine" and z in ["3", "4"]:
                continue
            for w in treatments:
                counter += 1
                print(counter)
                print(x,y,z)
                update_config(y, z, x,w,runnning_mode,analysis_mode,60,60)

                importlib.reload(core_scenarios)
                importlib.reload(partial_functions)
                importlib.reload(meta_parameters)

                ### import data before reloading the Posterior models

                importlib.reload(Model_ode)
                importlib.reload(Model_init)
                importlib.reload(meta_parameters)

                

                ##############################################
                ##### Modelling results from ODE compartmental model
                sols_actual = []
                # for i in range(len(theta_samples)):
                for i in range(len(theta_samples)):

                    params_series = theta_samples[i], beta_u_samples[i], beta_t_samples[i], beta_f_samples[i], delta_U_samples[i], delta_T_samples[i], delta_F_samples[i],delta_B_samples[i], h1_samples[i], h2_samples[i], eta_1_samples[i], eta_2_samples[i], eta_3_samples[i], eta_4_samples[i], eta_5_samples[i], eta_6_samples[i], b_asterisk_samples[i], b_k_samples[i], c_asterisk_samples[i], c_k_samples[i],  rho_asterisk_samples[i], counsel_samples[i], transfer_2ndline_samples[i], f1_rate_samples[i], f2_rate_samples[i], f3_rate_samples[i], f4_rate_samples[i], mu1_samples[i], qt_samples[i],year_OffInfs

                    result = ode_solver(params_series, initial_conditions, tspan,w)

                    df = pd. DataFrame(result, columns=['S','I_W', 'I_DR', 'D_W', 'D_DR', 'T_W1','T_W2', 'T_DR1','T_DR2', 'F_W1','F_W2', 'F_TDR1','F_ADR1','F_DR2','incidence','deaths','total_PLHIV','incidence_Resist','preDRtest','postDRtest_routine','postDRtest_POC','VLtest_routine','VLtest_POC','newAcquiredDR','newTransmittedDR','newMisTreatedDR','diagnoses','treat_inits','all_mortality'])

                    # Create the new column data
                    simulation_data = ["Simulation " + str(i) for _ in range(df.shape[0])]
                    
                    # Add the new column to the end of the DataFrame
                    df['Simulation'] = simulation_data

                    # proportion_df = copy.deepcopy(df)
                    actual_population = copy.deepcopy(df)


                    # ### Calculate the proportion Compartmenta model
                    # proportion_df['summation']  = proportion_df.iloc[:,0:14].sum(axis=1)
                    # proportion_df = proportion_df.iloc[:,0:].div(proportion_df.summation,axis=0)

                    # proportion_df = compute_results(proportion_df)
                    actual_population = compute_results(actual_population)

                    # proportion_df[proportion_df.select_dtypes(include=['number']).columns] *= 100
                    actual_population[actual_population.select_dtypes(include=['number']).columns] *= corrected_population['Total_Pop'][starting_epidemic_year]


                    # Store the results
                    sols_actual.append(actual_population)
                

                
                # Flatten the first dimension and reshape
                df = pd.concat(sols_actual, ignore_index=True)

                # filename = f"output/Posterior_results_B/All_Posteriors/{x}_{y}_{z}_{w}.xlsx"
                filename = os.path.join(base_path, 'output', 'Chp2_scenarios', f'{x}_{y}_{z}_{w}.xlsx')
                df.to_excel(filename, index=False)

                




                print(Monitoring_Type,scenario,DR_testing_scenario,treatment)

                Monitoring_Type = y
                scenario = z
                DR_testing_scenario = x
                treatment = w



end_time = datetime.now()
print(f"Code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")

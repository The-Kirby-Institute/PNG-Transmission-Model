### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

import json

def update_config(Monitoring_Type, scenario, DR_testing_scenario,treatment,running_mode,analysis_mode,year_OffInfs,year_OffInterv):
    config = {"Monitoring_Type": Monitoring_Type, "scenario": scenario, "DR_testing_scenario": DR_testing_scenario,"treatment":treatment,'running_mode':running_mode,'analysis_mode':analysis_mode,"year_OffInfs":year_OffInfs,"year_OffInterv":year_OffInterv}
    with open(os.path.join(base_path, 'model', 'config.json'), 'w') as file:
        json.dump(config, file)

starting_epidemic_year  = 1994
year_end_simulation = 2030
reference_year = 2000
adj =  reference_year - starting_epidemic_year
year_end_results = 2100
evaluation_period = year_end_simulation - starting_epidemic_year 

##### Set up the Baseline reference condition for Bayesian model - Continuous running mode for Bayesian time-varying functions
##### These time-varying functions implement the various interventions/ trends that happened in HIV epidemics in PNG
update_config("Routine", "2", "Baseline", "WithDolutegravir", "Study_3", "Inference",evaluation_period,evaluation_period)


from model.meta_parameters import *

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

import importlib

import random
from datetime import datetime


start_time = datetime.now()
print(f"Starting  code execution at {start_time}...")


random.seed(334)

from model.read_posterior_samples import*

from model import partial_functions
from model import meta_parameters
from model import core_scenarios
from model import Model_ode
from model import Model_init

from model.Model_ode import *
from model.Model_init import *
from model.core_scenarios import *

from EconEval.cea_functions import *
from EconEval.prb_data_matrix import*


#### determine the range of time the model being simulated for 
tspan = np.arange(0, year_end_results - starting_epidemic_year +1 , 1)

def setCorrectxAxis(frame, frequency_ticks =5,starting_position=0):
    plt.figure()
    default_x_ticks = list(range(len(frame.I_W)))
    new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+len(frame.I_W) ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, len(frame.I_W)+1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)



################################################################################################
################################################################################################
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
plt.savefig(os.path.join(base_path, 'output','diagnostics','Chp3_diag', 'Baseline_Newly infected HIV cases.png'), dpi=500)
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
scatter_points = [item for item in scatter_points if isinstance(item,int)] + [int(item.replace(" ",""))  for item in scatter_points if isinstance(item,str)]
lower_CI = lower_CI["Missing19"].tolist()
lower_CI = [item for item in lower_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in lower_CI if isinstance(item,str)]
higher_CI = higher_CI["Missing20"].tolist()
higher_CI = [item for item in higher_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in higher_CI if isinstance(item,str)]

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
plt.savefig(os.path.join(base_path, 'output','diagnostics','Chp3_diag', 'Baseline_Adults living with HIV.png'), dpi=500)
plt.show()









VL_percent = np.array([(sol["T"] * 100 / (sol["T"] + sol["F"])).fillna(0) for sol in sols])
median_Total = np.median(VL_percent, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(VL_percent, 25, axis=0)
Q3 = np.percentile(VL_percent, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(VL_percent, 2.5, axis=0)
P97_5 = np.percentile(VL_percent, 97.5, axis=0)


#### Viral Load Monitoring in PNG - combined Viral loads  graph of all Scenarios 
setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled Viral Suppression (%)",linewidth = 0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)



scatter_points = dat_testntreat.loc[dat_testntreat.index[dat_testntreat.index >= starting_epidemic_year],['Among people living with HIV, the percent with suppressed viral load']]
scatter_points[scatter_points == "..."] = np.nan
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'Reported viral suppression',c="red")



scatter_points = dat_testntreat.loc[dat_testntreat.index[dat_testntreat.index >= starting_epidemic_year],['Clinical Viral Suppresion in PNG']]
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points,c="#DC143C")
plt.grid(True)

plt.yticks(np.arange(0, 110, 10))
plt.legend(loc="lower right",prop={'size': 10})
plt.ylim(ymin=0)  # this line

plt.ylim(ymax=100)
plt.ylabel('Viral Suppression levels (%)')

plt.subplots_adjust(bottom=0.15)
# plt.savefig('Bayesian Predictive Pictures/Viral Suppression inPNG_IQR range_all(production).png', dpi=500)
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Baseline_Viral Suppression levels.png'), dpi=500)
plt.show()













TransmitDR = np.array([(sol.D_DR+sol.I_DR)/(sol["D"] + sol["I"])*100 for sol in sols])
median_Total = np.median(TransmitDR, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(TransmitDR, 25, axis=0)
Q3 = np.percentile(TransmitDR, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(TransmitDR, 2.5, axis=0)
P97_5 = np.percentile(TransmitDR, 97.5, axis=0)


setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled Drug Resistance in ART-naive",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)

scatter_points = dat_resistance.loc[dat_resistance.index[dat_resistance.index >= starting_epidemic_year],['Overall prevalence of PDR amongst those who are ART naïve (%)']]
lower_CI = dat_resistance.loc[dat_resistance.index[dat_resistance.index >= starting_epidemic_year],['Missing38']]
higher_CI = dat_resistance.loc[dat_resistance.index[dat_resistance.index >= starting_epidemic_year],['Missing39']]

scatter_points[scatter_points == "..."] = np.nan
x_values_scatters = list(range(0,len(scatter_points)))


#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='o',  capsize=2, c = "#DC143C",ecolor ="#DC143C",elinewidth = 0.3, alpha=0.8, label = "Reported prevalence of DR in ART-naive HIV people")

plt.ylabel("DR prevalence among Treatment-naive PLHIV (%)")
plt.legend(loc="upper left",prop={'size': 10})
plt.ylim(ymin=0,ymax=100)
plt.grid(True)
# plt.savefig('Bayesian Predictive Pictures/DR prevalence in Treatment-naive PLHIV_including IQR range_all(production).png', dpi=500)
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Baseline_DR prevalence transmitted.png'), dpi=500)
plt.show()







Deaths = np.array([np.gradient(sol.deaths,tspan) for sol in sols])
median_Total = np.median(Deaths, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(Deaths, 25, axis=0)
Q3 = np.percentile(Deaths, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(Deaths, 2.5, axis=0)
P97_5 = np.percentile(Deaths, 97.5, axis=0)

setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled deaths",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)



#### Updated UNAIDS/ SPECTRUM estimates in the 2022 
scatter_points = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['AIDS-related deaths among adults (15+)']]
lower_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing13']]
higher_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing14']]

scatter_points[scatter_points =="<100"] = np.nan
lower_CI[lower_CI =="<100"] = np.nan
higher_CI[higher_CI =="<100"] = np.nan

scatter_points[scatter_points =="<200"] = np.nan
lower_CI[lower_CI =="<200"] = np.nan
higher_CI[higher_CI =="<200"] = np.nan

scatter_points[scatter_points =="<500"] = np.nan
lower_CI[lower_CI =="<500"] = np.nan
higher_CI[higher_CI =="<500"] = np.nan



#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))


x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Adults (15-49) prevalence',c="yellow")
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='o',  capsize=2,c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates")




scatter_points = dat_estimates.loc[dat_estimates.index[dat_estimates.index >= starting_epidemic_year],['Re-estimated Mortality of HIV in GBD']]
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'Deaths estimated by Global Burden Disease (GBD)',c="cyan")

plt.grid(True)
plt.ylabel("Persons")
plt.legend(loc="upper left",prop={'size': 7})
plt.ylabel('Number of HIV-attributable deaths among\n adults (15+) in PNG')
# plt.savefig('Bayesian Predictive Pictures/Deaths in HIV_including IQR range_all(production).png', dpi=500)
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Baseline_New deaths in PNG.png'), dpi=500)
plt.show()







onART = np.array([sol.TreatTotal for sol in sols])
median_Total = np.median(onART, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(onART, 25, axis=0)
Q3 = np.percentile(onART, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(onART, 2.5, axis=0)
P97_5 = np.percentile(onART, 97.5, axis=0)


setCorrectxAxis(actual_population,5,1)
##### Total size of population 
plt.plot(median_Total,c="black", label="Modelled on ART",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)

# scatter_points = dat_testntreat.loc[dat_testntreat.index[dat_testntreat.index >= starting_epidemic_year],['Missing84']]
##### Updated figures in 2022 numbers
scatter_points = dat_testntreat_2023.loc[dat_testntreat_2023.index[dat_testntreat_2023.index >= starting_epidemic_year],['Missing80']]
scatter_points[scatter_points == "..."] = np.nan
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'UNAIDS/Spectrum Estimates',c="#0F52BA")


plt.grid(True)

plt.ylabel('Number of Adults with HIV on ART')
plt.subplots_adjust(bottom=0.15)

plt.legend(loc="upper left",prop={'size': 8})
# plt.savefig('Bayesian Predictive Pictures/on ART (adults)_IQR all(production).png', dpi=500)
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Baseline_on ART.png'), dpi=500)
plt.show()



#### trends in Total population in 


Total_PNG = np.array([sol.Total+sol.S for sol in sols])
median_Total = np.median(Total_PNG, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(Total_PNG, 25, axis=0)
Q3 = np.percentile(Total_PNG, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(Total_PNG, 2.5, axis=0)
P97_5 = np.percentile(Total_PNG, 97.5, axis=0)



setCorrectxAxis(actual_population,5,1)
##### Total size of population 
plt.plot(median_Total,c="red", label="Modelled Total PNG population ",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3, label='Uncertainty range around Total PNG population (~95%)')
plt.fill_between(tspan, Q1, Q3, color='cornflowerblue', alpha=0.5, label='IQR range around Total PNG population')

scatter_points = corrected_population.loc[corrected_population.index[corrected_population.index >= starting_epidemic_year],['Total_Pop']]
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'Literature',c="red")
plt.legend(loc="lower left",prop={'size': 5})
plt.ylabel("Persons (Millions)")
plt.xlabel("Year")
plt.title('Total population in PNG')
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Total pop in PNG_trends.png'), dpi=500)
plt.show()

















discount_rate = 5/100 

prob_discounted_cost = probabilistic_all_populations_aslist(sols, p_costs, discount_rate,0,year_end_results - starting_epidemic_year,tspan)

total_costs = prob_discounted_cost

median_Total = np.median(total_costs, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(total_costs, 25, axis=0)
Q3 = np.percentile(total_costs, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(total_costs, 2.5, axis=0)
P97_5 = np.percentile(total_costs, 97.5, axis=0)

setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled costs",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)

# plt.ticklabel_format(style='plain')
plt.gca().get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x / 1000), ',')))
plt.ylabel('Cost (in thousands USD)')


plt.grid(True)
plt.legend(loc="upper left",prop={'size': 7})
# plt.savefig('CEA_results/Baseline_Trends of Costs.png', dpi=500)
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Baseline_Trends of Costs.png'), dpi=500)
plt.show()




prob_discounted_daly = probabilistic_all_populations_aslist(sols, p_daly, discount_rate,0,year_end_results - starting_epidemic_year,tspan)

totals = prob_discounted_daly

median_Total = np.median(totals, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(totals, 25, axis=0)
Q3 = np.percentile(totals, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(totals, 2.5, axis=0)
P97_5 = np.percentile(totals, 97.5, axis=0)

setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled DALYs",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)

# plt.ticklabel_format(style='plain')
# plt.gca().get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x / 1000), ',')))
plt.ylabel('DALYs')


plt.grid(True)
plt.legend(loc="upper left",prop={'size': 7})
# plt.savefig('CEA_results/Baseline_Trends of DALYs.png', dpi=500)
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Baseline_Trends of DALYs.png'), dpi=500)
plt.show()






ratios = total_costs/totals

median_Total = np.median(ratios, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(ratios, 25, axis=0)
Q3 = np.percentile(ratios, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(ratios, 2.5, axis=0)
P97_5 = np.percentile(ratios, 97.5, axis=0)

setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled Cost/DALYs",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)

# plt.ticklabel_format(style='plain')
# plt.gca().get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x / 1000), ',')))
plt.ylabel('ICERs')


plt.grid(True)
plt.legend(loc="upper left",prop={'size': 7})
# plt.savefig('CEA_results/Baseline_Trends of Ratios.png', dpi=500)
plt.savefig(os.path.join(base_path, 'output', 'diagnostics', 'Chp3_diag', 'Baseline_Trends of Ratios.png'), dpi=500)
plt.show()



# Flatten the first dimension and reshape - write the baseline combinations in the model
df = pd.concat(sols, ignore_index=True)
# filename = f"CEA_results/Scenarios_results/{Monitoring_Type}_{scenario}_{DR_testing_scenario}_{treatment}_{year_OffInfs}_{year_OffInterv}.xlsx"
filename = os.path.join(base_path, 'output', 'Chp3_scenarios', f'{Monitoring_Type}_{scenario}_{DR_testing_scenario}_{treatment}_{year_OffInfs}_{year_OffInterv}.xlsx')
df.to_excel(filename, index=False)



runnning_mode = "Continuous"
analysis_mode = "Inference"

counter = 0

policy_evl = year_end_simulation - starting_epidemic_year 
direct_evl = 2025 - starting_epidemic_year 
designated_policy_eval = 2035 - starting_epidemic_year 
new_policy_evl = 2050 - starting_epidemic_year 
# Define the specific combinations explicitly
combinations = [
    #### Time evaluation of 5 year horizon only
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},

    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},

    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},

    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},


    ### Time evaluation of 10 year horizon only
    {"preartDR": "Baseline", "Monitoring_Type": "Routine", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},

    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},

    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},

    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":designated_policy_eval,"year_OffInterv":designated_policy_eval},


    ### Time evaluation of 25 year horizon only
    {"preartDR": "Baseline", "Monitoring_Type": "Routine", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},

    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},

    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Baseline", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},

    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "3", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},
    {"preartDR": "Timevarying", "Monitoring_Type": "POC_AcquiredDR", "Scenario": "4", "Treatment": "WithDolutegravir","yearOffInfs":new_policy_evl,"year_OffInterv":new_policy_evl},

    
    ### No Dolutegravir? 
    {"preartDR": "Baseline", "Monitoring_Type": "Routine", "Scenario": "2", "Treatment": "NoDolutegravir","yearOffInfs":policy_evl,"year_OffInterv":policy_evl},
    ### Turn off infections after 2025 and turn off evaluation 2030
    {"preartDR": "Baseline", "Monitoring_Type": "Routine", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":direct_evl,"year_OffInterv":policy_evl},
        ### Turn off infections after 2025 and turn off evaluation 2035
    {"preartDR": "Baseline", "Monitoring_Type": "Routine", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":direct_evl,"year_OffInterv":designated_policy_eval},
    ### Turn off infections after 2025 and turn off evaluation 2050
    {"preartDR": "Baseline", "Monitoring_Type": "Routine", "Scenario": "2", "Treatment": "WithDolutegravir","yearOffInfs":direct_evl,"year_OffInterv":new_policy_evl}


]
running_mode = "Study_3"

for comb in combinations:
    counter += 1
    print(counter)

    print (comb['preartDR'],comb['Monitoring_Type'],comb['Scenario'])

    ##### Update the config to rerun the models 
    update_config(comb['Monitoring_Type'],comb['Scenario'],comb['preartDR'],comb['Treatment'], running_mode,analysis_mode,comb['yearOffInfs'],comb['year_OffInterv'])


    # Reload necessary modules
    importlib.reload(core_scenarios)
    importlib.reload(partial_functions)
    importlib.reload(meta_parameters)
    importlib.reload(Model_ode)
    importlib.reload(Model_init)
    importlib.reload(meta_parameters)

    sols = []
    for i in range(len(theta_samples)):

        params_series = theta_samples[i], beta_u_samples[i], beta_t_samples[i], beta_f_samples[i], delta_U_samples[i], delta_T_samples[i], delta_F_samples[i],delta_B_samples[i], h1_samples[i], h2_samples[i], eta_1_samples[i], eta_2_samples[i], eta_3_samples[i],  eta_4_samples[i], eta_5_samples[i], eta_6_samples[i], b_asterisk_samples[i], b_k_samples[i],  c_asterisk_samples[i], c_k_samples[i],  rho_asterisk_samples[i], counsel_samples[i], transfer_2ndline_samples[i], f1_rate_samples[i], f2_rate_samples[i], f3_rate_samples[i], f4_rate_samples[i], mu1_samples[i], qt_samples[i],comb['yearOffInfs']


        result = ode_solver(params_series, initial_conditions, tspan,treatment)


        df = pd. DataFrame(result, columns=['S','I_W', 'I_DR', 'D_W', 'D_DR', 'T_W1','T_W2', 'T_DR1','T_DR2', 'F_W1','F_W2', 'F_TDR1','F_ADR1','F_DR2','incidence','deaths','total_PLHIV','incidence_Resist','preDRtest','postDRtest_routine','postDRtest_POC','VLtest_routine','VLtest_POC','newAcquiredDR','newTransmittedDR','newMisTreatedDR','diagnoses','treat_inits','all_mortality'])
        actual_population = compute_results(df)

        actual_population[actual_population.select_dtypes(include=['number']).columns] *= corrected_population['Total_Pop'][starting_epidemic_year]


        # Store the results
        sols.append(actual_population)

    # Flatten the first dimension and reshape - write the baseline combinations in the model
    df = pd.concat(sols, ignore_index=True)
    ### If this is the last element in a list, then it should be the scenario without infections since 2025
    # filename = f"CEA_results/Scenarios_results/{comb['Monitoring_Type']}_{comb['Scenario']}_{comb['preartDR']}_{comb['Treatment']}_{comb['yearOffInfs']}_{comb['year_OffInterv']}.xlsx"
    filename = os.path.join(base_path, 'output', 'Chp3_scenarios', f"{comb['Monitoring_Type']}_{comb['Scenario']}_{comb['preartDR']}_{comb['Treatment']}_{comb['yearOffInfs']}_{comb['year_OffInterv']}.xlsx")
    df.to_excel(filename, index=False)


end_time = datetime.now()
print(f"Code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")

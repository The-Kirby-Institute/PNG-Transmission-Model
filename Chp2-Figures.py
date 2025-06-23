# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 21:37:33 2023

@author: nguye
"""
import json

def update_config(Monitoring_Type, scenario, DR_testing_scenario,treatment,running_mode,analysis_mode,year_OffInfs,year_OffInterv):
    config = {"Monitoring_Type": Monitoring_Type, "scenario": scenario, "DR_testing_scenario": DR_testing_scenario,"treatment":treatment,'running_mode':running_mode,'analysis_mode':analysis_mode,"year_OffInfs":year_OffInfs,"year_OffInterv":year_OffInterv}
    with open('model/config.json', 'w') as file:
        json.dump(config, file)


##### Set up the Baseline reference condition for Bayesian model - Continuous running mode for Bayesian time-varying functions
##### These time-varying functions implement the various interventions/ trends that happened in HIV epidemics in PNG
##### Continuous functions gave smooth/continuous implement for these trends
##### While discrete functions are binary/boolean functions that have sharper - but more precise- change
update_config("Routine", "2", "Baseline", "WithDolutegravir", "Continuous", "Inference",60,60)


from model.meta_parameters import *

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats 

from model.Model_init import *
from scipy.stats import gaussian_kde



from scipy.interpolate import UnivariateSpline, CubicSpline
import random

from matplotlib.ticker import MaxNLocator


random.seed(334)

mode_cal = "Posteriors"
from model.read_posterior_samples import*


starting_epidemic_year  = 1994
# plt.style.use('ggplot')

reference_year = 2000
adj =  reference_year - starting_epidemic_year


def setCorrectxAxis(frame, frequency_ticks =5,starting_position=0):
    plt.figure()
    default_x_ticks = list(range(len(frame.I_W)))
    new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+len(frame.I_W) ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, len(frame.I_W)+1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)




excel_path = 'Bayesian Pictures/Bayesian_posterior.xlsx' 
parameter_name = 'theta_samples' 



################################################################################################
################################################################################################
#### Initialise all the appropriate parameters and conditions for the UNCERTAINTY analyses

S,I_W, I_DR, D_W, D_DR, T_W1,T_W2, T_DR1,T_DR2, F_W1,F_W2, F_TDR1,F_ADR1,F_DR2,incidence,deaths,total_PLHIV,incidence_Resist,preDRtest,postDRtest_routine,postDRtest_POC,VLtest_routine,VLtest_POC,newAcquiredDR,newTransmittedDR,newMisTreatedDR,diagnoses,treat_inits,all_mortality = initial_conditions    



sols = []
for i in range(len(theta_samples)):

    params_series = theta_samples[i], beta_u_samples[i], beta_t_samples[i], beta_f_samples[i], delta_U_samples[i], delta_T_samples[i], delta_F_samples[i], delta_B_samples[i],h1_samples[i], h2_samples[i], eta_1_samples[i], eta_2_samples[i], eta_3_samples[i],  eta_4_samples[i], eta_5_samples[i], eta_6_samples[i], b_asterisk_samples[i], b_k_samples[i],  c_asterisk_samples[i], c_k_samples[i],  rho_asterisk_samples[i], counsel_samples[i], transfer_2ndline_samples[i], f1_rate_samples[i], f2_rate_samples[i], f3_rate_samples[i], f4_rate_samples[i], mu1_samples[i], qt_samples[i],year_OffInfs


    result = ode_solver(params_series, initial_conditions, tspan,treatment)


    df = pd. DataFrame(result, columns=['S','I_W', 'I_DR', 'D_W', 'D_DR', 'T_W1','T_W2', 'T_DR1','T_DR2', 'F_W1','F_W2', 'F_TDR1','F_ADR1','F_DR2','incidence','deaths','total_PLHIV','incidence_Resist','preDRtest','postDRtest_routine','postDRtest_POC','VLtest_routine','VLtest_POC','newAcquiredDR','newTransmittedDR','newMisTreatedDR','diagnoses','treat_inits','all_mortality'])
    actual_population = compute_results(df)

    actual_population[actual_population.select_dtypes(include=['number']).columns] *= corrected_population['Total_Pop'][starting_epidemic_year]


    # Store the results
    sols.append(actual_population)

writing_sols = sols




































##### All graphs for IQR analyses in the model 


Total = np.array([sol.Total for sol in writing_sols])
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


    
#### Updated UNAIDS/ SPECTRUM estimates in 2024
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

plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.xticks(fontsize=10)
plt.grid(True)

plt.ylabel('Number of Adults (15+) living with HIV in PNG')

# plt.legend(loc="upper left",prop={'size': 10})
plt.savefig('output/figures/Uncertainty_graph_PLHIV_(production).png', dpi=500)
plt.show()


















VL_percent = np.array([(sol["T"] * 100 / (sol["T"] + sol["F"])).fillna(0) for sol in writing_sols])
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




### Updated UNAIDS/ SPECTRUM estimates in 2024
scatter_points = dat_testntreat_2023.loc[dat_testntreat_2023.index[dat_testntreat_2023.index >= starting_epidemic_year],['Missing64']]
lower_CI = dat_testntreat_2023.loc[dat_testntreat_2023.index[dat_testntreat_2023.index >= starting_epidemic_year],['Missing65']]
higher_CI = dat_testntreat_2023.loc[dat_testntreat_2023.index[dat_testntreat_2023.index >= starting_epidemic_year],['Missing66']]

scatter_points[scatter_points =="..."] = np.nan
lower_CI[lower_CI =="..."] = np.nan
higher_CI[higher_CI =="..."] = np.nan
higher_CI[higher_CI ==">98"] = 98


#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))



x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Adults (15+) newly infected with HIV',c="red")
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='o',  capsize=2, c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates")

# plt.axhline(y=95, color='r', linestyle='-.',linewidth=0.8)

plt.yticks(np.arange(0, 110, 10))
# plt.legend(loc="lower right",prop={'size': 10})
plt.ylim(ymin=0)  # this line

plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
plt.xticks(fontsize=10)

plt.ylim(ymax=100)
plt.ylabel('Level of suppressed VL')

plt.subplots_adjust(bottom=0.15)

plt.savefig('output/figures/Third Goal UNAIDS(production).png', dpi=500)
plt.show()























HIV_incident = np.array([np.gradient(sol.incidence ,tspan) for sol in writing_sols])
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

### Updated UNAIDS/ SPECTRUM estimates in 2024
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
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.xticks(fontsize=10)


plt.grid(True)
# plt.legend(loc="upper left",prop={'size': 10})
plt.ylabel('Number of new HIV infections')
plt.subplots_adjust(bottom=0.15,left=0.12)
plt.savefig('output/figures/Incidence in PNG (adults)_with IQR range_all(production).png', dpi=500)
plt.show()






incidenceHIV = np.array([(np.gradient(sol.incidence ,tspan))/np.array(sol.S.tolist()) *1000 for sol in writing_sols])
median_Total = np.median(incidenceHIV, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(incidenceHIV, 25, axis=0)
Q3 = np.percentile(incidenceHIV, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(incidenceHIV, 2.5, axis=0)
P97_5 = np.percentile(incidenceHIV, 97.5, axis=0)


setCorrectxAxis(actual_population,5,1)
plt.plot(median_Total,c="black", label="Modelled Incidence of HIV",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)


## Updated from UNAIDS/ Spectrum 2022
scatter_points = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Adults (15-49) incidence (per 1000 uninfected population)']]
lower_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing23']]
higher_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing24']]


#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))


x_values_scatters = list(range(0,len(scatter_points)))
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='--o',  capsize=2, c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates")


# plt.ylim(0,0.17)
# plt.yticks(np.arange(0, 0.17, 0.02))
plt.grid(True)
plt.legend(loc="upper right",prop={'size': 6})
plt.xlabel('Year')
plt.ylabel('Incidence (per 1000 uninfected population)')

plt.title('Adults (15+) Incidence (%) of HIV in PNG')
plt.savefig('output/figures/HIV incidence_adults_IQRsall.png', dpi=500)
plt.show()












Deaths = np.array([np.gradient(sol.deaths,tspan) for sol in writing_sols])
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

plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.xticks(fontsize=10)

plt.grid(True)
plt.ylabel("Persons")
# plt.legend(loc="upper left",prop={'size': 7})
plt.ylabel('Number of HIV-attributable deaths among\n adults (15+) in PNG')
plt.savefig('output/figures/Deaths in HIV_including IQR range_all(production).png', dpi=500)
plt.show()





















TransmitDR = np.array([(sol.D_DR+sol.I_DR)/(sol["D"] + sol["I"])*100 for sol in writing_sols])
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

plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
plt.xticks(fontsize=10)

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='o',  capsize=2, c = "#DC143C",ecolor ="#DC143C",elinewidth = 0.3, alpha=0.8, label = "Reported prevalence of DR in ART-naive HIV people")

plt.ylabel("Prevalence of transmitted DR")
# plt.legend(loc="upper left",prop={'size': 10})
plt.ylim(ymin=0,ymax=100)
plt.grid(True)
plt.savefig('output/figures/DR prevalence in Treatment-naive PLHIV_including IQR range_all(production).png', dpi=500)
plt.show()











combined_graphs= []
beta_f = copy.deepcopy(beta)
beta_r_n = copy.deepcopy(beta_r)
# rho_f = copy.deepcopy(rho)
setCorrectxAxis(actual_population,5,1)



##### extraction of each diagram, and then plot those diagram here from the posteriors
for i in range(len(theta_samples)):
    uncertainty_beta = partial(Richardsfunc,a=eta_1_samples[i],b= eta_2_samples[i],k=0,q=eta_3_samples[i])
    uncertainty_beta_r = partial(Richardsfunc,b=eta_5_samples[i], k=eta_4_samples[i],a = 0, q = eta_6_samples[i])
    # alternative
    # uncertainty_beta = partial(invRichardsfunc_fixc,a=eta_1_samples[i],b= eta_2_samples[i],k=0,q=eta_3_samples[i])
    # uncertainty_beta_r = partial(Richardsfunc,b=beta_q_samples[i], k=eta_4_samples[i],a = 0, q = k_samples[i])

    combined_graphs.append(1 + uncertainty_beta(tspan) +uncertainty_beta_r(tspan))

min_combined = np.min(combined_graphs, axis=0)
max_combined = np.max(combined_graphs, axis=0)
P2_5_combined = np.percentile(combined_graphs, 2.5, axis=0)
P97_5_combined = np.percentile(combined_graphs, 97.5, axis=0)
median_combined = np.median(combined_graphs, axis=0)

plt.fill_between(tspan, min_combined, max_combined, color='gray', alpha=0.3, label='Total Force of infection combined')
plt.fill_between(tspan, P2_5_combined, P97_5_combined, color='mediumseagreen', alpha=0.5, label='IQR Force of infection combined')
plt.plot(median_combined,c="black",label="Median Force of Infection combined",linewidth =1)



plt.ylim(ymin=0)  # this line for ymin = 0
plt.axhline(y=1, color='black', linestyle='-.',linewidth =0.5)
plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=12))  
plt.xlabel('Year')
# plt.legend(loc="upper right",prop={'size': 6})
# plt.title('Combined Force of Infection coefficient in PNG')
plt.savefig('output/figures/Combined Force of Infection_Appendix.png', dpi=500)
plt.show()











##### Plotting the incidence from UNAIDS previous fit with itself 




setCorrectxAxis(actual_population,5,1)


### Updated UNAIDS/ SPECTRUM estimates in 2024
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
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='-o',  capsize=2, c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates 2024")


#### Incidence reported by UNAIDS in 2023

scatter_points = dat_depraciated_2023.loc[dat_depraciated_2023.index[dat_depraciated_2023.index >= starting_epidemic_year],['Adults (15+) newly infected with HIV']]
lower_CI = dat_depraciated_2023.loc[dat_depraciated_2023.index[dat_depraciated_2023.index >= starting_epidemic_year],['Missing31']]
higher_CI = dat_depraciated_2023.loc[dat_depraciated_2023.index[dat_depraciated_2023.index >= starting_epidemic_year],['Missing32']]

scatter_points[scatter_points =="<500"] = np.nan
lower_CI[lower_CI =="<500"] = np.nan
higher_CI[higher_CI =="<500"] = np.nan
lower_CI[lower_CI =="<200"] = np.nan
higher_CI[higher_CI =="<200"] = np.nan

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))



x_values_scatters = list(range(0,len(scatter_points)))
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='-o',  capsize=2, c = "#ff7f0e",ecolor ="#ff7f0e",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates 2023")


#### Incidence reported by UNAIDS in 2022
scatter_points = dat_depraciated_2022.loc[dat_depraciated_2022.index[dat_depraciated_2022.index >= starting_epidemic_year],['Adults (15+) newly infected with HIV']]
lower_CI = dat_depraciated_2022.loc[dat_depraciated_2022.index[dat_depraciated_2022.index >= starting_epidemic_year],['Missing31']]
higher_CI = dat_depraciated_2022.loc[dat_depraciated_2022.index[dat_depraciated_2022.index >= starting_epidemic_year],['Missing32']]

scatter_points[scatter_points =="<500"] = np.nan
lower_CI[lower_CI =="<500"] = np.nan
higher_CI[higher_CI =="<500"] = np.nan
lower_CI[lower_CI =="<200"] = np.nan
higher_CI[higher_CI =="<200"] = np.nan

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))



x_values_scatters = list(range(0,len(scatter_points)))
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='-o',  capsize=2, c = "#2ca02c",ecolor ="#2ca02c",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates 2022")



plt.grid(True)
plt.legend(loc="upper left",prop={'size': 8})
plt.ylabel('Number of newly infected HIV cases in PNG')
plt.subplots_adjust(bottom=0.15,left=0.12)
plt.savefig('output/figures/UNAIDS HIV Incidence to itself.png', dpi=500)
plt.show()





















#### Comparing UNAIDS estimates to itself 


setCorrectxAxis(actual_population,5,1)

#### Updated UNAIDS/ SPECTRUM estimates in the 2022 
scatter_points = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Adults (15-49) prevalence (%)']]
lower_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing3']]
higher_CI = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Missing4']]

scatter_points[scatter_points =="<0.1"] = np.nan
lower_CI[lower_CI =="<0.1"] = np.nan
higher_CI[higher_CI =="<0.1"] = np.nan

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))


x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Adults (15-49) prevalence',c="yellow")
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='-o',  capsize=2,c = "#0F52BA",ecolor ="#0F52BA",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates 2024")
########

scatter_points = dat_estimates.loc[dat_estimates.index[dat_estimates.index >= starting_epidemic_year],['ACTUAL prevalence of HIV in reports']]
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'PNG National  estimate',c="#DC143C")




#### using the depractiated model from UNAIDS 2023

scatter_points = dat_depraciated_2023.loc[dat_depraciated_2023.index[dat_depraciated_2023.index >= starting_epidemic_year],['Adults (15-49) prevalence (%)']]
lower_CI = dat_depraciated_2023.loc[dat_depraciated_2023.index[dat_depraciated_2023.index >= starting_epidemic_year],['Missing3']]
higher_CI = dat_depraciated_2023.loc[dat_depraciated_2023.index[dat_depraciated_2023.index >= starting_epidemic_year],['Missing4']]

scatter_points[scatter_points =="<0.1"] = np.nan
lower_CI[lower_CI =="<0.1"] = np.nan
higher_CI[higher_CI =="<0.1"] = np.nan

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))


x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Adults (15-49) prevalence',c="yellow")
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='-o',  capsize=2,c = "#ff7f0e",ecolor ="#ff7f0e",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates from 2023")

#######

scatter_points = dat_depraciated_2022.loc[dat_depraciated_2022.index[dat_depraciated_2022.index >= starting_epidemic_year],['Adults (15-49) prevalence (%)']]
lower_CI = dat_depraciated_2022.loc[dat_depraciated_2022.index[dat_depraciated_2022.index >= starting_epidemic_year],['Missing3']]
higher_CI = dat_depraciated_2022.loc[dat_depraciated_2022.index[dat_depraciated_2022.index >= starting_epidemic_year],['Missing4']]

scatter_points[scatter_points =="<0.1"] = np.nan
lower_CI[lower_CI =="<0.1"] = np.nan
higher_CI[higher_CI =="<0.1"] = np.nan

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI))
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points))


x_values_scatters = list(range(0,len(scatter_points)))
# plt.scatter(x_values_scatters,scatter_points, label = 'Adults (15-49) prevalence',c="yellow")
plt.errorbar(x_values_scatters,scatter_points.to_numpy().flatten(),yerr=[lower_CI.to_numpy().flatten(), higher_CI.to_numpy().flatten()], fmt='-o',  capsize=2,c = "#2ca02c",ecolor ="#2ca02c",elinewidth = 0.3, alpha=0.8, label = "UNAIDS/Spectrum Estimates from 2022")




plt.ylabel("Prevalence of HIV in PNG amongst \nAdults (15-49) - %")
plt.grid(True)
# plt.ylim(0,1.5)
plt.yticks(np.arange(0, 2.2, 0.2))
plt.legend(loc="upper left",prop={'size': 8})
plt.savefig('output/figures/UNAIDS HIV prevalence to itself.png', dpi=500)
plt.show()





##### Prevalence of HIV among adults older than 15 years old 

prevalenceHIV = np.array([(sol.Total)/(sol["Total"] + sol["S"] )*100 for sol in writing_sols])
# prevalenceHIV = np.array([(sol.total_PLHIV)/(sol["total_PLHIV"] + sol["S"] )*100 for sol in writing_sols])
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
scatter_points = [item for item in scatter_points if isinstance(item,int)] + [int(item.replace(" ",""))  for item in scatter_points if isinstance(item,str)]
lower_CI = lower_CI["Missing19"].tolist()
lower_CI = [item for item in lower_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in lower_CI if isinstance(item,str)]
higher_CI = higher_CI["Missing20"].tolist()
higher_CI = [item for item in higher_CI if isinstance(item,int)] + [int(item.replace(" ",""))  for item in higher_CI if isinstance(item,str)]

#### calculation_lower_upper_limits
lower_CI = np.subtract(scatter_points, np.asarray(lower_CI)).tolist()
higher_CI =  np.subtract(higher_CI,np.asarray(scatter_points)).tolist()

### divide by the number of adults greater than 15 years old 
df2 = corrected_population.loc[1990:]['Total_Pop'].tolist()
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
plt.savefig('output/figures/HIV prevalence (with CIs)_adults15+_IQR_all(production).png', dpi=500)
plt.show()































onART = np.array([sol.TreatTotal for sol in writing_sols])
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

plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.xticks(fontsize=10)

# plt.legend(loc="upper left",prop={'size': 8})
plt.savefig('output/figures/on ART (adults)_IQR all(production).png', dpi=500)
plt.show()





















Diagnosed = np.array([sol.Aware for sol in writing_sols])
median_Total = np.median(Diagnosed, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(Diagnosed, 25, axis=0)
Q3 = np.percentile(Diagnosed, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(Diagnosed, 2.5, axis=0)
P97_5 = np.percentile(Diagnosed, 97.5, axis=0)

setCorrectxAxis(actual_population,5,1)
##### Total size of population 
plt.plot(median_Total,c="black", label="Modelled Detected with HIV ",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5)

# scatter_points = dat_testntreat.loc[dat_testntreat.index[dat_testntreat.index >= starting_epidemic_year],['Missing76']]
###### Updated figures in 2022 numbers 
scatter_points = dat_testntreat_2023.loc[dat_testntreat_2023.index[dat_testntreat_2023.index >= starting_epidemic_year],['Missing76']]
scatter_points[scatter_points == "..."] = np.nan
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points, label = 'UNAIDS/Spectrum Estimates',c="#0F52BA")

plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.xticks(fontsize=10)
plt.grid(True)

plt.ylabel('Number of Adults Diagnosed with HIV in PNG')
plt.subplots_adjust(bottom=0.15,left=0.14)

# plt.legend(loc="upper left",prop={'size': 8})
plt.savefig('output/figures/Aware of HIV (adults)_IRQ all(production).png', dpi=500)
plt.show()







Total_PNG = np.array([sol.Total+sol.S for sol in writing_sols])
median_Total = np.median(Total_PNG, axis=0)
# Calculate the first and third quartiles
Q1 = np.percentile(Total_PNG, 25, axis=0)
Q3 = np.percentile(Total_PNG, 75, axis=0)
# Calculate the 2.5th and 97.5th percentiles
P2_5 = np.percentile(Total_PNG, 2.5, axis=0)
P97_5 = np.percentile(Total_PNG, 97.5, axis=0)



setCorrectxAxis(actual_population,5,1)
##### Total size of population 
plt.plot(median_Total,c="black", label="Modelled Total PNG population ",linewidth=0.8)
plt.fill_between(tspan, P2_5, P97_5, color='gray', alpha=0.3)
plt.fill_between(tspan, Q1, Q3, color='mediumseagreen', alpha=0.5 )

scatter_points = corrected_population.loc[corrected_population.index[corrected_population.index >= starting_epidemic_year],['Total_Pop']]
x_values_scatters = list(range(0,len(scatter_points)))
plt.scatter(x_values_scatters,scatter_points,c="#0F52BA")
plt.ylabel("Persons (Millions)")
plt.xlabel("Year")
plt.title('Total population in PNG')
plt.savefig('output/figures/Total pop in PNG.png', dpi=500)
plt.show()








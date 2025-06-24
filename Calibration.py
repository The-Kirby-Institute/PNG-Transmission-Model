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


update_config("Routine", "2", "Baseline", "WithDolutegravir", "Continuous","Fitting",60,60)

from model.meta_parameters import *

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pymc as pm
from pymc.ode import DifferentialEquation
import arviz as az


from model.Model_init import *


from scipy.interpolate import UnivariateSpline, CubicSpline

import random

import pytensor
import pytensor.tensor as pt
from scipy.stats import gaussian_kde
from datetime import datetime



### setting up the optimiser for pytensor
### this line of code is used in a Mac. Windows do not need this 
pytensor.config.exception_verbosity = 'high'


np.random.seed(334)

random.seed(123)


starting_epidemic_year  = 1994

reference_year = 2000
adj =  reference_year - starting_epidemic_year


#### Load the priors distribution parameters from the Excel file 
with open('model/prior_path.json') as json_file:
    prior_path = json.load(json_file)

excel_path = prior_path['excel_path']
priors = pd.read_excel(excel_path, sheet_name='Sheet1')

priors = priors.rename(columns={'Unnamed: 0': 'Index'})


def get_bounded_normal(name):
    lower_bound = priors.loc[priors['Index'] == name, 'hdi_3%'].values[0]
    upper_bound = priors.loc[priors['Index'] == name, 'hdi_97%'].values[0]
    mu = priors.loc[priors['Index'] == name, 'mean'].values[0]
    sigma = priors.loc[priors['Index'] == name, 'mcse_mean'].values[0]
    TruncatedNormal = pm.TruncatedNormal
    return TruncatedNormal.dist(mu=mu, sigma=sigma, lower=lower_bound, upper=upper_bound)


def initialize_PMYC_samples(df, variable):
    std = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]
    lower = df.loc[df['Index'] == variable, 'hdi_3%'].values[0]
    upper = df.loc[df['Index'] == variable, 'hdi_97%'].values[0]
    mean = df.loc[df['Index'] == variable, 'mean'].values[0]
    return pm.TruncatedNormal(variable, mu=mean, sigma=std, lower=lower, upper=upper)


def initialize_PMYC_uniform_samples(df, variable):
    lower = df.loc[df['Index'] == variable, 'mean'].values[0]
    upper = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]
    return pm.Uniform(variable, lower=lower, upper=upper)

def get_uniform_distribution(name):
    lower = priors.loc[priors['Index'] == name, 'mean'].values[0]
    upper = priors.loc[priors['Index'] == name, 'mcse_mean'].values[0]    
    return pm.Uniform.dist(lower=lower, upper=upper)



##### Functions used in the Gamma samples distribution 
def initialize_PYMC_gamma_samples(df, variable):
    alpha = df.loc[df['Index'] == variable, 'mean'].values[0]
    beta = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]
    return pm.Gamma(variable, alpha=alpha, beta=beta)

def get_gamma_distribution(name):
    alpha = priors.loc[priors['Index'] == name, 'mean'].values[0]
    beta = priors.loc[priors['Index'] == name, 'mcse_mean'].values[0]    
    return pm.Gamma.dist(alpha=alpha, beta=beta)


def initialize_PYMC_poisson_samples(df, variable):
    lam = df.loc[df['Index'] == variable, 'mean'].values[0]  # Lambda parameter for Poisson
    return pm.Poisson(variable, mu=lam)
def get_poisson_distribution(name):
    lam = priors.loc[priors['Index'] == name, 'mean'].values[0]
    return pm.Poisson.dist(mu=lam)


#### For log-normal distribution 
def initialize_PYMC_lognormal_samples(df, variable):
    mu = df.loc[df['Index'] == variable, 'mean'].values[0]
    sigma = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]
    return pm.Lognormal(variable, mu=mu, sigma=sigma)

def get_lognormal_distribution(name):
    mu = priors.loc[priors['Index'] == name, 'mean'].values[0]
    sigma = priors.loc[priors['Index'] == name, 'mcse_mean'].values[0]    
    return pm.Lognormal.dist(mu=mu, sigma=sigma)


#### For beta distribution 
def initialize_PYMC_beta_samples(df, variable):
    alpha = df.loc[df['Index'] == variable, 'mean'].values[0]
    beta = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]
    return pm.Beta(variable, alpha=alpha, beta=beta)

def get_beta_distribution(name):
    alpha = priors.loc[priors['Index'] == name, 'mean'].values[0]
    beta = priors.loc[priors['Index'] == name, 'mcse_mean'].values[0]    
    return pm.Beta.dist(alpha=alpha, beta=beta)

def setCorrectxAxis(frame, frequency_ticks =5,starting_position=0):
    plt.figure()
    default_x_ticks = list(range(len(frame.I_W)))
    new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+len(frame.I_W) ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, len(frame.I_W)+1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)




###### Initiate variables for Bayesian MCMC analysis
t0=0
##### Year 2025
years = adj+25
tspan = np.arange(0, years, 1)

mcmc_ode = DifferentialEquation(
    func=ode_model,
    times=tspan,
    n_states=29,
    n_theta=30,
    t0=t0,
)
############ 
############ In this block of codes, the relevant modelling and epidemiological data are extracted to 
############ prepare for the Bayesian model calibration.
#### Extract out data from Spectrum, prepare for processing of the Bayesian model data points. 
TotalPLHIV_2023 = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Estimated adults (15+) living with HIV']]
TotalPLHIV_2023['Estimated adults (15+) living with HIV'] = TotalPLHIV_2023['Estimated adults (15+) living with HIV'].astype(str).str.replace(" ", "").astype(int)

onTreatment_2023 = dat_testntreat_2023.loc[dat_testntreat_2023.index[dat_testntreat_2023.index >= starting_epidemic_year],['Missing80']]
onTreatment_2023 = onTreatment_2023.rename(columns={'Missing80': 'Number of PLHIV on ART (Adults, ages 15+)'})

Diagnoses_2023 = dat_testntreat_2023.loc[dat_testntreat_2023.index[dat_testntreat_2023.index >= starting_epidemic_year],['Missing76']]
Diagnoses_2023 = Diagnoses_2023.rename(columns={'Missing76': 'Number of PLHIV knew status (Adults, ages 15+)'})

newinfections_2023 = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Adults (15+) newly infected with HIV']]
newinfections_2023 = newinfections_2023['Adults (15+) newly infected with HIV'].tolist()
#### new infections - selecting from 1995 onwards, because the limitation of PyMC calibration for new infections/ deaths
newinfections_2023 = newinfections_2023[1:]

deaths_2023 = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['AIDS-related deaths among adults (15+)']]
deaths_2023 =  deaths_2023['AIDS-related deaths among adults (15+)'].tolist()
## only pick the last 6 data points of HIV-related deaths for calibration
deaths_2023 = deaths_2023[-6:]

#### Updated UNAIDS data for people living with HIV in PNG
observed_data = dat_estimates_2023.loc[dat_estimates_2023.index[dat_estimates_2023.index >= starting_epidemic_year],['Estimated adults (15+) living with HIV']]
observed_data = observed_data["Estimated adults (15+) living with HIV"].tolist()
observed_data = [item for item in observed_data if isinstance(item,int)] + [int(item.replace(" ",""))  for item in observed_data if isinstance(item,str)]
len_prevalence = len(observed_data)


#### Update data for the viral suppression levels --> 
df1 = dat_testntreat.loc[dat_testntreat.index[dat_testntreat.index >= starting_epidemic_year],['Among people living with HIV, the percent with suppressed viral load']]
df2 = dat_testntreat.loc[dat_testntreat.index[dat_testntreat.index >= starting_epidemic_year],['Clinical Viral Suppresion in PNG']]

df1.columns.values[0] = 'Percent'
df2.columns.values[0] = 'Percent'

combined_vf_df = df1.combine_first(df2)


####### Given the pre-treatment DR prevalence of 12.3% in 2017, we will added the approximate number of pre-treatment Drug resistance 
DRtreatment_naive2017 = 12.3/100*(TotalPLHIV_2023['Estimated adults (15+) living with HIV'][2017] - onTreatment_2023['Number of PLHIV on ART (Adults, ages 15+)'][2017])
observed_data = observed_data +[DRtreatment_naive2017]
####### Given the current level of viral suppression of 86.7% in 2021, we have the additional data point and additional weights 
# VLsuppression_2021 = 86.7/100*onTreatment_2023['Number of PLHIV on ART (Adults, ages 15+)'][2021]
# observed_data = observed_data + [VLsuppression_2021]
#### All viral suppression data from 2018 to 2021 are added here
VLSuppression = combined_vf_df.loc[2018:2021,'Percent'].tolist() 
VLSuppression = [i/100* onTreatment_2023['Number of PLHIV on ART (Adults, ages 15+)'][2021] for i in VLSuppression]
observed_data = observed_data + VLSuppression

####### Given the pre-treatment DR prevalence of 15.7% in 2022, we will added the approximate number of pre-treatment Drug resistance
DRtreatment_naive2022 = 15.7/100*(TotalPLHIV_2023['Estimated adults (15+) living with HIV'][2022] - onTreatment_2023['Number of PLHIV on ART (Adults, ages 15+)'][2022])
observed_data = observed_data +[DRtreatment_naive2022]


###### Add new data from 2018 to 2021 - For Diagnoses info 
All_Diagnoses = Diagnoses_2023[Diagnoses_2023.index >= 2019]['Number of PLHIV knew status (Adults, ages 15+)'].tolist()
observed_data = observed_data + All_Diagnoses
##### Add new data from 2010 to 2021 - For Treatment info
All_Treatments = onTreatment_2023[onTreatment_2023.index >= 2010]['Number of PLHIV on ART (Adults, ages 15+)'].tolist()
observed_data = observed_data + All_Treatments

### adding new infections calibration targets into the modelling 
observed_data = observed_data + newinfections_2023

### adding the mortality targets for calibration in the modelling 
observed_data = observed_data + deaths_2023

#### adding the total population in PNG in 2022 as calibration target
observed_data = observed_data + [corrected_population['Total_Pop'][2022]]

### new prevalence updates, based on the reported prevalence (literature, reports) in PNG
plhiv94_23 = TotalPLHIV_2023['Estimated adults (15+) living with HIV']
prevalence_fix = plhiv94_23 / corrected_population['Total_Pop'][starting_epidemic_year]

#### calculate the number of people living with HIV based on the reported prevalence data, to prepare data for the calibration
report_prevalence = dat_estimates.loc[dat_estimates.index[dat_estimates.index >= starting_epidemic_year],['ACTUAL prevalence of HIV in reports']]['ACTUAL prevalence of HIV in reports']
pop_size = corrected_population.loc[1994:2022, 'Total_Pop']
report_plhiv = report_prevalence.astype(float)/100 * pop_size.astype(float)
report_prev_fix = report_plhiv / corrected_population['Total_Pop'][starting_epidemic_year]

### fill in values that are non-empty from report_prevalence
prevalence_fix_updated = prevalence_fix.combine_first(report_prev_fix)
prevalence_fix_updated = prevalence_fix_updated.tolist()
n = len(prevalence_fix_updated)  # Or choose a smaller number if you only want to replace a part of it


observed_data = [i/corrected_population['Total_Pop'][starting_epidemic_year] for i in observed_data]

### replace the prevalence (fixed), after adding the actual data into it
observed_data[:n] = prevalence_fix_updated[:n]

#################################################################################### 


####### These are the weights of the data points, which are used in the Bayesian model calibration
####### The weights are used to indicate the importance of each data point in the calibration process.
## new weights for July 13th 2023

revised_weights_prevalence = [1 if i < len_prevalence - 6 else  4**(i - len_prevalence + 7) for i in range(len_prevalence)] 
revised_weights_prevalence[adj+2] = 4000
revised_weights_prevalence[-1] = 8800
# ## Modelling weights used in the model 
### weights for viral suppression levels 
VLweights = [1 if i < len(VLSuppression) - 1 else 2500 for i in range(len(VLSuppression) )]
VLweights[0] = 2500
Treatmentweights = [1 if i < len(All_Treatments) - 1 else 2900 for i in range(len(All_Treatments) )] 
Treatmentweights[0] = 2900
newinfectionsweights = [1 if i < len(newinfections_2023) - 1 else 4200 for i in range(len(newinfections_2023) )] 
mortalityweights = [1200]*len(deaths_2023)
mortalityweights[-1] = 4200

weights = revised_weights_prevalence+ [1400] +  VLweights + [1400] +  [1 if i < len(All_Diagnoses) - 1 else 2400 for i in range(len(All_Diagnoses) )]  + Treatmentweights + newinfectionsweights + mortalityweights + [2900]


# Print the current time indicating that the PyMC code is starting
start_time = datetime.now()
print(f"Starting PyMC code execution at {start_time}...")



with pm.Model() as model:
    y = pm.Data("y", observed_data)
   
    theta = initialize_PYMC_lognormal_samples(priors,  'Infectivity_cofficient_DR' )

    c_asterisk = initialize_PYMC_gamma_samples(priors,  'Treatment_rate_b' )
    c_k = initialize_PYMC_gamma_samples(priors,  'Treatment_rate_k' )

        ##### New VN model parameters
    beta_u = initialize_PYMC_gamma_samples(priors,  'Infectivity_infected' )
    beta_t = initialize_PYMC_beta_samples(priors,  'RR_Infection_treated' )
    beta_f = initialize_PYMC_lognormal_samples(priors,  'RR_Infection_fail' )


    delta_U = initialize_PYMC_gamma_samples(priors,  'Mortality_Undiagnosed' )
    delta_T = initialize_PYMC_beta_samples(priors,  'Mortality_Treated' )
    delta_F = initialize_PYMC_beta_samples(priors,  'Mortality_Fail' )
    delta_B = initialize_PYMC_gamma_samples(priors, 'Mortality_Background')
        ##### VN model 
    h1 = initialize_PYMC_gamma_samples(priors,  'Rise_of_DR_1' )
    h2 = initialize_PYMC_beta_samples(priors,  'Rise_of_DR_2' )

    eta_1 = initialize_PYMC_gamma_samples(priors,  'Force_of_Infection_eta_1' )
    eta_2 = initialize_PYMC_gamma_samples(priors,  'Force_of_Infection_eta_2' )
    eta_3 = initialize_PMYC_uniform_samples(priors,  'Force_of_Infection_eta_3' )
    eta_4 = initialize_PYMC_gamma_samples(priors,  'Force_of_Infection_eta_4' )
    eta_5 = initialize_PYMC_gamma_samples(priors,  'Force_of_Infection_eta_5' )
    eta_6 = initialize_PMYC_uniform_samples(priors,  'Force_of_Infection_eta_6' )

    b_asterisk = initialize_PYMC_gamma_samples(priors,  'Diagnosis_rate_b' )
    b_k = initialize_PYMC_gamma_samples(priors,  'Diagnosis_rate_k' )


    rho_asterisk = initialize_PYMC_gamma_samples(priors,  'Rate_Population_Increase' )



    counsel = initialize_PYMC_beta_samples(priors,  'Counselling_effectiveness' )
    transfer_2ndline = initialize_PMYC_samples(priors,  'Transfer_to_second_lineART' )

    f1_rate = initialize_PYMC_beta_samples(priors,  'f1_VirologicalFailure' )
    f2_rate = initialize_PYMC_beta_samples(priors,  'f2_VirologicalFailure' )
    f3_rate = initialize_PYMC_lognormal_samples(priors,  'f3_VirologicalFailure' )
    f4_rate = initialize_PYMC_beta_samples(priors,  'f4_VirologicalFailure' )

    mu1 = initialize_PYMC_beta_samples(priors,  'Rate_LTFU_mu1' )

    qt = initialize_PYMC_lognormal_samples(priors,  'Calibration_parameter_VL_population' )


    phi_inv = pm.Exponential("phi_inv", lam=10)
     
    phi = pm.Deterministic('phi', 1. / phi_inv)



    mcmc_curves = mcmc_ode(y0=initial_conditions, theta=[theta,beta_u,beta_t,beta_f, delta_U, delta_T, delta_F,delta_B, h1, h2, eta_1, eta_2, eta_3, eta_4, eta_5, eta_6, b_asterisk,b_k, c_asterisk,c_k,rho_asterisk,counsel,transfer_2ndline,f1_rate,f2_rate,f3_rate,f4_rate,mu1,qt,year_OffInfs])


    ##### Declare tensor variables, according to correct years of DR prevalence. From beginning start_epideimc_year to 2022
    tensor_a = (mcmc_curves[:,16])[0:(adj+23+1)]
    ##### Levels of pre-treatment DR in 2017
    tensor_b = (mcmc_curves[:,2]+mcmc_curves[:,4])[adj+17]
    ##### Levels of viral suppression from 2018 to 2021
    tensor_c = (mcmc_curves[:,5]+mcmc_curves[:,6]+mcmc_curves[:,7]+mcmc_curves[:,8])[(adj+18):(adj+21+1)]
    ##### Levels of pre-treatment DR in 2022
    tensor_d = (mcmc_curves[:,2]+mcmc_curves[:,4])[adj+22]

    ##### Add data from MCMC model for comparing data points of Diagnoses and treatments
    treatments = mcmc_curves[:,5] + mcmc_curves[:,6] + mcmc_curves[:,7] + mcmc_curves[:,8] + mcmc_curves[:,9] + mcmc_curves[:,10] + mcmc_curves[:,11] + mcmc_curves[:,12] + mcmc_curves[:,13]
    diagnosis = mcmc_curves[:,3] + mcmc_curves[:,4] + treatments

    ##### Declare tensor variables, according to correct years of Diagnoses and treatments. From 2019 for treatments and 2010 for diagnoses - UNAIDS data
    tensor_e = diagnosis[(adj+19):(adj+23+1)]
    tensor_f = treatments[(adj+10):(adj+23+1)]

    ####declare tensor gradient for new infections in the modelling now
    modelled_cum_infs = mcmc_curves[:,14]
    modelled_infections =  modelled_cum_infs[1:] - modelled_cum_infs[:-1]
    ### the calibration target was removed of 1 element at the start so no need to add 1 to adj+23
    tensor_g = modelled_infections[0:(adj+23)]


    ####declare tensor gradient for mortality of AIDS-related deaths in the modelling now 
    modelled_cum_dths = mcmc_curves[:,15]
    modelled_deaths = modelled_cum_dths[1:] - modelled_cum_dths[:-1]
    #### Last 6 years of mortality estimates will be from 2018 (need to subtract 1- because of the tensor issue) up to 2023
    tensor_h = modelled_deaths[(adj+18 - 1): (adj+23)]
    
    ##### Add the tensor for population calibration in the model 
    tensor_i = (mcmc_curves[:,0]+mcmc_curves[:,16])[adj+22]

    tensor_a_flattened = pt.flatten(tensor_a)
    tensor_b_flattened = pt.flatten(tensor_b)
    tensor_c_flattened = pt.flatten(tensor_c)
    tensor_d_flattened = pt.flatten(tensor_d)
    tensor_e_flattened = pt.flatten(tensor_e)
    tensor_f_flattened = pt.flatten(tensor_f)
    tensor_g_flattened = pt.flatten(tensor_g)
    tensor_h_flattened = pt.flatten(tensor_h)
    tensor_i_flattened = pt.flatten(tensor_i)

    concatenated_tensor = pt.concatenate([tensor_a_flattened, tensor_b_flattened,tensor_c_flattened,tensor_d_flattened,tensor_e_flattened,tensor_f_flattened,tensor_g_flattened,tensor_h_flattened,tensor_i_flattened])

    ### This is for weighted model, DR point = 100
    # logp = pm.NegativeBinomial.dist(mu=concatenated_tensor, alpha=phi).logp(y)
    rv = pm.NegativeBinomial.dist(mu=concatenated_tensor, alpha=phi)
    logp = pm.logp(rv, y)
    weighted_logp = logp * weights
    potential = pm.Potential("weighted_logp", weighted_logp.sum())


with model:   
    #### multicore processing  
    trace = pm.sample(draws=250, tune=250, step=pm.NUTS(), chains=4, cores=4)
    # trace = pm.sample(draws=250, tune=250, step=pm.NUTS(), chains=2, cores=4, progressbar=True)

# # Plot the trace plots
# az.plot_trace(trace)
# plt.savefig('output/figures/Trace_plot.png', dpi=500)
# plt.show()


results = az.summary(trace)
results.to_excel("output/diagnostics/Bayesian_parameter(results).xlsx")




posterior_samples = {
    'theta_samples': trace.posterior['Infectivity_cofficient_DR'].values.flatten(), 
    'beta_u_samples': trace.posterior['Infectivity_infected'].values.flatten(),
     'beta_t_samples' : trace.posterior['RR_Infection_treated'].values.flatten(),
    'beta_f_samples':  trace.posterior['RR_Infection_fail'].values.flatten(), 
    'delta_U_samples': trace.posterior['Mortality_Undiagnosed'].values.flatten(), 
    'delta_T_samples': trace.posterior['Mortality_Treated'].values.flatten(), 
    'delta_F_samples': trace.posterior['Mortality_Fail'].values.flatten(), 
    'delta_B_samples': trace.posterior['Mortality_Background'].values.flatten(),
    'h1_samples': trace.posterior['Rise_of_DR_1'].values.flatten(), 
    'h2_samples': trace.posterior['Rise_of_DR_2'].values.flatten() ,
    'eta_1_samples':trace.posterior['Force_of_Infection_eta_1'].values.flatten(), 
    'eta_2_samples': trace.posterior['Force_of_Infection_eta_2'].values.flatten(),
    'eta_3_samples':trace.posterior['Force_of_Infection_eta_3'].values.flatten(),
    'eta_4_samples': trace.posterior['Force_of_Infection_eta_4'].values.flatten(),
    'eta_5_samples': trace.posterior['Force_of_Infection_eta_5'].values.flatten(), 
    'eta_6_samples': trace.posterior['Force_of_Infection_eta_6'].values.flatten() , 
    'b_asterisk_samples': trace.posterior['Diagnosis_rate_b'].values.flatten(),
    'b_k_samples': trace.posterior['Diagnosis_rate_k'].values.flatten(), 
    'c_asterisk_samples': trace.posterior['Treatment_rate_b'].values.flatten(),
    'c_k_samples': trace.posterior['Treatment_rate_k'].values.flatten(),
    'rho_asterisk_samples': trace.posterior['Rate_Population_Increase'].values.flatten(), 
    'counsel_samples': trace.posterior['Counselling_effectiveness'].values.flatten(), 
    'transfer_2ndline':trace.posterior['Transfer_to_second_lineART'].values.flatten(),
    'f1_rate':trace.posterior['f1_VirologicalFailure'].values.flatten(),
    'f2_rate':trace.posterior['f2_VirologicalFailure'].values.flatten(),
    'f3_rate':trace.posterior['f3_VirologicalFailure'].values.flatten(),
    'f4_rate':trace.posterior['f4_VirologicalFailure'].values.flatten(), 
    'mu1':trace.posterior['Rate_LTFU_mu1'].values.flatten(),
    'qt':trace.posterior['Calibration_parameter_VL_population'].values.flatten()}


df = pd.DataFrame(posterior_samples)

with pd.ExcelWriter('data/posteriors/Bayesian_posterior-recalibrate.xlsx') as writer:
    df.to_excel(writer, sheet_name='Sheet1',index="False")
    
    

az.plot_forest(trace.posterior, r_hat=True)

# plot the posterior distributions of the model parameters using Arviz:
az.plot_posterior(trace.posterior)
plt.savefig('output/diagnostics/Posterior_plot.png', dpi=500)

plt.show()






original_var_names = [name for name in trace.posterior.data_vars if not name.endswith('__')]
modified_var_names = original_var_names[:-2]




# Print the current time indicating that the PyMC code has finished
end_time = datetime.now()
print(f"PyMC code execution finished at {end_time}.")
print(f"Total execution time: {end_time - start_time}")



#### Reporting status for divergence analyses 

divergences = trace.sample_stats['diverging'].values.nonzero()[0]
print(f"Number of Divergences: {len(divergences)}")

#Analyzing Divergences
pm.plot_pair(trace, divergences=True)
plt.savefig('output/diagnostics/Pair Plot for Divergence.png')
plt.show()


az.plot_trace(trace, divergences="top")
plt.savefig('output/diagnostics/Trace plot (top) - divergence.png')
plt.show()

az.plot_energy(trace)
plt.savefig('output/diagnostics/Energy transition diagram.png')
plt.show()

# Calculate ESS
ess = pm.summary(trace)  # This gives ESS for the mean of each parameter
print(ess.ess_bulk)
print(ess.tail)
# Ideal Value: An R-hat value close to 1.0 indicates that the variance within each chain is about the same as the variance across the chains, suggesting convergence. Typically, values of R-hat less than 1.1 (or sometimes 1.05 in more stringent cases) are considered indicative of good convergence.
#https://mc-stan.org/rstan/reference/Rhat.html
# Gelman, A. and D. B. Rubin (1992) Inference from iterative simulation using multiple sequences (with discussion). Statistical Science, 7:457-511.
# Gelman, A., Carlin, J. B., Stern, H. S., and D. B. Rubin (2003) Bayesian Data Analysis, 2nd edition. Chapman and Hall/CRC.
# Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and Paul-Christian Bürkner (2019). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
print(ess.r_hat.values)
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 12:08:15 2023

@author: nguye
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import gamma, lognorm, uniform
from scipy.stats import beta as beta_dist
random_seed = 123
from scipy.stats import norm, truncnorm 


no_of_samples = 400

import random

np.random.seed(334)
random.seed(334)

#### Load the priors distribution parameters from the Excel file 
with open('model/prior_path.json') as json_file:
    prior_path = json.load(json_file)

excel_path = prior_path['excel_path']
df = pd.read_excel(excel_path, sheet_name='Sheet1')

# ##### Read the summary pandas dataframessss
df = df.rename(columns={'Unnamed: 0': 'Index'})

df['Index'] = df['Index'].replace({
    'Infectivity_cofficient_DR': 'theta_samples',
    'Infectivity_infected': 'beta_u_samples',
    'RR_Infection_treated': 'beta_t_samples',
    'RR_Infection_fail': 'beta_f_samples',
    'Mortality_Undiagnosed': 'delta_U_samples',
    'Mortality_Treated': 'delta_T_samples',
    'Mortality_Fail': 'delta_F_samples',
    'Mortality_Background': 'delta_B_samples',

    'Rise_of_DR_1': 'h1_samples',
    'Rise_of_DR_2': 'h2_samples',
    'Force_of_Infection_eta_1': 'eta_1_samples',
    'Force_of_Infection_eta_2': 'eta_2_samples',
    'Force_of_Infection_eta_3': 'eta_3_samples',
    'Force_of_Infection_eta_4': 'eta_4_samples',
    'Force_of_Infection_eta_5': 'eta_5_samples',
    'Force_of_Infection_eta_6': 'eta_6_samples',
    'Diagnosis_rate_b': 'b_asterisk_samples',
    'Diagnosis_rate_k': 'b_k_samples',
    'Diagnosis_rate_q': 'b_q_samples',
    'Diagnosis_rate_v': 'b_v_samples',

    'Treatment_rate_b': 'c_asterisk_samples',
    'Treatment_rate_k': 'c_k_samples',
    'Initial_HIV_cases': 'epsilon_samples',
    'Treatment_rate_v': 'c_v_samples',

    'RR_FailureVirl_churn': 'f_asterisk_samples',
    'Rate_Population_Increase': 'rho_asterisk_samples',
    'Counselling_effectiveness': 'counsel_samples',
    'Transfer_to_second_lineART': 'transfer_2ndline_samples',
    'Rate_LTFU_mu1': 'mu1_samples',
    'Rate_LTFU_mu2': 'mu2_samples',
    'Rate_churn_mu3':'mu3_samples',

    'LTFU_1st_line': 'first_LTFU_samples',
    'LTFU_2nd_line': 'second_LTFU_samples',

    'f1_VirologicalFailure': 'f1_rate_samples',
    'f2_VirologicalFailure': 'f2_rate_samples',
    'f3_VirologicalFailure': 'f3_rate_samples',
    'f4_VirologicalFailure': 'f4_rate_samples',

    'Calibration_parameter_VL_population': 'qt_samples',

    'ACTUP_VL_Effectiveness': 'ACTUP_VL_samples'
})


#### Function to generate random samples from the posterior distribution
def acquire_basevalues(df,variable):
    mean = df.loc[df['Index'] == variable, 'mean'].values[0]
    return mean

def acquire_gamma(df, variable, no_of_samples=1, random_seed=None):
    alpha = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha (shape parameter)
    beta = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta (rate parameter)
    theta = 1 / beta  # converting rate to scale

    return gamma.rvs(alpha, scale=theta, size=no_of_samples, random_state=random_seed)[0]


def acquire_lognormal(df, variable, no_of_samples=1, random_seed=None):
    mu = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha (shape parameter)
    sigma = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta (rate parameter)

    return lognorm.rvs(s=sigma, scale=np.exp(mu), size=no_of_samples, random_state=random_seed)[0]


def acquire_beta(df, variable, no_of_samples=1, random_seed=None):
    alpha = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha (shape parameter)
    beta = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta (rate parameter)

    return beta_dist.rvs(a =alpha, b=beta, size=no_of_samples, random_state=random_seed)[0]


def acquire_uniform(df, variable, no_of_samples=1, random_seed=None):
    lower = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha (shape parameter)
    upper = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta (rate parameter)

    return uniform.rvs(loc =lower, scale=upper - lower, size=no_of_samples, random_state=random_seed)[0]

# Generate random samples of beta and gamma
theta = acquire_lognormal(df,'theta_samples')
######
beta_u = acquire_gamma(df,'beta_u_samples')
beta_t= acquire_beta(df,'beta_t_samples')
beta_f = acquire_lognormal(df,'beta_f_samples')


delta_U = acquire_gamma(df,'delta_U_samples')
delta_T = acquire_beta(df,'delta_T_samples')
delta_F = acquire_beta(df,'delta_F_samples')
delta_B = acquire_gamma(df,'delta_B_samples')




h1 = acquire_gamma(df,'h1_samples')
h2 = acquire_beta(df,'h2_samples')




##### All time-varying parameters to be generated 
# eta_2_samples =  acquire_basevalues(df,'eta_2_samples',overwrite=1.09)
eta_1 =  acquire_gamma(df,'eta_1_samples')
######
eta_2 =  acquire_gamma(df,'eta_2_samples')
######
eta_3 = acquire_uniform(df,'eta_3_samples')
######
eta_4 =  acquire_gamma(df,'eta_4_samples')
######
eta_5 =  acquire_gamma(df,'eta_5_samples')
######
eta_6 = acquire_uniform(df,'eta_6_samples')
######
b_asterisk =  acquire_gamma(df,'b_asterisk_samples')
######
b_k =  acquire_gamma(df,'b_k_samples')
######
b_q =  acquire_uniform(df,'b_q_samples')
######
b_v =  acquire_basevalues(df,'b_v_samples')
######

c_asterisk =  acquire_gamma(df,'c_asterisk_samples')
######
c_k =  acquire_gamma(df,'c_k_samples')
######
epsilon =  acquire_basevalues(df,'epsilon_samples')
######
c_v =  acquire_uniform(df,'c_v_samples')
######
f_asterisk =  acquire_basevalues(df,'f_asterisk_samples')
######
rho_asterisk =  acquire_gamma(df,'rho_asterisk_samples')
######
counsel = acquire_beta(df,'counsel_samples')
######
transfer_2ndline= acquire_basevalues(df,'transfer_2ndline_samples')
######
first_LTFU = acquire_basevalues(df,'first_LTFU_samples')
######
second_LTFU= acquire_basevalues(df,'second_LTFU_samples')
######
mu1 = acquire_beta(df,'mu1_samples')
######
mu2 = acquire_basevalues(df,'mu2_samples')
###### 
mu3 = acquire_basevalues(df,'mu3_samples')
######
f1 = acquire_beta(df,'f1_rate_samples')
######
f2 = acquire_beta(df,'f2_rate_samples')
######
f3 = acquire_lognormal(df,'f3_rate_samples')
######
f4 = acquire_beta(df,'f4_rate_samples')

##### Calibration parameter for viral load levels in population
qt = acquire_lognormal(df,'qt_samples')


#### Effectiveness of ACTUP study in Viral suppression (theoretically)
actup = acquire_basevalues(df,'ACTUP_VL_samples')
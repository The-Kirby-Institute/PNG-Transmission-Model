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




import pandas as pd
from model.core_scenarios import *
import os 

# Get the path to the main directory
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Read in the Excel file
if running_mode == "Study_3":
    file_path = os.path.join(base_dir, 'data', 'posteriors', 'Bayesian_posterior_Chap3.xlsx')
else:
    file_path = os.path.join(base_dir, 'data', 'posteriors', 'Bayesian_posterior.xlsx')

# Read the Excel file
df = pd.read_excel(file_path, sheet_name='Sheet1')

theta_samples = df['theta_samples']
beta_u_samples =df['beta_u_samples']
beta_t_samples = df['beta_t_samples']
beta_f_samples = df['beta_f_samples']


delta_U_samples = df['delta_U_samples']
delta_T_samples = df['delta_T_samples']
delta_F_samples = df['delta_F_samples']
delta_B_samples = df['delta_B_samples']


h1_samples = df['h1_samples']
h2_samples = df['h2_samples']

eta_1_samples = df['eta_1_samples']
eta_2_samples =  df['eta_2_samples']
eta_3_samples = df['eta_3_samples']
eta_4_samples = df['eta_4_samples']
eta_5_samples = df['eta_5_samples']
eta_6_samples = df['eta_6_samples']
b_asterisk_samples =  df['b_asterisk_samples']
b_k_samples = df['b_k_samples']

c_asterisk_samples =  df['c_asterisk_samples']
c_k_samples =  df['c_k_samples']

rho_asterisk_samples =  df['rho_asterisk_samples']
counsel_samples =  df['counsel_samples']
transfer_2ndline_samples = df['transfer_2ndline']
#### New variables to include in July 11th version
f1_rate_samples = df['f1_rate']
f2_rate_samples = df['f2_rate']
f3_rate_samples = df['f3_rate']
f4_rate_samples = df['f4_rate']
mu1_samples = df['mu1']
qt_samples = df['qt']
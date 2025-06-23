import pandas as pd
from model.core_scenarios import *
# Read in the Excel file
if running_mode == "Study_3":
    df = pd.read_excel('data/posteriors/Bayesian_posterior_Chap3.xlsx', sheet_name='Sheet1')
else:
    df = pd.read_excel('data/posteriors/Bayesian_posterior.xlsx', sheet_name='Sheet1')

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
# b_q_samples =  df['b_q_samples']
# b_v_samples =  df['b_v_samples']

c_asterisk_samples =  df['c_asterisk_samples']
c_k_samples =  df['c_k_samples']
# epsilon_samples = df['epsilon_samples']
# c_v_samples =  df['c_v_samples']

# f_asterisk_samples =  df['f_asterisk_samples']
rho_asterisk_samples =  df['rho_asterisk_samples']
counsel_samples =  df['counsel_samples']
transfer_2ndline_samples = df['transfer_2ndline']
#### New variables to include in July 11th version
# first_LTFU_samples = df['first_LTFU']
# second_LTFU_samples = df['second_LTFU']
f1_rate_samples = df['f1_rate']
f2_rate_samples = df['f2_rate']
f3_rate_samples = df['f3_rate']
f4_rate_samples = df['f4_rate']
mu1_samples = df['mu1']
# mu2_samples = df['mu2']
# mu3_samples = df['mu3']
qt_samples = df['qt']
### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

from model.read_posterior_samples import*
import json
import pandas as pd
import numpy as np
import scipy.stats as stats

import pandas as pd
import scipy.stats as stats
import os
import matplotlib.pyplot as plt
import seaborn as sns

### new functions to generate distributions, instead of values
from model.acquire_distribution import *

from scipy.stats import gaussian_kde

np.random.seed(123)


# Function to calculate and normalize KDE
def plot_normalized_kde(samples, label, color,linestyle="-"):
    # Calculate the KDE
    kde = gaussian_kde(samples)
    x_range = np.linspace(min(samples), max(samples), 1000)
    density = kde(x_range)
    
    # Normalize the density values
    density /= density.max()
    
    # Plot the normalized KDE
    plt.plot(x_range, density, label=label, color=color, linestyle=linestyle)

random_seed = 123
number_of_samples = 5000
result_path = os.path.join(base_path, 'output', 'Chp2_tables')




def calculate_stats(samples, as_percentage=False):
    lower_95_ci = np.percentile(samples, 2.5)
    upper_95_ci = np.percentile(samples, 97.5)
    median = np.median(samples)
    
    if as_percentage:
        lower_95_ci *= 100
        upper_95_ci *= 100
        median *= 100
        return f"{median:.1f}% ({lower_95_ci:.1f}% - {upper_95_ci:.1f}%)"
    else:
        return f"{median:.2f} ({lower_95_ci:.2f} - {upper_95_ci:.2f})"

def calculate_range(samples, as_percentage=False):
    lower_95_ci = np.percentile(samples, 2.5)
    upper_95_ci = np.percentile(samples, 97.5)
    
    if as_percentage:
        lower_95_ci *= 100
        upper_95_ci *= 100
        return f"{lower_95_ci:.1f}% to {upper_95_ci:.1f}%"
    else:
        return f"{lower_95_ci:.2f} to {upper_95_ci:.2f}"

# Load data
# Build full path to the Excel file
posterior_path = os.path.join(base_path, 'data', 'posteriors', 'Bayesian_posterior.xlsx')

df = pd.read_excel(posterior_path, sheet_name='Sheet1')

# List of sample columns to process and their corresponding percentage reporting options
sample_columns = {
    'beta_u_samples': True,
    'beta_t_samples': False,
    'beta_f_samples': False,
    'eta_1_samples': False,
    'rho_asterisk_samples':True,
    'delta_U_samples':True,
    'delta_T_samples':False,
    'delta_F_samples':False,
    'delta_B_samples':True,
    'theta_samples':True,
    'b_k_samples':True,
    'c_k_samples':True,
    'h1_samples':True,
    'h2_samples':True,
    'f1_rate':True,
    'f2_rate':True,
    'f3_rate':False,
    'f4_rate':True,
    'counsel_samples':True,
    'mu1':True,
    'qt':False, 
    'eta_2_samples':False, 
    'eta_5_samples':False,
    'b_asterisk_samples':False,
    'c_asterisk_samples':False,

    'eta_3_samples':False,
    'eta_6_samples':False,
    'eta_4_samples':False
    }

# Calculate stats for each sample column based on the specified percentage reporting option
summary_stats = {col: calculate_stats(df[col], as_percentage=as_percentage) for col, as_percentage in sample_columns.items()}

# Create a summary DataFrame with a single column
summary_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['Summary'])

# Print the summary DataFrame
# print(summary_df)





# with open('model/prior_path.json') as json_file:
#     prior_path = json.load(json_file)

# Build the path to model/prior_path.json
prior_path_file = os.path.join(base_path, 'model', 'prior_path.json')

# Open and load the file
with open(prior_path_file, 'r') as json_file:
    prior_path = json.load(json_file)

# Combine base path with the relative path from JSON
excel_rel_path = prior_path['excel_path']
excel_path = os.path.join(base_path, excel_rel_path)
distribution = pd.read_excel(excel_path, sheet_name='Sheet1')
distribution = distribution.rename(columns={'Unnamed: 0': 'Index'})


distribution['Index'] = distribution['Index'].replace({
    'Infectivity_cofficient_DR': 'theta_samples',
    'Infectivity_infected': 'beta_u_samples',
    'RR_Infection_treated': 'beta_t_samples',
    'RR_Infection_fail': 'beta_f_samples',
    'Mortality_Undiagnosed': 'delta_U_samples',
    'Mortality_Treated': 'delta_T_samples',
    'Mortality_Fail': 'delta_F_samples',
    'Mortality_Background':'delta_B_samples',

    'Rise_of_DR_1': 'h1_samples',
    'Rise_of_DR_2': 'h2_samples',
    'Force_of_Infection_eta_1': 'eta_1_samples',
    'Force_of_Infection_eta_2': 'eta_2_samples',
    'Force_of_Infection_eta_5':'eta_5_samples',
    'Diagnosis_rate_b': 'b_asterisk_samples',
    'Diagnosis_rate_k': 'b_k_samples',

    'Treatment_rate_b': 'c_asterisk_samples',
    'Treatment_rate_k': 'c_k_samples',
    'Treatment_rate_v': 'c_v_samples',

    'Rate_Population_Increase': 'rho_asterisk_samples',
    'Counselling_effectiveness': 'counsel_samples',
    'Rate_LTFU_mu1': 'mu1',


    'f1_VirologicalFailure': 'f1_rate',
    'f2_VirologicalFailure': 'f2_rate',
    'f3_VirologicalFailure': 'f3_rate',
    'f4_VirologicalFailure': 'f4_rate',

    'Calibration_parameter_VL_population': 'qt',

    'Force_of_Infection_eta_3': 'eta_3_samples',
    'Force_of_Infection_eta_6':'eta_6_samples',
    'Force_of_Infection_eta_4': 'eta_4_samples'
})


# distribution = pd.read_excel('Bayesian Pictures/Bayesian_posterior.xlsx', sheet_name='Sheet1')

# List of variables with distribution type and percentage reporting options
variables = {
    'beta_u_samples': {'distribution': 'gamma', 'as_percentage': True},
    'beta_t_samples': {'distribution': 'beta', 'as_percentage': False},
    'beta_f_samples': {'distribution': 'lognormal', 'as_percentage': False},
    'eta_1_samples': {'distribution': 'gamma', 'as_percentage': False},
    'rho_asterisk_samples': {'distribution': 'gamma', 'as_percentage': True},
    'delta_B_samples': {'distribution': 'gamma', 'as_percentage': True},
    'delta_U_samples': {'distribution': 'gamma', 'as_percentage': True},
    'delta_T_samples': {'distribution': 'beta', 'as_percentage': False},
    'delta_F_samples': {'distribution': 'beta', 'as_percentage': False},
    'theta_samples': {'distribution': 'lognormal', 'as_percentage': True},
    'b_k_samples': {'distribution': 'gamma', 'as_percentage': True},
    'c_k_samples': {'distribution': 'gamma', 'as_percentage': True},
    'h1_samples': {'distribution': 'gamma', 'as_percentage': True},
    'h2_samples': {'distribution': 'beta', 'as_percentage': True},
    'f1_rate': {'distribution': 'beta', 'as_percentage': True},
    'f2_rate': {'distribution': 'beta', 'as_percentage': True},
    'f3_rate': {'distribution': 'lognormal', 'as_percentage': False},
    'f4_rate': {'distribution': 'beta', 'as_percentage': True},
    'counsel_samples': {'distribution': 'beta', 'as_percentage': True},
    'mu1': {'distribution': 'beta', 'as_percentage': True},
    'qt': {'distribution': 'lognormal', 'as_percentage': False},
    'eta_2_samples': {'distribution': 'gamma', 'as_percentage': False},
    'eta_5_samples': {'distribution': 'gamma', 'as_percentage': False},
    'b_asterisk_samples':{'distribution': 'gamma', 'as_percentage': False},
    'c_asterisk_samples':{'distribution': 'gamma', 'as_percentage': False},

    'eta_3_samples': {'distribution': 'uniform', 'as_percentage': False},
    'eta_6_samples': {'distribution': 'uniform', 'as_percentage': False},
    'eta_4_samples':{'distribution': 'gamma', 'as_percentage': False}
}

# Generate samples, calculate stats, and store results in a dictionary
summary_stats = {}
for variable, options in variables.items():
    dist_type = options['distribution']
    as_percentage = options['as_percentage']
    
    if dist_type == 'gamma':
        samples = acquire_gamma(distribution, variable, no_of_samples=number_of_samples)
    elif dist_type == 'beta':
        samples = acquire_beta(distribution, variable, no_of_samples=number_of_samples)
    elif dist_type == 'lognormal':
        samples = acquire_lognormal(distribution, variable, no_of_samples=number_of_samples)
    elif dist_type == 'uniform':
        samples = acquire_uniform(distribution, variable, no_of_samples=number_of_samples)
    else:
        raise ValueError(f"Unsupported distribution type: {dist_type}")
    
    summary_stats[variable] = calculate_range(samples, as_percentage=as_percentage)

# Create a summary DataFrame with a single column
summary_df2 = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['Summary'])

# Print the summary DataFrame
# print(summary_df2)

combined_summary_df = pd.concat([ summary_df2, summary_df], axis=1)


# output_file_path = os.path.join(result_path, "Part of Table 1-Manuscript Processed.xlsx")
# combined_summary_df.to_excel(output_file_path, index=True)
# combined_summary_df




summary_df.index = summary_df.index.str.strip()
summary_df2.index = summary_df2.index.str.strip()

label_list = {
    'beta_u_samples': {'legend': 'Beta U Samples', 'title': 'Transmission rate of infections from untreated HIV'},
    'beta_t_samples': {'legend': 'Beta T Samples', 'title': 'RR of Transmission of Unsuppressed VL vs those Untreated Individuals'},
    'beta_f_samples': {'legend': 'Beta F Samples', 'title': 'RR of HIV Transmission of Unsuppressed VL vs those Untreated Individuals'},
    'eta_1_samples': {'legend': 'Beta A Samples', 'title': 'Multiplier for the effect of early/acute infections on the early phase of HIV epidemic in PNG'},
    'rho_asterisk_samples': {'legend': 'Rho Samples', 'title': 'Net growth rate of adult population (>15 years old) in PNG'},
    'delta_B_samples': {'legend': 'Delta B Samples', 'title': 'Background mortality rate'},
    'delta_U_samples': {'legend': 'Delta U Samples', 'title': 'Annual mortality rate attributable to HIV/AIDS for those untreated PLHIV'},
    'delta_T_samples': {'legend': 'Delta T Samples', 'title': 'RR of Mortality Rate for Treated, suppressed VL vs. Undiagnosed HIV'},
    'delta_F_samples': {'legend': 'Delta F Samples', 'title': 'Relative Risk of Mortality Rate for Treated, unsuppressed VL vs. Undiagnosed HIV'},
    'theta_samples': {'legend': 'Theta Samples', 'title': 'Fitness of drug-resistant HIV viruses compared to wild-type HIV'},
    'b_k_samples': {'legend': 'B K Samples', 'title': 'Average rate of diagnoses of HIV for those infected and undiagnosed'},
    'c_k_samples': {'legend': 'C K Samples', 'title': 'Rate of treatment initiation of those who are diagnosed (%)'},
    'h1_samples': {'legend': 'H1 Samples', 'title': 'Rate of emergence of drug-resistant mutations for those on NNRTIs'},
    'h2_samples': {'legend': 'H2 Samples', 'title': 'Relative risk of resistance emergence when on dolutegravir compared to on NNRTIs'},
    'f1_rate': {'legend': 'F1 Rate', 'title': 'Rate of virological failure of the old ARTs first-line (NNRTIs) (year-1)'},
    'f2_rate': {'legend': 'F2 Rate', 'title': 'Rate of virological failure of the new ART first-line (dolutegravir)'},
    'f3_rate': {'legend': 'F3 Rate', 'title': 'Relative risk virological failure of dolutegravir \n with background drug resistance compared to failure rate in wild-type HIV'},
    'f4_rate': {'legend': 'F4 Rate', 'title': 'Rate of virological failure in second-line ART'},
    'counsel_samples': {'legend': 'Counsel Samples', 'title': 'Effectiveness of counselling following the result of viral load test'},
    'mu1': {'legend': 'Mu', 'title': 'Proportion lost-to-follow-up (LTFU)'},
    'qt': {'legend': 'lambda', 'title': 'Resuppression rate of centralised laboratory VL testing, if counselling is 100% effective'},
    'eta_2_samples': {'legend': 'beta_b', 'title': 'Rate of change in the time-varying effect of early/acute infections'},
    'eta_5_samples': {'legend': 'beta_d', 'title': 'Rate of change in the time-varying effect of later infections'},
    'b_asterisk_samples': {'legend': 'b_b', 'title': 'Rate of change in rate of diagnoses'},
    'c_asterisk_samples': {'legend': 'c_asterisk', 'title': 'Rate of change in rate of treatments'},

    'eta_3_samples': {'legend': 'beta_c', 'title': 'The inflection point in the first logistic curve'},
    'eta_6_samples': {'legend': 'beta_f', 'title': 'The inflection point in the second logistic curve'},
    'eta_4_samples': {'legend': 'eta_4', 'title': 'The right asymptope in the second logistic curve'}

}

# Generate density plots
for variable in variables.keys():
    plt.figure(figsize=(10, 6))
    
    sns.kdeplot(df[variable], label=f'{label_list[variable]["legend"]} Posterior Distribution', color='red')
    # plot_normalized_kde(df[variable],label=f'{label_list[variable]["legend"]} Prior Distribution', color='red')

    dist_type = variables[variable]['distribution']
    if dist_type == 'gamma':
        samples = acquire_gamma(distribution, variable, no_of_samples=number_of_samples)
    elif dist_type == 'beta':
        samples = acquire_beta(distribution, variable, no_of_samples=number_of_samples)
    elif dist_type == 'lognormal':
        samples = acquire_lognormal(distribution, variable, no_of_samples=number_of_samples)
    elif dist_type == 'uniform':
        samples = acquire_uniform(distribution, variable, no_of_samples=number_of_samples)
    else:
        raise ValueError(f"Unsupported distribution type: {dist_type}")
    
    sns.kdeplot(samples, label=f'{label_list[variable]["legend"]} Prior Distribution', color='#0F52BA',linestyle ="--")
    # plot_normalized_kde(samples, label=f'{label_list[variable]["legend"]} Posterior Distribution', color='#0F52BA',linestyle ="--")

    plt.title(f'Density Plot for {label_list[variable]["title"]}')
    plt.xlabel(label_list[variable]['legend'])
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True)
    
    # Define full output path
    plot_dir = os.path.join(base_path, 'output', 'priors-posteriors')
    os.makedirs(plot_dir, exist_ok=True)  

    # Build full file path
    plot_path = os.path.join(plot_dir, f'{label_list[variable]["legend"]}_density_plot.png')
    plt.savefig(plot_path, dpi=500)
    # plt.savefig(f'Bayesian Predictive Pictures/Posterior Graphs/{label_list[variable]["legend"]}_density_plot.png', dpi=500)

    plt.show()


### calculate the rates to fill in the Table 
### Prior distribution
fixed_f3_rate = acquire_lognormal(distribution, "f3_rate", no_of_samples=number_of_samples)
base_rate = acquire_beta(distribution, "f1_rate", no_of_samples=number_of_samples)
fixed_f3_rate = fixed_f3_rate * base_rate
### Posterior distribution
posterior_f3 = df['f3_rate']*df['f1_rate']

### plot the new graph 
plt.figure(figsize=(10, 6))
sns.kdeplot(posterior_f3, label=f'F3 rate Posterior Distribution', color='red')
sns.kdeplot(fixed_f3_rate, label=f'F3 rate Prior Distribution', color='#0F52BA',linestyle ="--")
plt.title(f'Density Plot for Relative risk virological failure of dolutegravir \n with background drug resistance compared to failure rate in wild-type HIV')
plt.xlabel('F3 Rate')
plt.ylabel('Density')
plt.legend()
plt.grid(True)
# plt.savefig(f'Bayesian Predictive Pictures/Posterior Graphs/f3_rate_density_plot_fixed.png', dpi=500)
plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'output', 'priors-posteriors', 'f3_rate_density_plot_fixed.png'), dpi=500)


fixed_h2_samples = acquire_beta(distribution, "h2_samples", no_of_samples=number_of_samples)
base_rate = acquire_gamma(distribution, "h1_samples", no_of_samples=number_of_samples)
fixed_h2_samples = fixed_h2_samples * base_rate
### Posterior distribution
posterior_h2 = df['h2_samples']*df['h1_samples'] 

### plot the new graph 
plt.figure(figsize=(10, 6))
sns.kdeplot(posterior_h2, label=f'h2 Posterior Distribution', color='red')
sns.kdeplot(fixed_h2_samples, label=f'h2 Prior Distribution', color='#0F52BA',linestyle ="--")
plt.title(f'Density Plot for Relative risk of resistance emergence when on dolutegravir compared to on NNRTIs')
plt.xlabel('h2 Rate')
plt.ylabel('Density')
plt.legend()
plt.grid(True)
# plt.savefig(f'Bayesian Predictive Pictures/Posterior Graphs/h2_rate_density_plot_fixed.png', dpi=500)
plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'output', 'priors-posteriors', 'h2_rate_density_plot_fixed.png'), dpi=500)

combined_summary_df.columns = ['Range', 'Dist']
combined_summary_df.loc['f3_rate'] = [calculate_range(fixed_f3_rate, as_percentage=True),
                                      calculate_stats(posterior_f3, as_percentage=True)]
combined_summary_df.loc['h2_samples'] = [calculate_range(fixed_h2_samples, as_percentage=True),
                                      calculate_stats(posterior_h2, as_percentage=True)]

### remove eta1 from Table 1 - there is no direct link
combined_summary_df = combined_summary_df.drop('eta_1_samples')

output_file_path = os.path.join(result_path, "Table 1-Manuscript FINAL.xlsx")
combined_summary_df.to_excel(output_file_path, index=True)
combined_summary_df

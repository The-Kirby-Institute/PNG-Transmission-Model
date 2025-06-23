from scipy.stats import gamma, lognorm, uniform
from scipy.stats import beta as beta_dist
import numpy as np

number_of_samples = 5000

def acquire_gamma(df, variable, no_of_samples=number_of_samples, random_seed=None):
    alpha = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha (shape parameter)
    beta_value = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta (rate parameter)
    theta = 1 / beta_value  # converting rate to scale
    return gamma.rvs(alpha, scale=theta, size=no_of_samples, random_state=random_seed)

def acquire_beta(df, variable, no_of_samples=number_of_samples, random_seed=None):
    alpha = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha parameter
    beta = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta parameter
    return beta_dist.rvs(alpha, beta, size=no_of_samples, random_state=random_seed)

def acquire_lognormal(df, variable, no_of_samples=number_of_samples, random_seed=None):
    mu = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha (shape parameter)
    sigma = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta (rate parameter)

    return lognorm.rvs(s=sigma, scale=np.exp(mu), size=no_of_samples, random_state=random_seed)

def acquire_uniform(df, variable, no_of_samples=number_of_samples, random_seed=None):
    lower = df.loc[df['Index'] == variable, 'mean'].values[0]  # alpha (shape parameter)
    upper = df.loc[df['Index'] == variable, 'mcse_mean'].values[0]  # beta (rate parameter)

    return uniform.rvs(loc =lower, scale=upper - lower, size=no_of_samples, random_state=random_seed)

def acquire_basevalues(df,variable):
    mean = df.loc[df['Index'] == variable, 'mean'].values[0]
    return mean

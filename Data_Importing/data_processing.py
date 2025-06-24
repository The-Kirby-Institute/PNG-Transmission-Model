#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:39:28 2022

@author: ducquangnguyen
"""


import pandas as pd
import math
import numpy as np
import re
import os 
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

pd.options.display.float_format = '{:.0f}'.format

# xls = pd.ExcelFile('Data/All_clinics_data.xlsx')
xls = pd.ExcelFile(os.path.join(base_path, 'Data', 'All_clinics_data.xlsx'))



dat_estimates = pd.read_excel(xls, 'HIV Estimates (records)')
dat_estimates.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900},inplace=True)

dat_estimates['Variable'] = dat_estimates['Variable'].fillna('Name')
n = dat_estimates['Variable'][(dat_estimates['Variable'] == 'Name')].count()

l1 = ['Missing']*n
l2 = list( range(1,n+1))
l2 =  [str(x) for x in l2]

new_column_names = [i + j for i, j in zip(l1, l2)]

dat_estimates.loc[dat_estimates['Variable'] == 'Name', 'Variable'] = new_column_names

dat_estimates = dat_estimates.transpose()
dat_estimates.columns = dat_estimates.iloc[0] 
dat_estimates = dat_estimates[1:]


dat_testntreat = pd.read_excel(xls, 'HIV Test&Treat (records)')
dat_testntreat.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900,'Year ':1990},inplace=True)


dat_testntreat['Variable'] = dat_testntreat['Variable'].fillna('Name')
n = dat_testntreat['Variable'][(dat_testntreat['Variable'] == 'Name')].count()

l1 = ['Missing']*n
l2 = list( range(1,n+1))
l2 =  [str(x) for x in l2]

new_column_names = [i + j for i, j in zip(l1, l2)]

# dat_testntreat['Variable'][(dat_testntreat['Variable'] == 'Name')] = new_column_names
dat_testntreat.loc[(dat_testntreat['Variable'] == 'Name'),'Variable'] = new_column_names


dat_testntreat = dat_testntreat.transpose()
dat_testntreat.columns = dat_testntreat.iloc[0] 
dat_testntreat = dat_testntreat[1:]


dat_resistance = pd.read_excel(xls, 'HIV Drug Resistance')
dat_resistance.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900,'Year ':1990},inplace=True)

dat_resistance['Variable'] = dat_resistance['Variable'].fillna('Name')
n = dat_resistance['Variable'][(dat_resistance['Variable'] == 'Name')].count()

l1 = ['Missing']*n
l2 = list( range(1,n+1))
l2 =  [str(x) for x in l2]

new_column_names = [i + j for i, j in zip(l1, l2)]

dat_resistance.loc[(dat_resistance['Variable'] == 'Name'),'Variable'] = new_column_names


dat_resistance = dat_resistance.transpose()
dat_resistance.columns = dat_resistance.iloc[0] 
dat_resistance = dat_resistance[3:]



# dat_key = pd.read_excel(xls, 'Key populations statistics')
# dat_key.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900,'Year ':1990},inplace=True)

# dat_key['Variable'] = dat_key['Variable'].fillna('Name')
# n = dat_key['Variable'][(dat_key['Variable'] == 'Name')].count()

# dat_key = dat_key.transpose()
# dat_key.columns = dat_key.iloc[0] 
# dat_key = dat_key[1:]







# dat_program = pd.read_excel(xls, 'Program-related parameters')
# dat_program.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900,'Year ':1990},inplace=True)


# dat_program['Variable'] = dat_program['Variable'].fillna('Name')
# n = dat_program['Variable'][(dat_program['Variable'] == 'Name')].count()

# l1 = ['Missing']*n
# l2 = list( range(1,n+1))
# l2 =  [str(x) for x in l2]

# new_column_names = [i + j for i, j in zip(l1, l2)]

# dat_program.loc[(dat_program['Variable'] == 'Name'),'Variable'] = new_column_names


# dat_program = dat_program.transpose()
# dat_program.columns = dat_program.iloc[0] 
# dat_program = dat_program[1:]








##### New UNAIDS SPECTRUM updated in 2022 

# xls = pd.ExcelFile('Data/Spectrum 2024/Extraction data_July22nd2024.xlsx')
xls = pd.ExcelFile(os.path.join(base_path, 'Data','Spectrum 2024', 'Extraction data_July22nd2024.xlsx'))


dat_estimates_2023 = pd.read_excel(xls, 'Estimates')
dat_estimates_2023.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900},inplace=True)

dat_estimates_2023['Variable'] = dat_estimates_2023['Variable'].fillna('Name')
n = dat_estimates_2023['Variable'][(dat_estimates_2023['Variable'] == 'Name')].count()

l1 = ['Missing']*n
l2 = list( range(1,n+1))
l2 =  [str(x) for x in l2]

new_column_names = [i + j for i, j in zip(l1, l2)]

dat_estimates_2023.loc[dat_estimates_2023['Variable'] == 'Name', 'Variable'] = new_column_names

dat_estimates_2023 = dat_estimates_2023.transpose()
dat_estimates_2023.columns = dat_estimates_2023.iloc[0] 
dat_estimates_2023 = dat_estimates_2023[1:]




dat_testntreat_2023 = pd.read_excel(xls, 'TestnTreat')
dat_testntreat_2023.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900,'Year ':1990,'Unnamed: 2':1900},inplace=True)


dat_testntreat_2023['Variable'] = dat_testntreat_2023['Variable'].fillna('Name')
n = dat_testntreat_2023['Variable'][(dat_testntreat_2023['Variable'] == 'Name')].count()

l1 = ['Missing']*n
l2 = list( range(1,n+1))
l2 =  [str(x) for x in l2]

new_column_names = [i + j for i, j in zip(l1, l2)]

# dat_testntreat_2023['Variable'][(dat_testntreat_2023['Variable'] == 'Name')] = new_column_names
dat_testntreat_2023.loc[(dat_testntreat_2023['Variable'] == 'Name'),'Variable'] = new_column_names


dat_testntreat_2023 = dat_testntreat_2023.transpose()
dat_testntreat_2023.columns = dat_testntreat_2023.iloc[0] 
dat_testntreat_2023 = dat_testntreat_2023[1:]




######## old data extraction from previous years of UNAIDS estimates

dat_depraciated_2023 = pd.read_excel(xls, 'Estimates-2023')
dat_depraciated_2023.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900},inplace=True)

dat_depraciated_2023['Variable'] = dat_depraciated_2023['Variable'].fillna('Name')
n = dat_depraciated_2023['Variable'][(dat_depraciated_2023['Variable'] == 'Name')].count()

l1 = ['Missing']*n
l2 = list( range(1,n+1))
l2 =  [str(x) for x in l2]

new_column_names = [i + j for i, j in zip(l1, l2)]

dat_depraciated_2023.loc[dat_depraciated_2023['Variable'] == 'Name', 'Variable'] = new_column_names

dat_depraciated_2023 = dat_depraciated_2023.transpose()
dat_depraciated_2023.columns = dat_depraciated_2023.iloc[0] 
dat_depraciated_2023 = dat_depraciated_2023[1:]





dat_depraciated_2022 = pd.read_excel(xls, 'Estimates-2022')
dat_depraciated_2022.rename(columns={'Unnamed: 0':'Variable','Unnamed: 1':1900},inplace=True)

dat_depraciated_2022['Variable'] = dat_depraciated_2022['Variable'].fillna('Name')
n = dat_depraciated_2022['Variable'][(dat_depraciated_2022['Variable'] == 'Name')].count()

l1 = ['Missing']*n
l2 = list( range(1,n+1))
l2 =  [str(x) for x in l2]

new_column_names = [i + j for i, j in zip(l1, l2)]

dat_depraciated_2022.loc[dat_depraciated_2022['Variable'] == 'Name', 'Variable'] = new_column_names

dat_depraciated_2022 = dat_depraciated_2022.transpose()
dat_depraciated_2022.columns = dat_depraciated_2022.iloc[0] 
dat_depraciated_2022 = dat_depraciated_2022[1:]






### Corrected population of children 0-14 years old 
# xls = pd.ExcelFile('Data/World Bank 2024_/API_SP.POP.0014.TO_DS2_en_excel_v2_1594023.xls')
xls = pd.ExcelFile(os.path.join(base_path, 'data','World Bank 2024_', 'API_SP.POP.0014.TO_DS2_en_excel_v2_1594023.xls'))


dat_population_14 = pd.read_excel(xls, 'Data')
dat_population_14 = dat_population_14.loc[dat_population_14['Country Name'] == 'Papua New Guinea',]
dat_population_14 = dat_population_14.transpose()
dat_population_14.columns = ['Total_Pop']
dat_population_14 = dat_population_14.iloc[4:,]
dat_population_14.index = dat_population_14.index.map(int)



### Updated population dataset from 2023 - corrected 
# xls = pd.ExcelFile('Data/World Bank 2024_/API_SP.POP.TOTL_DS2_en_excel_v2_1584408.xls')
xls = pd.ExcelFile(os.path.join(base_path, 'data','World Bank 2024_', 'API_SP.POP.TOTL_DS2_en_excel_v2_1584408.xls'))


dat_population_2023 = pd.read_excel(xls, 'Data')
dat_population_2023 = dat_population_2023.loc[dat_population_2023['Country Name'] == 'Papua New Guinea',]
dat_population_2023 = dat_population_2023.transpose()
dat_population_2023.columns = ['Total_Pop']
dat_population_2023 = dat_population_2023.iloc[4:,]
dat_population_2023.index = dat_population_2023.index.map(int)



#### Adults older than 15 years old 


# Subtracting the Total_pop columns
corrected_population = dat_population_2023['Total_Pop'] - dat_population_14['Total_Pop']

corrected_population = pd.DataFrame({'Total_Pop': corrected_population})


def extract_values(value):
    if not isinstance(value, str):
        return np.nan, np.nan, np.nan

    match = re.search(r'(\d+\.\d+) \[(\d+\.\d+) - (\d+\.\d+)\]', value)
    if match:
        estimate = float(match.group(1))
        lower_bound = float(match.group(2))
        upper_bound = float(match.group(3))
        return estimate, lower_bound, upper_bound
    else:
        return np.nan, np.nan, np.nan
### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
from model.acquire_distribution import *
# Set the seed for reproducibility
np.random.seed(334) 

#### Get the number of data from Posteriors samples
pd.options.display.float_format = '{:.0f}'.format


#### This function will get the Dictionary of data items and translate it into 
#### vectors/ list of items that will match the number of simulations in the Posterior list
def generate_p_values(p_values_specified, n):
    p_values = []
    
    for i in range(n):
        new_dict = {}
        for key, (value, flag) in p_values_specified.items():
            ### If this is an array of values, as in the randomly generated costs in my model
            if isinstance(value, np.ndarray):  
                new_dict[key] = (value[i], flag)
            else:
                new_dict[key] = (value, flag)
        p_values.append(new_dict)
    return p_values

#### just get the total number of posteriors that were sampled for this modelling 
# xls = pd.ExcelFile('Bayesian Pictures/Bayesian_posterior.xlsx')
xls = pd.ExcelFile(os.path.join(base_path, 'data',"posteriors", 'Bayesian_posterior.xlsx'))
posteriors = pd.read_excel(xls, 'Sheet1')
n = posteriors.shape[0]


#### acquiring the distributions for the modelling of costs
xls = pd.ExcelFile(os.path.join(base_path, 'data', 'cea', 'Cost_data_November29th2024.xlsx'))
all_costs = pd.read_excel(xls, 'Sheet1')


p_ART = acquire_gamma(all_costs,"c_ARTCost_Adults",no_of_samples=n)
p_ART_2ndline = acquire_gamma(all_costs,"c_ARTCost_Adults_Secondline",no_of_samples=n)
p_VLtest = acquire_gamma(all_costs,"c_VLTest",no_of_samples=n)
p_VLtest_labour = acquire_gamma(all_costs,"c_VLTest_Labour",no_of_samples=n)
p_AdherenceCounsel = acquire_gamma(all_costs,"c_HIVCare_AdherenceCounselling",no_of_samples=n)
p_DRTesting = acquire_gamma(all_costs,"c_DRTesting",no_of_samples=n)
p_DRTesting_labour = acquire_gamma(all_costs,"c_DRTesting_Labour",no_of_samples=n)
p_cDiagnosis = acquire_gamma(all_costs,"c_HIV_Diagnosis",no_of_samples=n)
p_ClinicVisit = acquire_gamma(all_costs,"c_ARTVisit",no_of_samples=n)

p_ARTinitiation = acquire_gamma(all_costs,"c_ART_initiation",no_of_samples=n)
p_probOldART = acquire_beta(all_costs,"prob_OldART",no_of_samples=n)
p_probNewART = acquire_beta(all_costs,"prob_NewART",no_of_samples=n)

p_costs_specified = {
    'I_W': (0, False),        
    'I_DR': (0, False),      
    'D_W': (0, False),
    'D_DR': (0,False),
    'T_W1': (p_ART + p_ClinicVisit,False),
    'T_W2': (p_ART + p_ClinicVisit,False),
    'T_DR1': (p_ART + p_ClinicVisit,False),
    'T_DR2': (p_ART_2ndline + p_ClinicVisit,False), 
    'F_W1': ((p_ART + p_ClinicVisit)*p_probOldART,False),
    'F_W2': ((p_ART + p_ClinicVisit)*p_probNewART,False), 
    'F_TDR1': ((p_ART + p_ClinicVisit)*p_probOldART,False),
    'F_ADR1': ((p_ART + p_ClinicVisit)*p_probOldART,False),
    # 'F_DR2': ((p_ART + p_ARTVisit_LTFU)*p_probNewART,False),
    ### Those in F_DR2 likely LTFU
    'F_DR2': ((p_ART_2ndline + p_ClinicVisit)*p_probNewART,False),
    ### This is to calculate VL test in the routine settings - it is not applicable here in the POC Comparison
    # 'VLtest_routine': (p_VLtest + p_VLtest_labour + p_AdherenceCounsel,True),
    #### Costs done in VL visit
    ### VL Test and Labour-related costs
    ### Adherence counselling
    ### Cost of visit to the clinics 
    'VLtest_POC': (p_VLtest + p_VLtest_labour + p_AdherenceCounsel + p_ClinicVisit,True),
    ##### Cost of DR testing pre-treatment
    ### DR Test and Labour-related costs
    ### VL Test and Labour-related costs
    ### Cost of visit to the clinics 
    ### Cost of re-initiating on a New line of ART for the PLHIV
    'preDRtest':(p_DRTesting + p_DRTesting_labour + p_ARTinitiation + 2*(p_VLtest + p_VLtest_labour + p_AdherenceCounsel  + p_ClinicVisit),True),
    ### Probably should ignore the costs associated with DR tests routine
    # 'postDRtest_routine':(p_DRTesting + p_DRTesting_labour + p_VLtest + p_VLtest_labour,True),
    ##### Cost of DR testing pre-treatment
    ### DR Test and Labour-related costs
    ### VL Test and Labour-related costs
    ### Cost of visit to the clinics 
    ### ART re-initiation
    'postDRtest_POC':(p_DRTesting + p_DRTesting_labour + p_ARTinitiation + 2*(p_VLtest + p_VLtest_labour + p_AdherenceCounsel  + p_ClinicVisit),True)
    
    ,'deaths': (0,True)
    ,'diagnoses':(p_cDiagnosis,True)
    ,'treat_inits':(p_ARTinitiation,True)
}

p_costs = generate_p_values(p_costs_specified,n)



#### acquiring the distributions for the modelling of costs
xls = pd.ExcelFile(os.path.join(base_path, 'data', 'cea', 'Disability_weights_data_November27th.xlsx'))
all_dalys = pd.read_excel(xls, 'Sheet1')


p_Suppressed = acquire_beta(all_dalys,"HIV/AIDS with antiretroviral treatment without anemia",no_of_samples=n)
### Untreated with symptoms and AIDS
p_Untreated_AIDS = acquire_beta(all_dalys,"AIDS with mild anemia",no_of_samples=n)
### Untreated without any symptoms yet
p_Untreated_Early = acquire_beta(all_dalys,"Early HIV with severe anemia",no_of_samples=n)
prob_AIDS = acquire_beta(all_dalys,"Proportion of AIDS",no_of_samples=n)
p_Untreated = prob_AIDS*p_Untreated_AIDS + (1-prob_AIDS)*p_Untreated_Early

u_mortality = acquire_basevalues(all_dalys,"Mortality")

p_daly_specified = {
    'I_W': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early, False),        
    'I_DR': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early, False),      
    'D_W': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early, False),
    'D_DR': (prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,False),
    #### main analyses all suppressed
    'T_W1': (p_Suppressed,False),
    'T_W2': (p_Suppressed,False),
    'T_DR1': (p_Suppressed,False),
    'T_DR2': (p_Suppressed,False), 
    'F_W1': (p_Suppressed,False),
    'F_W2': (p_Suppressed,False), 
    'F_TDR1': (p_Suppressed,False),
    'F_ADR1': (p_Suppressed,False),
    'F_DR2': (p_Suppressed,False),
    #### deaths will be have to be accumulated, in the DALYs calculations
    'all_mortality': (u_mortality,False)
    # 'deaths': (u_mortality,True)

}


p_daly = generate_p_values(p_daly_specified,n)

death_matrix = {
    'all_mortality': (1,True)
}

p_deaths = generate_p_values(death_matrix,n)

incidence_matrix = {
    'incidence': (1,True)
}
p_incidences = generate_p_values(incidence_matrix,n)

incidence_DR_matrix = {
    'incidence_Resist': (1,True)
}
p_incidences_DR = generate_p_values(incidence_DR_matrix,n)

##### accumulations of TF due ADR
TF_ADR_matrix = {
    'newAcquiredDR': (1,True)
}
p_ADR = generate_p_values(TF_ADR_matrix,n)

TF_TDR_matrix = {
    'newTransmittedDR': (1,True)
}
p_TDR = generate_p_values(TF_TDR_matrix,n)

secondline_matrix = {
    'T_DR2': (1,False)
}
p_secondline = generate_p_values(secondline_matrix,n)

mistreated_matrix = {
    'newMisTreatedDR': (1, True)
}
p_mistreated = generate_p_values(secondline_matrix,n)
#### accumulations of Total compartments
Total_ADR_matrix = {
    'F_ADR1': (1,False)
}
p_ADR_Total = generate_p_values(Total_ADR_matrix,n)

Total_TDR_matrix = {
    'F_TDR1': (1,False)
}
p_TDR_Total = generate_p_values(Total_TDR_matrix,n)


Total_treatmentF_matrix = {
    'F_TDR1': (1,False),
    'F_ADR1': (1,False)
}
p_Total_TFailure = generate_p_values(Total_treatmentF_matrix,n)

Viral_unsuppressed_matrix = {
    'F_W1': (1,False),
    'F_W2': (1,False),
    'F_TDR1': (1,False),
    'F_ADR1': (1,False),
    'F_DR2': (1,False)
}
p_Viral_unsuppressed = generate_p_values(Viral_unsuppressed_matrix,n)

new_TFailure_matrix = {
    'newAcquiredDR': (1,True),
    'newTransmittedDR': (1,True)
}
p_new_TFailure = generate_p_values(new_TFailure_matrix,n)


#### Exporting key variables for CEA Table range in PSA
cost_items = {
    "p_cDiagnosis": p_cDiagnosis,
    "p_ARTinitiation": p_ARTinitiation,
    "p_ClinicVisit": p_ClinicVisit,
    "p_ART": p_ART,
    "p_ART_2ndline": p_ART_2ndline,
    "p_VLtest": p_VLtest,
    "p_VLtest_labour": p_VLtest_labour,
    "p_AdherenceCounsel": p_AdherenceCounsel,
    "p_DRTesting": p_DRTesting,
    "p_DRTesting_labour": p_DRTesting_labour,
}

cost_ranges = {
    name: f"${round(values.min(), 1)} to ${round(values.max(), 1)}"
    for name, values in cost_items.items()
}


perc_items = {
    "p_probOldART": p_probOldART,
    "p_probNewART": p_probNewART,
    "prob_AIDS":prob_AIDS
}

perc_ranges = {
    name: f"{round(values.min()*100, 1)}% to {round(values.max()*100, 1)}%"
    for name, values in perc_items.items()
}

daly_items = {
    "d_weight_untreated": prob_AIDS * p_Untreated_AIDS + (1-prob_AIDS) * p_Untreated_Early,
    "d_weight_untreated_early":p_Untreated_Early,
    "d_weight_AIDS":p_Untreated_AIDS,
    "d_weight_treated":p_Suppressed
}

daly_ranges = {
    name: f"{round(values.min(), 3)} to {round(values.max(), 3)}"
    for name, values in daly_items.items()
}
cost_df = pd.DataFrame(list(cost_ranges.items()), columns=["Variable", "Range"])

perc_df = pd.DataFrame(list(perc_ranges.items()), columns=["Variable", "Range"])

daly_df = pd.DataFrame(list(daly_ranges.items()), columns=["Variable", "Range"])

combined_df = pd.concat([cost_df, perc_df, daly_df], ignore_index=True)

# output_file = "CEA_results/Output_Results/Final_CEA_Table1.xlsx"
output_file = os.path.join(base_path, 'output', 'Chp3_tables', 'Chp3 - Table1.xlsx')
combined_df.to_excel(output_file, index=False)

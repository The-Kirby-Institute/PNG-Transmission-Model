# from Data_Importing.importing import *
from Data_Importing.data_processing import *
import copy

from model.Model_ode import *

mode_cal = "Posteriors"

I_W = dat_estimates_2023['Estimated adults (15+) living with HIV'][starting_epidemic_year]  
### check here for population >15 years old
S=corrected_population['Total_Pop'][starting_epidemic_year] - I_W
I_DR = 0
D_W = 0
D_DR=0
T_W1=0
T_DR1=0
T_W2=0
T_DR2=0
F_W1=0
F_TDR1=0
F_ADR1=0
F_W2=0
F_DR2=0
incidence = 0
deaths = 0
total_PLHIV = I_W + I_DR + D_W + D_DR + T_W1 + T_DR1 + T_W2 + T_DR2+F_W1+ F_TDR1+ F_ADR1+ F_W2+ F_DR2
incidence_Resist = 0
preDRtest = 0
postDRtest_routine = 0
postDRtest_POC = 0
VLtest_routine = 0
VLtest_POC = 0
newAcquiredDR = 0
newTransmittedDR = 0
newMisTreatedDR = 0
diagnoses = 0
treat_inits = 0
all_mortality = 0

initial_conditions = [S,I_W, I_DR, D_W, D_DR, T_W1,T_W2, T_DR1,T_DR2, F_W1,F_W2, F_TDR1,F_ADR1,F_DR2,incidence,deaths,total_PLHIV,incidence_Resist,preDRtest,postDRtest_routine,postDRtest_POC,VLtest_routine,VLtest_POC,newAcquiredDR,newTransmittedDR,newMisTreatedDR,diagnoses,treat_inits,all_mortality]
initial_conditions = [x/sum(initial_conditions) for x in initial_conditions]

#### add in one more year to make 2030 results
years = 2030 - 1994 + 1
tspan = np.arange(0, years, 1)

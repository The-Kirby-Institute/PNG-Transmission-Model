from model.partial_functions import *
from functools import partial
import copy

from model.core_scenarios import *
from model.acquire_basevalues import *

switch_Infs = partial(reverse_expit,k=year_OffInfs)
###### New Force of infection 2023 April for beta_u beta_f model  
#### Force of infection 
beta = partial(Richardsfunc,k = 0)
beta_r = partial(Richardsfunc,a = 0)

###### New diagnosis rate 
#### UNAIDS 2024 decided to remove the historial 2019 data point. 
#### Diagnoses are thought to be started much later now 
b = partial(Richards_scaleupfunc,a=-0.01,q=10)



 ##### New treatment rates
c1 = partial(Richards_scaleupfunc,a=-0.05,q=13)
c_1 = partial(firstlineswitch,treatment_rate=c1, newfirstline=False)
c_1prime = partial(firstlineswitch,treatment_rate=c1, newfirstline=True)

c2 = partial(Richards_scaleupfunc,a=-0.05,q=13)
### standard rate of DR at the beginning = 0
c_2 = partial(drugresistancescale,treatment_rate=c2, secondline=False,drug_resistance_rate = 0)
c_2prime = partial(drugresistancescale,treatment_rate=c2, secondline=True,drug_resistance_rate = 0)



### Treatment transfer rate (starting from 2020)
k = partial(firstlinetransfer)








f1_rate=f1
f2_rate=f2
f3_rate=f3
f4_rate=f2
#### Functions related to rates of treatment failures 

################################################################################################
################################################################################################
################################################################################################


f1 = partial(virologicalfailure_oldfirstline,failure_rate = f1)
f2 = partial(virologicalfailure_firstline,failure_rate=f2)
f3 = partial(virologicalfailure_drugresistance,newdolutegravir_rate = f3)
f4 = partial(virologicalfailure_firstline,failure_rate=f4)






#### rate of being tested in any year up till 2021 from 2016 (for 5 years)
routine_coverage = partial(annualized_routineVL,coverage = 11/100)
#### Level of access to DR testing, if there is no POC VL in place = 0.0%
noaccess_coverage = partial(annualized_routineVL,coverage = 0.000)


if scenario == "2":
    routine_coverage = partial(annualized_routineVL,coverage = 80/100 * 11/100)
    POC_coverage = partial(annualized_POCVL,POC_rate = 20/100 * 40/100 )
    POC_access  = partial(annualized_POCVL, POC_rate = 20/100)

elif scenario =="3":

    routine_coverage = partial(annualized_routineVL,coverage = 60/100 * 11/100)
    POC_coverage = partial(annualized_POCVL,POC_rate = 40/100 * 40/100  )
    POC_access  = partial(annualized_POCVL, POC_rate = 40/100)

elif scenario =="4":

    routine_coverage = partial(annualized_routineVL,coverage = 40/100 * 11/100)
    POC_coverage = partial(annualized_POCVL,POC_rate = 60/100 * 40/100)
    POC_access  = partial(annualized_POCVL, POC_rate = 60/100)





##############################################################################################################################################################################################
##############################################################################################################################################################################################
############################################################################################################################################################################################################################################################################################################################################################################################










#### Time function to calibrate the levels of viral load suppression in the country 
time_function = partial(logisticfunc,k=1.85)



##### Effectiveness of ACT-UP Viral load (POC) is in the actup variable parameters (being read from Excel sheet)
g1 = partial(ScaleupVLfunc,time_T_F = time_function,scale_uprate = POC_access,type_VL = Monitoring_Type,ACTUP_e = actup)
g_1 = partial(firstlineswitch_VLtesting,VLtesting_rate=g1, newfirstline=False)
g_1prime = partial(firstlineswitch_VLtesting,VLtesting_rate=g1, newfirstline=True)
g2 = partial(ScaleupVLfunc,time_T_F = time_function,scale_uprate = POC_access,type_VL = Monitoring_Type,ACTUP_e = actup)
g3 = partial(ScaleupDRtestingfunc,time_T_F = time_function,scale_uprate = POC_coverage,type_testing = Monitoring_Type)
g4 = partial(ScaleupVLfunc,time_T_F = time_function,scale_uprate = POC_access,type_VL = Monitoring_Type,ACTUP_e = actup)


##### For the pre-treatment DR testing and the scenario is time-varying  
##### The reason for this piece of code is because in Time-varying conditions, the drug_resistance_rate is a function instead of a number --> require a different way to ahndle
if DR_testing_scenario == "Timevarying":
    if Monitoring_Type == "Routine":
        #### There is no drug resistance testing without ACT-UP interventions 
        c_2 = partial(drugresistancescale,treatment_rate=c2, secondline=False,drug_resistance_rate = noaccess_coverage)
        c_2prime = partial(drugresistancescale,treatment_rate=c2, secondline=True,drug_resistance_rate = noaccess_coverage)
    elif Monitoring_Type == "POC" or Monitoring_Type == "POC_AcquiredDR":
        c_2 = partial(drugresistancescale,treatment_rate=c2, secondline=False,drug_resistance_rate = POC_access)
        c_2prime = partial(drugresistancescale,treatment_rate=c2, secondline=True,drug_resistance_rate = POC_access)




##### Initialize time-varying rates that track the actual testing (DR and VL) rates of Viral Load and Drug Resistance based on our different scenarios 

g1Testingrate = partial(actualTestingfunc,routine_rate = routine_coverage,scale_uprate = POC_coverage,type_VL = Monitoring_Type)
g2Testingrate = partial(actualTestingfunc,routine_rate = routine_coverage,scale_uprate = POC_coverage,type_VL = Monitoring_Type)
g3Testingrate = partial(actualTestingfunc,routine_rate = routine_coverage,scale_uprate = POC_coverage,type_VL = Monitoring_Type)
g4Testingrate = partial(actualTestingfunc,routine_rate = routine_coverage,scale_uprate = POC_coverage,type_VL = Monitoring_Type)

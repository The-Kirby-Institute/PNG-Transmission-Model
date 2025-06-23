from model.function_response import *
from model.core_scenarios import *
from functools import partial
# from  pymc.math  import sigmoid

# from scipy.special import expit  # expit function provides the logistic function
import numpy as np

starting_epidemic_year  = 1994
reference_year = 2000
ending_forecast_year = 2030
dolutegravir_year = 2020
consistent_treatment_year = 2005
adj =  reference_year - starting_epidemic_year

yearofRoutine_VL = 2016
yearofPOC_VL = 2021
yearofDR_testing = 2023

def expit(x,k):
    if analysis_mode == "Fitting":
        return 1 / (1 + np.exp(-10*(x-k)))
    elif analysis_mode == "Inference":
        return np.heaviside(x - k, 0)
    else:
        raise ValueError("analysis_mode must be either 'Inference' or 'Fitting'")


def reverse_expit(x,k):
    if analysis_mode == "Fitting":
        return 1 - expit(x,k)
    elif analysis_mode == "Inference":
        return np.heaviside(k - x, 0)
    else:
        raise ValueError("analysis_mode must be either 'Inference' or 'Fitting'")



##### Using year_OffInterv, turn off interventions (test,treatment and POC VL and DR testing) after year_OffInterv (year)
switch_Intervns = partial(reverse_expit,k=year_OffInterv)




def Richardsfunc (t,b,a=1,k=1,q=1):
    return richardscurve(t,a,b,k,q)

def Richards_scaleupfunc (t,q,a=1,b=1,k=1):
    #### implement turning off tests and treatment by multiplying the swithch intervention 
    return richardscurve_scalingup(t,a,b,k,q) * switch_Intervns(t)


def logisticfunc(t,L=1, k=1, t0=yearofRoutine_VL-starting_epidemic_year):
    return logistic_curve(t,L, k, t0)



def annualized_routineVL(t, coverage = 0.11):
    value = coverage*expit(t,yearofRoutine_VL - starting_epidemic_year)
    return value




def annualized_POCVL(t,POC_rate):
    value = POC_rate * expit(t,yearofPOC_VL - starting_epidemic_year)
    #### implement the functions, while also multiplied by the switch to turn off interventions after year_OffInterv
    return value * switch_Intervns(t)

#### New implementation of viral load monitoring - with adjustment to the viral load monitoring rates
def ScaleupVLfunc(t,time_T_F,scale_uprate, baseline_CD4monitoring= -np.log(1-5/100)/5,type_VL = "CD4", counsel_f = 0.70,qt = 1.0,ACTUP_e = 2.69):
    if type_VL == "CD4":
        scaleupVL =  baseline_CD4monitoring
    elif type_VL == "Routine":
        scaleupVL = qt*counsel_f * time_T_F(t)
    elif type_VL == "POC" or type_VL == "POC_AcquiredDR":
        scaleupVL = ( qt*(1-scale_uprate(t)) + ACTUP_e*qt*scale_uprate(t) )*counsel_f  * time_T_F(t)
    else: 
        raise ValueError('Type of test not specified')
    return scaleupVL


#### New implementation of Drug Resistance testing DR in the model 
def ScaleupDRtestingfunc(t,time_T_F,scale_uprate,bg_drtest = 0.05,type_testing  = "CD4", qt = 1.0):
    if type_testing == "POC_AcquiredDR":
        if t < yearofDR_testing - starting_epidemic_year:
            ## background (bg) rate of switch to second-line ART (resuppression rate) 
            DRresuppressionrate =  bg_drtest * time_T_F(t)
        else:
            ## once DR testing is in place, the rate of viral load resuppression is the same as if VL test is available 
            ## resuppression rate due to DR testing 
            ## counselling is not in play here --> should be 100% effective as long as the result of the 
            ## DR test is known
            DRresuppressionrate = ( bg_drtest*(1-scale_uprate(t)) + qt*scale_uprate(t)  ) * time_T_F(t)
    else:
        DRresuppressionrate = bg_drtest*time_T_F(t)
    return DRresuppressionrate

####### Implementation the actual Point-of-care testing rates across different scenarios
####### The above function is for the effect of VL testing on the viral load suppression rates, not the actual testing rates
####### The difference is that we need a calibration factor max_failureate*time_T_F(t) to match the actual viral suppression levels 

def actualTestingfunc(t,routine_rate,scale_uprate,type_VL = "CD4"):
    ### routine_rate and scale_uprate are the coverage of VL tests for Routine VL and POC VL testing, respectively
    ### Coverage of POC VL also includes some that are done in central-lab based testing (which POCVL aims to gradually replace)
    ### This function measures the rate of each tests that are being conducted 
    if type_VL == "CD4":
        actualVLrate_routine = 0
        actualVLrate_POC = 0
        actualDRrate_routine = 0
        actualDRrate_POC = 0
    elif type_VL == "Routine":
        actualVLrate_routine = routine_rate(t)
        actualVLrate_POC = 0
        actualDRrate_routine = 0
        actualDRrate_POC = 0
    elif type_VL == "POC":
        actualVLrate_routine = routine_rate(t)
        actualVLrate_POC = scale_uprate(t) 
        actualDRrate_routine = 0
        actualDRrate_POC = 0   
    elif type_VL == "POC_AcquiredDR":
        #### This is because DR testing is only available after yearof DR testing
        #### but routine_rate(t) starts from 2016 - which is too early for DR testing 
        #### so the function needs to be adjusted
        #### the rate of suppression in the model is according to scale_uprate(t), not routine_rate(t)
        if t < yearofDR_testing - starting_epidemic_year:
            actualVLrate_routine = routine_rate(t)
            actualVLrate_POC = scale_uprate(t) 
            actualDRrate_routine = 0
            actualDRrate_POC = 0
        else:
            actualVLrate_routine = routine_rate(t)
            actualVLrate_POC = scale_uprate(t) 
            actualDRrate_routine = routine_rate(t)
            actualDRrate_POC = scale_uprate(t) 
    return actualVLrate_routine , actualVLrate_POC,actualDRrate_routine,actualDRrate_POC




def firstlineswitch(t,b,k, treatment_rate,newfirstline=False):
    if newfirstline == False:
        value = (1- expit(t,dolutegravir_year-starting_epidemic_year))*treatment_rate(t,b=b,k=k)
    elif newfirstline == True:
        value = expit(t,dolutegravir_year-starting_epidemic_year)*treatment_rate(t,b=b,k=k)

    return value




def firstlinetransfer(t):
    # value = logisticfunc(t,k=rate,t0=dolutegravir_year-starting_epidemic_year)
    value = expit(t,dolutegravir_year-starting_epidemic_year)
    return value







def _drugresistancescale(t,b,k,treatment_rate, secondline = False,drug_resistance_rate=0):
    if secondline == False:
        value = (1 - drug_resistance_rate*expit(t,yearofDR_testing - starting_epidemic_year)) * treatment_rate(t,b=b,k=k)
    elif secondline == True:
        value = drug_resistance_rate*expit(t,yearofDR_testing - starting_epidemic_year)  * treatment_rate(t,b=b,k=k)

    return value

def _drugresistancescale_timevarying(t,b,k,treatment_rate,drug_resistance_rate, secondline = False):
    if secondline == False:
        value = (1 - drug_resistance_rate(t)*expit(t,yearofDR_testing - starting_epidemic_year)) * treatment_rate(t,b=b,k=k)
    elif secondline == True:
        value = drug_resistance_rate(t)*expit(t,yearofDR_testing - starting_epidemic_year)  * treatment_rate(t,b=b,k=k)

    return value






def firstlineswitch_VLtesting(t,VLtesting_rate,counsel_f,qt, newfirstline=False):
    if newfirstline == False:
        value = (1- expit(t,dolutegravir_year - starting_epidemic_year))*VLtesting_rate(t,counsel_f=counsel_f,qt=qt)
    elif newfirstline == True:
        value = expit(t,dolutegravir_year - starting_epidemic_year)*VLtesting_rate(t,counsel_f=counsel_f,qt=qt)

    return value






def virologicalfailure_firstline(t,failure_rate):
    return failure_rate 


def virologicalfailure_drugresistance(t,newdolutegravir_rate):
    return 1.0 + (newdolutegravir_rate
                   - 1.0) * expit(t,dolutegravir_year- starting_epidemic_year)



def virologicalfailure_oldfirstline(t,failure_rate):
    return 1.0 + (failure_rate
                   - 1.0) * expit(t,consistent_treatment_year - starting_epidemic_year)






if DR_testing_scenario == "Timevarying":
    drugresistancescale = _drugresistancescale_timevarying
elif DR_testing_scenario in ("Baseline","Medium","High"):
    drugresistancescale = _drugresistancescale

import numpy as np
from model.core_scenarios import *



def richardscurve(t,a,b,k,q):
    return a + (k-a)/(1 +  np.exp(-b*(t-q))) 

def _richardscurve_scalingup_v1(t,a,b,k,q):
    value = a + (k-a)/(1 + np.exp(-b*(t-q))) 
    value = np.where(value < 0.00,0,value)
    return value

def _richardscurve_scalingup_v2(t,a,b,k,q):
    value = a + (k-a)/(1 + np.exp(-b*(t-q))) 
    return value


def logistic_curve(t, L, k, t0):
    return L / (1 + np.exp(-k * (t - t0)))


if analysis_mode == "Inference":
    richardscurve_scalingup = _richardscurve_scalingup_v1
elif analysis_mode == "Fitting":
    richardscurve_scalingup = _richardscurve_scalingup_v2
else:
    raise ValueError("analysis_mode must be either 'Inference' or 'Fitting'")
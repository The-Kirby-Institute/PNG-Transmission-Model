from scipy.integrate import odeint
from model.meta_parameters import *






def ode_model(state,t,params): 
    #### State all the compartments in the model    
    #### this line of code needs to be rechecked again, in case the order of infections is wrong 
    # N,S,I_W, I_DR, D_W, D_DR, T_W1,T_W2, T_DR1,T_DR2, F_W1,F_W2, F_DR1,F_DR2,incidence,deaths,total_PLHIV = state

    # S,I_W, I_DR, D_W, D_DR, T_W1,T_W2, T_DR1,T_DR2, F_W1,F_W2, F_TDR1,F_ADR1,F_DR2,incidence,deaths,incidence_Resist = state
    S = state[0]
    I_W= state[1]
    I_DR= state[2]
    D_W= state[3]
    D_DR= state[4]
    T_W1= state[5]
    T_W2= state[6]
    T_DR1= state[7]
    T_DR2= state[8]
    F_W1= state[9]
    F_W2= state[10]
    F_TDR1= state[11]
    F_ADR1= state[12]
    F_DR2= state[13]
    incidence= state[14]
    deaths= state[15]
    total_PLHIV = state[16]
    incidence_Resist = state[17]
    preDRtest = state[18]
    postDRtest_routine = state[19]
    postDRtest_POC = state[20]
    VLtest_routine = state[21]
    VLtest_POC = state[22]


    ### new compartments to be added to the model 
    newAcquiredDR = state[23]
    newTransmittedDR = state[24]
    newMisTreatedDR = state[25]


    ### compartments for new diagnoses and new treatment initiations
    diagnoses = state[26]
    treat_inits = state[27]


    ### all mortality including susceptible deaths
    all_mortality = state[28]


    theta = params[0]
    beta_u = params[1]
    beta_t = params[2]
    beta_f = params[3]
    delta_U = params[4]  
    delta_T = params[5]
    delta_F = params[6]
    delta_B = params[7]
    h1 = params[8]
    h2 = params[9]
    eta_1 = params[10]
    eta_2 = params[11]
    eta_3 = params[12]
    eta_4 = params[13]
    eta_5 = params[14]
    eta_6 = params[15]
    b_asterisk = params[16]
    b_k = params[17]
    c_asterisk = params[18]
    c_k = params[19]
    rho_asterisk = params[20]
    counsel = params[21]
    ##### New parameters to be added to posteriors model 
    transfer_2ndline = params[22]
    f1_rate = params[23]
    f2_rate = params[24]
    f3_rate = params[25]
    f4_rate = params[26]
    ##### Last two parameters to be added for ltfu
    mu1 = params[27]
    ##### Provide the calibration factor into the model
    qt = params[28]
    year_OffInfs = params[29]


    ######### Reinitalize the mortality ratios in according to the RR instead of actual rates 
    delta_U = delta_U*1.0
    ##### there is no different death rate for Diagnosed group
    delta_T = delta_T*delta_U
    delta_F = delta_F*delta_U

    #### h2 - which is the RR of Dolutegravir rate of developing Drug Resistance, is redefined here as RR instead
    h2 = h1*h2
    #### treatment failure rate when there is background drug resistance (after the introduction of dolutegravir only). based on f2_rate
    f3_rate = f3_rate*f2_rate


    N = S + I_W + I_DR + D_W + D_DR + T_W1 + T_W2 + T_DR1 + T_DR2 + F_W1 + F_W2 + F_TDR1 +F_ADR1 + F_DR2

    ######## Update the force of infections expression  
    ### Force of infection for Wild-type variants 
    I_Wa =   (1 + beta(t,a =eta_1,b=eta_2,q = eta_3)+ beta_r(t,k = eta_4,b=eta_5,q = eta_6) )*beta_u*(  (I_W+ D_W) + 
                        beta_t*(T_W1+T_W2) +   beta_f*(F_W1+F_W2))/N*S   * switch_Infs(t,k=year_OffInfs)
        ### Force of infection for Drug-resistant type variants
    I_Ra =   theta*  (1 + beta(t,a =eta_1,b=eta_2,q = eta_3)+ beta_r(t,k = eta_4,b=eta_5,q = eta_6))*beta_u*(  (I_DR+D_DR) + 
                                        beta_t*(T_DR1+T_DR2)+   beta_f*(F_TDR1+F_ADR1+F_DR2))/N*S    * switch_Infs(t,k=year_OffInfs)

    #### This was the adults population increase rates
    ### background mortality in the modelling 
    ### when the turn-off infections is in place, we also need to turn off Population increase and let people die out from the system
    S_prime = rho_asterisk*S  * switch_Infs(t,k=year_OffInfs) - (I_Wa+I_Ra) - delta_B*S


       
    ### infected people
    I_W_prime = I_Wa -   b(t,b=b_asterisk,k=b_k)*I_W -   (delta_U+delta_B)*I_W
    I_DR_prime= I_Ra -   b(t,b=b_asterisk,k=b_k)*I_DR -   (delta_U+delta_B)*I_DR
    

    ### diagnosed people
    D_W_prime =   b(t,b=b_asterisk,k=b_k)*I_W  -   c1(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W -   (delta_U+delta_B)*D_W
    D_DR_prime =   b(t,b=b_asterisk,k=b_k)*I_DR  -   c2(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR -   (delta_U+delta_B)*D_DR
    


    ### Viral Load monitoring is in effect
    ### the rate of LTFU here - mu2 but it should be the same as mu at the k(t)
    T_W1_prime =   c_1(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W +   g_1(t,counsel_f = counsel,qt=qt)*(1-  mu1)*F_W1 - (  f1(t,failure_rate=f1_rate)+  (delta_T+delta_B))*T_W1 -   k(t)*(1-  mu1)*T_W1
    T_W2_prime =   c_1prime(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W + (1-  mu1)*(  g_1prime(t,counsel_f = counsel,qt=qt)*F_W1 +   g2(t,counsel_f = counsel,qt=qt)*F_W2) - (  f2(t,failure_rate=f2_rate)+  (delta_T+delta_B))*T_W2 +   k(t)*(1-  mu1)*T_W1


    T_DR1_prime =   c_2(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR - (  f3(t,newdolutegravir_rate=f3_rate)+  (delta_T+delta_B))*T_DR1
    T_DR2_prime =   c_2prime(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR + (1-  mu1)*(  g3(t,bg_drtest = transfer_2ndline,qt=qt)*(F_TDR1 + F_ADR1)+  g4(t,counsel_f = counsel,qt=qt)*F_DR2) - (  f4(t,failure_rate=f4_rate)+  (delta_T+delta_B))*T_DR2

    

    F_W1_prime =   f1(t,failure_rate=f1_rate)*T_W1 -   g1(t,counsel_f = counsel,qt=qt)*(1-  mu1)*F_W1 - (  (delta_F+delta_B)+  h1)*F_W1
    F_W2_prime =   f2(t,failure_rate=f2_rate)*T_W2 -   g2(t,counsel_f = counsel,qt=qt)*(1-  mu1)*F_W2 - (  (delta_F+delta_B)+  h2)*F_W2

    # F_DR1_prime =   f3(t,newdolutegravir_rate=f3_rate)*T_DR1 + (  h1*F_W1 +   h2*F_W2)-   g3(t,bg_drtest = transfer_2ndline,qt=qt)*(1-  mu1)*F_DR1 -   delta_F*F_DR1
    F_TDR1_prime =   f3(t,newdolutegravir_rate=f3_rate)*T_DR1 -   g3(t,bg_drtest = transfer_2ndline,qt=qt)*(1-  mu1)*F_TDR1 -   (delta_F+delta_B)*F_TDR1
    F_ADR1_prime =    (  h1*F_W1 +   h2*F_W2)-   g3(t,bg_drtest = transfer_2ndline,qt=qt)*(1-  mu1)*F_ADR1 -   (delta_F+delta_B)*F_ADR1


    F_DR2_prime =   f4(t,failure_rate=f4_rate)*T_DR2 -   g4(t,counsel_f = counsel,qt=qt)*(1-  mu1)*F_DR2 -   (delta_F+delta_B)*F_DR2




    incidence =  I_Wa + I_Ra
    #### only model mortality related to HIV/AIDS here
    deaths =  delta_U*(I_W+I_DR) + delta_U*(D_W+D_DR) + delta_T*(T_W1+T_DR1+T_W2+T_DR2) + delta_F*(F_W1+F_TDR1+F_ADR1+F_W2+F_DR2)


    total_PLHIV_prime = I_Wa + I_Ra +c_1(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W + c_1prime(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W +c_2(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR + c_2prime(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR + g_1(t,counsel_f = counsel,qt=qt)*(1-  mu1)*F_W1 + g_1prime(t,counsel_f = counsel,qt=qt)*(1-  mu1)*F_W1 -    (delta_U+delta_B)*(I_W+I_DR) - (delta_U+delta_B)*(D_W+D_DR) - (delta_T+delta_B)*(T_W1+T_W2+T_DR1+T_DR2)  - (delta_F+delta_B)*(F_W1+F_W2+F_TDR1+F_ADR1+F_DR2)- c1(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W -   c2(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR -   g1(t,counsel_f = counsel,qt=qt)*(1-  mu1)*F_W1


    incidence_Resist =   I_Ra


    ###### calculate the number of tests (DR and VL) being done in a year
    ##### pre-treatment DR test
    preDRtest = (1-  mu1) * c_2prime(t,b=c_asterisk,k=c_k)*(D_DR + D_W)
    ##### VL actual testing rates --> counselling effectiveness wouold decrease this rate, however since we want the actual number of tests, this is the actual rate. 
    g1VLrate_routine,g1VLrate_POC ,g1DRrate_routine, g1DRrate_POC = g1Testingrate(t)
    g2VLrate_routine,g2VLrate_POC ,g2DRrate_routine, g2DRrate_POC = g2Testingrate(t)
    g3VLrate_routine,g3VLrate_POC ,g3DRrate_routine, g3DRrate_POC = g3Testingrate(t)
    g4VLrate_routine,g4VLrate_POC ,g4DRrate_routine, g4DRrate_POC = g4Testingrate(t)

    ###### post-treatment DR test
    ###### Remove the counselling, in the counting of number of VF tests
    postDRtest_routine =  (1-  mu1) * (g1DRrate_routine*F_W1 + g2DRrate_routine*F_W2 + g3DRrate_routine *(F_ADR1 + F_TDR1)) 
    postDRtest_POC =  (1-  mu1) * (g1DRrate_POC*F_W1 + g2DRrate_POC*F_W2 + g3DRrate_POC *(F_ADR1 + F_TDR1))

    
     
    ##### VL test
    ##### VL test in routine is only done once for those who are VF 
    VLtest_routine =  (1-  mu1) * (g1VLrate_routine*(T_W1+ F_W1) + g2VLrate_routine*(T_W2+ F_W2) + g3VLrate_routine*(T_DR1+ F_TDR1 + F_ADR1) + g4VLrate_routine*(T_DR2+ F_DR2))
    ##### VL test POC is done twice for those who are virological failure in a year 
    VLtest_POC =  (1-  mu1) * (g1VLrate_POC*(T_W1+ 2*F_W1) + g2VLrate_POC*(T_W2+ 2*F_W2) + g3VLrate_POC*(T_DR1+ 2*F_TDR1 + 2*F_ADR1) + g4VLrate_POC*(T_DR2+ 2*F_DR2))




    ##### Compute the negative impacts of Drug resistance (Transmitted, Acquired and Mis-treated)
    newAcquiredDR =  (  h1*F_W1 +   h2*F_W2)
    newTransmittedDR =  (  f3(t,newdolutegravir_rate=f3_rate)*T_DR1)
    newMisTreatedDR  = c_2(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR 


    ##### Compute the new diagnoses and treatment initiations in the model 
    diagnoses = b(t,b=b_asterisk,k=b_k)* (I_W + I_DR)
    treat_inits = c_1(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W + c_1prime(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_W + c_2(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR + c_2prime(t,b=c_asterisk,k=c_k)*(1-  mu1)*D_DR

    #### Compute all mortality, including background/ all-cause mortality in the modelling
    #### which also include mortality from susceptible

    all_mortality = (delta_U + delta_B)*(I_W+I_DR + D_W+D_DR) + (delta_T + delta_B)*(T_W1+T_DR1+T_W2+T_DR2) + (delta_F + delta_B)*(F_W1+F_TDR1+F_ADR1+F_W2+F_DR2) + delta_B*S


    return [S_prime,I_W_prime,I_DR_prime,D_W_prime,D_DR_prime,T_W1_prime,T_W2_prime,T_DR1_prime,T_DR2_prime,F_W1_prime,F_W2_prime,F_TDR1_prime,F_ADR1_prime,F_DR2_prime,incidence,deaths,total_PLHIV_prime,incidence_Resist,preDRtest,postDRtest_routine,postDRtest_POC,VLtest_routine,VLtest_POC,newAcquiredDR,newTransmittedDR,newMisTreatedDR,diagnoses,treat_inits, all_mortality]




def compute_results(d):
    ##### This is the original compartment that was splitted into two compartments
    d["F_DR1"] = d["F_TDR1"] + d["F_ADR1"]

    d["I"] = (d["I_W"]+d["I_DR"])
    d["D"] = (d["D_W"]+d["D_DR"])
    d["T"] = (d["T_W1"]+d["T_DR1"]+d["T_W2"]+d["T_DR2"])
    d["F"] = (d["F_W1"]+d["F_DR1"]+d["F_W2"]+d["F_DR2"])
    d["Total"] = d["I"] + d["D"] + d["T"] + d["F"]
    d["Aware"] = d["D"] + d["T"] + d["F"]
    d["TreatTotal"] = d["T"] + d["F"]
    d["DResistance"] = d["I_DR"] + d["D_DR"] + d["T_DR1"] + d["T_DR2"] + d["F_DR1"] + d["F_DR2"]


    d["adult_incidence"] = d["incidence"]/(d["S"])
    d["adult_prevalence"] = d["Total"]/(d["S"]+d["Total"])


    return(d)



def ode_solver(params,initial_conditions,tspan,treatment="Normal"):
    S,I_W, I_DR, D_W, D_DR, T_W1,T_W2, T_DR1,T_DR2, F_W1,F_W2, F_TDR1,F_ADR1,F_DR2,incidence,deaths,total_PLHIV,incidence_Resist,preDRtest,postDRtest_routine,postDRtest_POC,VLtest_routine,VLtest_POC,newAcquiredDR,newTransmittedDR,newMisTreatedDR,diagnoses,treat_inits,all_mortality = initial_conditions    

    theta,beta_u,beta_t,beta_f, delta_U,delta_T,delta_F,delta_B,h1, h2,eta_1, eta_2,eta_3,eta_4,eta_5,eta_6,  b_asterisk,b_k, c_asterisk,c_k,rho_asterisk, counsel,transfer_2ndline,f1_rate,f2_rate,f3_rate,f4_rate, mu1,qt,year_OffInfs = params

    if treatment =="NoDolutegravir":
        #### Failure rates are resetted to be the same as status quo 
        f3_rate = 1.0
        f2_rate = f1_rate
        f4_rate = f1_rate
        h2 = h1

    result = odeint(ode_model, [S,I_W, I_DR, D_W, D_DR, T_W1,T_W2, T_DR1,T_DR2, F_W1,F_W2, F_TDR1,F_ADR1,F_DR2,incidence,deaths,total_PLHIV,incidence_Resist,preDRtest,postDRtest_routine,postDRtest_POC,VLtest_routine,VLtest_POC,newAcquiredDR,newTransmittedDR,newMisTreatedDR,diagnoses,treat_inits,all_mortality], tspan, args=tuple([[theta,beta_u,beta_t,beta_f, delta_U,delta_T,delta_F,delta_B,h1, h2,eta_1, eta_2,eta_3,eta_4,eta_5, eta_6, b_asterisk,b_k, c_asterisk,c_k, rho_asterisk, counsel,transfer_2ndline,f1_rate,f2_rate,f3_rate,f4_rate,mu1,qt,year_OffInfs]]))

    return result

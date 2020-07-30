import numpy as np
import pandas as pd
import math
from replace import replace


'''
    sensitiv_anal_vect takes an arrya of 0s and 1s in the order [amount, CD, 
    S_price, H2SO4_price, elec_cost, elec_sell_price, PV_output, shift]
    where amount is kg H2, S_price is price of sulfur per kg, H2SO4_price is
    price of sulfuric acid per g, elec_cost is the cost of buying grid
    electricity, elec_sell_price, is the value of selling unused electricity,
    PV_output is peak power of PV in kW/m2, and shift, if the shift applied
    to John Weidner's CD data. This function then returns a cell vector
    of values where the data to cpmare are arrays and everything else is a
    double. The graph also returns the standard values for the comparison
    array as the last two values in the cell
'''

def Def_anal(switches, Skarn, Ave_basalt, tornado):

    #   default current density of the electrolysers
    DPkW = 100
    CO2_Tax = 0.0000000000000000000000001
    mineCost = 5

    heatCost = 3
    maxSCM = 6 # 0.3 + 0.1286;
    eCost = 0.02
    r = 0.08
    CH = 1.25
    kWhr_kg = 15

    # Skarns are all over https://www.researchgate.net/profile/Ewan_Pelleter/publication/233801803_World_Skarn_Deposits_-_Skarns_of_Western_Europe/links/00b4952b00d0be7304000000/World-Skarn-Deposits-Skarns-of-Western-Europe.pdf
    # https://journals.lib.unb.ca/index.php/GC/article/view/3773
    if Skarn == 1:
        CaO_Frac_Rock = 0.3 # 0.30 fraction of CaO in the source rock
        Al2O3_Frac_Rock = 0.05 # 0.05 fraction of Al2O3 in the source rock
        SiO2_Frac_Rock = 0.35 # 0.35 fraction of SiO2 in the source rock
        FeO_Frac_Rock = 0.00000000000001 # 0.00000000000001
        Fe2O3_Frac_Rock = 0.25 # 0.25 fraction of Fe2O3 in the source rock
        MgO_Frac_Rock = 0.05 # 0.05 fraction of MgO in the source rock
    elif (Ave_basalt == 1):
        CaO_Frac_Rock = 0.26 # 0.30 fraction of CaO in the source rock      0.26
        Al2O3_Frac_Rock = 0.14 # 0.05 fraction of Al2O3 in the source rock  0.14
        SiO2_Frac_Rock = 0.55 # 0.35 fraction of SiO2 in the source rock    0.55
        FeO_Frac_Rock = 0.05 # 0.00000000000001                             0.05
        Fe2O3_Frac_Rock = 0.00 # 0.25 fraction of Fe2O3 in the source rock  0.00
        MgO_Frac_Rock = 0.00 # 0.05 fraction of MgO in the source rock      0.00

    else:
        CaO_Frac_Rock = 0.305 #  0.30 fraction of CaO in the source rock 0.25 0.155 0.105
        Al2O3_Frac_Rock = 0.05 # 0.05 fraction of Al2O3 in the source rock 0.05 0.1 0.14
        SiO2_Frac_Rock = 0.595 # 0.35 fraction of SiO2 in the source rock 0.65 0.6 0.6
        FeO_Frac_Rock = 0.05 # 0.00000000000001  0.05 .1 .11
        Fe2O3_Frac_Rock = 0.000000000000000000000001# 0.25 fraction of Fe2O3 in the source rock 0 0 0
        MgO_Frac_Rock = 0.00000000000000000000000001 # 0.05 fraction of MgO in the source rock o 0.045 0.045

    CO2int = 56.1 / 1000 # T CO2/GJ heat
    eCO2int = 3.8440e-04 # T CO2/kWhr
    S_Cost = 20 # $/T buying sulfur

    CF = 0.97 # Kiln Capacity Factor
    TPY = 1000000 * CF # Kiln Output clinker in tonnes per year
    SA_ratio = 0.65 # concentration H2SO4 from the electrolyzer

    Eff = 0.08 # heat to power efficiency for organic rankine cycle, if used
    Rev = -0.000000000000000000000000001
    PPT_SCM = 60 # $/T SCM

    PPT_F = 90
    PPT_Al = 300
    PPT_Agg = 15
    V = 0.56
    W = 2.305
    CD = 0.000000000000000000000000000001

    Al_eff = 0.9
    Fe_eff = 0.9
    SCM_eff = 0.9
    OPC_eff = 0.9
    Agg_eff = 0.9
    CapEx_Fac = 1

    switches = replace(switches)
    constants = np.array([DPkW, CO2_Tax, CH, mineCost, heatCost, maxSCM, eCost, kWhr_kg, r, CaO_Frac_Rock, Al2O3_Frac_Rock,
                       SiO2_Frac_Rock, FeO_Frac_Rock, Fe2O3_Frac_Rock, MgO_Frac_Rock, CO2int, eCO2int, PPT_SCM, S_Cost,
                       CF, TPY, SA_ratio, Eff, Rev, PPT_F, PPT_Al, PPT_Agg, V, W, CD, Al_eff, Fe_eff, SCM_eff, OPC_eff,
                       Agg_eff, CapEx_Fac])

##############################################################################################
    # turn constants that are to be compared into zeros, and make the array a cell so it can have both doubles and arrays in the array

    input_a = constants*switches
    sensitiv_anal_vect = input_a

    # Python Syntax: np.arange(start, stop, step)
    #                   vs
    # MATLAB Syntax: start: step: stop

    DPkW = np.arange(50, 1050, 50) # 50:50: 1000
    CO2_Tax = np.arange(-190, 110, 10) #-190:10: 100
    CH =  np.arange( 1, 4.05, 0.05) # 1:.05: 4
    mineCost =  np.arange(1, 11, 1) #1:1: 10
    heatCost = np.arange(1, 11, 1) # 1:1: 10
    maxSCM =  np.arange(0.4, 5.4, 0.1) # 0.4:0.1: 5.3
    eCost = np.arange(0.01, 0.11, 0.01) # 0.01:0.01: 0.1
    kWhr_kg = np.arange(5, 31, 1) # 5:1: 30
    r = np.arange(0.04, 0.31, 0.01) #0.04:0.01: 0.3


    CaO_Frac_Rock = np.arange(0.05, 0.41, 0.01) # 0.05:0.01: 0.4
    Al2O3_Frac_Rock = np.arange(0.05, 0.41, 0.01) # 0.05:0.01: 0.4
    SiO2_Frac_Rock = np.arange(0.05, 0.41, 0.01) # 0.05:0.01: 0.4
    FeO_Frac_Rock = np.arange(0.05, 0.41, 0.01) # 0.05:0.01: 0.4
    Fe2O3_Frac_Rock = np.arange(0.05, 0.41, 0.01) # 0.05:0.01: 0.4
    MgO_Frac_Rock = np.arange(0.05, 0.41, 0.01) # 0.05:0.01: 0.4

    CO2int = np.arange(0, ((200/1000) + 0.005), (5/1000)) # 0:(5 / 1000): (200 / 1000);
    eCO2int = np.arange(0, 10e-04, 1e-04) # 0:0.1e-04: 10e-04;
    S_Cost = np.arange(0, 210, 10) # 0:10: 200;
    CF = np.arange(0.15, 1.01, 0.01) # 0.15:0.01: 1;
    TPY = np.arange(1000000, 5400000, ((5000000 - 1000000) / 10)) # 1000000:(5000000 - 1000000) / 10: 5000000;
    SA_ratio = np.arange(0.1, 0.75, 0.05) #0.1:0.05: 0.7;
    Eff = np.arange(0.0, 0.31, 0.01) # 0.0:0.01: 0.3;
    Rev = np.arange(0, 210, 10) # 0:10: 200;
    PPT_SCM = np.arange(0, 110, 10) # 0:10: 100;
    PPT_F = np.arange(0, 101, 1) # 0:1: 100;
    PPT_Al = np.arange(0, 410, 10) # 0:10: 400;
    PPT_Agg = np.arange(0, 21, 1) # 0:1: 20;
    V = np.arange(0.2, 1.25, 0.05) # 0.2:0.05: 1.2;
    W = np.arange(0.1, 5.1, 0.1) # 0.1:0.1: 5;
    CD = np.arange(-0.4, 0.5, 0.1) # -0.4:0.1: 0.4;

    Al_eff = np.arange(0.5, 1.1, 0.1) # 0.5:0.1: 1;
    Fe_eff = np.arange(0.5, 1.1, 0.1) # 0.5:0.1: 1;
    SCM_eff = np.arange(0.5, 1.1, 0.1) # 0.5:0.1: 1;
    OPC_eff = np.arange(0.5, 1.0, 0.1) # 0.5:0.1: 0.9;
    Agg_eff = np.arange(0.5, 1.1, 0.1) # 0.5:0.1: 1;
    CapEx_Fac = np.arange(0.5, 2.1, 0.1) # 0.5:0.1: 2;

    variables = np.array([DPkW, CO2_Tax, CH, mineCost, heatCost, maxSCM, eCost, kWhr_kg, r, CaO_Frac_Rock,
                         Al2O3_Frac_Rock, SiO2_Frac_Rock, FeO_Frac_Rock, Fe2O3_Frac_Rock, MgO_Frac_Rock,
                         CO2int, eCO2int, PPT_SCM, S_Cost, CF, TPY, SA_ratio, Eff, Rev, PPT_F, PPT_Al,
                         PPT_Agg, V, W, CD, Al_eff, Fe_eff, SCM_eff, OPC_eff, Agg_eff, CapEx_Fac], dtype = object)

 

    
# put the arrays that are to be compared into the sensitive_anal_vec
    sens_analysis = [[i] for i in sensitiv_anal_vect]
    sensitiv_anal_vect = [i for i in sensitiv_anal_vect]
    check = sensitiv_anal_vect[0]



    for i in range(len(sensitiv_anal_vect)):




        if sensitiv_anal_vect[i] == 0:

    # add the standard values of the array to be compared to the end of the cell
            sensitiv_anal_vect = np.append(sensitiv_anal_vect, constants[i])


            sensitiv_anal_vect = [i for i in sensitiv_anal_vect]


            sensitiv_anal_vect[i] = variables[i]
            sens_analysis[i] = variables[i]




    if tornado == 1:
        sensitiv_anal_vect = variables

    return sensitiv_anal_vect, constants, sens_analysis

#Aircraft Design Calculator
#Weight Estimation
import scipy as sc
import math as m
import matplotlib.pyplot as plt
from sympy import symbols,Eq,exp,solve
from tabulate import tabulate

#================================================================================================================#
#Dictionary of Values
parameters = {
    "R"             :         1640420,                       
    "E"             :            3600,
    "V_Emax"        :          118.22,
    "n_p"           :            0.85,
    "sfc_cruise"    : 2.20*(10**(-7)),
    "sfc_loiter"    : 2.26*(10**(-7)),
    "aspect_ratio"  :               7.5,
    "swet_sref"     :             4.5,
    "crew"          :           180*2,
    "payload"       :             100
}

#==================================================================================================================#
#Calculation Functions

def weight_estimation():
    #Requriements-----------------------------------------------------------#
    R = parameters["R"]                                
    E = parameters["E"]	                            
    V_Emax = parameters["V_Emax"]                       

    #Emperical Data
    n_p = parameters["n_p"]                             
    c_cruise = parameters["sfc_cruise"]                                  
    c_loiter = parameters["sfc_loiter"]                         
    Aspect_Ratio = parameters["aspect_ratio"]
    Swet_Sref = parameters["swet_sref"]
    K_LD = 10.5
    AVGAS100LL= 6.01

    #Segment Assumptions

    #Mission Profile
    crew = parameters["crew"]
    payload = parameters["payload"]

    #-----------------------------------------------------------------------#

    #Weight Estimation: Returns Empty WEight, and Fuel Weight
    #Define Symbols
    W_0, W_1, W_2, W_3, W_4, W_5, W_6 = symbols('W_0 W_1 W_2 W_3 W_4 W_5 W_6')
    W_empty, W_fuel, W_payload, W_crew = symbols('W_empty W_fuel W_crew W_payload')

    #Lift Estimation
    L_D = K_LD*((Aspect_Ratio/(Swet_Sref))**(1/2))

    #Weight Summation Equation
    eq1 = Eq(W_0,W_crew + W_payload + W_empty + W_fuel)

    #Empty Weight Equation
    eq2 = Eq(W_empty/W_0,(1.543e-5)*W_0 + 0.57)

    #Fuel Weight Equation
    eq3 = Eq(W_6/W_0,(W_1/W_0)*(W_2/W_1)*(W_3/W_2)*(W_4/W_3)*(W_5/W_4)*(W_6/W_5)) #Fuel General Equation
    eq4 = Eq(W_fuel/W_0, 1.05*(1 - W_6/W_0))
    eq5 = Eq(W_1/W_0,0.98)
    eq6 = Eq(W_2/W_1,0.97)
    eq7 = Eq(W_3/W_2, exp((-R * c_cruise)/(n_p*L_D)))
    eq9 = Eq(W_4/W_3, exp((-E * c_loiter * V_Emax)/(0.866 * L_D * n_p)))
    eq8 = Eq(W_5/W_4,0.99)
    eq10 = Eq(W_6/W_5, 0.997)

    #Crew and Payload Weight Estimation
    eq11 = Eq(W_crew, crew)
    eq12 = Eq(W_payload, payload)
        
    #Solve
    solution = solve([
        eq1, eq2, eq3, eq4, eq5, eq6,
        eq7, eq8, eq9, eq10, eq11, eq12
    ], [W_0, W_empty, W_fuel, W_crew, W_payload, W_1, W_2, W_3, W_4, W_5, W_6],dict = True )[0]

    # second_solution = solution [1]
    solution = {str(k): float(v.evalf()) for k, v in solution.items()}
    # Disctionary for solution
    sol_dict = solution
    sol_dict["L_D-maxCruise"] = L_D
    sol_dict["L_D-maxLoiter"] = L_D*0.866
    
    ratio = {
        "W_empty/W_0":sol_dict["W_empty"] / sol_dict["W_0"],
        "W_fuel/W_0":  sol_dict["W_fuel"] / sol_dict["W_0"],
        "W_1/W_0":        sol_dict["W_1"] / sol_dict["W_0"],
        "W_2/W_1":        sol_dict["W_2"] / sol_dict["W_1"],
        "W_3/W_2":        sol_dict["W_3"] / sol_dict["W_2"],
        "W_4/W_3":        sol_dict["W_4"] / sol_dict["W_3"],
        "W_5/W_4":        sol_dict["W_5"] / sol_dict["W_4"],
        "W_6/W_5":        sol_dict["W_6"] / sol_dict["W_5"],
        "Fuel_Density_AVGAS100LL_lb/gal":        AVGAS100LL,
        "Fuel_Volume_gal":    sol_dict["W_fuel"]/AVGAS100LL,
        "W_6/W_0":        sol_dict["W_6"] / sol_dict["W_0"]
    }


    return sol_dict,ratio
def trade_study():
    #This function is for trade study
    #It will be used to see how the weight changes with different parameters
    print()

#===================================================================================================================#
#Print Final Findings
print("""
Parameters
      """)
parm_list=[]
units = ["ft","s","","1/ft","1/ft","","","lb","lb"] 
for [key,val],units in zip(parameters.items(),units):
    parm_list.append([key,val,units])
print(tabulate(parm_list, headers = ["Parameter","Value","Unit"]))

print(f"""
Weight Estimation
      """)
weight_list = []
weight_units = ["lb","lb","lb","lb","lb","lb","lb","lb","lb","lb","lb","",""]
sol_dict,ratio= weight_estimation()
for [key,val],weight_units in zip(sol_dict.items(),weight_units):
    weight_list.append([key,val,weight_units])
print(tabulate(weight_list, headers = ("Parameters","Values","Units")))

print(f"""
Weight Ratios
      """)
ratio_list = []
for key,val in ratio.items():
    ratio_list.append([key,val])
print(tabulate(ratio_list, headers = ("Parameters","Values")))


#Final Dictionary===================================================================================================#
weight_estimation_dict = {**sol_dict,**parameters,**ratio}
#===================================================================================================================#


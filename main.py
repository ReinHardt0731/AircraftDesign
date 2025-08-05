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
    "R": 1640420,                       
    "E": 3600,
    "V_Emax": 118.22,
    "n_p" : 0.85,
    "sfc_cruise": 2.20*(10**(-7)),
    "sfc_loiter": 2.26*(10**(-7)),
    "aspect_ratio": 7,
    "s_wet_s_reference": 4.5,
    "crew":180*2,
    "payload":100
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
    Swet_Sref = parameters["s_wet_s_reference"]
    K_LD = 9
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
    ], [W_0, W_empty, W_fuel, W_crew, W_payload, W_1, W_2, W_3, W_4, W_5, W_6])

    #Show Solution
    symbol_list = W_0, W_empty, W_fuel, W_crew, W_payload, W_1, W_2, W_3, W_4, W_5, W_6
    first_solution = solution [0]
    second_solution = solution [1]

    for var,val, in zip(symbol_list,first_solution,):
        print(f"""  {var}: {val.evalf():.4f} lb""")
    
    weight_answer = input(f"""Would you like to see the second solution? yes/no:
                          """)

    if weight_answer == str("yes"):
        print("Second Solution")
        for var,val in zip(symbol_list,second_solution):
            print(f"   {var}: {val.evalf():.4f} lb")
    
    else:
        print("")
    # Disctionary for solution
    sol_dict = {var:val.evalf()
                for var,val in zip(symbol_list,first_solution)}
    sol_dict["L_D-maxCruise"] = L_D
    sol_dict["L_D-maxLoiter"] = L_D*0.866
    
    print(f""" Other Parameters
L_Dmax Cruise : {(L_D):.4f}
L_Dmax Loiter : {(sol_dict["L_D-maxLoiter"]):.4f}
W_empty/W_0   : {(sol_dict[W_empty]/sol_dict[W_0]):.4f}
W_fuel/W_0    : {(sol_dict[W_fuel]/sol_dict[W_0]):.4f}
W_1/W_0       : {(sol_dict[W_1]/sol_dict[W_0]):.4f}
W_2/W_1       : {(sol_dict[W_2]/sol_dict[W_1]):.4f}
W_3/W_2       : {(sol_dict[W_3]/sol_dict[W_2]):.4f}
W_4/W_3       : {(sol_dict[W_4]/sol_dict[W_3]):.4f}
W_5/W_4       : {(sol_dict[W_5]/sol_dict[W_4]):.4f}
W_6/W_5       : {(sol_dict[W_6]/sol_dict[W_5]):.4f}
Fuel AVGAS 100LL = 6.01 lb/gal
Fuel_Volume   : {(sol_dict[W_fuel]/AVGAS100LL):.4f}
W_6/W_0       : {(sol_dict[W_6]/sol_dict[W_0]):.4f}
          
          """)

    return sol_dict



def trade_study():
    #This function is for trade study
    #It will be used to see how the weight changes with different parameters
    print()

#===================================================================================================================#
#Print Final Findings
print("""Parameters
      """) 
for key,val in parameters.items():
    print(f"   {key}: {round(val,9)}")


#===================================================================================================================#
#Final Output

print(f"""
Weight Estimation
      """)
sol_dict= weight_estimation()


print(sol_dict[symbols('W_0')]) #sample syntax of calling the variables


#===================================================================================================================#


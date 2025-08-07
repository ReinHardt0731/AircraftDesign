import math as m
import matplotlib.pyplot as plt
from sympy import symbols,Eq,exp,solve
import numpy as np
from sympy import pi
from tabulate import tabulate
from main import weight_estimation_dict

#================================================================#
weight_parm = weight_estimation_dict
#================================================================#
#Aifoils
airfoil_root = {
    'Airfoil': 'NACA3413', 
    'Cl': 1.4, 
    'Cd': 0.068, 
    'L/Dmax': 0.11551777400717757, 
    'Cm': 0.16, 
    't/c': 0.16663377446214842, 
    'AoA_margin': 0.17, 
    'Total Score': 0.8402, 
    'Rank': 1}

airfoil_mid = {
    'Airfoil': 'NACA3413', 
    'Cl': 1.5, 
    'Cd': 0.068, 
    'L/Dmax': 0.11551777400717757, 
    'Cm': 0.16, 
    't/c': 0.16663377446214842, 
    'AoA_margin': 0.17, 
    'Total Score': 0.8402, 
    'Rank': 1}

airfoil_tip = {
    'Airfoil': 'NACA3413', 
    'Cl': 1.7, 
    'Cd': 0.068, 
    'L/Dmax': 0.11551777400717757, 
    'Cm': 0.16, 
    't/c': 0.16663377446214842, 
    'AoA_margin': 0.17, 
    'Total Score': 0.8402, 
    'Rank': 1}

#-----------------
#Functions
def initial_drag():
    L_D = weight_parm["L_D-maxCruise"]
    Aspect_Ratio = weight_parm["aspect_ratio"]
    
    #------------------------------------------
    K_clean, K_to, K_landing,e_clean,eTo_flaps, eLanding_flaps= symbols(
        ['K_clean', 'K_to', 'K_landing','e_clean', 'eTo_flaps', 'eLanding_flaps'])
    
    c_fe, cd0, cd0_LG, cd0_HLD_TO, cd0_TO,cd0_flaps_L, cd0_landing = symbols([
        'c_fe', 'cd0', 'cd0_LG', 'cd0_HLD_TO', 'cd0_TO','cd0_flaps_L', 'cd0_landing'])
    
    cl_maxave, cl_maxhld, cl_max, cl_c, cl_flapTO, cl_rotation = symbols([
        'cl_maxave', 'cl_maxhld', 'cl_max', 'cl_c', 'cl_flapTO', 'cl_rotation'])  

    #Cd -------------------------------------------------
    eq1  = Eq(c_fe, 0.0055)
    eq2  = Eq(cd0,c_fe*(weight_parm["swet_sref"]))
    eq3  = Eq(cd0_LG,0.008)
    eq4  = Eq(cd0_HLD_TO, 0.005)
    eq5  = Eq(cd0_TO , cd0 + cd0_LG + cd0_HLD_TO)
    eq6  = Eq(cd0_flaps_L, 0.06)
    eq7  = Eq(cd0_landing, cd0 + cd0_LG + cd0_flaps_L)
    #Cl -------------------------------------------------
    eq8  = Eq(cl_maxave,np.average([airfoil_root['Cl'],airfoil_mid['Cl'],airfoil_tip['Cl']]))
    eq9  = Eq(cl_maxhld, 0.9)   #Utilize Plain Flap
    eq10 = Eq(cl_max, 0.9*(cl_maxave + cl_maxhld))
    eq11 = Eq(cl_c, 0.3)
    eq12 = Eq(cl_flapTO, 0.8)
    eq13 = Eq(cl_rotation, cl_max/(1.1**2)) 
    #k -------------------------------------------------
    eq14 = Eq(L_D,(1/(4*K_clean*cd0))**(1/2) )
    eq15 = Eq(K_clean,1/(pi*e_clean*Aspect_Ratio))
    eq16 = Eq(eTo_flaps, 0.75)
    eq17 = Eq(eLanding_flaps,0.70)
    eq18 = Eq(K_to,1/(pi*eTo_flaps*Aspect_Ratio))
    eq19 = Eq(K_landing,1/(pi*eLanding_flaps*Aspect_Ratio))

    solution = solve([eq1,eq2,eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, 
                      eq13, eq14, eq15, eq16, eq17, eq18, eq19],[
            K_clean, K_to, K_landing,e_clean,eTo_flaps, eLanding_flaps,
             c_fe, cd0, cd0_LG, cd0_HLD_TO, cd0_TO,cd0_flaps_L, cd0_landing,
               cl_maxave, cl_maxhld, cl_max, cl_c, cl_flapTO, cl_rotation], dict = True)[0]
    
    solution = {str(key): float(value.evalf()) for key, value in solution.items()}

    return solution

def list(parasitic_drag,k_constant):
    cl_list = []
    cd_list = []
    start = -2
    end = 2 
    step_size = 0.1
    for cl in np.arange(start,end + step_size,step_size):
        cd = parasitic_drag + k_constant*(cl**2)
        cl_list.append(cl)
        cd_list.append(cd)
    
    return cl_list,cd_list


#======================================================================
solution= initial_drag()

steady = f"Cd_clean = {solution['cd0']} + {solution['K_clean']}Cl^2"
take_off = f"Cd_takeoff = {solution['cd0_TO']} + {solution['K_to']}Cl^2"
landing = f"Cd_clean = {solution['cd0_landing']} + {solution['K_landing']}Cl^2"

print(f"""
Coefficients of Cd,Cl and k values
""")
coeff_list = []
for [key, val] in solution.items():
    coeff_list.append([key,val])
print(tabulate(coeff_list, headers = ("Parameters", "Value"),))

print(f"""
Drag Polars 
      """)
print(steady)
print(take_off)
print(landing)

cl_steady, cd_steady = list(solution['cd0'],solution['K_clean'])
cl_takeoff , cd_takeoff = list(solution['cd0_TO'],solution['K_to'])
cl_landing, cd_landing = list(solution['cd0_landing'],solution['K_landing'])

#=====================================================================
#Warnings
if solution['eLanding_flaps'] < solution['eTo_flaps'] < solution['e_clean']:
    print()
else:
    print("")
    print("ERROR: Invalid Osswald Estimation")

#=====================================================================
#Plot
plt.style.use('classic')
plt.figure(figsize=(12, 8))
plt.plot(cd_steady, cl_steady, label = "Cd_clean", marker ="o")
plt.plot(cd_takeoff,cl_takeoff,label = "Cd_takeoff",marker = "o")
plt.plot(cd_landing,cl_landing, label =  "Cd_landing", marker = "o")

plt.grid(
    which='major',          
    color='gray',           
    linestyle='--',         
    linewidth=0.5,         
    alpha=0.7              
)

plt.xlabel("Coefficient of Drag")
plt.ylabel("Coefficient of Lift")
plt.title("Drag Polar",fontdict={'family':'rockwell', 'weight': 'bold',}, fontsize = 20)
plt.legend(loc = 'best')
plt.show()

#Final Dictionary=============================================================#
initial_drag_parms = solution
#=============================================================================#
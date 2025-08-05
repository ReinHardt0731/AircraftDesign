import math as m
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd
import numpy as np
#==============================================================================================#
#Dictionary of values
param = {
    "altitude" : 15000,
    "Wing Area" : 37.8,
    "Gross Weight" : 960,
    "C_d0" : 0.02,
    "k" : 0.062
}

#==============================================================================================#
#Functions

#---Altitude---
def rho(altitude):
    t_0 = 519
    a = -0.003566
    rho_0 = 0.002378
    p_0 = 2116.2

    if altitude <= 36000:
        temperature = t_0 + (altitude * a)
        pressure = p_0 * ((1 + (a * altitude) / t_0) ** 5.26)
        density = rho_0 * ((1 + (a * altitude) / t_0) ** 4.26)
    elif 36000 < altitude <= 82000:
        temperature = t_0 + (36000 * a)
        pressure = p_0 * (1.26 / (m.e ** ((4.805*10**-5) * altitude)))
        density = rho_0 * (1.68 / (m.e ** ((4.805*10**-5) * altitude)))
    else:
        print("36000 to 82000 ft only")

    return temperature, pressure, density

#---Dynamic Pressure---
def dynamic_pressure(density,velocity):
    return 0.5*(density)*(velocity**2)

#---Lift Coeffecient---
def c_l(q, S, W):
    return W/(q*S)

#--Induced Drag
def c_di(k,c_lift):
    return k*(c_lift**2)

#---Drag Coeffecient---
def c_d(c_d0,c_di):
    return c_d0 + c_di

#==============================================================================================#
#Define Variables

altitude = param["altitude"]
S = param["Wing Area"]
W = param["Gross Weight"]
c_d0 = param ["C_d0"]
k = param["k"]
_,_,density= rho(altitude)
param["Density"] = density
#============================================================================================================#
#Solving Commands
c_lift_list = []
c_di_list = []
c_parasitic_list = []
thrust_required_list = []
velocity_list = []

start = 100
end = 700
time_step = 5

for velocity in range(start,end,time_step):
    q = dynamic_pressure(density,velocity)
    c_lift = c_l(q,S,W)
    c_induced = (c_di(k,c_lift)) * q * S
    c_parasitic = c_d0 * q * S
    thrust_required = c_parasitic + c_induced

    c_lift_list.append(c_lift)
    c_di_list.append(c_induced)
    c_parasitic_list.append(c_parasitic)
    thrust_required_list.append(thrust_required)
    velocity_list.append(velocity)

#==============================================================================================================#
#Print and CSV generation
data = list(zip(velocity_list,c_lift_list,c_di_list,c_parasitic_list,thrust_required_list))
header = ["Velocity ft/s", "C_L","Induced Thrust Required","Parasitic Thrust Required","Thrust Required lbf"]
dict = pd.DataFrame(data,columns = header)
print(tabulate(data, headers = header, floatfmt = ".4f", tablefmt = "fancy_grid"))

file_name = "Thrust_required.csv"

with open(file_name,"w") as f:
    f.write("Parameter, Value\n")
    for key,val in param.items():
        f.write(f"{key},{val} \n")
    f.write(f"C_d,{c_d0}+{k}c_l \n")
    f.write("\n")

dict.to_csv(file_name, mode = 'a', index = False, )

#===============================================================================================================#
#Plot

min_index = thrust_required_list.index(min(thrust_required_list))
V_min = velocity_list[min_index]
Tr_min = thrust_required_list[min_index]

# Split for piecewise plotting
left_V = velocity_list[:min_index+1]
left_Tr = thrust_required_list[:min_index+1]
right_V = velocity_list[min_index:]
right_Tr = thrust_required_list[min_index:]

#-------------------------------------------------------------------------------
# Plot the left and right branches
plt.figure(figsize=(10,6))
plt.plot(left_V, left_Tr, label='Region of Velocity\n Instability', color='blue', marker='+')
plt.plot(right_V, right_Tr, label='Region of Velocity\n Stability', color='red', marker='+')
plt.plot(velocity_list,c_di_list, label = 'Induced Thrust Required', marker = '.')
plt.plot(velocity_list,c_parasitic_list, label = 'Parasitic Thrust Required', marker = '.')

# Min point marker
plt.scatter([V_min], [Tr_min], color='green', label='Minimum Thrust Required', zorder=5)
plt.axvline(x=V_min, color='green', linestyle='--', linewidth=1)

# Annotate regions
plt.text(left_V[1], left_Tr[1] + 10, 'Region of Velocity\n Instability', color='blue')
plt.text(right_V[-15], right_Tr[-5] + 5, 'Region of Velocity\n Stability', color='red')

# Plot formatting
plt.xlabel("Velocity (ft/s)")
plt.ylabel("Thrust Required (lb)")
plt.title("Thrust Required Curve @ 15,000 ft", fontdict={'family':'rockwell', 'weight': 'bold'})
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#This only shows thrust required at only 1 alttitude, To get the values for other altitudes,
#make the Velocity loop a function then make another function for altitude

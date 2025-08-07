import subprocess
import csv
import math as m
from turtle import color
import matplotlib.pyplot as plt
import numpy
from tabulate import tabulate

#Flow Condition

#Sim Controls-----------------------------------------------------------------
panel = int(300)
reynolds_number = int(5153748)
mach_number = 0.3
number_itterations = int(200)
ncrit = 3
file_name = "airfoil_data_optimized.csv"
file_name_for_normalized = "Normalized_Data.csv"
file_name_for_scored = "Airfoil_Scores.csv"

#Airfoil Configuration Range-------------------------------------------------
camber_min = 1
camber_max = 3
camber_location_min= 3
camber_location_max = 5
thickness_min = 12
thickness_max = 15

#Objective Weights------------------------------------------------------------
cl_W = 0.16
cd_W = 0.17
LD_W = 0.17
cm_W = 0.16
tc_W = 0.17
AoAmarg_W = 0.17
W_Total = int(cl_W + cd_W + LD_W + cm_W + tc_W + AoAmarg_W)

#Ideal Values---------------------------------------------------------------
ideal_tc = 13
deviation_tc = 5

# Functions #====================================================================================#
def airfoil_simulation():
    for m in range(camber_min,camber_max + 1):
        for p in range(camber_location_min,camber_location_max + 1):
            for t in range(thickness_min,thickness_max + 1):
                airfoil_code= f"{m}{p}{t:02d}"
                airfoil_name = f"NACA{airfoil_code}"
                #Command Sequence-----------------------------------------------------------#               
                xfoil_commands = f"""
NACA {airfoil_code}
PPAR
N {panel}


OPER
VPAR
N {ncrit}

ITER {number_itterations}
VISC
{reynolds_number}
MACH {mach_number}
PACC
{airfoil_name}.txt

ASEQ 0 20 0.5
PACC

QUIT
"""
                #initiates the commands
                proc = subprocess.Popen(['xfoil.exe'],
                        stdin = subprocess.PIPE,
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE,
                        text = True)
                proc.communicate(xfoil_commands)
                #Stop the simulation if dCl/dalpha changes its sign
                #Reads the {airfoil nadme} and outputs a list of tuples
                with open(f"{airfoil_name}.txt",'r') as f:
                    lines = f.readlines()
                    data = []
                    stall_trigger = False
                    post_stall_error = False

                    for i in range(len(lines) - 1):
                        tokens_current = lines[i].strip().split()
                        tokens_next = lines[i+1].strip().split()
                        
                        if len(tokens_current) == 7 and len(tokens_next) == 7:
                            try:
                                alpha = float(tokens_current[0])
                                cl = float(tokens_current[1])
                                cd = float(tokens_current[2])
                                cdp = float(tokens_current[3])
                                cm = float(tokens_current[4])
                                L_D = cl/cd

                                next_cl = float(tokens_next[1])

                                #Stall and Post Stall Detection Logic to clean data
                                if stall_trigger == False and post_stall_error == False:
                                    data.append((alpha , cl , cd , cdp , cm, L_D))
                                    if cl >= next_cl:
                                        stall_trigger = True

                                elif stall_trigger == True and post_stall_error == False:
                                    data.append((alpha , cl , cd , cdp , cm, L_D))
                                    if cl < next_cl:
                                        post_stall_error = True
                        
                                else:
                                    break
                                            
                            except ValueError:
                                continue
                            
                #Getting the performance paramters           
                cl_max = max(data, key = lambda x: x[1])
                cd_min = min(data, key = lambda x: x[2])
                L_Dmax = max(data, key = lambda x: x[5])
                cm_min = min(data, key = lambda x: x[4])
                t_ratio = t
                stall_angle = cl_max[0]
                L_Dmaxangle = L_Dmax[0]
                delta_AoA = stall_angle - L_Dmaxangle

                final_data = [
                cl_max[1],   # maximum Cl
                cd_min[2],   # minimum Cd
                L_Dmax[5],   # maximum L/D
                cm_min[4],   # minimum Cm
                t_ratio,     # thickness
                delta_AoA    # AoA margin
                ]

                print(f"{airfoil_name}: {final_data}")
                #Store Data as csv
                store_data(airfoil_name,final_data)

                with open("Sim_Configuration.txt","a") as config:
                    config.write(f"{airfoil_name}: {final_data}\n")

    return airfoil_name,final_data

#---Store Data
def store_data(airfoil,parameters):
    with open(file_name,"a") as f:
        f.write(f"{airfoil},{parameters[0]},{parameters[1]},{parameters[2]},{parameters[3]},{parameters[4]},{parameters[5]}\n")

#---Read Data to find min and max
def find_min_max(header):
    with open(file_name,"r") as f:
        var_list  = []
        reader = csv.DictReader(f)
        for row in reader:
            try:
                var = float(row[header])
                var_list.append(var)
            except (ValueError,KeyError):
                continue

    var_min = min(var_list)
    var_max = max(var_list)

    data = [[header,var_min,var_max]]
    fieldname = ["Parameter","minimum value","maximum value"]
    print(tabulate(data,headers = fieldname, floatfmt = ".4f", tablefmt = "fancy_grid"))

    return var_min,var_max

#---Normalize

def normalize(head_min,head_max,var):
        var_norm = (var - head_min)/(head_max-head_min)
        return var_norm

def inverse_norm(head_min,head_max,var):
    var_norm = 1 - ((var - head_min)/(head_max-head_min))
    return var_norm

def ideal_norm(ideal_value,deviation,var):
    var_norm = m.e**-(((var-ideal_value)**2)/(2*(deviation**2)))
    return var_norm

#---Normalize
def read_and_write_normalized_file():
    with open(file_name,"r") as infile, open(file_name_for_normalized, "w", newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames = fieldnames)

        writer.writeheader()
        for row in reader:
            try:
                row["Airfoil"] = row["Airfoil"]
                row["Cl"] = normalize(cl_min, cl_max,float(row["Cl"]))
                row["Cd"] = inverse_norm(cd_min, cd_max,float(row["Cd"]))
                row["L/Dmax"] = normalize(L_Dmin, L_Dmax,float(row["L/Dmax"]))
                row["Cm"] = inverse_norm(cm_min, cm_max, float(row["Cm"]))
                row["t/c"] = ideal_norm(ideal_tc, deviation_tc,float(row["t/c"]))
                row["AoA_margin"] = normalize(AoAmarg_min, AoAmarg_max, float(row["AoA_margin"]))

                print(row)
                writer.writerow(row)
            except Exception as e:
                print(f"ERROR: {e}")
#--Score Calculate

#Print Sim Controls and Save #==========================================================================================#
if W_Total == 1:
    print()
else:
    print("ERROR: W_Total is not equal to 1")
    exit()

with open("Sim_Configuration.txt","w") as f:
    f.write(f"""
==============================================================================================
Airfoil Optmizer Created By:Rexus Bryan L. Gan          
==============================================================================================
Airfoil Optimization Data
-------------------------------------------------            
Sim Controls             
    Panels:                {panel}                         
    Reynolds No:           {reynolds_number}
    Mach No:               {mach_number}
    Ncrit:                 {ncrit}
    Itterations:           {number_itterations}
-------------------------------------------------
 Objective Weights
    Cl :                   {cl_W*100}%            
    CD :                   {cd_W*100}%
    Cm :                   {cm_W*100}%
    L/D:                   {LD_W*100}%
    t/c:                   {tc_W*100}%
    AoAmarg:               {AoAmarg_W*100}%
    W_Total:               {W_Total*100}%
-------------------------------------------------
Airfoil Configuration Range (NACA 4_digits [camber, camber_loc, thickness])
    camber_min :           {camber_min}
    camber_max :           {camber_max}
    camber_location_min :  {camber_location_min}
    camber_location_max :  {camber_location_max}
    thickness_min :        {thickness_min}
    thickness_max :        {thickness_max}
    Total Airfoils:        {(camber_max-camber_min + 1)*(camber_location_max-camber_location_min + 1)*(thickness_max-thickness_min + 1 )}
-------------------------------------------------
Ideal Thickness: {ideal_tc}
Deviation:       {deviation_tc}

==============================================================================================
Airfoil Data: [Airfoil_Name]:[Cl,Cd,L/Dmax,Cm,t/c,AoA_margin]
==============================================================================================
""")

with open("Sim_Configuration.txt","r") as f:
    for line in f:
        lines = line.strip()
        print(lines)
#========================================================================================#
#File Generation and Process
with open(file_name,"w") as f:
    f.write(f"Airfoil,Cl,Cd,L/Dmax,Cm,t/c,AoA_margin\n")

#Function Calls #====================================================================================#
#Initiate simulation ; return those values as (airfoil_name: parameter lists); 
#Create an optmized parameter summary of the entire simulation range
airfoil,parameters = airfoil_simulation()
##Find Min and max for normalization--------------------

print(f"""
==============================================================================================
Find Min and Max
==============================================================================================
 """)
cl_min , cl_max = find_min_max("Cl")
cd_min, cd_max = find_min_max("Cd")
L_Dmin , L_Dmax = find_min_max("L/Dmax")
cm_min , cm_max = find_min_max("Cm")
AoAmarg_min, AoAmarg_max = find_min_max("AoA_margin")
print(f"""
==============================================================================================
Normalize Data
==============================================================================================
 """)

#Normalize the Datas and make it a new File--------------------
read_and_write_normalized_file()
#Make a new file for the weight application and Summation and sort---------------------

print(f"""
==============================================================================================
Done
==============================================================================================
 """)
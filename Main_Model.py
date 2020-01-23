"""
IDP 3 Group 9 - Tidal Lagoon Mathematical Model
Module: Main Program
Gianluca Cantone 2019-2020

Program Modules

 1. Main_Model.py - This file, runs simulation on the system to calculate how fast lagoon drains/fills and how much power will be generated.
 2. 3D_Model - 3D model of the bay (with visuals) plus integration function across whole model, function integrate and plot graph at all height levels.
 3. Grid_Element_Integration.py - Single element integration stand alone program (functions copied into 3D_Model.py)
 3. Optimizer - Runs simulator iterativly looking for best configuration
 
"""

#Libraries ====================================================================
import math as m
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#Variables ==================================================================== 
  
#Main Program =================================================================

def Run_Simulation():
   
    M = 1453217
    G = 9.81
    V_0 = 18891820
    Area = 40*5
    Time = [0]
    Time[0] = 0
    Volume = [0]
    Volume[0] = V_0
    
    V = M*np.power(( (np.sqrt((V_0)/(M)+5)) - ((Area*np.sqrt(2*G))/(2*M)) ), 2)-5*M
    print(V)
    
    for i in range(21600):
        Time.append(i)
        Volume.append(np.power(np.sqrt(V_0)-0.5*Area*i*np.sqrt((2*G)/(M)),2))
        #Volume.append(M*np.power(( (np.sqrt((V_0)/(M)+5)) - ((Area*i*np.sqrt(2*G))/(2*M)) ), 2)-5*M)
        
    plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes() #proj_type = 'ortho'
    plt.title("Volume Vs Time")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Volume (m^3)")
    ax.plot(Time, Volume)
    plt.minorticks_on()
    ax.grid(which='major', color='black', linestyle='-', linewidth=1)
    ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
    
    #consider implications of altering turbine diameter (alters height!)
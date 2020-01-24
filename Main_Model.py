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

    
def Run_Simulation(**kwargs):

    #Parameters passed:  
    
    Turbines = kwargs["turbines"]
    Turbine_Diameter = kwargs["diameter"]
    Step_Size = kwargs["step"]
    Sluices = kwargs["slucies"]
    Sluice_Size = kwargs["sluice_size"]
    
    Operation_Mode = kwargs["mode"]
    Run_Time = kwargs["time"]
    
    #Economic Parameters
    
    ###################
    
    V_0 = 18891820
    Euler_Volume = [0]
    Euler_Volume_Tide = [0]
    Euler_Volume[0] = V_0
    Euler_Volume_Tide[0] = V_0
    Time = [0]
    M = 1453217
    G = 9.807
    
    Area = np.pi*np.power((Turbine_Diameter/2),2)*Turbines
        
    for i in range(30000):
        
        Time.append(i)
        Euler_Volume.append(Euler_Volume[i]-Step_Size*Area*np.sqrt((2*G*Euler_Volume[i])/(M)))
        Euler_Volume_Tide.append(Euler_Volume_Tide[i]-Step_Size*Area*np.sqrt(2*G)*np.sqrt((np.abs(Euler_Volume_Tide[i]))/(M)-(2*np.sin(2*np.pi*0.0000463*i+1)+5)))
            
    #print(Euler_Volume_Tide)
    
    plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes() #proj_type = 'ortho'
    plt.title("Forward Euler")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Volume (m^3)")
    ax.plot(Time, Euler_Volume)
    ax.plot(Time, Euler_Volume_Tide)
    ax.legend(["No tide", "With tide"])
    plt.minorticks_on()
    ax.grid(which='major', color='black', linestyle='-', linewidth=1)
    ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
    
    
def Heads_Graph():          #Shows a graph of lagoon height and tidal height against time parameter. Tidal height is built into program, volume comes from analyitic function call.
    
    Volume, Time = Analyitic_Simulation_Simple()
    print(Volume)
    Tide_Heights = [0]
    M = 1453217
    
    for i in range(len(Volume)):
        Volume[i] = Volume[i]/M
        Tide_Heights.append(5*np.sin(2*np.pi*0.0000463*(i+1)+1)+5)
    
    print("Time: " + str(len(Time)))
    print("Volume: " + str(len(Volume)))
    print("Tide height: " + str(len(Tide_Heights)))
    plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes() #proj_type = 'ortho'
    plt.title("Lagoon Elevation Head And Tidal Height")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Head (m)")
    ax.plot(Time, Volume)
    Time.append(Time[len(Time)-1]+1)
    ax.plot(Time, Tide_Heights)
    plt.minorticks_on()
    ax.grid(which='major', color='black', linestyle='-', linewidth=1)
    ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
   
     
def Analyitic_Simulation_Simple(Duration=21600):        #Calculates the volume and graphs it over a time parameter analyitically.
   
    M = 1453217
    G = 9.807
    V_0 = 18891820
    Area = 40*5
    Time = [0]
    Volume = [0]
    Volume[0] = V_0
 
    for i in range(Duration):
        Time.append(i)
        Volume.append(np.power(np.sqrt(V_0)-0.5*Area*i*np.sqrt((2*G)/(M)),2))
    
    plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes() #proj_type = 'ortho'
    plt.title("Volume Vs Time")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Volume (m^3)")
    ax.plot(Time, Volume)
    plt.minorticks_on()
    ax.grid(which='major', color='black', linestyle='-', linewidth=1)
    ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)     
    
    return Volume, Time
        
        
def test(**kwargs):
    print(kwargs)
    print("\n")
    print(type(kwargs))
    name = kwargs["name"]
    print(name)    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
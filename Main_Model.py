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
import datetime

#Global Variables =============================================================

Global_Time = [0]
Global_Volume = [0]  
Global_Tide = [0]
Global_Head = [0] 
Profile_List = []

Civil_Price = 0
Blade_Price = 0
Turbine_Price = 0
Gearbox_Price = 0
Generator_Price = 0
Sluice_Price = 0

#Main Program =================================================================
print("\n")
print("======================================================================")
print("IDP 3 Group 9: Tidal Lagoon Mathematical Model and Simulation Software")
print("======================================================================")
print("This software simulates the performance of our tidal lagoon system under different configurations and operating conditions. Type help() to see a list of commands and how to use them")
    
def Run_Simulation(**kwargs):

    #Parameters passed:  
    
    Step_Size = kwargs["step"]
    Turbines = kwargs["turbines"]
    Turbine_Diameter = kwargs["diameter"]
    Sluices = kwargs["slucies"]
    Sluice_Size = kwargs["sluice_size"]
    Profile_Number = kwargs["profile"]
    Run_Time = kwargs["time"]
    Econ = kwargs["econ"]
    Output = kwargs["output"]
    Graphs = kwargs["graphs"]
    
    #Efficiency terms:
    
    Eff_Turbine = 0
    Eff_Gearbox = 0
    Eff_Generator = 0
    
    #Global Variables
    
    global Civil_Price, Blade_Price, Turbine_Price, Gearbox_Price, Generator_Price, Sluice_Price
    global Global_Time, Global_Volume, Global_Tide, Global_Head
    Global_Time = [0]
    Global_Volume = [0]  
    Global_Tide = [0]
    Global_Head = [0]
    
    #Local Variables
    
    V_0 = 18891820
    Euler_Volume = [0]
    Euler_Volume_Tide = [0]
    Euler_Volume[0] = V_0
    Euler_Volume_Tide[0] = V_0
    Time = [0]
    M = 1453217
    G = 9.807
    Tidal_Function = 0
    Discharge_Coefficient= 0.65
    Sluicing_Discharge_Coefficient = 0.95
    State = 0                   #0 = Waiting, 1 = filling sluice, 2 = draining generation, 3 = filling generation, 4 = draining sluice
    Current_Time = 1
    Operational_Profile = Profile_List[Profile_Number-1]
    Profile_Stage = 0
    Total_Stages = len(Operational_Profile)
    #Initial Calculations
    
    Area = np.pi*np.power((Turbine_Diameter/2),2)*Turbines
    Sluice_Area =  (np.pi*np.power((Turbine_Diameter/2),2)*Turbines)+(Sluice_Size*Sluices)
    Pipe_Loss = 0.0
    Turbine_Loss = 0
    
    '''
    Read operational algorithm profile, then switch between 5 different states, checking stopping conditions each time. Put Euler approximation in different function
    Operational states:
        -Filling_Sluice
        -Filling_Generation
        -Draining_Sluice
        -Draining_Generation
        -Transition (waiting)
    Save data to global arrays
    Plot graph using given parameters
    If value in square root turns negative, return an error. (This should only happen if sluice gates are not shut at the right time, it means the tide is higher than the lagoon)
    '''
    
    #State = operational profile.startstate
    State = Operational_Profile[0][0]
    
    while Current_Time < Run_Time:
        
        if State == 0:
            
            print("In state 0: Waiting for tidal shift")
            
            while State == 0:
                
                if Current_Time > Run_Time:
                    State = 5
                    print("Runtime elapsed")
                elif Global_Tide[Current_Time-1] < Operational_Profile[Profile_Stage][1]: 
                    Profile_Stage = (Profile_Stage+1)%Total_Stages
                    State = Operational_Profile[Profile_Stage][0]
                else:
                    Global_Time.append(Current_Time)
                    Global_Tide.append(6*np.cos(2*np.pi*0.0000231*Current_Time-np.pi)+6)
                    Global_Volume.append(Global_Volume[-1])
                    Global_Head.append((Global_Volume[-1])/(M))
                    Current_Time += 1
            
        if State == 1:
            
            print("In state 1: Filling lagoon via sluicing")
            
            while State == 1:
                
                if Current_Time > Run_Time:
                    State = 5
                    print("Runtime elapsed")
                elif Global_Head[Current_Time-1] > Operational_Profile[Profile_Stage][1]:
                    Profile_Stage = (Profile_Stage+1)%Total_Stages
                    State = Operational_Profile[Profile_Stage][0]  
                else:
                
                    Global_Time.append(Current_Time)
                    Global_Tide.append(6*np.cos(2*np.pi*0.0000231*Current_Time-np.pi)+6)
                    #Global_Volume.append(Global_Volume[Current_Time-1]+Step_Size*Sluice_Area*Sluicing_Discharge_Coefficient*np.sqrt(2*G)*np.sqrt(Global_Tide[Current_Time]-((Global_Volume[Current_Time-1])/(M))-Pipe_Loss))
                    
                    if (Global_Tide[Current_Time])-((Global_Volume[Current_Time-1])/(M)) < 0:
                        print("WARNING: Target value of " + str(Operational_Profile[Profile_Stage][1]) + "m in state 1 could not be met!")
                        Global_Volume.append(Global_Volume[-1])
                        Profile_Stage = (Profile_Stage+1)%Total_Stages
                        State = Operational_Profile[Profile_Stage][0]
                    else:    
                        Global_Volume.append(Global_Volume[Current_Time-1]+Step_Size*Sluice_Area*Sluicing_Discharge_Coefficient*np.sqrt(2*G)*np.sqrt((Global_Tide[Current_Time])-((Global_Volume[Current_Time-1])/(M))-Pipe_Loss))
                    
                    Global_Head.append((Global_Volume[Current_Time])/(M))
                    Current_Time += 1
                    #print("the time is: " + str(Current_Time))
        
        if State == 2:
            
            print("In state 2: Draining lagoon via energy generation")
            
            while State == 2:
                
                if Current_Time > Run_Time:
                    State = 5
                    print("Runtime elapsed")
                elif Global_Head[Current_Time-1] < Operational_Profile[Profile_Stage][1]: 
                    Profile_Stage = (Profile_Stage+1)%Total_Stages
                    State = Operational_Profile[Profile_Stage][0]
                    #Current_Time = Run_Time+1
                else:
                    
                    Global_Time.append(Current_Time)
                    Global_Tide.append(6*np.cos(2*np.pi*0.0000231*Current_Time-np.pi)+6)
                    
                    if (Global_Volume[Current_Time-1])/(M)-(Global_Tide[Current_Time]) < 0:
                        print("WARNING: Target value of " + str(Operational_Profile[Profile_Stage][1]) + "m in state 2 could not be met!")
                        Global_Volume.append(Global_Volume[-1])
                        Profile_Stage = (Profile_Stage+1)%Total_Stages
                        State = Operational_Profile[Profile_Stage][0]
                    else:
                        Global_Volume.append(Global_Volume[Current_Time-1]-Step_Size*Area*Discharge_Coefficient*np.sqrt(2*G)*np.sqrt((Global_Volume[Current_Time-1])/(M)-(Global_Tide[Current_Time])-Turbine_Loss-Pipe_Loss))
                    
                    Global_Head.append((Global_Volume[Current_Time])/(M))
                    Current_Time += 1
            
        if State == 3:
            
            print("In state 3: Filling lagoon via energy generation")
            #Check exit condition
            #Break or increment tide and update lagoon volume
            #NEED EQUATION FOR THIS
            
        if State == 4:
            
            print("In state 4: Draining lagoon via sluicing")
            #Check exit condition
            #Break or increment tide and update lagoon volume
            #NEED EQUATION FOR THIS
        print("Current time: " + str(Current_Time))    
    
    #Energy generation calculations
    
    #Log headloss across turbine at every step and store in an array
    #Multiply each value by thing to get energy/power
    #Make graph of power over time and total cumlative energy generated
    #Print peak power and total energy
    #Pass energy through inefficiencies of electrical system
    #Print final power and energy and percentage lost
    
    
    if Econ == True:                                                            #Section for calculating economic assessment.
        
        print("Running economic assessment of configuration...")
    #Fixed costs - Start up costs, running costs, etc
    #Variable costs
    #Take into account maintencance downtimes
    
    if Graphs == True:
    
        plt.figure(figsize=plt.figaspect(1)*2)
        
        ax = plt.axes()
        second_ax = ax.twinx()
        plt.title("Forward Euler")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Volume (m^3)")
        second_ax.set_ylabel('Head (m)')
       
        Filling_Plot = ax.plot(Global_Time, Global_Volume, label="Lagoon Volume", color="deepskyblue", linewidth=2)
        Lagoon_Head_Plot = second_ax.plot(Global_Time, Global_Head, "--", label="Lagoon head", color="green")
        Tide_Height_Plot = second_ax.plot(Global_Time, Global_Tide, "--", label="Tide height", color ="blue")
        
        Lines = Filling_Plot+Lagoon_Head_Plot+Tide_Height_Plot
        Labels =[l.get_label() for l in Lines]
        ax.legend(Lines, Labels)
        
        plt.minorticks_on()
        ax.grid(which='major', color='black', linestyle='-', linewidth=1)
        ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
        ax.set_ylim(0,2e7)
    

def Print_Costs():
    print("Here are some costs")
    
def Set_Costs(**kwargs):
    print("Setting costs")
    
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
   
def Tidal_Function_Testing(Interval=86400):     #Shows a graph of a desired tidal function hardcoded in. The time interval is taken as the parameter.
     
    xlim = np.arange(0, 60* 60 *26, 60*60*2)
    Time = [0]
    Tide_Height = [0]
    
    for i in range(Interval):
        
        Time.append(i)
        Tide_Height.append(6*np.cos(2*np.pi*0.0000231*i-np.pi)+6)
    
    plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes() #proj_type = 'ortho'
    plt.title("Ideal Semidiurnal Tide Over 24 Hours")
    ax.set_xlabel("Time (H)")
    ax.set_ylabel("Head (m)")
    ax.plot(Time, Tide_Height)
    plt.minorticks_on()
    plt.xticks(xlim, [str(n).zfill(2) + ':00' for n in np.arange(0, 26,2)])
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
        
        
def Setup_Profile(states=[]):
    
    #Format [[State number, triggering lagoon head], ...]
    print("Adding profile to list of saved profiles. This is profile number " + str(len(Profile_List)+1))
    Profile_List.append(states)

def help(*args):
    
    if len(args) == 0:
    
        print("Help:")
        print("\n1. To run a simulation, you must first setup an operating profile to tell the lagoon how it should operate. To learn how to do this, type help('profile'). \
               Alternativly, you can use the defaul profile already built into the program.")
        print("\n2. To learn how to run a simulation, type help('sim').")
    
    elif args[0] == "profile":
        
        print("Help profile:")
        print("\nHere is some help info on profiles")
        
    elif args[0] == "sim":
        
        print("Help simulation:")
        print("\nHere is some help info on running a simulation")
        
    else:
        
        print("Help command unknown, are you sure you typed it correctly?")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
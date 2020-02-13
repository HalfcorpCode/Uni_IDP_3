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
import logging

#Global Variables =============================================================

Global_Time = [0]
Global_Volume = [0]  
Global_Tide = [0]
Global_Head = [0]
Global_Head_Difference = [0]
Global_Velocity = [0]
Global_Discharge = [0] 
Profile_List = []
Global_Head_Loss = [0]
Global_Power = [0]
Global_Power_Elec = []

Civil_Price = 60000000        #Guess
Blade_Price = 1500000        #Good    
Turbine_Price = 2000000        #(Just guide vanes atm)
Gearbox_Price = 150000       #Good
Generator_Price = 150000     #Check but should be good
Electrical_Price = 100000    #Guess
Sluice_Price = 200000        #Good

#Main Program =================================================================

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger()
log.setLevel(logging.WARNING)

print("\n")
print("======================================================================")
print("IDP 3 Group 9: Tidal Lagoon Mathematical Model and Simulation Software")
print("======================================================================")
print("This software simulates the performance of our tidal lagoon system under different configurations and operating conditions. Type help() to see a list of commands and how to use them")
    
def Run_Simulation(**kwargs):

    #Parameters passed:  
    
    Step_Size = kwargs["step"]
    Tidal_Function = kwargs["tidal_function"]
    Turbines = kwargs["turbines"]
    Turbine_Diameter = kwargs["diameter"]
    Sluices = kwargs["slucies"]
    Sluice_Size = kwargs["sluice_size"]
    Profile_Number = kwargs["profile"]
    Run_Time = kwargs["time"]
    Econ = kwargs["econ"]
    Output = kwargs["output"]
    Graphs = kwargs["graphs"]
    Graph_Head = kwargs["graph_head"]
    Graph_QV = kwargs["graph_QV"]
    Graph_P = kwargs["graph_P"]
    
    #Efficiency terms:
    
    Eff_Turbine = 0.85
    Eff_Gearbox = 0.85
    Eff_Generator = 0.97
    
    #Global Variables
    
    global Civil_Price, Blade_Price, Turbine_Price, Gearbox_Price, Generator_Price, Sluice_Price
    global Global_Time, Global_Volume, Global_Tide, Global_Head, Global_Head_Difference, Global_Velocity, Global_Discharge, Global_Head_Loss, Global_Power, Global_Power_Elec
    Global_Time = [0]
    Global_Volume = [0]  
    Global_Tide = [0]
    Global_Head = [0]
    Global_Head_Difference = [0]
    Global_Velocity = [0]
    Global_Discharge = [0]
    Global_Head_Loss = [0]
    Global_Power = [0]
    Global_Power_Elec = []
    
    #Local Variables
    
    V_0 = 18891820
    Euler_Volume = [0]
    Euler_Volume_Tide = [0]
    Euler_Volume[0] = V_0
    Euler_Volume_Tide[0] = V_0
    Time = [0]
    M = 1453217
    G = 9.807
    rho = 1029
    Discharge_Coefficient= 0.8
    Sluicing_Discharge_Coefficient = 0.98
    Friction_Factor = 0.00035
    Draft_Length = 10
    Draft_Diameter = 10
    State = 0                   #0 = Waiting, 1 = filling sluice, 2 = draining generation, 3 = filling generation, 4 = draining sluice
    Current_Time = Step_Size
    Operational_Profile = Profile_List[Profile_Number-1]
    Profile_Stage = 0
    Total_Stages = len(Operational_Profile)
    Total_Mechanical_Energy = 0
    Total_Electrical_Energy = 0
    
    #Initial Calculations
    
    Area = np.pi*np.power((Turbine_Diameter/4),2)*Turbines
    Sluice_Area =  (np.pi*np.power((Turbine_Diameter/4),2)*Turbines)+(Sluice_Size*Sluices)
    
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
    
    if Output == True:
        log.setLevel(logging.INFO)
    else:
        log.setLevel(logging.WARNING)
    
    print("\nRunning simulation...\n")
    
    State = Operational_Profile[0][0]
    
    while Current_Time < Run_Time:
        
        if State == 0:
            
            log.info("In state 0: Waiting for tidal shift")
            
            while State == 0:
                
                if Current_Time > Run_Time:
                    State = 5
                    log.info("Runtime elapsed")
                elif Global_Tide[-1] < Operational_Profile[Profile_Stage][1]: 
                    Profile_Stage = (Profile_Stage+1)%Total_Stages
                    State = Operational_Profile[Profile_Stage][0]
                else:
                    Global_Time.append(Current_Time)
                    Global_Tide.append(Evaluate_Tidal_Function(Tidal_Function, Current_Time))
                    Global_Volume.append(Global_Volume[-1])
                    Global_Head.append((Global_Volume[-1])/(M))
                    Global_Head_Difference.append(Global_Head[-1]-Global_Tide[-1])
                    Global_Velocity.append(0)
                    Global_Discharge.append(0)
                    Global_Head_Loss.append(0)
                    Global_Power.append(0)
                    Current_Time += Step_Size
            
        if State == 1:
            
            log.info("In state 1: Filling lagoon via sluicing")
            
            while State == 1:
                
                if Current_Time > Run_Time:
                    State = 5
                    log.info("Runtime elapsed")
                elif Global_Head[-1] > Operational_Profile[Profile_Stage][1]:
                    Profile_Stage = (Profile_Stage+1)%Total_Stages
                    State = Operational_Profile[Profile_Stage][0]  
                else:
                
                    Global_Time.append(Current_Time)
                    Global_Tide.append(Evaluate_Tidal_Function(Tidal_Function, Current_Time))
                    #Global_Volume.append(Global_Volume[Current_Time-1]+Step_Size*Sluice_Area*Sluicing_Discharge_Coefficient*np.sqrt(2*G)*np.sqrt(Global_Tide[Current_Time]-((Global_Volume[Current_Time-1])/(M))-Pipe_Loss))
                    
                    if (Global_Tide[-1])-((Global_Volume[-1])/(M)) < 0:
                        log.info("WARNING: Target value of " + str(Operational_Profile[Profile_Stage][1]) + "m in state 1 could not be met!")
                        Global_Volume.append(Global_Volume[-1])
                        Profile_Stage = (Profile_Stage+1)%Total_Stages
                        State = Operational_Profile[Profile_Stage][0]
                        Global_Velocity.append(Global_Velocity[-1])
                        Global_Discharge.append(Global_Discharge[-1])
                    else:    
                        #Global_Volume.append(Global_Volume[-1]+Step_Size*Sluice_Area*Sluicing_Discharge_Coefficient*np.sqrt(2*G)*np.sqrt((Global_Tide[-1])-((Global_Volume[-1])/(M))-Pipe_Loss))
                        #Global_Velocity.append(np.sqrt(2*G)*np.sqrt((Global_Tide[-1])-((Global_Volume[-1])/(M))-Pipe_Loss))
                        
                        Global_Velocity.append(np.sqrt(2*G*((Global_Tide[-1])-((Global_Volume[-1])/(M)))/(1+Friction_Factor*(Draft_Length/Draft_Diameter))))
                        Global_Volume.append(Global_Volume[-1]+Step_Size*Sluice_Area*Sluicing_Discharge_Coefficient*Global_Velocity[-1])
                        
                        Global_Discharge.append(Global_Velocity[-1]*Sluice_Area*Sluicing_Discharge_Coefficient)
                        
                    Global_Head.append((Global_Volume[-1])/(M))
                    Global_Head_Difference.append(Global_Head[-1]-Global_Tide[-1])
                    Global_Head_Loss.append(0)
                    Global_Power.append(0)
                    Current_Time += Step_Size
                    #print("the time is: " + str(Current_Time))
        
        if State == 2:
            
            log.info("In state 2: Draining lagoon via energy generation")
            
            while State == 2:
                
                if Current_Time > Run_Time:
                    State = 5
                    log.info("Runtime elapsed")
                elif Global_Head[-1] < Operational_Profile[Profile_Stage][1]: 
                    Profile_Stage = (Profile_Stage+1)%Total_Stages
                    State = Operational_Profile[Profile_Stage][0]
                    #Current_Time = Run_Time+1
                else:
                    
                    Global_Time.append(Current_Time)
                    Global_Tide.append(Evaluate_Tidal_Function(Tidal_Function, Current_Time))
                    
                    if (Global_Volume[-1])/(M)-(Global_Tide[-1]) < 0:
                        log.info("WARNING: Target value of " + str(Operational_Profile[Profile_Stage][1]) + "m in state 2 could not be met!")
                        Global_Volume.append(Global_Volume[-1])
                        Profile_Stage = (Profile_Stage+1)%Total_Stages
                        State = Operational_Profile[Profile_Stage][0]
                        Global_Velocity.append(Global_Velocity[-1])
                        Global_Discharge.append(Global_Discharge[-1])
                        Global_Head_Loss.append(Global_Head_Loss[-1])
                        Global_Power.append(Global_Power[-1])
                    else:
                        #Global_Volume.append(Global_Volume[-1]-Step_Size*Area*Discharge_Coefficient*np.sqrt(2*G)*np.sqrt((Global_Volume[-1])/(M)-(Global_Tide[-1])-Turbine_Loss-Pipe_Loss))
                        #Global_Velocity.append(np.sqrt(2*G)*np.sqrt((Global_Volume[-1])/(M)-(Global_Tide[-1])-Turbine_Loss-Pipe_Loss))
                        
                        Global_Velocity.append(np.sqrt((2*G*((Global_Volume[-1])/(M)-(Global_Tide[-1]))*(1-Eff_Turbine))/(1-Eff_Turbine*Friction_Factor*(Draft_Length/Draft_Diameter))))
                        Global_Volume.append(Global_Volume[-1]-Step_Size*Area*Discharge_Coefficient*Global_Velocity[-1])
                        
                        Global_Discharge.append(Global_Velocity[-1]*Area*Discharge_Coefficient)
                        Global_Head_Loss.append(0.5*Global_Velocity[-1])
                        Global_Power.append((rho*G)*Global_Discharge[-1]*((Global_Volume[-1])/(M)-Global_Tide[-1]))
                    
                    Global_Head.append((Global_Volume[-1])/(M))
                    Global_Head_Difference.append(Global_Head[-1]-Global_Tide[-1])
                    Current_Time += Step_Size
            
        if State == 3:
            
            log.info("In state 3: Filling lagoon via energy generation")
            #Check exit condition
            #Break or increment tide and update lagoon volume
            #NEED EQUATION FOR THIS
            
        if State == 4:
            
            log.info("In state 4: Draining lagoon via sluicing")
            #Check exit condition
            #Break or increment tide and update lagoon volume
            #NEED EQUATION FOR THIS
        log.info("Current time: " + str(Current_Time))    
    
    print("\nSimulation complete\n")
    
    #Energy generation calculations
    print("\n================================================================")
    print("Running energy calculations...")
    
    for Power in Global_Power:
        
        if np.isnan(Power) == True:
            Power = 0
            
        Global_Power_Elec.append(Power*Eff_Turbine*Eff_Gearbox*Eff_Generator)
        Total_Mechanical_Energy += int(Power*Step_Size)
        Total_Electrical_Energy += int(Global_Power_Elec[-1]*Step_Size)

    Energy_Lost = Total_Mechanical_Energy-Total_Electrical_Energy
    Efficiency = (Eff_Turbine*Eff_Gearbox*Eff_Generator)
    print("Total mechanical energy generated: " + str(Total_Mechanical_Energy))
    print("Energy lost: " + str(Energy_Lost))
    print("System efficiency: " + str(Efficiency))
    print("Total electrical energy generated (J): " + str(Total_Electrical_Energy))
    print("Total electrical energy generated (kWh): " + str(Total_Electrical_Energy/3.6e+6))
    #Print peak power and total energy
    #Print final power and energy and percentage lost
    
    if Econ == True:                                                            #Section for calculating economic assessment.
        
        print("\n================================================================")
        print("Running economic assessment of configuration...")
        
        Startup_Costs = Civil_Price+(Turbines*(Blade_Price+Turbine_Price+Gearbox_Price+Generator_Price+Electrical_Price+(2*Sluice_Price)))+(Sluices*Sluice_Price)
        Running_Costs_Day = 823 #a day
        Total_Running_Costs = Running_Costs_Day*(Run_Time/86400)
        Energy_Price = 50/1000 #per kWh
        Turnover = (Total_Electrical_Energy/(3.6e+6))*Energy_Price
        Gross_Profit = Turnover-Total_Running_Costs
        Net_Profit = Gross_Profit-Startup_Costs
        Payback_Time = Startup_Costs/Gross_Profit
        
        Total_Costs_Gradient = (Total_Running_Costs)
        #Total_Costs_Equation = Total_Running_Costs*x + Startup_Costs
        Revenue_Gradient = Turnover
        #Revenue_Equation = Turnover*x
        
        print("Startup costs: " + str(Startup_Costs))
        print("Running costs: " + str(Total_Running_Costs))
        print("Turnover: " + str(Turnover))
        print("Gross profit: " + str(Gross_Profit))
        print("Net profit: " + str(Net_Profit))
        print("Pay back time in years: " + str(Payback_Time))
        
        plt.figure(figsize=plt.figaspect(1)*2)
        ax = plt.axes()
        plt.title("Break-Even Analysis")
        ax.set_xlabel("Time (Years)")
        ax.set_ylabel("Equity (£)")
            
        Fixed_Costs_Plot = ax.plot([0,Payback_Time,2*Payback_Time], [Startup_Costs, Startup_Costs, Startup_Costs], "--", label="Fixed costs", color="red", linewidth=3)
        Total_Costs_Plot = ax.plot([0,Payback_Time,2*Payback_Time], [Startup_Costs, (Total_Running_Costs*Payback_Time + Startup_Costs), (Total_Running_Costs*Payback_Time*2 + Startup_Costs)], "--", label="Total costs", color="darkred", linewidth=3)
        Revenue_Plot = ax.plot([0,Payback_Time,2*Payback_Time], [0, (Turnover*Payback_Time), (Turnover*Payback_Time*2)], label="Revenue", color="blue", linewidth=3)
        
        Lines = Fixed_Costs_Plot+Total_Costs_Plot+Revenue_Plot
        Labels =[l.get_label() for l in Lines]
        ax.legend(Lines, Labels)
        plt.minorticks_on()
        ax.grid(which='major', color='black', linestyle='-', linewidth=1)
        ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
        
    #Take into account maintencance downtimes
    
    if Graphs == True:
    
        plt.figure(figsize=plt.figaspect(1)*2)
        
        ax = plt.axes()
        second_ax = ax.twinx()
        plt.title("Lagoon Volume, Head and Tide Vs Time")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Volume (m^3)")
        second_ax.set_ylabel('Elevation (m)')
       
        Filling_Plot = ax.plot(Global_Time, Global_Volume, label="Lagoon Volume", color="deepskyblue", linewidth=3)
        Lagoon_Head_Plot = second_ax.plot(Global_Time, Global_Head, "--", label="Lagoon head", color="green", linewidth=2)
        Tide_Height_Plot = second_ax.plot(Global_Time, Global_Tide, "--", label="Tide height", color ="blue", linewidth=2)
        
        if Graph_Head == True:  
            Head_Difference_Plot = second_ax.plot(Global_Time, Global_Head_Difference, "--", label="Head difference", color ="orange", linewidth=2)
            Lines = Filling_Plot+Lagoon_Head_Plot+Tide_Height_Plot+Head_Difference_Plot
        else:
            Lines = Filling_Plot+Lagoon_Head_Plot+Tide_Height_Plot
            
        Labels =[l.get_label() for l in Lines]
        ax.legend(Lines, Labels)
        
        plt.minorticks_on()
        ax.grid(which='major', color='black', linestyle='-', linewidth=1)
        ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
        ax.set_ylim(0,2e7)
        
        if Graph_QV == True:
            
            #Red graphs
            
            plt.figure(figsize=plt.figaspect(1)*2)
            QVax = plt.axes()
            second_QVax = QVax.twinx()
            plt.title("Velocity & Discharge Vs Time")
            QVax.set_xlabel("Time (s)")
            QVax.set_ylabel("Velocity (ms^-1)")
            second_QVax.set_ylabel('Discharge (m^3s^-1)')
            
            Velocity_Plot = QVax.plot(Global_Time, Global_Velocity, label="Flow velocity", color="red", linewidth=2)
            Discharge_Plot = second_QVax.plot(Global_Time, Global_Discharge, "--", label="Total discharge", color="maroon", linewidth=2)
            Lines = Velocity_Plot+Discharge_Plot
            Labels =[l.get_label() for l in Lines]
            QVax.legend(Lines, Labels)
            plt.minorticks_on()
            QVax.grid(which='major', color='black', linestyle='-', linewidth=1)
            QVax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
            
        if Graph_P == True:
            
            plt.figure(figsize=plt.figaspect(1)*2)
            Pax = plt.axes()
            plt.title("Mechanical & Electrical Power Vs Time")
            Pax.set_xlabel("Time (s)")
            Pax.set_ylabel("Power (w)")
            
            Mech_Power_Plot = Pax.plot(Global_Time, Global_Power, label="Mechanical Power", color="gold", linewidth=2)
            Elec_Power_Plot = Pax.plot(Global_Time, Global_Power_Elec, "--", label="Electrical Power", color="y", linewidth=2)
            Lines = Mech_Power_Plot+Elec_Power_Plot
            Labels =[l.get_label() for l in Lines]
            Pax.legend(Lines, Labels)
            plt.minorticks_on()
            Pax.grid(which='major', color='black', linestyle='-', linewidth=1)
            Pax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
 
def Evaluate_Tidal_Function(Type, Time):
    
    if Type == "sine":
        return (6.5*np.cos(2*np.pi*0.0000231*Time-np.pi)+6.5)
    
    elif Type == "sine_average":
        return (4*np.cos(2*np.pi*0.0000231*Time-np.pi)+6.5)
    elif Type == "Newport_1":
        Time += 14000
        return ( (3.678*np.cos(2*np.pi*(2.236e-5)*Time+(48.06*(np.pi/180)))) + (1.4298*np.cos(2*np.pi*(2.239e-5)*Time-(130.6*(np.pi/180)))) + (1.4986*np.cos(2*np.pi*(2.315e-5)*Time+(176.7*(np.pi/180)))) + (0.7298*np.cos(2*np.pi*(2.194e-5)*Time-(53.75*(np.pi/180)))) + 6.278)
           
def Optimize(Item):

    global Global_Power, Profile_List
    
    if Item == "blade_size":
        
        Power = []
        Diameters = np.arange(1,12.5,0.5)
    
        print("Calculating blade diameter optimization")
        Setup_Profile([[1,12],[0,5],[2,0]])
        
        for Turbine in np.arange(1,21,1):

            Power.append([])
            
            for Diameter in np.arange(1,12.5,0.5):
                
                Total_Mechanical_Energy = 0
                Run_Simulation(step=10, tidal_function="sine", turbines=Turbine, diameter=Diameter, slucies=0, sluice_size=80, profile=1, time=60000, econ=False, output=False, graphs=False, graph_head=False, graph_QV=False, graph_P=False)
                
                for Energy in Global_Power:
                    #Global_Power_Elec.append(i*Eff_Turbine*Eff_Gearbox*Eff_Generator)
                    Total_Mechanical_Energy += Energy
                    
                Power[Turbine-1].append(Total_Mechanical_Energy)

        plt.figure(figsize=plt.figaspect(1)*2)
        ax = plt.axes()
        plt.title("Energy Output Vs Blade Diameter For Different Number of Turbines")
        ax.set_xlabel("Blade Diameter (m)")
        ax.set_ylabel("Energy Output (J)")
        
        Count = 0
        
        for i in Power:
            Temp = ax.plot(Diameters, Power[Count], label=("Turbines: " + str(Count+1)), linewidth=2)
            Count += 1
            
        #ax.legend("1 turbine","2 turbines","3 turbines","4 turbines","5 turbines","6 turbines","7 turbines","8 turbines","9 turbines","10 turbines","11 turbines","12 turbines","13 turbinse","14 turbines","15 turbines","16 turbines","17 turbines","18 turbines","19 turbines","20 turbines",)
        ax.legend()
        plt.minorticks_on()
        ax.grid(which='major', color='black', linestyle='-', linewidth=1)
        ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
        
    elif Item == "turbine_number":
        
        Power = []
        Turbines = np.arange(1,21,1)
        Diameters = np.arange(1,12.5,0.5)
    
        print("Calculating turbine number optimization")
        Setup_Profile([[1,12],[0,5],[2,0]])
        Count = 0
        
        for Diameter in np.arange(1,12.5,0.5):

            Power.append([])
            Count += 1
            
            for Turbine in np.arange(1,21,1):
                
                Total_Mechanical_Energy = 0
                Run_Simulation(step=10, tidal_function="sine", turbines=Turbine, diameter=Diameter, slucies=0, sluice_size=80, profile=1, time=60000, econ=False, output=False, graphs=False, graph_head=False, graph_QV=False, graph_P=False)
                
                for Energy in Global_Power:
                    #Global_Power_Elec.append(i*Eff_Turbine*Eff_Gearbox*Eff_Generator)
                    Total_Mechanical_Energy += Energy
                    
                Power[Count-1].append(Total_Mechanical_Energy)

        plt.figure(figsize=plt.figaspect(1)*2)
        ax = plt.axes()
        plt.title("Energy Output Vs Number of Turbines For Different Blasde Diameters")
        ax.set_xlabel("Number of Turbines")
        ax.set_ylabel("Energy Output (J)")
        
        Count = 0
        
        for i in Power:
            Temp = ax.plot(Turbines, Power[Count], label=("Diameter: " + str(Diameters[Count])), linewidth=2)
            Count += 1
            
        #ax.legend("1 turbine","2 turbines","3 turbines","4 turbines","5 turbines","6 turbines","7 turbines","8 turbines","9 turbines","10 turbines","11 turbines","12 turbines","13 turbinse","14 turbines","15 turbines","16 turbines","17 turbines","18 turbines","19 turbines","20 turbines",)
        ax.legend()
        plt.minorticks_on()
        ax.grid(which='major', color='black', linestyle='-', linewidth=1)
        ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
        
    elif Item == "algorithm":
        
        Profile_List = []
        Power = []
        Triggers = np.arange(12,-1,-0.5)
        Count = 1
        Turbine_Count = 0
        
        for Turbines in np.arange(1,21,1):
        
            Power.append([])
            Profile_List = []
            Turbine_Count += 1
            
            for i in np.arange(12,-1,-0.5):
                
                Total_Mechanical_Energy = 0
                Setup_Profile([[1,12],[0,i],[2,0]])
                Run_Simulation(step=10, tidal_function="sine", turbines=Turbine_Count, diameter=7, slucies=0, sluice_size=80, profile=Count, time=60000, econ=False, output=False, graphs=False, graph_head=False, graph_QV=False, graph_P=False)
            
                for Energy in Global_Power:
                        #Global_Power_Elec.append(i*Eff_Turbine*Eff_Gearbox*Eff_Generator)
                        Total_Mechanical_Energy += Energy
                        
                Power[Turbine_Count-1].append(Total_Mechanical_Energy)
                Count += 1

        plt.figure(figsize=plt.figaspect(1)*2)
        ax = plt.axes()
        plt.title("Energy Output Vs Sluice Gate Triggering Height")
        ax.set_xlabel("Lagoon Triggering Height (m)")
        ax.set_ylabel("Energy Output (J)")
        
        Count = 0
        
        for i in Power:
            Temp = ax.plot(Triggers, Power[Count], label=("Turbines: " + str(Count+1)), linewidth=2)
            Count += 1
        
        ax.legend()
        plt.minorticks_on()
        ax.grid(which='major', color='black', linestyle='-', linewidth=1)
        ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)

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
    FFT_Tidal_Function = [0]
    
    for i in range(Interval):
        
        Time.append(i)
        #Tide_Height.append(6*np.cos(2*np.pi*0.0000231*i-np.pi)+6)
        #FFT_Tidal_Function.append((3.68*np.sin(2*np.pi*(2.236e-5)*i*10))+(1.499*np.cos(2*np.pi*(2.315e-5)*i*10))+(0.73*np.cos(2*np.pi*(2.195e-5)*i*10))+(0.3581*np.sin(2*np.pi*(2.322e-5)*i*10)))
        FFT_Tidal_Function.append((3.678*np.cos(2*np.pi*(2.236e-5)*i+(48.06*(np.pi/180))))+(1.4298*np.cos(2*np.pi*(2.239e-5)*i-(130.6*(np.pi/180))))+(1.4986*np.cos(2*np.pi*(2.315e-5)*i+(176.7*(np.pi/180))))+(0.7298*np.cos(2*np.pi*(2.194e-5)*i-(53.75*(np.pi/180)))) + 6.278)
    
    plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes() #proj_type = 'ortho'
    plt.title("Ideal Semidiurnal Tide Over 24 Hours")
    ax.set_xlabel("Time (H)")
    ax.set_ylabel("Head (m)")
    ax.plot(Time, FFT_Tidal_Function)
    plt.minorticks_on()
    #plt.xticks(xlim, [str(n).zfill(2) + ':00' for n in np.arange(0, 26,2)])
    ax.grid(which='major', color='black', linestyle='-', linewidth=1)
    ax.grid(which='minor', color='black', linestyle='--', linewidth=0.5)
    
def Analyitic_Simulation_Simple(Duration=21600):        #Calculates the volume and graphs it over a time parameter analyitically.
   
    M = 1453217
    G = 9.807
    V_0 = 18891820
    Area = (np.power(7.5,2)/4)*np.pi
    Time = [0]
    Volume = [0]
    Volume[0] = V_0
    Euler_Volume = [V_0]
 
    for i in range(Duration):
        Time.append(i)
        Volume.append(np.power(np.sqrt(V_0)-0.5*Area*i*np.sqrt((2*G)/(M)),2))
        Euler_Volume.append(Euler_Volume[-1]-(Area*np.sqrt(2*G*(Euler_Volume[-1]/M))))
    
    plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes() #proj_type = 'ortho'
    plt.title("Volume Vs Time")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Volume (m^3)")
    ax.plot(Time, Volume)
    ax.plot(Time, Euler_Volume)
    ax.legend("Analytical", "Numerical")
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
        print("\n1. To run a simulation, you must first setup an operating profile to tell the lagoon how it should operate. To learn how to do this, type help('profile'). Alternativly, you can use the defaul profile already built into the program.")
        print("\n2. To learn how to run a simulation, type help('sim').")
    
    elif args[0] == "profile":
        
        print("Help profile:")
        print("\nThe lagoon can be in one of 5 states at any given time: \
State 0 - This is the waiting state, the lagoon is waiting for a head difference to build up before opening the sluice gates.\n\
State 1 - The lagoon is filling WITHOUT generating electricity. In this state, turbines are not moving, water is entering the lagoon through the turbine shafts and through sluice gates in the wall.\n\
State 2 - The lagoon is draining whilst generating electricity.\n\
State 3 - The lagoon is filling whilst generating electricity.\n\
State 4 - The lagoon is draining WITHOUT generating electricity.\n")
        print("The operational algorithm governs which states the lagoon switches between and when. The algorithm is built up of stages, each stage coresponds to one state and the condition that must be met to move to the next stage (state).")
    
        
    elif args[0] == "sim":
        
        print("Help simulation:")
        print("\nHere is some help info on running a simulation")
        
    else:
        
        print("Help command unknown, are you sure you typed it correctly?")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
"""
IDP 3 Group 9 - Tidal Lagoon Mathematical Model
Module: Integration of Grid Element
Gianluca Cantone 2019
"""

#Libraries ====================================================================
import math as m
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

#Variables ==================================================================== 
  


#Main Program =================================================================
"""
Components:
    
Surface - Triangular plot of 4 coords
Ocean - Surface plot (flat grid) computed from height input and Surface
Element 1 Side 1 (E1S1) - First side of first element, surface plot comnputed from Ocean and Surface
Element 1 Side 2 (E1S2) - Second side of first element, surface plot comnputed from Ocean and Surface
Element 2 Side 1 (E2S1) - First side of second element, surface plot comnputed from Ocean and Surface 
Element 2 Side 2 (E2S2) - Second side of second element, surface plot comnputed from Ocean and Surface    
Side 3 - Element 1 & 2 boundary, surface plot comnputed from Ocean and Surface

Boundary exists between coordinate 2 and 4 
"""

def view_element(Surface, Height):
    
    Surface_X = np.array([Surface[0][0], Surface[1][0], Surface[2][0], Surface[3][0]])
    Surface_Y = np.array([Surface[0][1], Surface[1][1], Surface[2][1], Surface[3][1]])
    Surface_Z = np.array([Surface[0][2], Surface[1][2], Surface[2][2], Surface[3][2]])
    
    Ocean_X = np.array([[Surface[0][0], Surface[1][0]], [Surface[3][0], Surface[2][0]]])
    Ocean_Y = np.array([[Surface[0][1], Surface[1][1]], [Surface[3][1], Surface[2][1]]])
    Ocean_Z = np.array([[Height, Height], [Height, Height]])
    
    E1S1_X = np.array([[Surface[0][0], Surface[1][0]], [Surface[0][0], Surface[1][0]]])
    E1S1_Y = np.array([[Surface[0][1], Surface[1][1]], [Surface[0][1], Surface[1][1]]])
    E1S1_Z = np.array([[Height, Height], [Surface[0][2], Surface[1][2]]])
    
    E1S2_X = np.array([[Surface[0][0], Surface[3][0]], [Surface[0][0], Surface[3][0]]])
    E1S2_Y = np.array([[Surface[0][1], Surface[3][1]], [Surface[0][1], Surface[3][1]]])
    E1S2_Z = np.array([[Height, Height], [Surface[0][2], Surface[3][2]]])
    
    E2S1_X = np.array([[Surface[3][0], Surface[2][0]], [Surface[3][0], Surface[2][0]]])
    E2S1_Y = np.array([[Surface[3][1], Surface[2][1]], [Surface[3][1], Surface[2][1]]])
    E2S1_Z = np.array([[Height, Height], [Surface[3][2], Surface[2][2]]])
    
    E2S2_X = np.array([[Surface[1][0], Surface[2][0]], [Surface[1][0], Surface[2][0]]])
    E2S2_Y = np.array([[Surface[1][1], Surface[2][1]], [Surface[1][1], Surface[2][1]]])
    E2S2_Z = np.array([[Height, Height], [Surface[1][2], Surface[2][2]]])
    
    Side_3_X = np.array([[Surface[1][0], Surface[3][0]], [Surface[1][0], Surface[3][0]]])
    Side_3_Y = np.array([[Surface[1][1], Surface[3][1]], [Surface[1][1], Surface[3][1]]])
    Side_3_Z = np.array([[Height, Height], [Surface[1][2], Surface[3][2]]])

    
    fig = plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes(projection='3d')
    
    plt.title("Integration Element (Two Triangular Pieces)")
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")
    
    ax.plot_trisurf(Surface_X, Surface_Y, Surface_Z, color="grey")
    ax.plot_surface(Ocean_X, Ocean_Y, Ocean_Z, color="blue")
    ax.plot_surface(E1S1_X, E1S1_Y, E1S1_Z, color="green", alpha=0.5, edgecolor="black")
    ax.plot_surface(E1S2_X, E1S2_Y, E1S2_Z, color="green", alpha=0.5, edgecolor="black")
    ax.plot_surface(E2S1_X, E2S1_Y, E2S1_Z, color="green", alpha=0.5, edgecolor="black")
    ax.plot_surface(E2S2_X, E2S2_Y, E2S2_Z, color="green", alpha=0.5, edgecolor="black")
    ax.plot_surface(Side_3_X, Side_3_Y, Side_3_Z, color="red", alpha=0.5, edgecolor="black")
    
def Integrate():
    
    '''
    Process:
        
    Step 1: Normalize coordinates so O1 is origin (O)
    Step 2: 
    Step 3:
    Step 3:Ocean 1 and 3 become your x axis, Ocean 1 is origin
    Step 3:Integrate x axis from the line from Ocean 1 to Ocean 2 and the line Ocean 2 to Ocean 3

def test():
    
    S1_X = np.array([0,1,1])
    S1_Y = np.array([0,0,1])
    S1_Z = np.array([1,0,0])
    
    S2_X = np.array([0,0,1])
    S2_Y = np.array([0,1,1])
    S2_Z = np.array([1,1,0])
    
    S3_X = np.array([0,0,1,1])
    S3_Y = np.array([0,1,1,0])
    S3_Z = np.array([2,2,2,2])
    
    S4_X = np.array([[0,1],
                     [0,1]])
    S4_Y = np.array([[0,0],
                     [0,0]])
    S4_Z = np.array([[2,2],
                     [1,0]])
    
    fig = plt.figure(figsize=plt.figaspect(1)*2)
    ax = plt.axes(projection='3d')
    
    plt.title("Integration Element (Two Triangular Pieces)")
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")
    
    ax.plot_trisurf(S1_X, S1_Y, S1_Z)
    ax.plot_trisurf(S2_X, S2_Y, S2_Z)
    ax.plot_trisurf(S3_X, S3_Y, S3_Z)
    ax.plot_surface(S4_X, S4_Y, S4_Z, alpha=0.5)
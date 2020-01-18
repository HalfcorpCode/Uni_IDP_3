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

def view_element(Surface, Height, Tet):
    
    #Calculate surface coordinates
    
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
    ax.plot_surface(Side_3_X, Side_3_Y, Side_3_Z, color="darkgreen", alpha=0.5, edgecolor="black")
    
    '''
    Check how many corners of the sea floor intersect the ocean plane:
        0 = Calculate integral of triangular prizm and two tetrahedrons
        1 = Calculate volume of two tetrahedrons if ocean was at max height, take away volume of prizm between full height and actual height, add volume of extra tetrahedron sticking out
        2 = 
        3 = 
        4 = Return 0 
            
    *Split square into two triangular sections. Order coordinates in terms of height.
    *Find two tetrahedrons by drawing straight lines out from highest coordinate
    *Switch case for adding/taking away together volume elements (Find volumes of tetrahedrons and prism, add together)
    *Return volume
    '''
    
    Sorted_Coords_A = sorted([Surface[0],Surface[1],Surface[3]], key=lambda x: x[2])
    Sorted_Coords_B = sorted([Surface[1],Surface[2],Surface[3]], key=lambda x: x[2])
    
    Coords_Above_A = 0
    Coords_Above_B = 0
    Element_A_Volume = 0
    Element_B_Volume = 0
    
    if Surface[0][2] > Height:
        Coords_Above_A += 1
    if Surface[1][2] > Height:
        Coords_Above_A += 1
    if Surface[3][2] > Height:
        Coords_Above_A += 1    
    
    print("Coords above A: " + str(Coords_Above_A))
    
    '''
    Get highest coordinate, generate tetrahedron 1 and 2 coordinates from it
    '''
    Tetrahedron_A_1_Coords = [Sorted_Coords_A[2],[Sorted_Coords_A[1][0],Sorted_Coords_A[1][1],Sorted_Coords_A[2][2]],[Sorted_Coords_A[0][0],Sorted_Coords_A[0][1],Sorted_Coords_A[2][2]],Sorted_Coords_A[1]]  #highest sea bed, other two sea beds at highest height, middle seabed coord
    Tetrahedron_A_2_Coords = [Sorted_Coords_A[2],[Sorted_Coords_A[0][0],Sorted_Coords_A[0][1],Sorted_Coords_A[2][2]],Sorted_Coords_A[1],Sorted_Coords_A[0]]  #highest sea bed, lowest sea bed at highest height, middle sea bed, lowest seabed coord
    print(Tetrahedron_A_1_Coords)
    if Coords_Above_A == 0:
        '''
        Element_A = Height*area of ocean triangle
        Element_A += tetrahedron 1 + tetrahedron 2
        '''  
        Delta_Height = Height - Sorted_Coords_A[2][2]   #Difference between ocean height and height of highest coordinate
        Tri_Area = Triangle_Area([[Surface[0][0],Surface[0][1],Height],
                                  [Surface[1][0],Surface[1][1],Height],
                                  [Surface[3][0],Surface[3][1],Height]])
        Element_A_Volume = (Delta_Height*Tri_Area) + Tetra_Volume(Tetrahedron_A_1_Coords) + Tetra_Volume(Tetrahedron_A_2_Coords)
        print("Volume is delta height " + str(Delta_Height) + " times triangle area " + str(Tri_Area) + " plus both tetrahedrons equals: " + str(Element_A_Volume))
        print("Tetrahedron volumes: " + str(Tetra_Volume(Tetrahedron_A_1_Coords)) + " and " + str(Tetra_Volume(Tetrahedron_A_2_Coords)))
        
    elif Coords_Above_A == 3:    
        
        Element_A_Volume = 0
    
    elif Coords_Above_A == 1:
        print("1 coord above")
    
    elif Coords_Above_A == 2:
        print("2 coord above")

    if Tet == "-t":
        
        #Calculate tetrahedrons coordinates
        
        Tet1_X = np.array([Tetrahedron_A_1_Coords[0][0],Tetrahedron_A_1_Coords[1][0],Tetrahedron_A_1_Coords[2][0]])
        Tet1_Y = np.array([Tetrahedron_A_1_Coords[0][1],Tetrahedron_A_1_Coords[1][1],Tetrahedron_A_1_Coords[2][1]])
        Tet1_Z = np.array([Tetrahedron_A_1_Coords[0][2],Tetrahedron_A_1_Coords[1][2],Tetrahedron_A_1_Coords[2][2]])
        ax.plot_trisurf(Tet1_X, Tet1_Y, Tet1_Z, color="red", alpha=0.5)
        
        Tet2_X = np.array([Tetrahedron_A_2_Coords[0][0],Tetrahedron_A_2_Coords[1][0],Tetrahedron_A_2_Coords[2][0]])
        Tet2_Y = np.array([Tetrahedron_A_2_Coords[0][1],Tetrahedron_A_2_Coords[1][1],Tetrahedron_A_2_Coords[2][1]])
        Tet2_Z = np.array([Tetrahedron_A_2_Coords[0][2],Tetrahedron_A_2_Coords[1][2],Tetrahedron_A_2_Coords[2][2]])       
        ax.plot_trisurf(Tet2_X, Tet2_Y, Tet2_Z, color="red", alpha=0.5)
        
    return Element_A_Volume + Element_B_Volume

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
    
def Tetra_Volume(Coords):
    
    '''Volume = 1/6 ((A to D) DOT ((A to B) X (A to C)))'''
    return np.abs(1/6*(np.dot(np.cross(Vector_Diff([Coords[0], Coords[1]]),Vector_Diff([Coords[0], Coords[2]])), Vector_Diff([Coords[0], Coords[3]]))))

def Triangle_Area(Coords):
    
    A = Vector_Diff([Coords[0],Coords[1]])
    B = Vector_Diff([Coords[1],Coords[2]])
    C = Vector_Diff([Coords[2],Coords[0]])
    S = (np.linalg.norm(A)+np.linalg.norm(B)+np.linalg.norm(C))*0.5
    print(S)
    A = np.sqrt(S*(S-np.linalg.norm(A))*(S-np.linalg.norm(B))*(S-np.linalg.norm(C)))    
    return round(A,4)
    
def Vector_Diff(Vectors):
    
    Result = [0,0,0]
    Result[0]= Vectors[1][0]-Vectors[0][0]
    Result[1]= Vectors[1][1]-Vectors[0][1]
    Result[2]= Vectors[1][2]-Vectors[0][2]
    
    return Result
        
        
        
        
        
        
        
        
        
        
        
        
        
        
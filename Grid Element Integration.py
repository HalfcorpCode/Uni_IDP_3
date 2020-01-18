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
    
    Sorted_Coords_A = sorted([Surface[0], Surface[1], Surface[3]], key=lambda x: x[2])
    Sorted_Coords_B = sorted([Surface[1], Surface[2], Surface[3]], key=lambda x: x[2])
    
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
        
    if Surface[1][2] > Height:
        Coords_Above_B += 1
    if Surface[2][2] > Height:
        Coords_Above_B += 1
    if Surface[3][2] > Height:
        Coords_Above_B += 1  
        
    print("Coords above A: " + str(Coords_Above_A) + " Coords above B: " + str(Coords_Above_B))
    
    if Coords_Above_A+Coords_Above_B == 6:
        
        return 0
    
    else:
    
        Tetrahedron_A_1_Coords = [Sorted_Coords_A[2], [Sorted_Coords_A[1][0], Sorted_Coords_A[1][1], Sorted_Coords_A[2][2]], [Sorted_Coords_A[0][0], Sorted_Coords_A[0][1], Sorted_Coords_A[2][2]], Sorted_Coords_A[1]]  #highest sea bed, other two sea beds at highest height, middle seabed coord
        Tetrahedron_A_2_Coords = [Sorted_Coords_A[2], [Sorted_Coords_A[0][0], Sorted_Coords_A[0][1], Sorted_Coords_A[2][2]], Sorted_Coords_A[1], Sorted_Coords_A[0]]  #highest sea bed, lowest sea bed at highest height, middle sea bed, lowest seabed coord
        
        Tetrahedron_B_1_Coords = [Sorted_Coords_B[2], [Sorted_Coords_B[1][0], Sorted_Coords_B[1][1], Sorted_Coords_B[2][2]], [Sorted_Coords_B[0][0], Sorted_Coords_B[0][1], Sorted_Coords_B[2][2]], Sorted_Coords_B[1]]  #highest sea bed, other two sea beds at highest height, middle seabed coord
        Tetrahedron_B_2_Coords = [Sorted_Coords_B[2], [Sorted_Coords_B[0][0], Sorted_Coords_B[0][1], Sorted_Coords_B[2][2]], Sorted_Coords_B[1], Sorted_Coords_B[0]]  #highest sea bed, lowest sea bed at highest height, middle sea bed, lowest seabed coord
        
        Triangle_Area_Ocean_A = Triangle_Area([[Surface[0][0],Surface[0][1],Height],
                                               [Surface[1][0],Surface[1][1],Height],
                                               [Surface[3][0],Surface[3][1],Height]])
    
        Triangle_Area_Ocean_B = Triangle_Area([[Surface[1][0],Surface[1][1],Height],
                                               [Surface[2][0],Surface[2][1],Height],
                                               [Surface[3][0],Surface[3][1],Height]])
    
        print("Triangle ocean area A: " + str(Triangle_Area_Ocean_A) + " Triangle ocean area B: " + str(Triangle_Area_Ocean_B))
    
        #Element A
    
        if Coords_Above_A == 3:
        
            print("3 coords above element A, volume is 0")
            Element_A_Volume = 0
            
        elif Coords_Above_A == 0:    
            
            print("0 coords above element A")
            Delta_Height = Height - Sorted_Coords_A[2][2]   #Difference between ocean height and height of highest coordinate
            Element_A_Volume = (Delta_Height*Triangle_Area_Ocean_A) + Tetra_Volume(Tetrahedron_A_1_Coords) + Tetra_Volume(Tetrahedron_A_2_Coords)
            print("Triangle area: " + str(Triangle_Area_Ocean_A) + " Delta height: " + str(Delta_Height) + " Tetra 1:  " + str(Tetra_Volume(Tetrahedron_A_1_Coords)) + " Tetra 2: " + str(Tetra_Volume(Tetrahedron_A_2_Coords)))
        
        elif Coords_Above_A == 1:
            
            print("1 coord above element A")
            Subtract_Prism = Triangle_Area_Ocean_A*(Sorted_Coords_A[2][2]-Height) 
            High_Middle_Intersection = Line_Intersect([Sorted_Coords_A[2], Sorted_Coords_A[1], [Sorted_Coords_A[2][0],Sorted_Coords_A[2][1], Height], [Sorted_Coords_A[1][0],Sorted_Coords_A[1][1], Height]])    #Highest point to middle point, ocean below highest point to ocean above middle point
            High_Low_Intersection = Line_Intersect([Sorted_Coords_A[2], Sorted_Coords_A[0], [Sorted_Coords_A[2][0],Sorted_Coords_A[2][1], Height], [Sorted_Coords_A[0][0],Sorted_Coords_A[0][1], Height]])    #Highest point to low point, ocean below highest point to ocean above low point
            Extra_Tetra_Coords = [Sorted_Coords_A[2], [Sorted_Coords_A[2][0],Sorted_Coords_A[2][1], Height], High_Middle_Intersection, High_Low_Intersection]  #Highest sea floor, highest sea floor at ocean height, line between highest and middle sea floor, line between highest and lowest sea floor.
            Element_A_Volume = Tetra_Volume(Tetrahedron_A_1_Coords) + Tetra_Volume(Tetrahedron_A_2_Coords) - Subtract_Prism + Tetra_Volume(Extra_Tetra_Coords)
            print("Both tetra volumes: " + str(Tetra_Volume(Tetrahedron_A_1_Coords) + Tetra_Volume(Tetrahedron_A_2_Coords)) + " Prism (to subtract) " + str(Subtract_Prism) + " Extra tetra volume (to add): " + str(Tetra_Volume(Extra_Tetra_Coords)))
        
        elif Coords_Above_A == 2:
            
            print("2 coord above elements A")
            Low_Middle_Intersection = Line_Intersect([Sorted_Coords_A[0], Sorted_Coords_A[1], [Sorted_Coords_A[0][0], Sorted_Coords_A[0][1], Height], [Sorted_Coords_A[1][0], Sorted_Coords_A[1][1], Height]]) #Low sea floor, middle sea floor, low sea floor at height, middle sea floor at height
            Low_High_Intersection = Line_Intersect([Sorted_Coords_A[0], Sorted_Coords_A[2], [Sorted_Coords_A[0][0], Sorted_Coords_A[0][1], Height], [Sorted_Coords_A[2][0], Sorted_Coords_A[2][1], Height]]) #Low sea floor, high sea floor, low sea floor at height, high sea floor at height
            Element_A_Volume = Tetra_Volume([Sorted_Coords_A[0], [Sorted_Coords_A[0][0], Sorted_Coords_A[0][1], Height], Low_Middle_Intersection, Low_High_Intersection])
            print("Only tetra volume: " + str(Element_A_Volume))
        
        print("Element A volume: " + str(Element_A_Volume))
        
        #Element B
        
        if Coords_Above_B == 3:
        
            print("3 coords above element B, volume is 0")
            Element_B_Volume = 0
            
        elif Coords_Above_B == 0:    
            
            print("0 coords above element B")
            Delta_Height = Height - Sorted_Coords_B[2][2]   #Difference between ocean height and height of highest coordinate
            Element_B_Volume = (Delta_Height*Triangle_Area_Ocean_B) + Tetra_Volume(Tetrahedron_B_1_Coords) + Tetra_Volume(Tetrahedron_B_2_Coords)
            print("Triangle area: " + str(Triangle_Area_Ocean_B) + " Delta height: " + str(Delta_Height) + " Tetra 1:  " + str(Tetra_Volume(Tetrahedron_B_1_Coords)) + " Tetra 2: " + str(Tetra_Volume(Tetrahedron_B_2_Coords)))
        
        elif Coords_Above_B == 1:
            
            print("1 coord above element B")
            Subtract_Prism = Triangle_Area_Ocean_B*(Sorted_Coords_B[2][2]-Height) 
            High_Middle_Intersection = Line_Intersect([Sorted_Coords_B[2], Sorted_Coords_B[1], [Sorted_Coords_B[2][0],Sorted_Coords_B[2][1], Height], [Sorted_Coords_B[1][0],Sorted_Coords_B[1][1], Height]])    #Highest point to middle point, ocean below highest point to ocean above middle point
            High_Low_Intersection = Line_Intersect([Sorted_Coords_B[2], Sorted_Coords_B[0], [Sorted_Coords_B[2][0],Sorted_Coords_B[2][1], Height], [Sorted_Coords_B[0][0],Sorted_Coords_B[0][1], Height]])    #Highest point to low point, ocean below highest point to ocean above low point
            Extra_Tetra_Coords = [Sorted_Coords_B[2], [Sorted_Coords_B[2][0],Sorted_Coords_B[2][1], Height], High_Middle_Intersection, High_Low_Intersection]  #Highest sea floor, highest sea floor at ocean height, line between highest and middle sea floor, line between highest and lowest sea floor.
            Element_B_Volume = Tetra_Volume(Tetrahedron_B_1_Coords) + Tetra_Volume(Tetrahedron_B_2_Coords) - Subtract_Prism + Tetra_Volume(Extra_Tetra_Coords)
            print("Both tetra volumes: " + str(Tetra_Volume(Tetrahedron_B_1_Coords) + Tetra_Volume(Tetrahedron_B_2_Coords)) + " Prism (to subtract) " + str(Subtract_Prism) + " Extra tetra volume (to add): " + str(Tetra_Volume(Extra_Tetra_Coords)))
        
        elif Coords_Above_B == 2:
            
            print("2 coord above elements B")
            Low_Middle_Intersection = Line_Intersect([Sorted_Coords_B[0], Sorted_Coords_B[1], [Sorted_Coords_B[0][0], Sorted_Coords_B[0][1], Height], [Sorted_Coords_B[1][0], Sorted_Coords_B[1][1], Height]]) #Low sea floor, middle sea floor, low sea floor at height, middle sea floor at height
            Low_High_Intersection = Line_Intersect([Sorted_Coords_B[0], Sorted_Coords_B[2], [Sorted_Coords_B[0][0], Sorted_Coords_B[0][1], Height], [Sorted_Coords_B[2][0], Sorted_Coords_B[2][1], Height]]) #Low sea floor, high sea floor, low sea floor at height, high sea floor at height
            Element_B_Volume = Tetra_Volume([Sorted_Coords_B[0], [Sorted_Coords_B[0][0], Sorted_Coords_B[0][1], Height], Low_Middle_Intersection, Low_High_Intersection])
            print("Only tetra volume: " + str(Element_B_Volume))
    
        print("Element B volume: " + str(Element_B_Volume))

        if Tet == "-t":
            
            #Calculate tetrahedrons coordinates
            
            TetA1_X = np.array([Tetrahedron_A_1_Coords[0][0],Tetrahedron_A_1_Coords[1][0],Tetrahedron_A_1_Coords[2][0]])
            TetA1_Y = np.array([Tetrahedron_A_1_Coords[0][1],Tetrahedron_A_1_Coords[1][1],Tetrahedron_A_1_Coords[2][1]])
            TetA1_Z = np.array([Tetrahedron_A_1_Coords[0][2],Tetrahedron_A_1_Coords[1][2],Tetrahedron_A_1_Coords[2][2]])
            ax.plot_trisurf(TetA1_X, TetA1_Y, TetA1_Z, color="red", alpha=0.5)
            
            TetA2_X = np.array([Tetrahedron_A_2_Coords[0][0],Tetrahedron_A_2_Coords[1][0],Tetrahedron_A_2_Coords[2][0]])
            TetA2_Y = np.array([Tetrahedron_A_2_Coords[0][1],Tetrahedron_A_2_Coords[1][1],Tetrahedron_A_2_Coords[2][1]])
            TetA2_Z = np.array([Tetrahedron_A_2_Coords[0][2],Tetrahedron_A_2_Coords[1][2],Tetrahedron_A_2_Coords[2][2]])       
            ax.plot_trisurf(TetA2_X, TetA2_Y, TetA2_Z, color="red", alpha=0.5)
            
            TetB1_X = np.array([Tetrahedron_B_1_Coords[0][0],Tetrahedron_B_1_Coords[1][0],Tetrahedron_B_1_Coords[2][0]])
            TetB1_Y = np.array([Tetrahedron_B_1_Coords[0][1],Tetrahedron_B_1_Coords[1][1],Tetrahedron_B_1_Coords[2][1]])
            TetB1_Z = np.array([Tetrahedron_B_1_Coords[0][2],Tetrahedron_B_1_Coords[1][2],Tetrahedron_B_1_Coords[2][2]])
            ax.plot_trisurf(TetB1_X, TetB1_Y, TetB1_Z, color="red", alpha=0.5)
            
            TetB2_X = np.array([Tetrahedron_B_2_Coords[0][0],Tetrahedron_B_2_Coords[1][0],Tetrahedron_B_2_Coords[2][0]])
            TetB2_Y = np.array([Tetrahedron_B_2_Coords[0][1],Tetrahedron_B_2_Coords[1][1],Tetrahedron_B_2_Coords[2][1]])
            TetB2_Z = np.array([Tetrahedron_B_2_Coords[0][2],Tetrahedron_B_2_Coords[1][2],Tetrahedron_B_2_Coords[2][2]])       
            ax.plot_trisurf(TetB2_X, TetB2_Y, TetB2_Z, color="red", alpha=0.5)

    return Element_A_Volume + Element_B_Volume
    
def Tetra_Volume(Coords):
    
    return round(np.abs(1/6*(np.dot(np.cross(Vector_Diff([Coords[0], Coords[1]]),Vector_Diff([Coords[0], Coords[2]])), Vector_Diff([Coords[0], Coords[3]])))),5)    #Volume = 1/6 ((A to D) DOT ((A to B) X (A to C)))

def Triangle_Area(Coords):
    
    A = Vector_Diff([Coords[0],Coords[1]])
    B = Vector_Diff([Coords[1],Coords[2]])
    C = Vector_Diff([Coords[2],Coords[0]])
    S = (np.linalg.norm(A)+np.linalg.norm(B)+np.linalg.norm(C))*0.5
    A = np.sqrt(S*(S-np.linalg.norm(A))*(S-np.linalg.norm(B))*(S-np.linalg.norm(C)))    
    return round(A,4)
    
def Vector_Diff(Vectors):
    
    Result = [0,0,0]
    Result[0]= Vectors[1][0]-Vectors[0][0]
    Result[1]= Vectors[1][1]-Vectors[0][1]
    Result[2]= Vectors[1][2]-Vectors[0][2]
    
    return Result

def Line_Intersect(Coords):
    
    AP = Coords[0]
    AD = Vector_Diff([Coords[0], Coords[1]])
    BP = Coords[2]
    BD = Vector_Diff([Coords[2], Coords[3]])
    
    Divisor = ((AD[0])*(-BD[2])-(-BD[0])*(AD[2]))
    
    if Divisor == 0:
        Divisor = ((AD[1])*(-BD[2])-(-BD[1])*(AD[2]))
        Parameter = ((BP[1]-AP[1])*(-BD[2])-(-BD[1])*(BP[2]-AP[2]))/((AD[1])*(-BD[2])-(-BD[1])*(AD[2]))
    else:
        Parameter = ((BP[0]-AP[0])*(-BD[2])-(-BD[0])*(BP[2]-AP[2]))/Divisor
    
    X = AP[0]+(AD[0]*Parameter)
    Y = AP[1]+(AD[1]*Parameter)
    Z = AP[2]+(AD[2]*Parameter)
    
    return [X,Y,Z]
        
        
        

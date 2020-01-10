"""
IDP 3 Group 9 - Tidal Lagoon Mathematical Model
Module: 3D Model of Lagoon
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

fig = plt.figure()
ax = plt.axes(projection='3d')

plt.title("My Plot")
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")

'''
1. Program initializes coordinate arrays
2. Known coordinates are laoded into arrays
3. Intermediate values are interpolated and added into arrays
4. Integration algorithm integrates across grid at a given height, some new coordiantes will have to be interpolated from known coordinates.

Dummy test: 10x10 grid with edges initialized at 0. Dip in middle that will be interpolated outwards.
'''

X_Steps = np.linspace(0, 10, 11)
Y_Steps = np.linspace(0, 10, 11)
X_Grid, Y = np.meshgrid(X_Steps, y)

X_Grid = np.array([[-1,1,2,3],
                   [0,1,2,3],
                   [0,1,2,3],
                   [0,1,2,3]])

Y_Grid = np.array([[-1,0,0,0],
                   [1,1,1,1],
                   [2,2,2,2],
                   [3,3,3,3]])

Z_Grid = np.array([[1,0,0,0],
                   [0,0,0,0],
                   [0,0,0,0],
                   [0,0,0,0]])

ax.plot_wireframe(X_Grid, Y_Grid, Z_Grid, color="red")

#x = np.linspace(0, 10, 11)
#y = np.linspace(0, 10, 11)
#
#X, Y = np.meshgrid(x, y)
# 
#Z = np.array([[1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
#              [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
#              [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
#              [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
#              [1, 1, 1, 1, 0.5, 0, 0, 0, 0, 0, 0],
#              [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
#              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])
#
#Z2 = np.array([[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
#              [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]])
#
#ax.plot_wireframe(X, Y, Z2, color="green")



























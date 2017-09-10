# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import matplotlib.pyplot as plt
from visualize import drawPatches

'''
element patches
'''
#coarse World
NWorldCoarse = np.array([11,11])
NpCoarse = np.prod(NWorldCoarse+1)

#ratio between Fine and Coarse
A = np.zeros(NWorldCoarse)
A[5][5] = 2
ABase = A.flatten()   

fig = plt.figure('1')                                                               
ax = fig.add_subplot(1,3,1)     
ax.set_title('$U_{0}(T)$')                                                 
drawPatches(NWorldCoarse, ABase, fig, ax, 11)

A[4][4:7] = 1.5
A[6][4:7] = 1.5
A[5,4:8:2] = 1.5
ABase = A.flatten()
ax = fig.add_subplot(1,3,2)  
ax.set_title('$U_{1}(T)$')                                                    
drawPatches(NWorldCoarse, ABase, fig, ax, 11)

A[3][3:8] = 1
A[7][3:8] = 1
A[3:7,3] = 1
A[3:7,7] = 1
ABase = A.flatten()
ax = fig.add_subplot(1,3,3) 
ax.set_title('$U_{2}(T)$')                                                     
drawPatches(NWorldCoarse, ABase, fig, ax, 11)

'''
nodal patches
'''
#coarse World
NWorldCoarse = np.array([10,10])
NpCoarse = np.prod(NWorldCoarse+1)

#ratio between Fine and Coarse
A = np.zeros(NWorldCoarse)
A[4:6,4:6] = 2
ABase = A.flatten()   

fig = plt.figure('2')                                                               
ax = fig.add_subplot(1,3,1)     
ax.set_title('$\omega_{x,1}$')                                                 
drawPatches(NWorldCoarse, ABase,fig, ax, 10)

A[3][3:7] = 1.5
A[6][3:7] = 1.5
A[4:6,3] = 1.5
A[4:6,6] = 1.5
ABase = A.flatten()
ax = fig.add_subplot(1,3,2)  
ax.set_title('$\omega_{x,2}$')                                                    
drawPatches(NWorldCoarse, ABase,fig,ax, 10)

A[2][2:8] = 1
A[7][2:8] = 1
A[3:6,2] = 1
A[3:7,7] = 1
A[6][2] = 1
ABase = A.flatten()
ax = fig.add_subplot(1,3,3) 
ax.set_title('$\omega_{x,3}$')                                                     
drawPatches(NWorldCoarse, ABase,fig,ax, 10)

plt.show()
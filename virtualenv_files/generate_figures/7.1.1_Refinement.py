# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

#coarse World
NWorldCoarse = np.array([11,11])
NpCoarse = np.prod(NWorldCoarse+1)

A = np.zeros(NWorldCoarse)
ABase = A.flatten()   

aCube = ABase.reshape(NWorldCoarse)

# coarse mesh
fig = plt.figure('Coarse Mesh')                                                               
ax = fig.add_subplot(1,1,1)     
major_ticks = np.arange(0, 11, 1)
ax.imshow(aCube, extent=[0, 11, 0, 11], cmap='Greys')
ax.axis([0, 11, 0, 11])
ax.set_xticks(major_ticks)                                                       
ax.set_yticks(major_ticks)  
ax.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
fig.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.98,wspace=0.2,hspace=0.2)
ax.grid(which='major', linestyle="-", color="black", alpha=0.8)                                                

# refinement
fig = plt.figure('Refinement')                                                               
ax = fig.add_subplot(1,1,1)     
major_ticks = np.arange(0, 11, 1)
minor_ticks = np.arange(0, 11, 0.2)                                                
ax.imshow(aCube, extent=[0, 11, 0, 11], cmap='Greys')
ax.axis([0, 11, 0, 11])
ax.set_xticks(minor_ticks)
ax.set_xticks(minor_ticks)
ax.set_xticks(major_ticks)                                                       
ax.set_yticks(major_ticks)  
ax.minorticks_on()                                                     
ax.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
fig.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.98,wspace=0.2,hspace=0.2)
ax.grid(which='major', linestyle="-", color="black", alpha=0.8)                                                
ax.grid(which='minor', linestyle="-",linewidth=0.5 , color="black", alpha=0.4)

plt.show()
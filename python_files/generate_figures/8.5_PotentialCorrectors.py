# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import matplotlib.pyplot as plt

from gridlod import interp, coef
import pg_pert, buildcoef2d
from gridlod.world import World

from visualize import drawCoefficient, drawCoefficientGrid, drawCoefficientwg
    
bg = 0.05      #background
val = 1        #values

#fine World
NWorldFine = np.array([256, 256])
NpFine = np.prod(NWorldFine+1)                                                                               

#coarse World
NWorldCoarse = np.array([16,16])
NpCoarse = np.prod(NWorldCoarse+1)

#ratio between Fine and Coarse
NCoarseElement = NWorldFine//NWorldCoarse

boundaryConditions = np.array([[0, 0],
                               [0, 0]])

world = World(NWorldCoarse, NCoarseElement, boundaryConditions)

#righthandside
f = np.ones(NpCoarse)

#coefficient
CoefClass = buildcoef2d.Coefficient2d(NWorldFine, 
                                    bg                  = 0.05, 
                                    val                 = 1, 
                                    length              = 8, 
                                    thick               = 8, 
                                    space               = 8, 
                                    probfactor          = 1, 
                                    right               = 1, 
                                    down                = 0, 
                                    diagr1              = 0, 
                                    diagr2              = 0, 
                                    diagl1              = 0, 
                                    diagl2              = 0, 
                                    LenSwitch           = None, 
                                    thickSwitch         = None, 
                                    equidistant         = True, 
                                    ChannelHorizontal   = None, 
                                    ChannelVertical     = None,
                                    BoundarySpace       = True)

A = CoefClass.BuildCoefficient()
ABase = A.flatten()

plt.figure("OriginalCoefficient")
drawCoefficient(NWorldFine, ABase)

#potentially updated
fig = plt.figure("Correcting")
ax = fig.add_subplot(2,3,1)
ax.set_title('Original', fontsize=7)
drawCoefficientwg(NWorldFine, ABase,fig,ax,Greys=True)

numbers = [2,70,97,153,205]

D = CoefClass.SpecificVanish(           Number              = numbers,
                                        probfactor          = 1,
                                        PartlyVanish        = None,
                                        Original            = True)
                                
value2 = 50
R2 = CoefClass.SpecificValueChange(     ratio               = value2,
                                        Number              = numbers,
                                        probfactor          = 1,
                                        randomvalue         = None,
                                        negative            = None,
                                        ShapeRestriction    = True,
                                        ShapeWave           = None,
                                        Original            = True,
                                        NewShapeChange      = True)

# precompute
NWorldFine = world.NWorldFine
NWorldCoarse = world.NWorldCoarse
NCoarseElement = world.NCoarseElement

boundaryConditions = world.boundaryConditions
NpFine = np.prod(NWorldFine+1)
NpCoarse = np.prod(NWorldCoarse+1)

#interpolant
IPatchGenerator = lambda i, N: interp.L2ProjectionPatchMatrix(i, N, NWorldCoarse, NCoarseElement, boundaryConditions)

#old Coefficient (need flatten form)
Aold = coef.coefficientFine(NWorldCoarse, NCoarseElement, ABase)

for k in range(1,6):
    print(('<<<<<<<<<<<<<<<< ' + str(k) + ' >>>>>>>>>>>>>>>>'))
    pglod = pg_pert.PerturbedPetrovGalerkinLOD(Aold, world, k, IPatchGenerator, 1)
    pglod.originCorrectors(clearFineQuantities=False)
    
    #new Coefficient
    ANew = R2.flatten()
    Anew = coef.coefficientFine(NWorldCoarse, NCoarseElement, ANew)
    
    # tolerance = 0
    vis, eps = pglod.updateCorrectors(Anew, 0,clearFineQuantities=False, Computing = False)
    
    elemente = np.arange(np.prod(NWorldCoarse))
    
    if k==5:
        plt.figure("Error indicator")
        plt.plot(elemente,eps,label="e_{u,T}")
        plt.ylabel('$e_{u,T}$')
        plt.xlabel('Element')
        plt.subplots_adjust(left=0.09,bottom=0.09,right=0.99,top=0.99,wspace=0.2,hspace=0.2)
        plt.grid()
        
        eps.sort()
        es = []
        for i in range(0,np.size(eps)):
            es.append(eps[np.size(eps)-i-1])
    
        plt.figure("Error indicator log")
        plt.semilogy(elemente,es,label='e_{u,T}')
        elemente = np.arange(np.prod(NWorldCoarse))
        plt.ylabel('$e_{u,T}$',fontsize=20)
        plt.xlabel('Elements',fontsize=20)
        plt.subplots_adjust(left=0.12,bottom=0.14,right=0.99,top=0.99,wspace=0.2,hspace=0.2)
        plt.grid()
        plt.tick_params(axis='both', which='major', labelsize=17)
        plt.tick_params(axis='both', which='minor', labelsize=17)
    
    #potentially updated
    fig = plt.figure("Correcting")
    if k == 1:
        ax = fig.add_subplot(2,3,k+2)
        ax.set_title('Updated for $k=$'+str(k),fontsize=7)
    elif k ==2:
        ax = fig.add_subplot(2,3,k+2)
        ax.set_title('Updated for $k=$'+str(k),fontsize=7)
    elif k ==3:
        ax = fig.add_subplot(2,3,k+2)
        ax.set_title('Updated for $k=$'+str(k),fontsize=7)
    elif k ==4:
        ax = fig.add_subplot(2,3,k+2)
        ax.set_title('Updated for $k=$'+str(k),fontsize=7)
    if k!=5:
        drawCoefficientGrid(NWorldCoarse, vis,fig,ax)

fig = plt.figure("Correcting")                                                               
ax = fig.add_subplot(2,3,2)     
ax.set_title('Defects', fontsize=7)                                                 
drawCoefficientwg(NWorldFine, R2,fig,ax, Greys=True)

plt.show()
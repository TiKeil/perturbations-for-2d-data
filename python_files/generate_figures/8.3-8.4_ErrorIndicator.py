# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import matplotlib.pyplot as plt

from gridlod import interp, coef, util, fem, world, linalg
from gridlod.world import World
import femsolverCoarse
import pg_rand
import buildcoef2d

from visualize import drawCoefficient

def result(pglod, world, A, R, f, k, String):
    print(("------------------------------------- " + String + " -------------------------------------------")) 
    NWorldFine = world.NWorldFine
    NWorldCoarse = world.NWorldCoarse
    NCoarseElement = world.NCoarseElement
    
    boundaryConditions = world.boundaryConditions
    NpFine = np.prod(NWorldFine+1)
    NpCoarse = np.prod(NWorldCoarse+1)
        
    #new Coefficient
    ANew = R.flatten()
    Anew = coef.coefficientFine(NWorldCoarse, NCoarseElement, ANew)
    
    # tolerance = 0
    vis, eps = pglod.updateCorrectors(Anew, 0, f, 1, Computing = False)
    
    elemente = np.arange(np.prod(NWorldCoarse))
    
    plt.figure("Error indicators")
    plt.plot(elemente,eps,label=String)
    plt.ylabel('$e_{u,T}$')
    plt.xlabel('Element')
    plt.subplots_adjust(left=0.09,bottom=0.09,right=0.99,top=0.99,wspace=0.2,hspace=0.2)
    plt.legend(loc='upper right') #Legende
    plt.grid()

bg = 0.05       #background
val = 1         #values

#fine World
NWorldFine = np.array([256, 256])
NpFine = np.prod(NWorldFine+1)                                                                               

#coarse World
NWorldCoarse = np.array([16,16])
NpCoarse = np.prod(NWorldCoarse+1)

#ratio between Fine and Coarse
NCoarseElement = NWorldFine/NWorldCoarse

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

numbers = [2,70,97,153,205]

value1 = 3
R1 = CoefClass.SpecificValueChange( ratio               = value1,
                                    Number              = numbers,
                                    probfactor          = 1,
                                    randomvalue         = None,
                                    negative            = None,
                                    ShapeRestriction    = True,
                                    ShapeWave           = None,
                                    Original            = True)

plt.figure("Change in value to 3")
drawCoefficient(NWorldFine, R1)

value2 = 50
R2 = CoefClass.SpecificValueChange( ratio               = value2,
                                    Number              = numbers,
                                    probfactor          = 1,
                                    randomvalue         = None,
                                    negative            = None,
                                    ShapeRestriction    = True,
                                    ShapeWave           = None,
                                    Original            = True)

plt.figure("Change in value to 50")
drawCoefficient(NWorldFine, R2)
        
D = CoefClass.SpecificVanish(       Number              = numbers,
                                    probfactor          = 1,
                                    PartlyVanish        = None,
                                    Original            = True)
                                    
plt.figure("Disappearance")
drawCoefficient(NWorldFine, D)

E2 = CoefClass.SpecificMove(        Number              = numbers,
                                    steps               = 3,
                                    randomstep          = None,
                                    randomDirection     = None,
                                    Right               = 0,
                                    BottomRight         = 0,
                                    Bottom              = 0,
                                    BottomLeft          = 0,
                                    Left                = 0,
                                    TopLeft             = 1,
                                    Top                 = 0,
                                    TopRight            = 0,
                                    Original            = True)

plt.figure("Shift one step")
drawCoefficient(NWorldFine, E2)

E3 = CoefClass.SpecificMove(        Number              = numbers,
                                    steps               = 7,
                                    randomstep          = None,
                                    randomDirection     = None,
                                    Right               = 0,
                                    BottomRight         = 0,
                                    Bottom              = 0,
                                    BottomLeft          = 0,
                                    Left                = 0,
                                    TopLeft             = 1,
                                    Top                 = 0,
                                    TopRight            = 0,
                                    Original            = True)

plt.figure("Shift two steps")
drawCoefficient(NWorldFine, E3)

# precomputations 
NWorldFine = world.NWorldFine
NWorldCoarse = world.NWorldCoarse
NCoarseElement = world.NCoarseElement

boundaryConditions = world.boundaryConditions
NpFine = np.prod(NWorldFine+1)
NpCoarse = np.prod(NWorldCoarse+1)

#interpolant
IPatchGenerator = lambda i, N: interp.L2ProjectionPatchMatrix(i, N, NWorldCoarse, NCoarseElement, boundaryConditions)

#old Coefficient
Aold = coef.coefficientFine(NWorldCoarse, NCoarseElement, ABase)

k = 5

pglod = pg_rand.VcPetrovGalerkinLOD(Aold, world, k, IPatchGenerator, 0)
pglod.originCorrectors(clearFineQuantities=False)

# Change in value 1
result(pglod ,world, A, R1, f, k, 'Change in value to' + str(value1))

# Change in value 2 
result(pglod ,world, A, R2, f, k, 'Change in value to' + str(value2))

# Disappearance 
result(pglod, world, A, D, f, k, 'Disappearance')

# Shift one step
result(pglod, world, A, E2, f, k, 'Shift one step')

# Shift two steps
result(pglod, world, A, E3, f, k, 'Shift two steps')

plt.show()
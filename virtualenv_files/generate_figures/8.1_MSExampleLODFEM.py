```
# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
```

import numpy as np
import scipy.sparse as sparse

import matplotlib.pyplot as plt
from gridlod import util, world, fem, coef, interp
from gridlod.world import World
import femsolverCoarse, buildcoef2d, pg_rand

from visualize import drawCoefficient

def PGsolver(world, ABase, f,k):
    NWorldFine = world.NWorldFine
    NWorldCoarse = world.NWorldCoarse
    NCoarseElement = world.NCoarseElement
    
    boundaryConditions = world.boundaryConditions
    NpFine = np.prod(NWorldFine+1)
    NpCoarse = np.prod(NWorldCoarse+1)
    
    #interpolant
    IPatchGenerator = lambda i, N: interp.L2ProjectionPatchMatrix(i, N, NWorldCoarse, NCoarseElement, boundaryConditions)
    
    #Coefficient (need flatten form)
    aCoef = coef.coefficientFine(NWorldCoarse, NCoarseElement, ABase)

    pglod = pg_rand.VcPetrovGalerkinLOD(aCoef, world, k, IPatchGenerator, 0)
    pglod.originCorrectors(clearFineQuantities=False)

    KFull = pglod.assembleMsStiffnessMatrix()                                    
    MFull = fem.assemblePatchMatrix(NWorldCoarse, world.MLocCoarse)
    free  = util.interiorpIndexMap(NWorldCoarse)                                 

    bFull = MFull*f

    KFree = KFull[free][:,free]
    bFree = bFull[free]

    xFree = sparse.linalg.spsolve(KFree, bFree)

    basis = fem.assembleProlongationMatrix(NWorldCoarse, NCoarseElement)
    basisCorrectors = pglod.assembleBasisCorrectors()
    modifiedBasis = basis - basisCorrectors
    
    xFull = np.zeros(NpCoarse)
    xFull[free] = xFree
    uLodCoarse = xFull
    uLodFine = modifiedBasis*xFull
    
    return uLodCoarse, uLodFine
    
# Example from Peterseim, Variational Multiscale Stabilization and the Exponential Decay of correctors, p. 2
# Two modifications: A with minus and u(here) = 1/4*u(paper).
fine = 1024
NFine = np.array([fine])
NpFine = np.prod(NFine+1)
NList = [2,4,8,16,32]

epsilon = 2**(-5)

pi = np.pi
xt = util.tCoordinates(NFine).flatten() #whats that mean values of the intervalls xp
xp = util.pCoordinates(NFine).flatten()

aFine = (2 - np.cos(2*pi*xt/epsilon))**(-1)

uSol  = 4*(xp - xp**2) - 4*epsilon*(1/(4*pi)*np.sin(2*pi*xp/epsilon) -
                                    1/(2*pi)*xp*np.sin(2*pi*xp/epsilon) -
                                    epsilon/(4*pi**2)*np.cos(2*pi*xp/epsilon) +
                                    epsilon/(4*pi**2))

uSol = uSol/4

plt.figure('Coefficient')
plt.plot(xt,aFine, label='$A_{\epsilon}(x)$')
plt.yticks((0,np.max(aFine)+np.min(aFine)),fontsize="small")
plt.ylabel('$y$', fontsize="small")
plt.xlabel('$x$', fontsize="small")
plt.legend(frameon=False,fontsize="large")

newErrorFine = []
x = []
y = []


for k in range(1,5):
    newErrorFine = []
    x = []
    y = []
    
    for N in NList:
        NWorldCoarse = np.array([N])
        boundaryConditions = np.array([[0, 0]])

        NCoarseElement = NFine/NWorldCoarse
        world = World(NWorldCoarse, NCoarseElement, boundaryConditions)
        AFine = fem.assemblePatchMatrix(NFine, world.ALocFine, aFine)
        #grid nodes
        xpCoarse = util.pCoordinates(NWorldCoarse).flatten()
        NpCoarse = np.prod(NWorldCoarse+1)
    
        f = np.ones(NpCoarse)
        uCoarseFull, uLodCoarse = PGsolver(world,aFine,f,k)
    
    
    
        newErrorFine.append(np.sqrt(np.dot(uSol - uLodCoarse, AFine*(uSol - uLodCoarse))))
        x.append(N)
        y.append(1./N)
    
    #extension of the plot
    x1 = []
    y = []
    for i in [2,4,8,16,32,64,128,256,512]:
        x1.append(i)
        y.append(1./i)
    
    plt.figure("Error")
    if k == 1:
       plt.loglog(x,newErrorFine,'--o', basex=2, basey=2, label = 'PG-LOD $k=1$')
       #plt.loglog(x1,y,'--k',basex=2, basey=2, linewidth=1, alpha=0.3)
    if k == 2:
        plt.loglog(x,newErrorFine,'-o', basex=2, basey=2, label = 'PG-LOD $k=2$')
        #plt.loglog(x1,y,'--k',basex=2, basey=2, linewidth=1, alpha=0.3)
    if k == 3:
       plt.loglog(x,newErrorFine,'-v', basex=2, basey=2, label = 'PG-LOD $k=3$')
       plt.loglog(x1,y,'--k',basex=2, basey=2, linewidth=1, alpha=0.3)
    if k == 4:
        plt.loglog(x,newErrorFine,'--<', basex=2, basey=2, label = 'PG-LOD $k=4$')
    if k == 5:
        plt.loglog(x,newErrorFine,'--*', basex=2, basey=2, label = '$k=5$')
    if k == 6:
        plt.loglog(x,newErrorFine,'-o', basex=2, basey=2, label = '$k=6$')

    plt.grid(True,which="both",ls="--")
    plt.ylabel('Energy error')
    plt.xlabel('$1/H$')
    plt.title('Energy error for FEM and PG-LOD')
    #plt.legend(frameon=False,fontsize="small") #Legende

# Example from Peterseim, Variational Multiscale Stabilization and the Exponential Decay of correctors, p. 2
# Two modifications: A with minus and u(here) = 1/4*u(paper).
fine = 4096
NFine = np.array([fine])
NpFine = np.prod(NFine+1)
NList = [2,4,8,16, 32, 64, 128, 256, 512]
#NList = [10]


pi = np.pi
xt = util.tCoordinates(NFine).flatten() #whats that mean values of the intervalls xp
xp = util.pCoordinates(NFine).flatten()

aFine = (2 - np.cos(2*pi*xt/epsilon))**(-1)

uSol  = 4*(xp - xp**2) - 4*epsilon*(1/(4*pi)*np.sin(2*pi*xp/epsilon) -
                                    1/(2*pi)*xp*np.sin(2*pi*xp/epsilon) -
                                    epsilon/(4*pi**2)*np.cos(2*pi*xp/epsilon) +
                                    epsilon/(4*pi**2))

uSol = uSol/4


newErrorFine = []
x = []
y = []
for N in NList:
    NWorldCoarse = np.array([N])
    boundaryConditions = np.array([[0, 0]])

    NCoarseElement = NFine/NWorldCoarse
    world = World(NWorldCoarse, NCoarseElement, boundaryConditions)
    AFine = fem.assemblePatchMatrix(NFine, world.ALocFine, aFine)
    
    #grid nodes
    xpCoarse = util.pCoordinates(NWorldCoarse).flatten()
    NpCoarse = np.prod(NWorldCoarse+1)
    f = np.ones(NpCoarse)
    uCoarseFull = femsolverCoarse.solveCoarse_fem(world, aFine, f, boundaryConditions)
    
    basis = fem.assembleProlongationMatrix(NWorldCoarse, NCoarseElement)
    uLodCoarse = basis*uCoarseFull
    newErrorFine.append(np.sqrt(np.dot(uSol - uLodCoarse, AFine*(uSol - uLodCoarse))))
    x.append(N)
    y.append(1./N)

plt.figure("Error")
plt.loglog(x,newErrorFine,'*-', basex=2, basey=2, label= 'FEM')
plt.loglog(x,y,'--k',basex=2, basey=2, linewidth=1, alpha=0.3)
plt.grid(True,which="both",ls="--")
plt.legend(fontsize="small") #Legende
    
plt.show()
```
# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
```

import numpy as np
import scipy.sparse as sparse

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MultipleLocator
from matplotlib import cm

from gridlod import util, world, fem, coef, interp
from gridlod.world import World
import pg_rand
import femsolverCoarse
import buildcoef2d

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
NList = [2,4,8,16,32,64]

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


for k in range(2,5):
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
        if k == 4:
            if np.size(x)==1:
                plt.figure('FEM-Solutions')
                plt.subplots_adjust(left=0.01,bottom=0.04,right=0.99,top=0.95,wspace=0,hspace=0.2)
                plt.subplot(231)
            elif np.size(x)==2:
                plt.subplot(232)
            elif np.size(x)==3:
                plt.subplot(233)
            elif np.size(x)==4:
                plt.subplot(234)
            elif np.size(x)==5:
                plt.subplot(235)
            elif np.size(x)==6:
                plt.subplot(236)
            
            plt.plot(xp,uSol,'k', label='$u_{\epsilon}(x)$')
            plt.plot(xpCoarse,uCoarseFull,'o--', label= '$u^{PG}(x)$')
            plt.title('1/H= ' + str(N),fontsize="small")
            plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
            plt.legend(frameon=False,fontsize="small")
            
    #extension of the plot
    x1 = []
    y = []
    for i in [2,4,8,16,32,64,128,256,512,1024]:
        x1.append(i)
        y.append(1./i)
    
    plt.figure("Error")
    if k == 1:
        plt.loglog(x,newErrorFine,'-o', basex=2, basey=2, label = '$k=1$')
        plt.loglog(x1,y,'--k',basex=2, basey=2, linewidth=1, alpha=0.3)
    if k == 2:
        plt.loglog(x,newErrorFine,'-o', basex=2, basey=2, label = '$k=2$')
    if k == 3:
        plt.loglog(x,newErrorFine,'-*', basex=2, basey=2, label = '$k=3$')
    if k == 4:
        plt.loglog(x,newErrorFine,'--o', basex=2, basey=2, label = '$k=4$')
        plt.loglog(x1,y,'--k',basex=2, basey=2, linewidth=1, alpha=0.3)
    if k == 5:
        plt.loglog(x,newErrorFine,'-o', basex=2, basey=2, label = '$k=5$')
    if k == 6:
        plt.loglog(x,newErrorFine,'-o', basex=2, basey=2, label = '$k=6$')

    plt.grid(True,which="both",ls="--")
    plt.ylabel('Error', fontsize="small")
    plt.xlabel('1/H', fontsize="small")
    plt.legend(frameon=False,fontsize="small") #Legende
    
     
plt.show()
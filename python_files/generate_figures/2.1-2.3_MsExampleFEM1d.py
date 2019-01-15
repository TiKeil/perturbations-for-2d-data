# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)


import numpy as np
import matplotlib.pyplot as plt

from gridlod import util, world, fem
from gridlod.world import World
import femsolverCoarse, buildcoef2d

# Example from Peterseim, Variational Multiscale Stabilization and the Exponential Decay of correctors, p. 2
# Two modifications: A with minus and u(here) = 1/4*u(paper).
fine = 4096
NFine = np.array([fine])
NpFine = np.prod(NFine+1)
NList = [2,4,8,16, 32, 64, 128, 256]

epsilon = 2**(-5)
epsilon1 = 2**(-6)

pi = np.pi
xt = util.tCoordinates(NFine).flatten()
xp = util.pCoordinates(NFine).flatten()

aFine = (2 - np.cos(2*pi*xt/epsilon))**(-1)
aFine1 = (2 - np.cos(2*pi*xt/epsilon1))**(-1)

uSol  = 4*(xp - xp**2) - 4*epsilon*(1/(4*pi)*np.sin(2*pi*xp/epsilon) -
                                    1/(2*pi)*xp*np.sin(2*pi*xp/epsilon) -
                                    epsilon/(4*pi**2)*np.cos(2*pi*xp/epsilon) +
                                    epsilon/(4*pi**2))

uSol = uSol/4

#plot1
plt.figure('Coefficient')
plt.plot(xt,aFine, label='$A_{\epsilon}(x)$')
plt.yticks((0,np.max(aFine)+np.min(aFine)),fontsize=16)
plt.xticks(fontsize=16)
plt.ylabel('$y$', fontsize=16)
plt.xlabel('$x$', fontsize=16)
plt.legend(frameon=False,fontsize=16)

#plot2
plt.figure('Coefficient1')
plt.plot(xt,aFine1, label='$A_{\epsilon}(x)$')
plt.yticks((0,np.max(aFine)+np.min(aFine1)),fontsize=16)
plt.xticks(fontsize=16)
plt.ylabel('$y$', fontsize=16)
plt.xlabel('$x$', fontsize=16)
plt.legend(frameon=False,fontsize=16)

newErrorFine = []
x = []
y = []
for N in NList:
    NWorldCoarse = np.array([N])
    boundaryConditions = np.array([[0, 0]])

    NCoarseElement = NFine//NWorldCoarse
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
    if np.size(x)==1:
        plt.figure('FEM-Solutions')
        plt.subplots_adjust(left=0.01,bottom=0.04,right=0.99,top=0.95,wspace=0,hspace=0.2)
        plt.subplot(241)
    elif np.size(x)==2:
        plt.subplot(242)
    elif np.size(x)==3:
        plt.subplot(243)
    elif np.size(x)==4:
        plt.subplot(244)
    elif np.size(x)==5:
        plt.subplot(245)
    elif np.size(x)==6:
        plt.subplot(246)
    elif np.size(x)==7:
        plt.subplot(247)
    elif np.size(x)==8:
        plt.subplot(248)
    
    if np.size(x)<5:
        plt.plot(xp,uSol,'k', label='$u_{\epsilon}(x)$')
        plt.plot(xpCoarse,uCoarseFull,'o--', label= 'u_h(x)')
        plt.title('1/h= ' + str(N),fontsize="small")
        plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
        plt.legend(frameon=False,fontsize="small")
    else:
        plt.plot(xp,uSol,'k', label='$u_{\epsilon}(x)$')
        plt.plot(xpCoarse,uCoarseFull,'--', label= 'u_h(x)')
        plt.title('1/h= ' + str(N),fontsize="small")
        plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
        plt.legend(frameon=False,fontsize="small")

plt.figure("Error")
plt.loglog(x,newErrorFine,'o-', basex=2, basey=2)
plt.loglog(x,y,'--k',basex=2, basey=2, linewidth=1, alpha=0.3)
plt.ylabel('Energy error')
plt.xlabel('$1/h$')
plt.subplots_adjust(left=0.1,bottom=0.1,right=0.98,top=0.95,wspace=0.2,hspace=0.2)
plt.title('Energy error for the standard FEM')
plt.grid(True,which="both",ls="--")

plt.show()

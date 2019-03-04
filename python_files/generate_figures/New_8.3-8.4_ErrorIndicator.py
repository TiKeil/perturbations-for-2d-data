# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import matplotlib.pyplot as plt

from gridlod import interp, coef, lod
from gridlod.world import World, Patch
import buildcoef2d
from visualize import drawCoefficientGrid

from visualize import drawCoefficient

bg = 0.05  # background
val = 1  # values

# fine World
NWorldFine = np.array([128, 128])
NpFine = np.prod(NWorldFine + 1)

# coarse World
NWorldCoarse = np.array([8, 8])
NpCoarse = np.prod(NWorldCoarse + 1)

# ratio between Fine and Coarse
NCoarseElement = NWorldFine // NWorldCoarse

boundaryConditions = np.array([[0, 0],
                               [0, 0]])

print("starting and precomputing ...")
world = World(NWorldCoarse, NCoarseElement, boundaryConditions)

# coefficient
CoefClass = buildcoef2d.Coefficient2d(NWorldFine,
                                      bg=bg,
                                      val=val,
                                      length=8,
                                      thick=8,
                                      space=8,
                                      probfactor=1,
                                      right=1,
                                      down=0,
                                      diagr1=0,
                                      diagr2=0,
                                      diagl1=0,
                                      diagl2=0,
                                      LenSwitch=None,
                                      thickSwitch=None,
                                      equidistant=True,
                                      ChannelHorizontal=None,
                                      ChannelVertical=None,
                                      BoundarySpace=True)

A = CoefClass.BuildCoefficient()
ABase = A.flatten()

plt.figure("OriginalCoefficient")
drawCoefficient(NWorldFine, ABase)

numbers = [8, 34, 37]

D = CoefClass.SpecificVanish(Number=numbers,
                             probfactor=1,
                             PartlyVanish=None,
                             Original=True)

Defect = D.flatten()
plt.figure("Disappearance")
drawCoefficient(NWorldFine, D)

# precomputations
NWorldFine = world.NWorldFine
NWorldCoarse = world.NWorldCoarse
NCoarseElement = world.NCoarseElement

boundaryConditions = world.boundaryConditions
NpFine = np.prod(NWorldFine + 1)
NpCoarse = np.prod(NWorldCoarse + 1)

k = 5


def computeKmsij(TInd):
    patch = Patch(world, k, TInd)
    IPatch = lambda: interp.L2ProjectionPatchMatrix(patch, boundaryConditions)
    aPatch = lambda: coef.localizeCoefficient(patch, ABase)

    correctorsList = lod.computeBasisCorrectors(patch, IPatch, aPatch)
    csi = lod.computeBasisCoarseQuantities(patch, correctorsList, aPatch)
    return patch, correctorsList, csi.Kmsij, csi

def computeIndicators(TInd):
    aPatch = lambda: coef.localizeCoefficient(patchT[TInd], ABase)
    rPatch = lambda: coef.localizeCoefficient(patchT[TInd], Defect)

    epsFine = lod.computeBasisErrorIndicatorFine(patchT[TInd], correctorsListT[TInd], aPatch, rPatch)
    epsCoarse = lod.computeErrorIndicatorCoarseFromCoefficients(patchT[TInd], csiT[TInd].muTPrime,  aPatch, rPatch)
    return epsFine, epsCoarse

# Use mapper to distribute computations (mapper could be the 'map' built-in or e.g. an ipyparallel map)
patchT, correctorsListT, KmsijT, csiT = zip(*map(computeKmsij, range(world.NtCoarse)))

print('compute error indicators')
epsFine, epsCoarse = zip(*map(computeIndicators, range(world.NtCoarse)))

elemente = np.arange(np.prod(NWorldCoarse))
plt.figure("Error indicators")
#plt.plot(elemente, epsFine, label="Fine")
plt.plot(elemente, epsCoarse, label="Coarse")
plt.ylabel('$e_{u,T}$')
plt.xlabel('Element')
plt.subplots_adjust(left=0.09, bottom=0.09, right=0.99, top=0.99, wspace=0.2, hspace=0.2)
plt.legend(loc='upper right')  # Legende
plt.grid()

fig = plt.figure("error indicator")
ax = fig.add_subplot(1,1,1)
np_eps = np.einsum('i,i -> i', np.ones(np.size(epsFine)), epsFine)
drawCoefficientGrid(NWorldCoarse, np_eps,fig,ax, original_style=True)

fig = plt.figure("error indicator coarse")
ax = fig.add_subplot(1,1,1)
np_eps = np.einsum('i,i -> i', np.ones(np.size(epsCoarse)), epsCoarse)
drawCoefficientGrid(NWorldCoarse, np_eps,fig,ax, original_style=True)


plt.show()
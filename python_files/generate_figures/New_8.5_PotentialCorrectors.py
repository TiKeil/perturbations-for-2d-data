# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import matplotlib.pyplot as plt

from gridlod import interp, coef, lod
import buildcoef2d
from gridlod.world import World, Patch

from visualize import drawCoefficient, drawCoefficientGrid, drawCoefficientwg
from gridlod_on_perturbations.visualization_tools import drawCoefficient_origin
bg = 0.05  # background
val = 1  # values

# fine World
NWorldFine = np.array([256, 256])
NpFine = np.prod(NWorldFine + 1)

# coarse World
NWorldCoarse = np.array([16, 16])
NpCoarse = np.prod(NWorldCoarse + 1)

# ratio between Fine and Coarse
NCoarseElement = NWorldFine // NWorldCoarse

boundaryConditions = np.array([[0, 0],
                               [0, 0]])

world = World(NWorldCoarse, NCoarseElement, boundaryConditions)

# righthandside
f = np.ones(NpCoarse)

# coefficient
CoefClass = buildcoef2d.Coefficient2d(NWorldFine,
                                      bg=0.05,
                                      val=1,
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

# potentially updated
fig = plt.figure("Correcting")
ax = fig.add_subplot(2, 3, 1)
ax.set_title('Original', fontsize=7)

drawCoefficientwg(NWorldFine, ABase, fig, ax)

numbers = [2, 70, 97, 153, 205]

D = CoefClass.SpecificVanish(Number=numbers,
                             probfactor=1,
                             PartlyVanish=None,
                             Original=True)

value2 = 50
R2 = CoefClass.SpecificValueChange(ratio=value2,
                                   Number=numbers,
                                   probfactor=1,
                                   randomvalue=None,
                                   negative=None,
                                   ShapeRestriction=True,
                                   ShapeWave=None,
                                   Original=True,
                                   NewShapeChange=True)

# precompute
NWorldFine = world.NWorldFine
NWorldCoarse = world.NWorldCoarse
NCoarseElement = world.NCoarseElement

boundaryConditions = world.boundaryConditions
NpFine = np.prod(NWorldFine + 1)
NpCoarse = np.prod(NWorldCoarse + 1)

ANew = D.flatten()
elemente = np.arange(np.prod(NWorldCoarse))

def computeKmsij(TInd):
    patch = Patch(world, TInd[1], TInd[0])
    IPatch = lambda: interp.L2ProjectionPatchMatrix(patch, boundaryConditions)
    aPatch = lambda: coef.localizeCoefficient(patch, ABase)

    correctorsList = lod.computeBasisCorrectors(patch, IPatch, aPatch)
    csi = lod.computeBasisCoarseQuantities(patch, correctorsList, aPatch)
    return patch, correctorsList, csi.Kmsij, csi

def computeIndicators(TInd):
    aPatch = lambda: coef.localizeCoefficient(patchT[TInd], ABase)
    rPatch = lambda: coef.localizeCoefficient(patchT[TInd], ANew)

    epsFine = lod.computeBasisErrorIndicatorFine(patchT[TInd], correctorsListT[TInd], aPatch, rPatch)
    epsCoarse = 0
    return epsFine, epsCoarse


for k in range(1, 6):
    print(('<<<<<<<<<<<<<<<< ' + str(k) + ' >>>>>>>>>>>>>>>>'))
    patchT, correctorsListT, KmsijT, csiT = zip(*map(computeKmsij, zip(range(world.NtCoarse),np.repeat(k,world.NtCoarse))))
    print('compute error indicators')
    epsFine, epsCoarse = zip(*map(computeIndicators, range(world.NtCoarse)))
    vis = np.copy(epsFine)
    for i in range(np.size(epsFine)):
        if epsFine[i] != 0:
            vis[i] = 1

    fig = plt.figure("Correcting")
    if k == 1:
        ax = fig.add_subplot(2, 3, k + 2)
        ax.set_title('Affected for $k=$' + str(k), fontsize=7)
    elif k == 2:
        ax = fig.add_subplot(2, 3, k + 2)
        ax.set_title('Affected for $k=$' + str(k), fontsize=7)
    elif k == 3:
        ax = fig.add_subplot(2, 3, k + 2)
        ax.set_title('Affected for $k=$' + str(k), fontsize=7)
    elif k == 4:
        ax = fig.add_subplot(2, 3, k + 2)
        ax.set_title('Affected for $k=$' + str(k), fontsize=7)
    if k != 5:
        drawCoefficientGrid(NWorldCoarse, vis, fig, ax, original_style=True, colorbar=False)

fig = plt.figure("Correcting")
ax = fig.add_subplot(2, 3, 2)
ax.set_title('Defects', fontsize=7)
drawCoefficientwg(NWorldFine, D, fig, ax)

plt.show()
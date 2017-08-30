# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np
import scipy.sparse as sparse

import matplotlib.pyplot as plt 
from gridlod import interp, coef, util, fem, world, linalg, femsolver
from gridlod.world import World
import femsolverCoarse
import pg_rand
import buildcoef2d

from visualize import drawCoefficient, d3plotter

bg = 0.05   #background
val = 1     #values

#fine World
NWorldFine = np.array([128, 128])
NpFine = np.prod(NWorldFine+1)                                                                               

#coarse World
NWorldCoarse = np.array([8,8])
NpCoarse = np.prod(NWorldCoarse+1)

#ratio between Fine and Coarse
NCoarseElement = NWorldFine/NWorldCoarse

boundaryConditions = np.array([[0, 0],
                               [0, 0]])

world = World(NWorldCoarse, NCoarseElement, boundaryConditions)

#righthandside
f = np.ones(NpCoarse)

#coefficient
WormsUltCoefClass = buildcoef2d.Coefficient2d(NWorldFine, 
                                        bg                  = bg, 
                                        val                 = val, 
                                        length              = 1, 
                                        thick               = 1, 
                                        space               = 1, 
                                        probfactor          = 1, 
                                        right               = 1, 
                                        down                = 1, 
                                        diagr1              = 1, 
                                        diagr2              = 1, 
                                        diagl1              = 1, 
                                        diagl2              = 1, 
                                        LenSwitch           = [1,2,3,4,5,6], 
                                        thickSwitch         = [1,2,3], 
                                        equidistant         = None, 
                                        ChannelHorizontal   = None, 
                                        ChannelVertical     = None,
                                        BoundarySpace       = True,
                                        Boxes2n             = None,
                                        Channels2n          = None)


A = WormsUltCoefClass.BuildCoefficient()
ABase = A.flatten()

plt.figure("OriginalCoefficient")
drawCoefficient(NWorldFine, ABase)

#fine-fem
f_fine = np.ones(NpFine)
uFineFem, AFine, MFine = femsolver.solveFine(world, ABase, f_fine, None, boundaryConditions)

#PGLOD
#interpolant
IPatchGenerator = lambda i, N: interp.L2ProjectionPatchMatrix(i, N, NWorldCoarse, NCoarseElement, boundaryConditions)

#Coefficient (need flatten form)
aCoef = coef.coefficientFine(NWorldCoarse, NCoarseElement, ABase)

k=4

pglod = pg_rand.VcPetrovGalerkinLOD(aCoef, world, k, IPatchGenerator, 2)
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
diff = basis*uLodCoarse-uLodFine

# plot multiscale splitting
d3plotter(NWorldFine, uFineFem,boundary=[np.min(uFineFem),np.max(uFineFem)])
d3plotter(NWorldCoarse, uLodCoarse, String='Coarse1', boundary=[np.min(uFineFem),np.max(uFineFem)])
d3plotter(NWorldFine, diff, String='uf2',boundary=[np.min(uFineFem),np.max(uFineFem)], Blues= True)

# plot for basis functions
schauen = np.zeros(NpCoarse)
schauen[40] = 1
schau = basis*schauen
schau1 = basisCorrectors*schauen
schau2 = modifiedBasis*schauen
zmax = np.max(schau2)
zmin = np.min(schau2)
d3plotter(NWorldFine, schau, '1', zmax=zmax, zmin=zmin, Blues=True)
d3plotter(NWorldFine, schau1, '2', zmax=zmax, zmin=zmin, Blues=True)
d3plotter(NWorldFine, schau2, '3', zmax=zmax, zmin=zmin)

plt.show()
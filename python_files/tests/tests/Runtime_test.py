# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import os
import sys
import numpy as np
import scipy.sparse as sparse
import random
import csv

import matplotlib.pyplot as plt
from visualize import drawCoefficient
from data import *

from gridlod import interp, coef, util, fem, world, linalg, femsolver
import pg_rand, femsolverCoarse, buildcoef2d
from gridlod.world import World

import timeit

def result(pglod, world, A, R, f, k, String):
    print "-------------- " + String + " ---------------"

    NWorldFine = world.NWorldFine
    NWorldCoarse = world.NWorldCoarse
    NCoarseElement = world.NCoarseElement

    boundaryConditions = world.boundaryConditions
    NpFine = np.prod(NWorldFine + 1)
    NpCoarse = np.prod(NWorldCoarse + 1)

    # new Coefficient
    ANew = R.flatten()
    Anew = coef.coefficientFine(NWorldCoarse, NCoarseElement, ANew)

    start_solver = timeit.default_timer()

    # reference solution
    f_fine = np.ones(NpFine)
    uFineFem, AFine, MFine = femsolver.solveFine(world, ANew, f_fine, None, boundaryConditions)
    print('Runtime: it took {} sec to compute the reference solution'.format(timeit.default_timer() - start_solver))

    start_solve_system = timeit.default_timer()
    # worst solution
    KFull = pglod.assembleMsStiffnessMatrix()
    MFull = fem.assemblePatchMatrix(NWorldCoarse, world.MLocCoarse)
    free = util.interiorpIndexMap(NWorldCoarse)

    bFull = MFull * f
    KFree = KFull[free][:, free]
    bFree = bFull[free]

    xFree = sparse.linalg.spsolve(KFree, bFree)

    basis = fem.assembleProlongationMatrix(NWorldCoarse, NCoarseElement)

    basisCorrectors = pglod.assembleBasisCorrectors()
    modifiedBasis = basis - basisCorrectors

    xFull = np.zeros(NpCoarse)
    xFull[free] = xFree
    uCoarse = xFull
    uLodFine = modifiedBasis * xFull

    print('Runtime: It took {} sec to apply LOD with given correctors'.format(timeit.default_timer()-start_solve_system))
    uLodFineWorst = uLodFine

    # energy error
    errorworst = np.sqrt(np.dot(uFineFem - uLodFineWorst, AFine * (uFineFem - uLodFineWorst)))

    start_runtime = timeit.default_timer()
    # tolerance = 0
    vis, eps = pglod.updateCorrectors(Anew, 0, f, 1, clearFineQuantities=False, Computing=False)

    print('Runtime: It took {} sec to compute the error indicators.'.format(timeit.default_timer()-start_runtime))

    PotentialCorrectors = np.sum(vis)
    elemente = np.arange(np.prod(NWorldCoarse))

    # identify tolerances
    epsnozero = filter(lambda x: x != 0, eps)

    assert (np.size(epsnozero) != 0)

    mini = np.min(epsnozero)
    minilog = int(round(np.log10(mini) - 0.49))
    epsnozero.append(10 ** (minilog))
    ToleranceListcomplete = []
    for i in range(0, int(np.size(epsnozero))):
        ToleranceListcomplete.append(epsnozero[i])

    ToleranceListcomplete.sort()
    ToleranceListcomplete = np.unique(ToleranceListcomplete)

    # with tolerance
    errorplotinfo = []
    tolerancesafe = []
    errorBest = []
    errorWorst = []
    recomputefractionsafe = []
    recomputefraction = 0
    Correctors = 0
    runtime = []
    total_time = 0
    leng = np.size(ToleranceListcomplete)
    for k in range(leng - 1, -1, -1):
        tol = ToleranceListcomplete[k]
        print " --- " + str(-k + leng) + "/" + str(leng) + " --- Tolerance: " + str(
            round(tol, 5)) + " in " + String + " ---- ",

        start_runtime = timeit.default_timer()
        vistol, time_to_compute = pglod.updateCorrectors(Anew, tol, f, clearFineQuantities=False, Testing=True, runtime=True)
        total_time += timeit.default_timer()-start_runtime- time_to_compute
        runtime.append(total_time)

        Correctors += np.sum(vistol)

        recomputefraction += float(np.sum(vistol)) / PotentialCorrectors * 100
        recomputefractionsafe.append(recomputefraction)

        KFull = pglod.assembleMsStiffnessMatrix()
        MFull = fem.assemblePatchMatrix(NWorldCoarse, world.MLocCoarse)
        free = util.interiorpIndexMap(NWorldCoarse)

        bFull = MFull * f
        KFree = KFull[free][:, free]
        bFree = bFull[free]

        xFree = sparse.linalg.spsolve(KFree, bFree)
        basis = fem.assembleProlongationMatrix(NWorldCoarse, NCoarseElement)

        basisCorrectors = pglod.assembleBasisCorrectors()

        modifiedBasis = basis - basisCorrectors

        xFull = np.zeros(NpCoarse)
        xFull[free] = xFree
        uCoarse = xFull
        uLodFine = modifiedBasis * xFull

        # energy error
        errortol = np.sqrt(np.dot(uFineFem - uLodFine, AFine * (uFineFem - uLodFine)))

        errorplotinfo.append(errortol)
        tolerancesafe.append(tol)

    # 100% updating
    uLodFinebest = uLodFine
    errorbest = np.sqrt(np.dot(uFineFem - uLodFinebest, AFine * (uFineFem - uLodFinebest)))

    print('Runtime: Total time of updating: {} sec.'.format(total_time))
    for k in range(leng - 1, -1, -1):
        errorBest.append(errorbest)
        errorWorst.append(errorworst)

    return vis, eps, PotentialCorrectors, recomputefractionsafe, errorplotinfo, errorWorst, errorBest, runtime


bg = 0.05  # background
val = 1  # values

# fine World
NWorldFine = np.array([256, 256])
NpFine = np.prod(NWorldFine + 1)

# coarse World
NWorldCoarse = np.array([16, 16])
NpCoarse = np.prod(NWorldCoarse + 1)

# ratio between Fine and Coarse
NCoarseElement = NWorldFine / NWorldCoarse

boundaryConditions = np.array([[0, 0],
                               [0, 0]])

world = World(NWorldCoarse, NCoarseElement, boundaryConditions)

# righthandside
f = np.ones(NpCoarse)

#Coefficient 1
CoefClass = buildcoef2d.Coefficient2d(NWorldFine,
                        bg                  = bg,
                        val                 = val,
                        length              = 2,
                        thick               = 2,
                        space               = 2,
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

ROOT = '../../../test_data/Coef1/'

# safe NworldFine
with open("%s/NWorldFine.txt" % ROOT, 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for val in NWorldFine:
        writer.writerow([val])

# safe NworldCoarse
with open("%s/NWorldCoarse.txt" % ROOT, 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for val in NWorldCoarse:
        writer.writerow([val])

# ABase
with open("%s/OriginalCoeff.txt" % ROOT, 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for val in ABase:
        writer.writerow([val])

# fine-fem
f_fine = np.ones(NpFine)

uFineFem, AFine, MFine = femsolver.solveFine(world, ABase, f_fine, None, boundaryConditions)

# fine solution
with open("%s/finescale.txt" % ROOT, 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for val in uFineFem:
        writer.writerow([val])

plt.figure("Original")
drawCoefficient(NWorldFine, ABase, greys=True)
plt.title("Original coefficient")
plt.show()

# random seed
random.seed(20)

# decision
valc = np.shape(CoefClass.ShapeRemember)[0]
numbers = []
decision = np.zeros(100)
decision[0] = 1

for i in range(0, valc):
    a = random.sample(decision, 1)[0]
    if a == 1:
        numbers.append(i)

value1 = 3
C1 = CoefClass.SpecificValueChange(ratio=value1,
                                   Number=numbers,
                                   probfactor=1,
                                   randomvalue=None,
                                   negative=None,
                                   ShapeRestriction=True,
                                   ShapeWave=None,
                                   ChangeRight=1,
                                   ChangeDown=1,
                                   ChangeDiagr1=1,
                                   ChangeDiagr2=1,
                                   ChangeDiagl1=1,
                                   ChangeDiagl2=1,
                                   Original=True,
                                   NewShapeChange=True)

V = CoefClass.SpecificVanish(Number=numbers,
                             probfactor=1,
                             PartlyVanish=None,
                             ChangeRight=1,
                             ChangeDown=1,
                             ChangeDiagr1=1,
                             ChangeDiagr2=1,
                             ChangeDiagl1=1,
                             ChangeDiagl2=1,
                             Original=True)

M1 = CoefClass.SpecificMove(probfactor=1,
                            Number=numbers,
                            steps=1,
                            randomstep=None,
                            randomDirection=None,
                            ChangeRight=1,
                            ChangeDown=1,
                            ChangeDiagr1=1,
                            ChangeDiagr2=1,
                            ChangeDiagl1=1,
                            ChangeDiagl2=1,
                            Right=1,
                            BottomRight=0,
                            Bottom=0,
                            BottomLeft=0,
                            Left=0,
                            TopLeft=0,
                            Top=0,
                            TopRight=0,
                            Original=True)

k = 4

NWorldFine = world.NWorldFine
NWorldCoarse = world.NWorldCoarse
NCoarseElement = world.NCoarseElement

boundaryConditions = world.boundaryConditions
NpFine = np.prod(NWorldFine + 1)
NpCoarse = np.prod(NWorldCoarse + 1)

# interpolant
IPatchGenerator = lambda i, N: interp.L2ProjectionPatchMatrix(i, N, NWorldCoarse, NCoarseElement, boundaryConditions)

# old Coefficient
ABase = A.flatten()
Aold = coef.coefficientFine(NWorldCoarse, NCoarseElement, ABase)

pglod = pg_rand.VcPetrovGalerkinLOD(Aold, world, k, IPatchGenerator, 0)

start = timeit.default_timer()
pglod.originCorrectors(clearFineQuantities=False)
stop = timeit.default_timer()

print('Runtime: It took {} seconds to compute the reference configuration.'.format(stop-start))

vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest, runtime = result(pglod, world, A, V, f,
                                                                                                            k, 'Vanish')

# runtime
with open("%s/runtime.txt" % ROOT, 'wb') as csvfile:
    writer = csv.writer(csvfile)
    for val in runtime:
        writer.writerow([val])
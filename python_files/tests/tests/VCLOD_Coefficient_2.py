# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import scipy.sparse as sparse
import random

import matplotlib.pyplot as plt
from visualize import drawCoefficient
from data import * 

from gridlod import interp, coef, util, fem, femsolver
import pg_rand, buildcoef2d
from gridlod.world import World

def result(pglod, world, A, R, f, k, String):
    print("-------------- " + String + " ---------------") 
    NWorldFine = world.NWorldFine
    NWorldCoarse = world.NWorldCoarse
    NCoarseElement = world.NCoarseElement
    
    boundaryConditions = world.boundaryConditions
    NpFine = np.prod(NWorldFine+1)
    NpCoarse = np.prod(NWorldCoarse+1)
        
    # new Coefficient
    ANew = R.flatten()
    Anew = coef.coefficientFine(NWorldCoarse, NCoarseElement, ANew)
    
    # reference solution
    f_fine = np.ones(NpFine)
    uFineFem, AFine, MFine = femsolver.solveFine(world, ANew, f_fine, None, boundaryConditions)
    
    # worst solution
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
    uCoarse = xFull
    uLodFine = modifiedBasis*xFull
    
    uLodFineWorst = uLodFine
    
    # energy error
    errorworst = np.sqrt(np.dot(uFineFem - uLodFineWorst, AFine*(uFineFem - uLodFineWorst)))
    
    # tolerance = 0 
    vis, eps = pglod.updateCorrectors(Anew, 0, f, 1, clearFineQuantities=False, Computing=False)
    
    PotentialCorrectors = np.sum(vis)
    elemente = np.arange(np.prod(NWorldCoarse))
            
    # identify tolerances
    epsnozero = [x for x in eps if x!=0]
    
    assert(np.size(epsnozero) != 0)
    
    mini = np.min(epsnozero)
    minilog = int(round(np.log10(mini)-0.49))
    epsnozero.append(10**(minilog))
    ToleranceListcomplete = []
    for i in range(0,int(np.size(epsnozero))):
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
    leng = np.size(ToleranceListcomplete)
    for k in range(leng-1,-1,-1):
        tol = ToleranceListcomplete[k]
        print(" --- "+ str(-k+leng) + "/" + str(leng)+ " --- Tolerance: " + str(round(tol,5)) + " in "+ String +" ---- ", end=' ') 
        vistol = pglod.updateCorrectors(Anew, tol, f, clearFineQuantities=False, Testing=True)
        
        Correctors += np.sum(vistol)
        
        recomputefraction += float(np.sum(vistol))/PotentialCorrectors * 100
        recomputefractionsafe.append(recomputefraction)
        
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
        uCoarse = xFull
        uLodFine = modifiedBasis*xFull
        
        #energy error
        errortol = np.sqrt(np.dot(uFineFem - uLodFine, AFine*(uFineFem - uLodFine)))
        
        errorplotinfo.append(errortol)
        tolerancesafe.append(tol)
    
    # 100% updating
    uLodFinebest = uLodFine
    errorbest = np.sqrt(np.dot(uFineFem - uLodFinebest, AFine*(uFineFem - uLodFinebest)))
    
    for k in range(leng-1,-1,-1):
        errorBest.append(errorbest)
        errorWorst.append(errorworst)

    return vis, eps, PotentialCorrectors, recomputefractionsafe, errorplotinfo, errorWorst, errorBest

bg = 0.05       #background
val = 1         #values

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

#Coefficient 2
CoefClass = buildcoef2d.Coefficient2d(NWorldFine, 
                        bg                  = bg,
                        val                 = val,
                        length              = 1,
                        thick               = 1,
                        space               = 2,
                        probfactor          = 1,
                        right               = 0,
                        down                = 0,
                        diagr1              = 0,
                        diagr2              = 0,
                        diagl1              = 0,
                        diagl2              = 0,
                        LenSwitch           = None,
                        thickSwitch         = None,
                        equidistant         = True,
                        ChannelHorizontal   = None,
                        ChannelVertical     = True,
                        BoundarySpace       = True)
                        
A = CoefClass.BuildCoefficient()
ABase = A.flatten()

ROOT = '../../../test_data/Coef2/'

#safe NworldFine
with open("%s/NWorldFine.txt" % ROOT, 'w') as csvfile:
    writer = csv.writer(csvfile)
    for val in NWorldFine:
        writer.writerow([val])

#safe NworldCoarse
with open("%s/NWorldCoarse.txt" % ROOT, 'w') as csvfile:
    writer = csv.writer(csvfile)
    for val in NWorldCoarse:
        writer.writerow([val])

#ABase
with open("%s/OriginalCoeff.txt" % ROOT, 'w') as csvfile:
    writer = csv.writer(csvfile)
    for val in ABase:
        writer.writerow([val])

#fine-fem
f_fine = np.ones(NpFine)
uFineFem, AFine, MFine = femsolver.solveFine(world, ABase, f_fine, None, boundaryConditions)

#fine solution
with open("%s/finescale.txt" % ROOT, 'w') as csvfile:
    writer = csv.writer(csvfile)
    for val in uFineFem:
        writer.writerow([val])
        
plt.figure("Original")
drawCoefficient(NWorldFine, ABase,greys=True)
plt.title("Original coefficient")
plt.show()

# random seed
random.seed(20)

# decision
valc = np.shape(CoefClass.ShapeRemember)[0]
numbers = []
decision = np.zeros(100)
decision[0] = 1


for i in range(0,valc):
    a = random.sample(list(decision),1)[0]
    if a == 1:
        numbers.append(i)

value1 = 3
C1 = CoefClass.SpecificValueChange(ratio=value1,
                                    Number = numbers,
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
                                    Original = True,
                                    NewShapeChange = True)

V = CoefClass.SpecificVanish(Number = numbers,
                                probfactor=1,
                                PartlyVanish=None,
                                ChangeRight=1,
                                ChangeDown=1,
                                ChangeDiagr1=1,
                                ChangeDiagr2=1,
                                ChangeDiagl1=1,
                                ChangeDiagl2=1,
                                Original = True)

M1 = CoefClass.SpecificMove(probfactor=1,
                            Number = numbers,
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
                            Original = True)
                            
k = 4

NWorldFine = world.NWorldFine
NWorldCoarse = world.NWorldCoarse
NCoarseElement = world.NCoarseElement

boundaryConditions = world.boundaryConditions
NpFine = np.prod(NWorldFine+1)
NpCoarse = np.prod(NWorldCoarse+1)

#interpolant
IPatchGenerator = lambda i, N: interp.L2ProjectionPatchMatrix(i, N, NWorldCoarse, NCoarseElement, boundaryConditions)

#old Coefficient
ABase = A.flatten()
Aold = coef.coefficientFine(NWorldCoarse, NCoarseElement, ABase)

pglod = pg_rand.VcPetrovGalerkinLOD(Aold, world, k, IPatchGenerator, 0)
pglod.originCorrectors(clearFineQuantities=False)

vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest = result(pglod ,world, A, C1, f, k, 'Specific value change' + str(value1))

safeChange(ROOT, C1, vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest)

vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest = result(pglod ,world, A, V, f, k, 'Vanish')

safeVanish(ROOT, V, vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest)

vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest = result(pglod ,world, A, M1, f, k, 'One Step Move')

safeShift(ROOT, M1, vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest)
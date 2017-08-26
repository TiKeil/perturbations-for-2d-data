```
# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
```

"""
DEPRECATED
"""

import numpy as np

import scipy.sparse as sparse
import scipy.stats as stats
import random
import csv

from gridlod import pg_rand, interp, coef, util, fem, world, linalg, femsolverCoarse, buildcoef2d
from gridlod.world import World

def result(pglod, world, A, R, f, k, String):
    print "------------------------------------- " + String + " -------------------------------------------" 
    NWorldFine = world.NWorldFine
    NWorldCoarse = world.NWorldCoarse
    NCoarseElement = world.NCoarseElement
    
    boundaryConditions = world.boundaryConditions
    NpFine = np.prod(NWorldFine+1)
    NpCoarse = np.prod(NWorldCoarse+1)
        
    #new Coefficient
    ANew = R.flatten()
    Anew = coef.coefficientFine(NWorldCoarse, NCoarseElement, ANew)
    
    ##################### worst solution ############################
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
    
    ################ tolerance = 0 #############
    vis, eps = pglod.updateCorrectors(Anew, 0, f, 1, clearFineQuantities=False)
    
    PotentialCorrectors = np.sum(vis)
    
    elemente = np.arange(np.prod(NWorldCoarse))
        
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
    
    uLodFinebest = uLodFine
    
    #search for tolerances
    epsnozero = filter(lambda x: x!=0, eps)
    
    assert(np.size(epsnozero) != 0)
    
    mini = np.min(epsnozero)
    minilog = int(round(np.log10(mini)-0.49))
    epsnozero.append(10**(minilog))
    ToleranceListcomplete = []
    for i in range(0,int(np.size(epsnozero))):
        ToleranceListcomplete.append(epsnozero[i])

    ToleranceListcomplete.sort()
    ToleranceListcomplete = np.unique(ToleranceListcomplete)

    ############ with tolerance ################
    errorplotinfo = []
    tolerancesafe = []
    errorbest = []
    errorworst = []
    recomputefractionsafe = []
    recomputefraction = 0
    Correctors = 0
    leng = np.size(ToleranceListcomplete)
    for k in range(leng-1,-1,-1):
        tol = ToleranceListcomplete[k]
        print " --------- "+ str(-k+leng) + "/" + str(leng)+ " ------------------- Tolerance: " + str(round(tol,5)) + " in "+ String +" ------------------------------------------- ", 
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
        
        ErrorFine, ErrorBest, ErrorDifference, ErrorFineWorst = Errors(world, R, uLodFinebest, uLodFine, uLodFineWorst)
        
        errorplotinfo.append(ErrorFine)
        tolerancesafe.append(tol)
        errorbest.append(ErrorBest)
        errorworst.append(ErrorFineWorst)
        
    # plt.figure("Errors for " + String)
    # plt.semilogy(recomputefractionsafe,errorplotinfo, label = 'WsLOD error')
    # plt.semilogy(recomputefractionsafe,errorworst, label = 'Worst LOD error')
    # plt.semilogy(recomputefractionsafe,errorbest, label = 'Best LOD error')
    # ymin, ymax = plt.ylim()
    # plt.subplots_adjust(left=0.20)
    # #plt.ticker.set_major_formatter(FormatStrFormatter('%.02f'))
    # plt.yticks((ymin,ymax),fontsize="small")
    #
    # plt.ylabel('Error', fontsize="small")
    # plt.xlabel('Updated correctors in %')
    # plt.legend(frameon=False,fontsize="small") #Legende
    
    # plt.figure("Errors")
    # plt.loglog(tolerancesafe,errorplotinfo, label = 'indicator ' + String)
    # plt.loglog(tolerancesafe,errorworst, label = 'worst ' + String)
    # plt.loglog(tolerancesafe,errorbest, label = 'best ' + String)
    # plt.ylabel('Error')
    # plt.xlabel('Tolerance')
    # plt.legend(loc='upper left', frameon=False,fontsize="small") #Legende
    
    return vis, eps, PotentialCorrectors, recomputefractionsafe, errorplotinfo, errorworst, errorbest
    
def Errors(world, R, uLodFinebest, uLodFinetol, uLodFineworst):
    NWorldFine = world.NWorldFine
    NpFine = np.prod(NWorldFine+1)
    
    ANew = R.flatten()
    boundaryConditions = world.boundaryConditions
    f_fine = np.ones(NpFine)
    uFineFem, AFine, MFine = femsolver.solveFine(world, ANew, f_fine, None, boundaryConditions)
    
    errorbest = np.sqrt(np.dot(uFineFem - uLodFinebest, AFine*(uFineFem - uLodFinebest)))
    errortol = np.sqrt(np.dot(uFineFem - uLodFinetol, AFine*(uFineFem - uLodFinetol)))
    errorworst = np.sqrt(np.dot(uFineFem - uLodFineworst, AFine*(uFineFem - uLodFineworst)))
    
    return errortol, errorbest, errortol-errorbest, errorworst

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

differentcoefs = 4

#path
pfade = ['/test_data/Coef1', '/test_data/Coef2', '/test_data/Coef3', '/test_data/Coef4']

#Coefficient 1
CoefClass1 = buildcoef2d.Coefficient2d(NWorldFine,
                        bg	    	        = bg,
                        val		            = val,
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

#Coefficient 2
CoefClass2 = buildcoef2d.Coefficient2d(NWorldFine, 
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

#Coefficient 3
CoefClass3 = buildcoef2d.Coefficient2d(NWorldFine, 
                        bg                  = bg, 
                        val                 = val, 
                        length              = 8, 
                        thick               = 2, 
                        space               = 4, 
                        probfactor          = 1, 
                        right               = 0, 
                        down                = 0, 
                        diagr1              = 0, 
                        diagr2              = 0, 
                        diagl1              = 0, 
                        diagl2              = 1, 
                        LenSwitch           = None, 
                        thickSwitch         = None, 
                        equidistant         = None, 
                        ChannelHorizontal   = None, 
                        ChannelVertical     = None,
                        BoundarySpace       = True)

#Coefficient 4
CoefClass4 = buildcoef2d.Coefficient2d(NWorldFine, 
                        bg                  = bg, 
                        val                 = val, 
                        thick               = 1, 
                        space               = 0, 
                        probfactor          = 1, 
                        right               = 0, 
                        down                = 0, 
                        diagr1              = 1, 
                        diagr2              = 0, 
                        diagl1              = 1, 
                        diagl2              = 0, 
                        LenSwitch           = [4,5,6,7,8], 
                        thickSwitch         = None, 
                        equidistant         = None, 
                        ChannelHorizontal   = None, 
                        ChannelVertical     = None,
                        BoundarySpace       = None) 


classes = [ CoefClass1,
            CoefClass2,
            CoefClass3,
            CoefClass4]

for s in range(0,4):
    print '<<<<<<<<<<<<<<<<< ' + str(s+1) + ' >>>>>>>>>>>>>>>>>>>>>>'
    ROOT = pfade[s]
    CoefClass = classes[s]
    
    A = CoefClass.BuildCoefficient()
    ABase = A.flatten()
    
    #safe NworldFine
    with open("%s/NWorldFine.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in NWorldFine:
            writer.writerow([val])
    
    #safe NworldCoarse
    with open("%s/NWorldCoarse.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in NWorldCoarse:
            writer.writerow([val])
    
            
    #ABase
    with open("%s/OriginalCoeff.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in ABase:
            writer.writerow([val])
    
        
    #fine-fem
    f_fine = np.ones(NpFine)
    uFineFem, AFine, MFine = femsolver.solveFine(world, ABase, f_fine, None, boundaryConditions)
    
    #fine solution
    with open("%s/finescale.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in uFineFem:
            writer.writerow([val])

    #decision whether randomize
    valc = np.shape(CoefClass.ShapeRemember)[0]
    numbers = []
    decision = np.zeros(100)
    decision[0] = 1
    
    # random seed
    random.seed(20)
    
    for i in range(0,valc):
        a = random.sample(decision,1)[0]
        if a == 1:
            numbers.append(i)

    print numbers
    #numbers = None

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

    ################################ precompute ###########################
    NWorldFine = world.NWorldFine
    NWorldCoarse = world.NWorldCoarse
    NCoarseElement = world.NCoarseElement

    boundaryConditions = world.boundaryConditions
    NpFine = np.prod(NWorldFine+1)
    NpCoarse = np.prod(NWorldCoarse+1)

    #interpolant
    IPatchGenerator = lambda i, N: interp.L2ProjectionPatchMatrix(i, N, NWorldCoarse, NCoarseElement, boundaryConditions)

    #old Coefficient (need flatten form)
    ABase = A.flatten()
    Aold = coef.coefficientFine(NWorldCoarse, NCoarseElement, ABase)

    pglod = pg_rand.VcPetrovGalerkinLOD(Aold, world, k, IPatchGenerator, 0)
    pglod.originCorrectors(clearFineQuantities=False)
    #pglod.originRhsCorrectors(clearFineQuantities=False)

    ################# Value Changes 1 ###########################
    vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest = result(pglod ,world, A, C1, f, k, 'Specific value change' + str(value1))
    #safe vis
    C1Base = C1.flatten()
    with open("%s/C1Base.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in C1Base:
            writer.writerow([val])
    
    
    with open("%s/C1vis.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in vis:
            writer.writerow([val])
    
    #safe eps
    with open("%s/C1eps.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in eps:
            writer.writerow([val])
    
    #safe PotentialUpdated
    with open("%s/C1PotentialUpdated.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in [PotentialUpdated]:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/C1recomputefractionsafe.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in recomputefractionsafe:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/C1errorplotinfo.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorplotinfo:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/C1errorworst.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorworst:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/C1errorbest.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorbest:
            writer.writerow([val])
                
    vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest = result(pglod ,world, A, V, f, k, 'Vanish')
    
    VBase = V.flatten()
    with open("%s/VBase.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in VBase:
            writer.writerow([val])
    
    
    #safe vis
    with open("%s/Vvis.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in vis:
            writer.writerow([val])
    
    #safe eps
    with open("%s/Veps.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in eps:
            writer.writerow([val])
    
    #safe PotentialUpdated
    with open("%s/VPotentialUpdated.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in [PotentialUpdated]:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/Vrecomputefractionsafe.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in recomputefractionsafe:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/Verrorplotinfo.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorplotinfo:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/Verrorworst.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorworst:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/Verrorbest.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorbest:
            writer.writerow([val])
    
    vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest = result(pglod ,world, A, M1, f, k, 'One Step Move')
    #safe vis
    
    M1Base = M1.flatten()
    with open("%s/M1Base.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in M1Base:
            writer.writerow([val])
    
    
    with open("%s/M1vis.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in vis:
            writer.writerow([val])
    
    #safe eps
    with open("%s/M1eps.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in eps:
            writer.writerow([val])
    
    #safe PotentialUpdated
    with open("%s/M1PotentialUpdated.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in [PotentialUpdated]:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/M1recomputefractionsafe.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in recomputefractionsafe:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/M1errorplotinfo.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorplotinfo:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/M1errorworst.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorworst:
            writer.writerow([val])
    
    #safe NworldFine
    with open("%s/M1errorbest.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorbest:
            writer.writerow([val])
    
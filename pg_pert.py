# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/MasterthesisLOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# This file is motivated by gridlod: https://github.com/TiKeil/gridlod.git

import numpy as np
from copy import deepcopy
import scipy.sparse as sparse

from gridlod import lod, util, fem, ecworker, eccontroller
import timeit

class PerturbedPetrovGalerkinLOD:
    def __init__(self, origincoef, world, k, IPatchGenerator, printLevel=0):
        self.world = world
        self.k = k
        self.IPatchGenerator = IPatchGenerator
        self.printLevel = printLevel
                
        self.epsilonList = None
        self.ageList = None
        
        self.ecList = None
        self.ecListtesting = None
        self.ecListOrigin = None
        self.Kms = None
        self.K = None
        self.basisCorrectors = None
        
        #for testing
        self.currentTestingCorrector = None
        #coefficient without defects
        self.origincoef = origincoef
        
        eccontroller.clearWorkers()
        
    def originCorrectors(self, clearFineQuantities=True):
        world = self.world
        k = self.k
        IPatchGenerator = self.IPatchGenerator
        coefficient = self.origincoef
        
        NtCoarse = np.prod(world.NWorldCoarse)

        # Reset all caches
        self.Kms = None
        self.K = None
        self.basisCorrectors = None
        
        self.ecListOrigin = [None]*NtCoarse
        
        if self.printLevel >= 2:
            print('Setting up workers for origin Correctors')
        eccontroller.setupWorker(world, coefficient, IPatchGenerator, k, clearFineQuantities, self.printLevel)
        if self.printLevel >= 2:
            print('Done')
        
        #element corrector list has coarse element size 
        ecListOrigin = self.ecListOrigin
        ecComputeList = []
                
        for TInd in range(NtCoarse):
            iElement = util.convertpLinearIndexToCoordIndex(world.NWorldCoarse-1, TInd)
                
            ecComputeList.append((TInd, iElement))    
        
        if self.printLevel >= 2:
            print('Waiting for results', len(ecComputeList))
                    
        ecResultList = eccontroller.mapComputations(ecComputeList, self.printLevel)
        for ecResult, ecCompute in zip(ecResultList, ecComputeList):
            ecListOrigin[ecCompute[0]] = ecResult
            
        self.ecList = deepcopy(ecListOrigin)
        self.ecListtesting = deepcopy(ecListOrigin)

    def CorrectorsToOrigin(self):
        self.ecListtesting = self.ecListOrigin
            
    def updateCorrectors(self, coefficient, epsilonTol, clearFineQuantities=True,
                         Testing = None, Computing= True, mc=0, new_correctors = False, runtime=False):
        assert(self.ecListOrigin is not None)

        # Note: most of the optional arguments of this function are used for purposes of the Masterthesis.
        # it became quite much and unsorted.

        start = timeit.default_timer()

        world = self.world
        k = self.k
        IPatchGenerator = self.IPatchGenerator
        
        NtCoarse = np.prod(world.NWorldCoarse)

        # Reset all caches
        self.Kms = None
        self.K = None
        self.basisCorrectors = None
                
        self.ageList = [0]*NtCoarse
        
        if self.epsilonList == None:
            self.epsilonList = [np.nan]*NtCoarse

        if Testing:
            ecListOrigin = self.ecListtesting
        else:
            ecListOrigin = self.ecListOrigin

        ecList = deepcopy(ecListOrigin)

        if new_correctors:
            ecListOrigin = self.ecList

        if self.printLevel >= 2:
            print('Setting up workers')
        eccontroller.setupWorker(world, coefficient, IPatchGenerator, k, clearFineQuantities, self.printLevel)
        if self.printLevel >= 2:
            print('Done')
        
        #For coarse coefficient
        if self.ecList is not None and hasattr(coefficient, 'rCoarse'):
            ANew = coefficient._aBase
            AOld = deepcopy(self.origincoef.aFine)
            if ANew.ndim == 1:
                delta = np.abs((AOld-ANew)/np.sqrt(AOld*ANew))
            else:
                delta = np.abs((AOld-ANew) * np.linalg.inv(np.sqrt(AOld)) )
                # This here is failing because of division by zero:
                # delta = np.abs((AOld - ANew) * np.linalg.inv(np.sqrt(AOld)) * np.linalg.inv(np.sqrt(ANew)))

        ageList = self.ageList
        
        if epsilonTol == 0:
            epsilonList = self.epsilonList
        else:
            epsilonList = deepcopy(self.epsilonList)
        
        recomputeCount = 0
        ecComputeList = []
        for TInd in range(NtCoarse):
            if self.printLevel >= 3:
                print(str(TInd) + ' / ' + str(NtCoarse), end=' ')
            
            ageList[TInd] += 1
            
            #mapper
            iElement = util.convertpLinearIndexToCoordIndex(world.NWorldCoarse-1, TInd)
            ecT = ecListOrigin[TInd]
            if Testing:
                epsilonT = epsilonList[TInd]
            else:    
                if hasattr(coefficient, 'aLagging'):
                    coefficientPatch = coefficient.localize(ecT.iPatchWorldCoarse, ecT.NPatchCoarse)
                    epsilonT = ecList[TInd].computeErrorIndicatorFineWithLagging(coefficientPatch.aFine, coefficientPatch.aLagging)
                if hasattr(coefficient, 'rCoarse'):
                    epsilonT = ecListOrigin[TInd].computeLocalCoarseErrorIndicator(delta)
                elif hasattr(ecT, 'fsi'):
                    coefficientPatch = coefficient.localize(ecT.iPatchWorldCoarse, ecT.NPatchCoarse)
                    epsilonT = ecListOrigin[TInd].computeErrorIndicatorFine(coefficientPatch)
                epsilonList[TInd] = epsilonT
            
            if self.printLevel >= 2:
                print('epsilonT = ' + str(epsilonT), end=' ') 
                
            if epsilonT > epsilonTol:
                if self.printLevel >= 2:
                    print('C')
                if Testing:
                    epsilonList[TInd] = 0
                    self.currentTestingCorrector = TInd
                ecComputeList.append((TInd, iElement))
                ecList[TInd] = None
                ageList[TInd] = 0
                recomputeCount += 1
            else:
                if self.printLevel > 1:
                    print('N')    

        time_to_compute = timeit.default_timer() - start

        if self.printLevel >= 2:
            print('Waiting for results', len(ecComputeList))
        
        if self.printLevel > 0 or Testing:
            if mc == 0:
                print("To be recomputed: ", float(recomputeCount)/NtCoarse*100, '%')
        
        self.printLevel = 0

        if Computing:
            ecResultList = eccontroller.mapComputations(ecComputeList, self.printLevel)
            for ecResult, ecCompute in zip(ecResultList, ecComputeList):
                ecList[ecCompute[0]] = ecResult
        else:
            if self.printLevel >= 1:
                print("Not Recomputed!")
                
        if Computing:
            self.ecList = ecList
        
        if epsilonTol != 0:
            self.ecListtesting = ecList
        
        if Testing:
            self.epsilonList = epsilonList    
            
        ageListinv = np.ones(np.size(ageList))
        ageListinv = ageListinv - ageList
        
        if runtime:
            return ageListinv, time_to_compute
        else:
            return ageListinv, epsilonList

    def clearCorrectors(self):
        NtCoarse = np.prod(self.world.NWorldCoarse)
        self.ecList = None
        self.coefficient = None

    def computeCorrection(self, ARhsFull=None, MRhsFull=None):
        assert(self.ecList is not None)
        assert(self.origincoef is not None)

        world = self.world
        NCoarseElement = world.NCoarseElement
        NWorldCoarse = world.NWorldCoarse
        NWorldFine = NWorldCoarse*NCoarseElement

        NpFine = np.prod(NWorldFine+1)
        
        coefficient = self.origincoef
        IPatchGenerator = self.IPatchGenerator

        localBasis = world.localBasis
        
        TpIndexMap = util.lowerLeftpIndexMap(NCoarseElement, NWorldFine)
        TpStartIndices = util.pIndexMap(NWorldCoarse-1, NWorldFine, NCoarseElement)
        
        uFine = np.zeros(NpFine)
        
        NtCoarse = np.prod(world.NWorldCoarse)
        for TInd in range(NtCoarse):
            if self.printLevel > 0:
                print(str(TInd) + ' / ' + str(NtCoarse))
                
            ecT = self.ecList[TInd]
            
            coefficientPatch = coefficient.localize(ecT.iPatchWorldCoarse, ecT.NPatchCoarse)
            IPatch = IPatchGenerator(ecT.iPatchWorldCoarse, ecT.NPatchCoarse)

            if ARhsFull is not None:
                ARhsList = [ARhsFull[TpStartIndices[TInd] + TpIndexMap]]
            else:
                ARhsList = None
                
            if MRhsFull is not None:
                MRhsList = [MRhsFull[TpStartIndices[TInd] + TpIndexMap]]
            else:
                MRhsList = None
                
            correctorT = ecT.computeElementCorrector(coefficientPatch, IPatch, ARhsList, MRhsList)[0]
            
            NPatchFine = ecT.NPatchCoarse*NCoarseElement
            iPatchWorldFine = ecT.iPatchWorldCoarse*NCoarseElement
            patchpIndexMap = util.lowerLeftpIndexMap(NPatchFine, NWorldFine)
            patchpStartIndex = util.convertpCoordinateToIndex(NWorldFine, iPatchWorldFine)

            uFine[patchpStartIndex + patchpIndexMap] += correctorT

        return uFine
    
    def assembleBasisCorrectors(self):
        if self.basisCorrectors is not None:
            return self.basisCorrectors

        assert(self.ecList is not None)

        world = self.world
        NWorldCoarse = world.NWorldCoarse
        NCoarseElement = world.NCoarseElement
        NWorldFine = NWorldCoarse*NCoarseElement
        
        NtCoarse = np.prod(NWorldCoarse)
        NpCoarse = np.prod(NWorldCoarse+1)
        NpFine = np.prod(NWorldFine+1)
        
        TpIndexMap = util.lowerLeftpIndexMap(np.ones_like(NWorldCoarse), NWorldCoarse)
        TpStartIndices = util.lowerLeftpIndexMap(NWorldCoarse-1, NWorldCoarse)
        
        cols = []
        rows = []
        data = []
        ecList = self.ecList
        for TInd in range(NtCoarse):
            ecT = ecList[TInd]
            assert(ecT is not None)
            assert(hasattr(ecT, 'fsi'))
            
            NPatchFine = ecT.NPatchCoarse*NCoarseElement
            iPatchWorldFine = ecT.iPatchWorldCoarse*NCoarseElement
            
            patchpIndexMap = util.lowerLeftpIndexMap(NPatchFine, NWorldFine)
            patchpStartIndex = util.convertpCoordIndexToLinearIndex(NWorldFine, iPatchWorldFine)
            
            colsT = TpStartIndices[TInd] + TpIndexMap
            rowsT = patchpStartIndex + patchpIndexMap
            dataT = np.hstack(ecT.fsi.correctorsList)
            
            cols.extend(np.repeat(colsT, np.size(rowsT)))
            rows.extend(np.tile(rowsT, np.size(colsT)))
            data.extend(dataT)
        
        basisCorrectors = sparse.csc_matrix((data, (rows, cols)), shape=(NpFine, NpCoarse))

        self.basisCorrectors = basisCorrectors
        return basisCorrectors
    
    def assembleBasisCorrectorsFast(self):
        ''' Is that even possible '''
        if self.basisCorrectors is not None:
            return self.basisCorrectors

        assert(self.ecList is not None)

        world = self.world
        NWorldCoarse = world.NWorldCoarse
        NCoarseElement = world.NCoarseElement
        NWorldFine = NWorldCoarse*NCoarseElement
        
        NtCoarse = np.prod(NWorldCoarse)
        NpCoarse = np.prod(NWorldCoarse+1)
        NpFine = np.prod(NWorldFine+1)
        
        TpIndexMap = util.lowerLeftpIndexMap(np.ones_like(NWorldCoarse), NWorldCoarse)
        TpStartIndices = util.lowerLeftpIndexMap(NWorldCoarse-1, NWorldCoarse)
        
        cols = []
        rows = []
        data = []
        ecList = self.ecList
        for TInd in range(NtCoarse):
            ecT = ecList[TInd]
            assert(ecT is not None)
            assert(hasattr(ecT, 'fsi'))

            NPatchFine = ecT.NPatchCoarse*NCoarseElement
            iPatchWorldFine = ecT.iPatchWorldCoarse*NCoarseElement
            
            patchpIndexMap = util.lowerLeftpIndexMap(NPatchFine, NWorldFine)
            patchpStartIndex = util.convertpCoordIndexToLinearIndex(NWorldFine, iPatchWorldFine)
            
            colsT = TpStartIndices[TInd] + TpIndexMap
            rowsT = patchpStartIndex + patchpIndexMap
            dataT = np.hstack(ecT.fsi.correctorsList)

            cols.extend(np.repeat(colsT, np.size(rowsT)))
            rows.extend(np.tile(rowsT, np.size(colsT)))
            data.extend(dataT)

        basisCorrectors = sparse.csc_matrix((data, (rows, cols)), shape=(NpFine, NpCoarse))

        self.basisCorrectors = basisCorrectors
        return basisCorrectors


        
    def assembleMsStiffnessMatrix(self):
        if self.Kms is not None:
            return self.Kms

        assert(self.ecList is not None)

        world = self.world
        NWorldCoarse = world.NWorldCoarse
        
        NtCoarse = np.prod(world.NWorldCoarse)
        NpCoarse = np.prod(world.NWorldCoarse+1)
        
        TpIndexMap = util.lowerLeftpIndexMap(np.ones_like(NWorldCoarse), NWorldCoarse)
        TpStartIndices = util.lowerLeftpIndexMap(NWorldCoarse-1, NWorldCoarse)

        cols = []
        rows = []
        data = []
        ecList = self.ecList
        for TInd in range(NtCoarse):
            ecT = ecList[TInd]
            assert(ecT is not None)

            NPatchCoarse = ecT.NPatchCoarse

            patchpIndexMap = util.lowerLeftpIndexMap(NPatchCoarse, NWorldCoarse)
            patchpStartIndex = util.convertpCoordIndexToLinearIndex(NWorldCoarse, ecT.iPatchWorldCoarse)
            
            colsT = TpStartIndices[TInd] + TpIndexMap
            rowsT = patchpStartIndex + patchpIndexMap
            dataT = ecT.csi.Kmsij.flatten()

            cols.extend(np.tile(colsT, np.size(rowsT)))
            rows.extend(np.repeat(rowsT, np.size(colsT)))
            data.extend(dataT)

        Kms = sparse.csc_matrix((data, (rows, cols)), shape=(NpCoarse, NpCoarse))

        self.Kms = Kms
        return Kms


    def assembleStiffnessMatrix(self):
        if self.K is not None:
            return self.K

        assert(self.ecList is not None)

        world = self.world
        NWorldCoarse = world.NWorldCoarse
        
        NtCoarse = np.prod(world.NWorldCoarse)
        NpCoarse = np.prod(world.NWorldCoarse+1)
        
        TpIndexMap = util.lowerLeftpIndexMap(np.ones_like(NWorldCoarse), NWorldCoarse)
        TpStartIndices = util.lowerLeftpIndexMap(NWorldCoarse-1, NWorldCoarse)

        cols = []
        rows = []
        data = []
        ecList = self.ecList
        for TInd in range(NtCoarse):
            ecT = ecList[TInd]
            assert(ecT is not None)

            NPatchCoarse = ecT.NPatchCoarse

            colsT = TpStartIndices[TInd] + TpIndexMap
            rowsT = TpStartIndices[TInd] + TpIndexMap
            dataT = ecT.csi.Kij.flatten()

            cols.extend(np.tile(colsT, np.size(rowsT)))
            rows.extend(np.repeat(rowsT, np.size(colsT)))
            data.extend(dataT)

        K = sparse.csc_matrix((data, (rows, cols)), shape=(NpCoarse, NpCoarse))

        self.K = K
        return K
    
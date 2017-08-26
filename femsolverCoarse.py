```
# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# This file is motivated by gridlod: https://github.com/TiKeil/gridlod.git
```

import numpy as np
import scipy.sparse as sparse

from gridlod.world import World
from gridlod import util, fem, linalg

#coarse solver for testing normal fem with coarse right hand side, fine stiffness but without boundary conditions    
def solveCoarse_fem(world, aFine, MbCoarse, boundaryConditions):
    NWorldCoarse = world.NWorldCoarse
    NWorldFine = world.NWorldCoarse*world.NCoarseElement
    NCoarseElement = world.NCoarseElement
    
    NpFine = np.prod(NWorldFine+1)
    NpCoarse = np.prod(NWorldCoarse+1)
    
    if MbCoarse is None:
        MbCoarse = np.zeros(NpCoarse)
        
    boundaryMap = boundaryConditions==0
    fixedCoarse = util.boundarypIndexMap(NWorldCoarse, boundaryMap=boundaryMap)
    freeCoarse  = np.setdiff1d(np.arange(NpCoarse), fixedCoarse)
    
    AFine = fem.assemblePatchMatrix(NWorldFine, world.ALocFine, aFine)
    MCoarse = fem.assemblePatchMatrix(NWorldCoarse, world.MLocCoarse)

    bCoarse = MCoarse*MbCoarse

    basis = fem.assembleProlongationMatrix(NWorldCoarse, NCoarseElement)
    ACoarse = basis.T*(AFine*basis)

    ACoarseFree = ACoarse[freeCoarse][:,freeCoarse]
    bCoarseFree = bCoarse[freeCoarse]

    uCoarseFree = linalg.linSolve(ACoarseFree, bCoarseFree)
    uCoarseFull = np.zeros(NpCoarse)
    uCoarseFull[freeCoarse] = uCoarseFree
    uCoarseFull = uCoarseFull
    
    return uCoarseFull

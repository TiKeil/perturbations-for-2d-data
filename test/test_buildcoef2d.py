import numpy as np
import matplotlib.pyplot as plt

from gridlod.world import World, Patch

from MasterthesisLOD import buildcoef2d
from gridlod_on_perturbations.visualization_tools import drawCoefficient_origin


potenz = 8
factor = 2**(potenz - 8)
fine = 2**potenz

N = 2**5
print('log H: ' ,np.abs(np.log(np.sqrt(2*(1./N**2)))))
k = 4  # goes like log H

NFine = np.array([fine,fine])
NpFine = np.prod(NFine + 1)
NWorldCoarse = np.array([N, N])

# boundary Conditions
boundaryConditions = np.array([[0, 0], [0, 0]])

NCoarseElement = NFine // NWorldCoarse
world = World(NWorldCoarse, NCoarseElement, boundaryConditions)

'''
Construct diffusion coefficient
'''

space = 3 * factor
thick = 6 * factor

bg = 0.1		#background
val = 1			#values

soilinput = np.array([[8, 6, 3],[8, 3, 6],[10, 3, 4]])
soilMatrix = buildcoef2d.soil_converter(soilinput,NFine, BoundarySpace=space)
print(soilMatrix)

CoefClass = buildcoef2d.Coefficient2d(NFine,
                        bg                  = bg,
                        val                 = val,
                        length              = thick,
                        thick               = thick,
                        space               = space,
                        probfactor          = 1,
                        right               = 1,
                        equidistant         = True,
                        BoundarySpace       = True,
                        soilMatrix          = soilMatrix)

aFine = CoefClass.BuildCoefficient().flatten()

plt.figure("Coefficient")
drawCoefficient_origin(NFine, aFine)

plt.show()
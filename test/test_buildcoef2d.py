import numpy as np
import matplotlib.pyplot as plt

from MasterthesisLOD import buildcoef2d
from visualize import drawCoefficient_origin

potenz = 8
factor = 2**(potenz - 8)
fine = 2**potenz

NFine = np.array([fine,fine])

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

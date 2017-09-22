# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import buildcoef2d

import matplotlib.pyplot as plt
from visualize import drawCoefficient, AllshapesSixdrawCoefficient

bg = 0.05   #background
val = 1     #values
        
NWorldFine = np.array([11, 11])
NWorldFine2 = np.array([38,42]) 
           
########################
# Print all shapes
#######################

CoefClass = buildcoef2d.Coefficient2d(NWorldFine, 
                                bg              = 0.05, 
                                val             = 1, 
                                length          = 7, 
                                thick           = 1, 
                                space           = 2, 
                                probfactor      = 1, 
                                right           = 1, 
                                down            = 0, 
                                diagr1          = 0, 
                                diagr2          = 0, 
                                diagl1          = 0, 
                                diagl2          = 0, 
                                BoundarySpace   = True)
                                                

CoefClass1 = buildcoef2d.Coefficient2d(NWorldFine, 
                                bg              = 0.05, 
                                val             = 1, 
                                length          = 7, 
                                thick           = 1, 
                                space           = 2, 
                                probfactor      = 1, 
                                right           = 0, 
                                down            = 1, 
                                diagr1          = 0, 
                                diagr2          = 0, 
                                diagl1          = 0, 
                                diagl2          = 0, 
                                BoundarySpace   = True)


CoefClass2 = buildcoef2d.Coefficient2d(NWorldFine, 
                                bg              = 0.05, 
                                val             = 1, 
                                length          = 6, 
                                thick           = 1, 
                                space           = 2, 
                                probfactor      = 1, 
                                right           = 0, 
                                down            = 0, 
                                diagr1          = 1, 
                                diagr2          = 0, 
                                diagl1          = 0, 
                                diagl2          = 0, 
                                BoundarySpace   = True)

CoefClass3 = buildcoef2d.Coefficient2d(NWorldFine, 
                                bg              = 0.05, 
                                val             = 1, 
                                length          = 6, 
                                thick           = 1, 
                                space           = 2, 
                                probfactor      = 1, 
                                right           = 0, 
                                down            = 0, 
                                diagr1          = 0, 
                                diagr2          = 1, 
                                diagl1          = 0, 
                                diagl2          = 0, 
                                BoundarySpace   = True)

CoefClass4 = buildcoef2d.Coefficient2d(NWorldFine, 
                                bg              = 0.05, 
                                val             = 1, 
                                length          = 6, 
                                thick           = 1, 
                                space           = 2, 
                                probfactor      = 1, 
                                right           = 0, 
                                down            = 0, 
                                diagr1          = 0, 
                                diagr2          = 0, 
                                diagl1          = 1, 
                                diagl2          = 0, 
                                BoundarySpace   = True)

CoefClass5 = buildcoef2d.Coefficient2d(NWorldFine, 
                                bg              = 0.05, 
                                val             = 1, 
                                length          = 6, 
                                thick           = 1, 
                                space           = 2, 
                                probfactor      = 1, 
                                right           = 0, 
                                down            = 0, 
                                diagr1          = 0, 
                                diagr2          = 0, 
                                diagl1          = 0, 
                                diagl2          = 1, 
                                BoundarySpace   = True)

CoefClass.BuildCoefficient()
CoefClass1.BuildCoefficient()
C = CoefClass2.BuildCoefficient()
D = CoefClass3.BuildCoefficient()
E = CoefClass4.BuildCoefficient()
F = CoefClass5.BuildCoefficient()

# modifications
A = CoefClass.SpecificVanish(Number=[0,2])
B = CoefClass1.SpecificVanish(Number=[0,2])

A1 = A.flatten()
A4 = B.flatten()
A2 = C.flatten()
A3 = D.flatten()
A5 = E.flatten()
A6 = F.flatten()

plt.figure("Shape 1")
drawCoefficient(NWorldFine, A1)
plt.figure("Shape 2")
drawCoefficient(NWorldFine, A2)
plt.figure("Shape 3")
drawCoefficient(NWorldFine, A3)
plt.figure("Shape 4")
drawCoefficient(NWorldFine, A4)
plt.figure("Shape 5")
drawCoefficient(NWorldFine, A5)
plt.figure("Shape 6")
drawCoefficient(NWorldFine, A6)

plt.figure("All shapes")
AllshapesSixdrawCoefficient(NWorldFine, A1, A2, A3, A4, A5, A6)

plt.show()
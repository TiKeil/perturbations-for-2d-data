# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import buildcoef2d

import matplotlib.pyplot as plt
from visualize import drawCoefficient

bg = 0.05   #background
val = 1     #values
        
############### new Shape ################
NWorldFine = np.array([11, 11])
CoefClassNew = buildcoef2d.Coefficient2d(NWorldFine,
                                         bg = bg, 
                                         val = val, 
                                         space = 2, 
                                         probfactor=1, 
                                         BoundarySpace = True)
                                         
newshape = np.array([[1,6,1,1,1,6,0], 
                     [4,6,1, 1,1,6,0], 
                     [1,6,1,-1,1,6,0],
                     [4,6,1,-1,0,0,0]])

CoefClassNew.NewShape(newshape)
New = CoefClassNew.BuildCoefficient()

N = New.flatten()

plt.figure("newShape")
drawCoefficient(NWorldFine, N)

CoefClassNew1 = buildcoef2d.Coefficient2d(NWorldFine,
                                         bg = bg, 
                                         val = val, 
                                         space = 2, 
                                         probfactor=1, 
                                         BoundarySpace = True)
                                         
newshape = np.array([[1,6,1,1,2,6,0,0,1], 
                     [4,6,1,1,1,6,0,0,0], 
                     [4,6,1,1,0,0,0,0,0],
                     [1,6,1,-1,0,0,0,0,0]])

CoefClassNew1.NewShape(newshape)
New = CoefClassNew1.BuildCoefficient()

N = New.flatten()

plt.figure("newShape1")
drawCoefficient(NWorldFine, N)

plt.show()

############### shape6 #############
CoefClassShape6 = buildcoef2d.Coefficient2d(NWorldFine,
                                         bg = bg, 
                                         val = val, 
                                         space = 2, 
                                         probfactor=1, 
                                         BoundarySpace = True)
                                         
shape6 = np.array([[1,2,1,1,1,0,1],
                      [1,2,1,-1,1,1,-1],
                      [1,2,1,-1,1,1,-1],
                      [1,2,1,-1,1,1,-1],
                      [1,2,1,-1,1,1,-1],
                      [1,2,1,-1,0,0,0]])

CoefClassShape6.NewShape(shape6)
Shape6 = CoefClassShape6.BuildCoefficient()

S6 = Shape6.flatten()

plt.figure("Shape 6")
drawCoefficient(NWorldFine, S6)

plt.show()
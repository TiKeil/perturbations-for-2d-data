# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np

import buildcoef2d

import matplotlib.pyplot as plt
from visualize import drawCoefficient, ExtradrawCoefficient

bg = 0.05     #background
val = 1       #values

NWorldFine = np.array([42, 42])

CoefClass = buildcoef2d.Coefficient2d(NWorldFine,
                bg                = bg,       # background
                val               = val,      # values
                length            = 2,        # length
                thick             = 2,        # thickness
                space             = 2,        # space between values
                probfactor        = 1,        # probability of an value
                right             = 1,        # shape 1
                down              = 0,        # shape 2
                diagr1            = 0,        # shape 3
                diagr2            = 0,        # shape 4
                diagl1            = 0,        # shape 5
                diagl2            = 0,        # shape 6
                LenSwitch         = None,     # various length
                thickSwitch       = None,     # various thickness
                ChannelHorizontal = None,     # horizontal Channels
                ChannelVertical   = None,     # vertical Channels
                BoundarySpace     = True      # additional space on the boundary
                )

A = CoefClass.BuildCoefficient()              # coefficient in a numpy array

numbers = [13,20,27,44,73]  #What entries will be changed

# Change in Value
B = CoefClass.SpecificValueChange( Number           = numbers,  
                                   ratio            = -0.4,
                                   randomvalue      = None,
                                   negative         = None,
                                   ShapeRestriction = True,
                                   ShapeWave        = None,
                                   probfactor       = 1,
                                   Original         = True)

# Disappearance
C = CoefClass.SpecificVanish(      Number           = numbers,
                                   PartlyVanish     = None,
                                   probfactor       = 1,
                                   Original          = True)

# Shift
D = CoefClass.SpecificMove(        Number           = numbers,
                                   steps            = 1,
                                   randomstep       = None,
                                   randomDirection  = None,
                                   Right            = 1,
                                   BottomRight      = 1,
                                   Bottom           = 1,
                                   BottomLeft       = 1,
                                   Left             = 1,
                                   TopLeft          = 1,
                                   Top              = 1,
                                   TopRight         = 1,
                                   Original         = True)
    
A = A.flatten()
B = B.flatten()
C = C.flatten()
D = D.flatten()

plt.figure("original")
drawCoefficient(NWorldFine, A)
        
plt.figure("1")
drawCoefficient(NWorldFine, B)
plt.figure("2")
drawCoefficient(NWorldFine, C)
plt.figure("3")
drawCoefficient(NWorldFine, D)

plt.figure('all')
ExtradrawCoefficient(NWorldFine, A, B, C, D)

plt.show()
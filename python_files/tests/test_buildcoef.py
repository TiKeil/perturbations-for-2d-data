# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import numpy as np
import random
import matplotlib.pyplot as plt
from visualize import drawCoefficient

import buildcoef2d

def test_1():
    bg = 0  # background
    val = 2  # values

    space = 10
    thick = 1

    NWorldFine = np.array([128, 128])

    CoefClass = buildcoef2d.Coefficient2d(NWorldFine,
                            bg                  = bg,
                            val                 = val,
                            length              = 1,
                            thick               = thick,
                            space               = space,
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
                            ChannelVertical     = True,
                            BoundarySpace       = True)

    A = CoefClass.BuildCoefficient().flatten()

    plt.figure("Coefficient_test_1")
    drawCoefficient(NWorldFine, A)

def test_2():
    bg = 0  # background
    val = 2  # values

    space = 10
    thick = 1
    fine = 64
    NWorldFine = np.array([fine, fine])

    ChoosingShapes = np.array([
        # shape , len, thick, space
        [ 1  , 1,    1,       1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [4, fine-2, 1, 5],
        [4, fine - 2, 1, 5],
        [4, fine - 2, 1, 5]
    ])

    CoefClass = buildcoef2d.Coefficient2d(NWorldFine,
                            bg                  = bg,
                            val                 = val,
                            right               = 1,
                            LenSwitch           = None,
                            thickSwitch         = None,
                            equidistant         = True,
                            ChannelHorizontal   = None,
                            ChannelVertical     = None,
                            BoundarySpace       = True,
                            probfactor=1,
                            ChoosingShapes=ChoosingShapes)

    A = CoefClass.BuildCoefficient().flatten()
    Correct = CoefClass.SpecificVanish(Number=[10,11])
    plt.figure("Coefficient_test_1")
    drawCoefficient(NWorldFine, Correct.flatten())

    Number = [10,11]
    random.seed(32)
    lis = np.zeros(80)
    lis[0] = 1
    for i in range(np.shape(CoefClass.ShapeRemember)[0]):
        Number.append(i * random.sample(list(lis),1)[0])
    plt.figure("Perturbation")
    Perturbed = CoefClass.SpecificVanish(Number=Number, Original=True).flatten()
    drawCoefficient(NWorldFine, Perturbed)

#test_1()
test_2()
plt.show()
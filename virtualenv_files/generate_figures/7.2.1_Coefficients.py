```
# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
```

import numpy as np

import matplotlib.pyplot as plt
from visualize import drawCoefficient

import buildcoef2d

bg = 0.05 		#background
val = 1			#values

NWorldFine = np.array([256, 256])

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

#Build Coefficients
A1 = CoefClass1.BuildCoefficient()
A2 = CoefClass2.BuildCoefficient()
A3 = CoefClass3.BuildCoefficient()
A4 = CoefClass4.BuildCoefficient()
        
A1 = A1.flatten()
A2 = A2.flatten()
A3 = A3.flatten()
A4 = A4.flatten()

plt.figure("Coefficient 1")
drawCoefficient(NWorldFine, A1, greys=True)
plt.figure("Coefficient 2")
drawCoefficient(NWorldFine, A2, greys=True)
plt.figure("Coefficient 3")
drawCoefficient(NWorldFine, A3, greys=True)
plt.figure("Coefficient 4")
drawCoefficient(NWorldFine, A4, greys=True)

plt.show()
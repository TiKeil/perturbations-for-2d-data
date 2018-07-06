# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import os
import sys
import numpy as np

import matplotlib.pyplot as plt
from visualize import *
from data import *
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MultipleLocator
import csv

NWorldFine = np.array([256, 256])
NpFine = np.prod(NWorldFine+1)

#coarse World
NWorldCoarse = np.array([16,16])
NpCoarse = np.prod(NWorldCoarse+1)

#ratio between Fine and Coarse
NCoarseElement = NWorldFine/NWorldCoarse

ROOT = '../../../test_data/Coef1'

OriginalCoeff = []
f = open("%s/OriginalCoeff.txt" % ROOT, 'rb')
reader = csv.reader(f)
for row in reader:
    OriginalCoeff.append(float(row[0]))
f.close()

Abase = np.array([])
Abase = np.append(Abase,OriginalCoeff)

finescale = []
f = open("%s/finescale.txt" % ROOT, 'rb')
reader = csv.reader(f)
for row in reader:
    finescale.append(float(row[0]))
f.close()

ufine = np.array([])
ufine = np.append(ufine,finescale)

plt.figure("Original")
drawCoefficientwt(NWorldFine, Abase,greys=True)
#plt.savefig('pic/Original.pdf')

plt.figure("Solution")
drawCoefficientwt(NWorldFine+1,ufine)
#plt.savefig('pic/Solution.pdf')

runtime = []
f = open("%s/runtime.txt" % ROOT, 'rb')
reader = csv.reader(f)
for row in reader:
    runtime.append(float(row[0]))
f.close()


# Disappearance
VBase, Veps, Verrorbest, Verrorplotinfo, Verrorworst, Vvis, Vrecomputefractionsafe = RegainVanish(ROOT)

updated = np.size(Vrecomputefractionsafe)-1
Vrecomputefractionsafe = [i * updated/256. for i in Vrecomputefractionsafe]
Vrecomputefractionsafe.append(100)
Verrorplotinfo.append(Verrorplotinfo[updated])
Verrorbest.append(Verrorbest[updated])


# Perturbation
Er = np.array([])
Er = np.append(Er,VBase)
plt.figure("Coefficient for Disappearance")
drawCoefficientwt(NWorldFine, Er,greys=True)
#plt.savefig('pic/V.pdf')

# Error indicator plot
plot_error_indicator(Veps,Vrecomputefractionsafe, NWorldCoarse, 'Disappearance')

# PGLOD error plot
plot_VCLOD_error(Verrorbest, Verrorworst, Verrorplotinfo, Vrecomputefractionsafe,'Disappearance')
runtime.append(runtime[256])
plt.figure('Runtime')
plt.plot(Vrecomputefractionsafe,runtime)
plt.title("Runtime of updating process")
plt.ylabel('runtime in sec', fontsize=14)
plt.xlabel('Updated correctors in %', fontsize=14)

plt.show()
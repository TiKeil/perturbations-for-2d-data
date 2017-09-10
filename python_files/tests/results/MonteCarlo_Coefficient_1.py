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

ROOT = '../../../test_data/MonteCarlo/Coef1/p20'
mum, a, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

Abase = np.array([])
Abase = np.append(Abase,ABase)

NWorldFine = np.array([256, 256])

plt.figure("Original")
drawCoefficient(NWorldFine, Abase,greys=True)
plt.title("Original coefficient")

plt.show()

plt.figure('confidence u_LoD')
plt.title('confidence u_LoD')
plt.semilogy(ems,plottingx,'--r', linewidth = 0.6)
plt.semilogy(ems,plottingy,'-', linewidth = 0.4)
plt.semilogy(ems,plottingz,'--r', linewidth = 0.6)

plt.figure('confidence u_VcLoD')
plt.title('confidence u_VcLoD')
plt.semilogy(ems,plotting2x,'--b', linewidth = 0.6)
plt.semilogy(ems,plotting2y,'-g', linewidth = 0.4)
plt.semilogy(ems,plotting2z,'--b', linewidth = 0.6)

plt.figure('comparison')
plt.title('comparison')
plt.semilogy(ems,plottingx,'--r', linewidth = 0.6)
plt.semilogy(ems,plottingy,'-', label='uLod', linewidth = 0.4)
plt.semilogy(ems,plottingz,'--r', linewidth = 0.6)

plt.semilogy(ems,plotting2x,'--b', linewidth = 0.6)
plt.semilogy(ems,plotting2y,'-g', label='VcLod', linewidth = 0.4)
plt.semilogy(ems,plotting2z,'--b', linewidth = 0.6)
ymin, ymax = plt.ylim()
plt.yticks((ymin,ymax),fontsize="small")

plt.legend()

print 'uLod  ' + str(mum[0]) + ' +- ' + str(a[0])
print 'uVcLod ' + str(mum[1]) + ' +- ' + str(a[1])
print 'uLodVcLod ' + str(mum[2]) + ' +- ' + str(a[2])

plt.show()

ROOT = '../../../test_data/MonteCarlo/Coef1/p40'
mum, a, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

Abase = np.array([])
Abase = np.append(Abase,ABase)

NWorldFine = np.array([256, 256])

plt.figure("Original")
drawCoefficient(NWorldFine, Abase,greys=True)
plt.title("Original coefficient")

plt.show()

plt.figure('confidence u_LoD1')
plt.title('confidence u_LoD')
plt.semilogy(ems,plottingx,'--r', linewidth = 0.6)
plt.semilogy(ems,plottingy,'-', linewidth = 0.4)
plt.semilogy(ems,plottingz,'--r', linewidth = 0.6)

plt.figure('confidence u_VcLoD1')
plt.title('confidence u_VcLoD')
plt.semilogy(ems,plotting2x,'--b', linewidth = 0.6)
plt.semilogy(ems,plotting2y,'-g', linewidth = 0.4)
plt.semilogy(ems,plotting2z,'--b', linewidth = 0.6)

plt.figure('comparison1')
plt.title('comparison')
plt.semilogy(ems,plottingx,'--r', linewidth = 0.6)
plt.semilogy(ems,plottingy,'-', label='uLod', linewidth = 0.4)
plt.semilogy(ems,plottingz,'--r', linewidth = 0.6)

plt.semilogy(ems,plotting2x,'--b', linewidth = 0.6)
plt.semilogy(ems,plotting2y,'-g', label='VcLod', linewidth = 0.4)
plt.semilogy(ems,plotting2z,'--b', linewidth = 0.6)
ymin, ymax = plt.ylim()
plt.yticks((ymin,ymax),fontsize="small")

plt.legend()

print 'uLod  ' + str(mum[0]) + ' +- ' + str(a[0])
print 'uVcLod ' + str(mum[1]) + ' +- ' + str(a[1])
print 'uLodVcLod ' + str(mum[2]) + ' +- ' + str(a[2])

plt.show()
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

ploti = []
ploti.append(0)
ploti.append(5)
ploti.append(10)
ploti.append(20)
ploti.append(30)
ploti.append(100)
ems = np.arange(99)
printer = []
printer2 = []
printer3 = []
printer4 = []

ROOT = '../../../test_data/MonteCarlo/Coef1/p1/0'
mum0, a0, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y0, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef1/p1/5'
mum5, a5, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y5, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef1/p1/10'
mum10, a10, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y10, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef1/p1/20'
mum20, a20, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y20, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef1/p1/30'
mum30, a30, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y30, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

mumi = []
mumi.append(mum0[1]/mum0[0])
mumi.append(mum5[1]/mum0[0])
mumi.append(mum10[1]/mum0[0])
mumi.append(mum20[1]/mum0[0])
mumi.append(mum30[1]/mum0[0])
mumi.append(mum0[0]/mum0[0])
mumi1 = mumi
printer.append(mum0[1])
printer.append(mum5[1])
printer.append(mum10[1])
printer.append(mum20[1])
printer.append(mum30[1])
printer.append(mum0[0])
printer2.append(mum0[2])
printer2.append(mum5[2])
printer2.append(mum10[2])
printer2.append(mum20[2])
printer2.append(mum30[2])
printer2.append(0)
printer3.append(a0[1])
printer4.append(a0[2])
printer3.append(a5[1])
printer4.append(a5[2])
printer3.append(a10[1])
printer4.append(a10[2])
printer3.append(a20[1])
printer4.append(a20[2])
printer3.append(a30[1])
printer4.append(a30[2])
printer3.append(a0[0])
printer4.append(0)

fig = plt.figure('comparison Coefficient 1')
ax = fig.add_subplot(111)
ax.set_title('Comparison Coefficient 1')
ax.semilogy(ems,plotting2y0,'-', label=str(0) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y5,'-', label=str(5) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y10,'-', label=str(10) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y20,'-', label=str(20) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y30,'-', label=str(30) + '%', linewidth = 0.4)
ax.semilogy(ems,plottingy,'-', label='100%', linewidth = 0.4)
ymin, ymax = ax.set_ylim()
ax.set_yticks((ymin,ymax))
ax.set_ylabel('$e_{L^2}(u_{vc},u_h)$',fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
ax.set_xlabel('$M$',fontsize=17)
ax.legend(fontsize=13) #Legende
fig.subplots_adjust(left=0.20,bottom=0.12,right=0.99,top=0.95,wspace=0,hspace=0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.e'))

ROOT = '../../../test_data/MonteCarlo/Coef2/p1/0'
mum0, a0, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y0, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef2/p1/5'
mum5, a5, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y5, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef2/p1/10'
mum10, a10, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y10, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef2/p1/20'
mum20, a20, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y20, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef2/p1/30'
mum30, a30, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y30, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

mumi = []
mumi.append(mum0[1]/mum0[0])
mumi.append(mum5[1]/mum0[0])
mumi.append(mum10[1]/mum0[0])
mumi.append(mum20[1]/mum0[0])
mumi.append(mum30[1]/mum0[0])
mumi.append(mum0[0]/mum0[0])
mumi2 = mumi
printer.append(mum0[1])
printer.append(mum5[1])
printer.append(mum10[1])
printer.append(mum20[1])
printer.append(mum30[1])
printer.append(mum0[0])
printer2.append(mum0[2])
printer2.append(mum5[2])
printer2.append(mum10[2])
printer2.append(mum20[2])
printer2.append(mum30[2])
printer2.append(0)
printer3.append(a0[1])
printer4.append(a0[2])
printer3.append(a5[1])
printer4.append(a5[2])
printer3.append(a10[1])
printer4.append(a10[2])
printer3.append(a20[1])
printer4.append(a20[2])
printer3.append(a30[1])
printer4.append(a30[2])
printer3.append(a0[0])
printer4.append(0)

fig = plt.figure('comparison Coefficient 2')
ax = fig.add_subplot(111)
ax.set_title('Comparison Coefficient 2')
ax.semilogy(ems,plotting2y0,'-', label=str(0) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y5,'-', label=str(5) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y10,'-', label=str(10) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y20,'-', label=str(20) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y30,'-', label=str(30) + '%', linewidth = 0.4)
ax.semilogy(ems,plottingy,'-', label='100%', linewidth = 0.4)
ymin, ymax = ax.set_ylim()
ax.set_yticks((ymin,ymax))
ax.set_ylabel('$e_{L^2}(u_{vc},u_h)$',fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
ax.set_xlabel('$M$',fontsize=17)
ax.legend(fontsize=13) #Legende
fig.subplots_adjust(left=0.20,bottom=0.12,right=0.99,top=0.95,wspace=0,hspace=0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.e'))

ROOT = '../../../test_data/MonteCarlo/Coef3/p1/0'
mum0, a0, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y0, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef3/p1/5'
mum5, a5, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y5, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef3/p1/10'
mum10, a10, ABase, ems, plottingx10, plottingy10, plottingz10, plotting2x10, plotting2y10, plotting2z10, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef3/p1/20'
mum20, a20, ABase, ems, plottingx20, plottingy20, plottingz20, plotting2x20, plotting2y20, plotting2z20, plotting3x20, plotting3y20, plotting3z20 = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef3/p1/30'
mum30, a30, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y30, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

fig = plt.figure('Comparison for Coefficient 3')
ax = fig.add_subplot(111)
ax.set_title('Comparison')
ax.semilogy(ems,plottingx,':r', linewidth = 0.6)
ax.semilogy(ems,plottingy,'-r', label='PGLOD', linewidth = 0.4)
ax.semilogy(ems,plottingz,':r', linewidth = 0.6)

ax.semilogy(ems,plotting2x20,':b', linewidth = 0.6)
ax.semilogy(ems,plotting2y20,'-b', label='VcLod 20%', linewidth = 0.4)
ax.semilogy(ems,plotting2z20,':b', linewidth = 0.6)
ymin, ymax = ax.set_ylim()
ax.set_yticks((ymin,ymax))
ax.set_ylabel('$e_{L^2}(u_{vc},u_h)$',fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
ax.set_xlabel('$M$',fontsize=17)
ax.legend(fontsize=13) #Legende
fig.subplots_adjust(left=0.20,bottom=0.12,right=0.99,top=0.95,wspace=0,hspace=0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.e'))

fig = plt.figure('Comparison2 for Coefficient 3')
ax = fig.add_subplot(111)
ax.set_title('Comparison for PGLOD and VcLod at 20%')
ax.semilogy(ems,plotting3x20,'--b', linewidth = 0.6)
ax.semilogy(ems,plotting3y20,'-g', label='PGLOD-VcLod 20%', linewidth = 0.4)
ax.semilogy(ems,plotting3z20,'--b', linewidth = 0.6)
ymin, ymax = ax.set_ylim()
ax.set_yticks((ymin,ymax))
ax.set_ylabel('$e_{L^2}(u_{vc},u_k)$',fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
ax.set_xlabel('$M$',fontsize=17)
ax.legend(fontsize=13) #Legende
fig.subplots_adjust(left=0.20,bottom=0.12,right=0.99,top=0.95,wspace=0,hspace=0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.e'))

mumi = []
mumi.append(mum0[1]/mum0[0])
mumi.append(mum5[1]/mum0[0])
mumi.append(mum10[1]/mum0[0])
mumi.append(mum20[1]/mum0[0])
mumi.append(mum30[1]/mum0[0])
mumi.append(mum0[0]/mum0[0])
mumi3 = mumi
printer.append(mum0[1])
printer.append(mum5[1])
printer.append(mum10[1])
printer.append(mum20[1])
printer.append(mum30[1])
printer.append(mum0[0])
printer2.append(mum0[2])
printer2.append(mum5[2])
printer2.append(mum10[2])
printer2.append(mum20[2])
printer2.append(mum30[2])
printer2.append(0)
printer3.append(a0[1])
printer4.append(a0[2])
printer3.append(a5[1])
printer4.append(a5[2])
printer3.append(a10[1])
printer4.append(a10[2])
printer3.append(a20[1])
printer4.append(a20[2])
printer3.append(a30[1])
printer4.append(a30[2])
printer3.append(a0[0])
printer4.append(0)

fig = plt.figure('comparison Coefficient 3')
ax = fig.add_subplot(111)
ax.set_title('Comparison Coefficient 3')
ax.semilogy(ems,plotting2y0,'-', label=str(0) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y5,'-', label=str(5) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y10,'-', label=str(10) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y20,'-', label=str(20) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y30,'-', label=str(30) + '%', linewidth = 0.4)
ax.semilogy(ems,plottingy,'-', label='100%', linewidth = 0.4)
ymin, ymax = ax.set_ylim()
ax.set_yticks((ymin,ymax))
ax.set_ylabel('$e_{L^2}(u_{vc},u_h)$',fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
ax.set_xlabel('$M$',fontsize=17)
ax.legend(fontsize=13) #Legende
fig.subplots_adjust(left=0.20,bottom=0.12,right=0.99,top=0.95,wspace=0,hspace=0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.e'))

ROOT = '../../../test_data/MonteCarlo/Coef4/p1/0'
mum0, a0, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y0, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef4/p1/5'
mum5, a5, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y5, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef4/p1/10'
mum10, a10, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y10, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef4/p1/20'
mum20, a20, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y20, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)

ROOT = '../../../test_data/MonteCarlo/Coef4/p1/30'
mum30, a30, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y30, plotting2z, plotting3x, plotting3y, plotting3z = regainer(ROOT)


mumi = []
mumi.append(mum0[1]/mum0[0])
mumi.append(mum5[1]/mum0[0])
mumi.append(mum10[1]/mum0[0])
mumi.append(mum20[1]/mum0[0])
mumi.append(mum30[1]/mum0[0])
mumi.append(mum0[0]/mum0[0])
mumi4 = mumi
printer.append(mum0[1])
printer.append(mum5[1])
printer.append(mum10[1])
printer.append(mum20[1])
printer.append(mum30[1])
printer.append(mum0[0])
printer2.append(mum0[2])
printer2.append(mum5[2])
printer2.append(mum10[2])
printer2.append(mum20[2])
printer2.append(mum30[2])
printer2.append(0)
printer3.append(a0[1])
printer4.append(a0[2])
printer3.append(a5[1])
printer4.append(a5[2])
printer3.append(a10[1])
printer4.append(a10[2])
printer3.append(a20[1])
printer4.append(a20[2])
printer3.append(a30[1])
printer4.append(a30[2])
printer3.append(a0[0])
printer4.append(0)

fig = plt.figure('comparison Coefficient 4')
ax = fig.add_subplot(111)
ax.set_title('Comparison Coefficient 4')
ax.semilogy(ems,plotting2y0,'-', label=str(0) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y5,'-', label=str(5) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y10,'-', label=str(10) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y20,'-', label=str(20) + '%', linewidth = 0.4)
ax.semilogy(ems,plotting2y30,'-', label=str(30) + '%', linewidth = 0.4)
ax.semilogy(ems,plottingy,'-', label='100%', linewidth = 0.4)
ymin, ymax = ax.set_ylim()
ax.set_yticks((ymin,ymax))
ax.set_ylabel('$e_{L^2}(u_{vc},u_h)$',fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)
ax.set_xlabel('$M$',fontsize=17)
ax.legend(fontsize=13) #Legende
fig.subplots_adjust(left=0.20,bottom=0.12,right=0.99,top=0.95,wspace=0,hspace=0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.e'))

plt.figure('Expectation value comparison')
plt.plot(ploti,mumi1, label='Coefficient 1')
plt.plot(ploti,mumi2, label='Coefficient 2')
plt.plot(ploti,mumi3, label='Coefficient 3')
plt.plot(ploti,mumi4, label='Coefficient 4')
plt.legend()
plt.show()

for i in range(0,np.size(printer)):
    print str(printer[i]) + ' +- ' + str(printer3[i]) 
    print str(printer2[i]) + ' +- ' + str(printer4[i])
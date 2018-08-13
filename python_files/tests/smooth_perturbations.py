# This file is part of the paper for "Localization of multiscale problems with random defects"
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil


import numpy as np
import matplotlib.pyplot as plt

from gridlod import util, world, fem, femsolver
from gridlod.world import World
import femsolverCoarse, buildcoef2d

def index_search(x, xpcoarse):
    # This function gets the correct index of the non transformed solution
    index = 0
    for k in range(0,np.shape(xpcoarse)[0]):
        if x > xpcoarse[k]:
            index += 1
    return index

# this fine resolution should be enough
fine = 4096
NFine = np.array([fine])
NpFine = np.prod(NFine + 1)
# list of coarse meshes
NList = [2, 4, 8, 16 ,32, 64, 128, 256]

# construction of the coefficients. 
# The plot is ecactly what we want and the perturbation is
#
# phi^-1(x) = x + (1-x) * x     (alpha is 1)
# phi(y) = 1- sqrt(1-y) 
#

# aFine is the reference coefficient
aFine = np.ones(fine)
aFine /= 10
aPert = np.copy(aFine)

for i in range(int(fine/2.)-1,int(fine*3/4.)-1):
    aFine[i] = 1
    
# aPert is the perturbed coefficient
for i in range(int(fine*3/4.)-1,int(fine*0.9375)-1):
    aPert[i] = 1

# this is the mesh size
delta_h = 1./fine

# jAj is the perturbed reference coefficient
jAj = np.copy(aFine)

for i in range(fine):
    point = i * delta_h + delta_h / 2
    transformed_point = 2*point - point**2
    # print('x point {} goes to y point {}'.format(point, transformed_point))
    detJ = 0.5 * 1/(np.sqrt(1-transformed_point))
    jAj[i] *=  detJ     # what value is supposed to be here ???

# interior nodes for plotting
xt = util.tCoordinates(NFine).flatten()

# This is the right hand side and the one occuring from the transformation
f = np.ones(fine+1)
f_pert = np.copy(f)
for i in range(0, fine):
    point = i * delta_h + delta_h / 2
    transformed_point = 2 * point - point ** 2
    detJ = 0.5 * 1. / (np.sqrt(1 - transformed_point))
    f_pert[i] *= 1/detJ

plt.figure('right hand side')
plt.plot(np.arange(fine+1), f, label='NT')
plt.plot(np.arange(fine+1), f_pert, label='TR')
plt.legend()

# plot coefficients and compare them
plt.figure('Coefficient_pert')
plt.plot(xt, aFine, label='$A$')
plt.grid(True)
plt.yticks((0, np.max(aFine) + np.min(aFine)), fontsize=16)
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],fontsize=16)
plt.ylabel('$y$', fontsize=16)
plt.xlabel('$x$', fontsize=16)
plt.legend(frameon=False, fontsize=16)

plt.plot(xt, aPert, label='$Apert$')
plt.yticks((0, np.max(aPert) + np.min(aPert)), fontsize=16)
plt.xticks(fontsize=16)
plt.ylabel('$y$', fontsize=16)
plt.xlabel('$x$', fontsize=16)
plt.legend(frameon=False, fontsize=16)

jAj_plotter = np.copy(jAj)
for i in range(fine):
    if jAj_plotter[i] > 5:
        jAj_plotter[i] = 5
plt.plot(xt, jAj_plotter, label='$JAJ$')
plt.yticks((0, np.max(jAj_plotter) + np.min(jAj_plotter)), fontsize=16)
plt.xticks(fontsize=16)
plt.ylabel('$y$', fontsize=16)
plt.xlabel('$x$', fontsize=16)
plt.legend(frameon=False, fontsize=16)

# Now we first compute the reference soulution (the perturbed solution with the perturbed domain)
exact_problem = []
newErrorFine = []
x = []
y = []
for N in NList:
    NWorldCoarse = np.array([N])
    boundaryConditions = np.array([[0, 0]])

    NCoarseElement = NFine / NWorldCoarse
    world = World(NWorldCoarse, NCoarseElement, boundaryConditions)
    AFine = fem.assemblePatchMatrix(NFine, world.ALocFine, aPert)

    # grid nodes
    xpCoarse = util.pCoordinates(NFine).flatten()
    NpCoarse = np.prod(NWorldCoarse + 1)
    uCoarseFull, nothing = femsolver.solveCoarse(world, aPert, f, None, boundaryConditions)

    basis = fem.assembleProlongationMatrix(NWorldCoarse, NCoarseElement)
    uCoarseFull = basis * uCoarseFull

    exact_problem.append(uCoarseFull)
    x.append(N)
    y.append(1. / N)
    if np.size(x) == 1:
        plt.figure('FEM-Solutions perturbed')
        plt.subplots_adjust(left=0.01, bottom=0.04, right=0.99, top=0.95, wspace=0, hspace=0.2)
        ax_1 = plt.subplot(241)
    elif np.size(x) == 2:
        ax_2 = plt.subplot(242)
    elif np.size(x) == 3:
        ax_3 = plt.subplot(243)
    elif np.size(x) == 4:
        ax_4 = plt.subplot(244)
    elif np.size(x) == 5:
        ax_5 = plt.subplot(245)
    elif np.size(x) == 6:
        ax_6 = plt.subplot(246)
    elif np.size(x) == 7:
        ax_7 = plt.subplot(247)
    elif np.size(x) == 8:
        ax_8 = plt.subplot(248)

    # EX means exact
    plt.plot(xpCoarse, uCoarseFull, '--', label='EX')
    plt.title('1/h= ' + str(N), fontsize="small")
    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False,
                        labelleft=False)
    plt.legend(frameon=False, fontsize="small")

# Now we compute the perturbed reference problem and transform it
transformed_problem = []
x = []
y = []
for N in NList:
    NWorldCoarse = np.array([N])
    boundaryConditions = np.array([[0, 0]])

    NCoarseElement = NFine / NWorldCoarse
    world = World(NWorldCoarse, NCoarseElement, boundaryConditions)
    AFine = fem.assemblePatchMatrix(NFine, world.ALocFine, jAj)

    # grid nodes
    xpCoarse = util.pCoordinates(NFine).flatten()
    NpCoarse = np.prod(NWorldCoarse + 1)

    uCoarseFull, nothing = femsolver.solveCoarse(world, jAj, f_pert, None, boundaryConditions)

    basis = fem.assembleProlongationMatrix(NWorldCoarse, NCoarseElement)
    uCoarseFull = basis * uCoarseFull
    uCoarseFull_transformed = np.copy(uCoarseFull)
    k = 0

    for k in range(0,np.shape(xpCoarse)[0]):
        transformed_x = 1 - np.sqrt(1-xpCoarse[k])   # this should be the correcto direction
        # print('point {} has been transformed to {} and the new index is {}'.format(xpCoarse[k],transformed_x,index_search(transformed_x, xpCoarse)))
        uCoarseFull_transformed[k] = uCoarseFull[index_search(transformed_x, xpCoarse)]

    transformed_problem.append(uCoarseFull_transformed)
    x.append(N)
    y.append(1. / N)
    if np.size(x) == 1:
        plt.figure('FEM-Solutions perturbed')
        plt.subplots_adjust(left=0.01, bottom=0.04, right=0.99, top=0.95, wspace=0, hspace=0.2)
        ax = ax_1
    elif np.size(x) == 2:
        ax = ax_2
    elif np.size(x) == 3:
        ax = ax_3
    elif np.size(x) == 4:
        ax = ax_4
    elif np.size(x) == 5:
        ax = ax_5
    elif np.size(x) == 6:
        ax = ax_6
    elif np.size(x) == 7:
        ax = ax_7
    elif np.size(x) == 8:
        ax = ax_8

    # T stands for transformed and NT for non transformed
    ax.plot(xpCoarse, uCoarseFull_transformed, '--', label='T')
    ax.plot(xpCoarse, uCoarseFull, '--', label='NT')
    ax.set_title('1/h= ' + str(N), fontsize="small")
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False,
                    labelleft=False)
    ax.legend(frameon=False, fontsize="small")

# here, we compare the solutions.
i = 0
for N in NList:
    NWorldCoarse = np.array([N])
    NCoarseElement = NFine / NWorldCoarse
    # grid nodes
    xpCoarse = util.pCoordinates(NFine).flatten()
    NpCoarse = np.prod(NWorldCoarse + 1)
    plt.figure('error')
    error = exact_problem[i]-transformed_problem[i]
    # plt.plot(xpCoarse, error, label=str(N))
    plt.legend(frameon=False, fontsize="small")
    i += 1

# print index_search(0.5,xpCoarse)

plt.show()

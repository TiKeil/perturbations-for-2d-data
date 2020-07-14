# ~~~
# This file is part of the project: perturbations-for-2d-data
#   https://github.com/TiKeil/perturbations-for-2d-data.git
# Copyright 2017-2020: Tim Keil. All rights reserved.
# License: Licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Tim Keil (2017 - 2020)
#
# This file is part of the module "Variational crimes in the Localized orthogonal decomposition method":
# ~~~

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MultipleLocator
from matplotlib import cm

def d3sol(N, s, String='FinescaleSolution'):
    '''
    3d solution
    '''
    fig = plt.figure(String)
    ax = fig.add_subplot(111, projection='3d')

    from gridlod import util
    xp = util.pCoordinates(N)
    X = xp[0:,1:].flatten()
    Y = xp[0:,:1].flatten()
    X = np.unique(X)
    Y = np.unique(Y)

    X, Y = np.meshgrid(X, Y)

    uLodFine = s.reshape(N+1)

    # Plot the surface.
    ymin, ymax = ax.set_zlim([0,0.021])
    ax.set_zticks((ymin, ymax))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.set_zlabel('$z$', size=16)
    ax.set_xlabel('$x$', size=16)
    ax.set_ylabel('$y$', size=16)
    # ax.view_init(50, 35)
    ax.zaxis.set_major_locator(LinearLocator(4))
    ax.tick_params(labelsize=12)
    #ax.axis('off')

    # Add a color bar which maps values to colors.
    surf = ax.plot_surface(X, Y, uLodFine, cmap=cm.hot, vmin=0, vmax=0.021)
    #fig.colorbar(surf, shrink=0.5)


def drawCoefficient_origin(N, a, transformed = False, lim=None, logNorm=True, colorbar=True):
    # This is drawCoefficient from test_pgtransport.py in gridlod
    if a.ndim == 3:
        a = np.linalg.norm(a, axis=(1, 2), ord=2)

    aCube = a.reshape(N, order='F')
    aCube = np.ascontiguousarray(aCube.T)

    plt.clf()

    cmap = plt.cm.hot_r
    if lim is not None:
        plt.imshow(aCube,
                   origin='lower_left',
                   interpolation='none', cmap=cmap,
                   norm=matplotlib.colors.LogNorm(), vmax=lim[0], vmin=lim[1])
    else:
        if logNorm:
            plt.imshow(aCube,
                       origin='lower_left',
                       interpolation='none', cmap=cmap,
                       norm=matplotlib.colors.LogNorm())
        else:
            plt.imshow(aCube,
                       origin='lower_left',
                       interpolation='none', cmap=cmap)

    plt.xticks([])
    plt.yticks([])
    if colorbar:
        plt.colorbar()
    if transformed:
        cb = plt.colorbar()
        font_size = 14  # Adjust as appropriate.
        cb.ax.tick_params(labelsize=font_size)
        cb.set_ticks([0.1,1,10])

def drawCoefficient(N, a, greys=False, normalize=None, cmap_default = False, colorbar=False):
    '''
    visualizing the 2d coefficient
    '''
    aCube = a.reshape(N)
    plt.clf()

    if greys:
        cmap = 'Greys'
    elif cmap_default:
        cmap = None
    else:
        cmap = cm.hot_r

    if normalize is None:
        plt.imshow(aCube,
                   origin='lower_left',
                   interpolation='none',
                   cmap=cmap)

    else:
        plt.imshow(aCube,
                   origin='lower_left',
                   interpolation='none',
                   cmap=cmap,
                   vmin=normalize[0],
                   vmax=normalize[1])

    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False,
                    labelleft=False)
    plt.subplots_adjust(left=0.00, bottom=0.02, right=1, top=0.95, wspace=0.2, hspace=0.2)
    if colorbar:
        plt.colorbar()

def drawCoefficientwt(N, a, greys=False):
    '''
    visualizing the 2d coefficient without a title
    '''
    aCube = a.reshape(N)
    plt.clf()

    if greys:
        cmap = 'Greys'
    else:
        cmap = cm.hot_r

    plt.imshow(aCube,
               origin='upper',
               interpolation='none',
               cmap=cmap)
    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    plt.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.98,wspace=0.2,hspace=0.2)


def ExtradrawCoefficient(N, a, b, c ,d):
    aCube = a.reshape(N)
    bCube = b.reshape(N)
    cCube = c.reshape(N)
    dCube = d.reshape(N)

    plt.clf()

    cmap = plt.cm.viridis

    plt.subplot(221)
    plt.title("Original", fontsize=10)
    plt.imshow(aCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    plt.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.95,wspace=0.00,hspace=0.15)

    plt.subplot(222)
    plt.title("Change in value", fontsize=10)
    plt.imshow(bCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    plt.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.95,wspace=0.00,hspace=0.15)

    plt.subplot(223)
    plt.title("Disappearance", fontsize=10)
    plt.imshow(cCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    plt.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.95,wspace=0.00,hspace=0.15)
    plt.subplot(224)
    plt.title("Shift", fontsize=10)
    plt.imshow(dCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    plt.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.95,wspace=0.00,hspace=0.15)

def drawPatches(N, a, fig, ax, te):
    aCube = a.reshape(N)

    if te == 10:
        ax.plot(5,5,'*w')

    major_ticks = np.arange(0, te, 1)
    ax.imshow(aCube, cmap=cm.Blues, extent=[0, te, 0, te])

    ax.axis([0, te, 0, te])
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    ax.grid(which='both')
    ax.grid(which='major', linestyle="-", color="black")

def AllshapesSixdrawCoefficient(N, a, b, c , d, e, f):
    aCube = a.reshape(N)
    bCube = b.reshape(N)
    cCube = c.reshape(N)
    dCube = d.reshape(N)
    eCube = e.reshape(N)
    fCube = f.reshape(N)

    plt.clf()

    plt.subplot(161)
    plt.title("Shape 1.", fontsize=10)
    plt.imshow(aCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.axis('off')
    plt.subplot(162)
    plt.title("Shape 2.", fontsize=10)
    plt.imshow(bCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.axis('off')
    plt.subplot(163)
    plt.title("Shape 3.", fontsize=10)
    plt.imshow(cCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.axis('off')
    plt.subplot(164)
    plt.title("Shape 4.", fontsize=10)
    plt.imshow(dCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.axis('off')
    plt.subplot(165)
    plt.title("Shape 5.", fontsize=10)
    plt.imshow(eCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.axis('off')
    plt.subplot(166)
    plt.title("Shape 6.", fontsize=10)
    plt.imshow(fCube,
               origin='upper',
               interpolation='none',
               cmap=cm.hot_r)
    plt.axis('off')

def drawCoefficientGrid(N, a, fig, ax, Greys=False, original_style = False, logplot=False, colorbar=True, Gridsize = 4):
    #now we only use log plots
    aCube = a.reshape(N, order ='F')
    aCube = np.ascontiguousarray(aCube.T)

    te = Gridsize
    major_ticks = np.arange(0, te, 1)
    if Greys:
        im  = ax.imshow(aCube, cmap='Greys', extent=[0, te, 0, te])
    else:
        if original_style:
            im = ax.imshow(aCube, cmap=cm.hot_r, origin='lower', extent=[0, te, 0, te], norm=matplotlib.colors.LogNorm())
        else:
            im = ax.imshow(aCube, cmap=cm.hot_r, extent=[0, te, 0, te], norm=matplotlib.colors.LogNorm())
    ax.axis([0, te, 0, te])
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    fig.subplots_adjust(left=0.00, bottom=0.02, right=1, top=0.95, wspace=0.2, hspace=0.2)
    ax.grid(which='both')
    ax.grid(which='major', linestyle="-", color="grey")
    if colorbar:
        fig.colorbar(im)



def drawCoefficientwg(N, a, fig, ax, Greys=False):
    aCube = a.reshape(N, order='F')
    aCube = np.ascontiguousarray(aCube)
    if Greys:
        ax.imshow(aCube, cmap='Greys', origin='lower')
    else:
        ax.imshow(aCube, cmap=cm.hot_r, origin='lower')
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    
    # or if you want differnet settings for the grids:                               
    fig.subplots_adjust(left=0.02,bottom=0.00,right=0.98,top=1,wspace=0.05,hspace=0.00)
    


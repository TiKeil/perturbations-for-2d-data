import unittest
import numpy as np

#plots
import matplotlib.pyplot as plt
from matplotlib import cm

def drawCoefficient(N, a, fig, ax):
    aCube = a.reshape(N) #why order f and transponation? 
    te = 11
    major_ticks = np.arange(0, te, 1)
    minor_ticks = np.arange(0, te, 0.2)                                                
    ax.imshow(aCube, cmap=cm.Blues, extent=[0, te, 0, te])
    ax.axis([0, te, 0, te])
    ax.set_xticks(minor_ticks)
    ax.set_xticks(minor_ticks)
    ax.set_xticks(major_ticks)                                                       
    ax.set_yticks(major_ticks)  
    ax.minorticks_on()                                                     
    ax.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
    fig.subplots_adjust(left=0.00,bottom=0.02,right=1,top=0.98,wspace=0.2,hspace=0.2)
    # or if you want differnet settings for the grids:                               
    ax.grid(which='major', linestyle="-", color="blue")                                                
    ax.grid(which='minor', linestyle="-",linewidth=0.5 , color="grey")
    
    
def Patchprint1():
    #try to run normal pg
    #background
    bg = 0.01
    #values
    val = 1
    #fine World

    #coarse World
    NWorldCoarse = np.array([11,11])
    NpCoarse = np.prod(NWorldCoarse+1)

    #ratio between Fine and Coarse
    A = np.zeros(NWorldCoarse)
    # A[5,5] = 1.5
    # A[4][4:7] = 1.5
    # A[6][4:7] = 1.5
    # A[5][4] = 1.5
    # A[5][6] = 1.5

    ABase = A.flatten()   


    fig = plt.figure('bla2')                                                               
    ax = fig.add_subplot(1,1,1)     
    drawCoefficient(NWorldCoarse, ABase,fig,ax)

    plt.show()
    
    
class PetrovGalerkinRandom(unittest.TestCase):
    def test_Patchprint1(self):
        Patchprint1()


if __name__ == '__main__':
    unittest.main()
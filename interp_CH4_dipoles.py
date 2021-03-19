import sys
import os
import math
import numpy as np
import shutil
import scipy
from subprocess import call
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

griddipole = np.load('griddipole.npy')

gridsize = 0.1

gridz = np.ceil(32.2/gridsize) + 1 #191 fast  381 more accurate
gridx = np.ceil(3.056/gridsize) + 1 #30  59
gridy = np.ceil(5.15/gridsize) + 1 #51  101
gridx = int(gridx)
gridy = int(gridy)
gridz = int(gridz)
print("grid")
print(gridx)
print(gridy)
print(gridz)
#gridcoords = np.zeros((gridx,gridy,gridz,3))
#griddipole = np.zeros((gridx,gridy,gridz,3))
zc = np.linspace(-4.1,28.1,gridz) #.1 A extra on each side
xc = np.linspace(-.1,2.956231451,gridx) #.1 A extra on each side
yc = np.linspace(-2.573458221,2.573458221,gridy) #.1 A extra on each side
spacez = 32.2/(gridz-1)
spacex = 3.056231451/(gridx-1)
spacey = 5.146916442/(gridy-1)


dipole_interp = scipy.interpolate.RegularGridInterpolator((xc,yc,zc), griddipole)

import copy

def dipole_CH4(xyz):
    xyzp = copy.deepcopy(xyz)
    xyzp[:,0] = np.mod(xyz[:,0],2.856231451)
    xyzp[:,1] = np.mod(xyz[:,1]+ 2.473458221,4.946916442) - 2.473458221
    dipole = dipole_interp(xyzp)
    return dipole

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import copy

data_source = 'interpolated_data/'
griddipole = np.load(data_source+'griddipole.npy')

gridsize = 0.1

gridz = int(np.ceil(32.2/gridsize)) + 1 #191 fast  381 more accurate
gridx = int(np.ceil(3.056/gridsize)) + 1 #30  59
gridy = int(np.ceil(5.15/gridsize)) + 1 #51  101
zc = np.linspace(-4.1,28.1,gridz) #.1 A extra on each side
xc = np.linspace(-.1,2.956231451,gridx) #.1 A extra on each side
yc = np.linspace(-2.573458221,2.573458221,gridy) #.1 A extra on each side

dipole_interp = RegularGridInterpolator((xc,yc,zc), griddipole)

def dipole_CH4(xyz):
    xyzp = copy.deepcopy(xyz)
    xyzp[:,0] = np.mod(xyz[:,0],2.856231451)
    xyzp[:,1] = np.mod(xyz[:,1]+ 2.473458221,4.946916442) - 2.473458221
    return dipole_interp(xyzp)

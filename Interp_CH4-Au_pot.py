#!/usr/bin/python

import numpy as np
from scipy.interpolate import RegularGridInterpolator

data_source = 'DFT_interpolated_data/'
potgrad = np.load(data_source+'potgrad_v4.npy')
gridpot = np.load(data_source+'gridpot.npy')

gridsize = .1

gridx = int(np.ceil(29.2/gridsize) + 1)
gridy = int(np.ceil(3.056/gridsize) + 1)
gridz = int(np.ceil(5.15/gridsize) + 1)
xc = np.linspace(-1.1,28.1,gridx) #.1 A extra on each side
yc = np.linspace(-.1,2.956231451,gridy) #.1 A extra on each side
zc = np.linspace(-2.573458221,2.573458221,gridz) #.1 A extra on each side

potgrad_xz_fix = np.zeros(np.shape(potgrad))
potgrad_xz_fix[0,:,:,:] = potgrad[1,:,:,:]
potgrad_xz_fix[1,:,:,:] = potgrad[2,:,:,:]
potgrad_xz_fix[2,:,:,:] = potgrad[0,:,:,:]

force_interp = RegularGridInterpolator((yc,zc,xc), np.moveaxis(-potgrad_xz_fix,(0,1),(-1,-2)))
energy_interp = RegularGridInterpolator((yc,zc,xc), np.moveaxis(gridpot,0,-1))

def force_CH4_v4(xyz):
    xyzp = 1.0*xyz
    xyzp[:,0] = np.mod(xyz[:,0],2.856231451)
    xyzp[:,1] = np.mod(xyz[:,1],4.946916442) - 2.473458221
    return force_interp(xyzp)

def energyf_v4(xyz):
    xyzp = 1.0*xyz
    xyzp[:,0] = np.mod(xyz[:,0],2.856231451)
    xyzp[:,1] = np.mod(xyz[:,1],4.946916442) - 2.473458221
    return energy_interp(xyzp)

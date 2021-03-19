#!/usr/bin/python

import numpy as np
from scipy.interpolate import RegularGridInterpolator

potgrad = np.load('potgrad_v4.npy')
gridpot = np.load('gridpot.npy')

#Why do we need these?:
xyzpotgrad = np.zeros(np.shape(np.moveaxis(potgrad,(0,1),(-1,-2)))) #rearranged for z out of plane
xyzpotgrad[:,:,:,0] = np.moveaxis(potgrad,(0,1),(-1,-2))[:,:,:,1]
xyzpotgrad[:,:,:,1] = np.moveaxis(potgrad,(0,1),(-1,-2))[:,:,:,2]
xyzpotgrad[:,:,:,2] = np.moveaxis(potgrad,(0,1),(-1,-2))[:,:,:,0]

gridsize = .1

gridx = np.ceil(29.2/gridsize) + 1 #191 fast  381 more accurate (291 now??)
gridy = np.ceil(3.056/gridsize) + 1 #30  59
gridz = np.ceil(5.15/gridsize) + 1 #51  101
gridx = int(gridx)
gridy = int(gridy)
gridz = int(gridz)
print("grid")
print(gridx)
print(gridy)
print(gridz)
xc = np.linspace(-1.1,28.1,gridx) #.1 A extra on each side
yc = np.linspace(-.1,2.956231451,gridy) #.1 A extra on each side
zc = np.linspace(-2.573458221,2.573458221,gridz) #.1 A extra on each side
spacex = 29.2/(gridx-1)
spacey = 3.056231451/(gridy-1)
spacez = 5.146916442/(gridz-1)
#19 (29??)  2.856231451  4.946916442

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
    force = force_interp(xyzp)
    return force

def energyf_v4(xyz):
    xyzp = 1.0*xyz
    xyzp[:,0] = np.mod(xyz[:,0],2.856231451)
    xyzp[:,1] = np.mod(xyz[:,1],4.946916442) - 2.473458221
    E = energy_interp(xyzp)
    return E

"""
thing0 = np.array([0.,2.,1.])
thing1 = np.array([1.,1.,2.])
thing2 = np.array([2.,0.,1.])
thing3 = np.array([0.,1.,2.])

things = np.vstack((thing0,thing1,thing2,thing3))

print(things)

#print(energyf_v4(thing0))
#print(energyf_v4(thing1))
#print(energyf_v4(thing2))
#print(energyf_v4(thing3))

print(energyf_v4(things))
"""



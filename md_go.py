#!/usr/bin/env python3
#SBATCH -n 1
#SBATCH -J MD_run
#SBATCH -t 200:00:00
##SBATCH --account=brubenst-condo
##SBATCH --mem=8G

from argparse import ArgumentParser
from numpy import mod, zeros
from ase import Atoms, units
from ase.io import read
from ase.md.langevin import Langevin

# Another way to connect to ASE's path
#from sys import path
#from os import getcwd
#path.append(getcwd())

from interp_pot_custom_md import InterpPot_CH4

""" Arg Parsing """
ap = ArgumentParser()
ap.add_argument("-N", "--num", help="# of particles", type=int, required=True)
ap.add_argument("--mdstep", help="MD timestep (fs)", required=True)
ap.add_argument("-T", "--temp", help="Temperature (K)", required=True)
ap.add_argument("-L", "--Lang", help="Langevin friction coefficient (atomic units)", required=True)
ap.add_argument("--ns", help="Number of ns to run for", type=int, default=1000)
ap.add_argument("--name_extra", help="Name tag", default='main', type=str)
ap.add_argument("--short_time", help="Run for short time?", type=int, default=0)
ap.add_argument("--fixcm", help="fixcm tag for Lang", type=int, default=1)
ap.add_argument("--zero-vels-init", help="Initialize all velocities to 0 if no vellast or vel file?", type=int, default=0)
ap.add_argument("--xcellfactor", help="x6 or x5 for xcell", type=int, default=6)
ap.add_argument("--startpath", help="Startpath", type=str, default='')
ap.add_argument("--lastpath", help="Path to 'last' folders", type=str, default='')
ap.add_argument("--write-freq", help="Write to file every how many rec_amts?", type=int, default=10000)
ap.add_argument("--rec-amt", help="Record positions and velocities every how many fs?", type=int, default=1000)
args = vars(ap.parse_args())
print(args)

""" Particle setup """
N = args['num']
x_spacing = 2.856231451
y_spacing = 4.946916442
xcell = x_spacing*args['xcellfactor']*2
ycell = y_spacing*3*2
zcell = 15.0

""" File and Folders Info """
startpath = args['startpath']
lastpath = args['lastpath']
xyzout = 'xyzs_and_posfiles/'
xyzlast = 'xyzlasts/'
velpath = 'vels/'
vellast = 'vellasts/'

""" MD Simulation setup """
T = float(args['temp']) # Kelvin
MD_frame = float(args['mdstep']) # Time between each MD frame (and thus each energy and force calculation) in femtoseconds
ns = args['ns'] # Total simulation time (in nanoseconds)
Lang_coeff = float(args['Lang']) # Langevin friction coeff

""" Unit cell and starting molecule arrangement setup """
atoms = Atoms()
recstr = ''
if args['rec_amt'] != 1000:
	recstr = 'rec'+str(args['rec_amt'])+'fs_'
name = "%iX_%iK_%ifs_Lang%.4f_%s%s" % (N, int(T), int(MD_frame), Lang_coeff, recstr, args['name_extra'])

""" Reconstruct atoms object """
try:
    atoms = read(lastpath + xyzlast + name + '_pos_last.xyz')
except:
    atoms = read(startpath + xyzout + name + '_pos.xyz')
try:
    getvels = read(lastpath + vellast + name + '_vels_last.xyz')
except:
    try:
        getvels = read(startpath + velpath + name + '_vels.xyz')
    except:
        if bool(args['zero_vels_init']):
            print('No vellast or vel. So vels to zero to start')
            getvels = None
try:
    momenta = getvels.get_positions()*16.0
except:
    momenta = zeros((N,3))
try:
    del getvels
except:
    pass
masses = [16.0]*len(atoms)
symbols = ['X']*len(atoms)
atoms.set_masses(masses)
atoms.set_chemical_symbols(symbols)
atoms.set_momenta(momenta)
atoms.set_cell([xcell, ycell, zcell])
atoms.set_pbc([True,True,False])
atoms.set_calculator(InterpPot_CH4())

""" Quality Control of Atoms Object Before MD Sim """
print('-'*50)
print(name)
print('-'*50)
print(atoms)
print(len(atoms))
print(sum(atoms.get_masses())/len(atoms))
print(atoms.get_chemical_symbols()[0])
print(atoms.get_temperature())
if (sum(atoms.get_masses())/len(atoms) != 16.0):
    print('MASS ERROR')
    exit()

""" MD Simulation """
dyn = Langevin(atoms, MD_frame*units.fs, T*units.kB, Lang_coeff, fixcm=bool(args['fixcm']))
print(dyn.fixcm)
nsteps = ns/MD_frame # Total number simulation timesteps (in MILLIONS of timesteps)

def atom_wrap():
    pos = atoms.positions*1.0
    atoms.positions[:,0:2] = mod(pos[:,0:2], [xcell, ycell])
dyn.attach(atom_wrap, interval=1)

# Write to .xyz file every ps (or other timeframe)
#ps = int(1000/MD_frame)
record = int(args['rec_amt']/MD_frame)
X = str(atoms[0].symbol)
mass = atoms[0].mass

def xyz_acc(pos, xyz_to_write, vel_to_write):
    xyz_to_write += '%i\n%iK_%ifs_Lang%.4f\n' % (N, int(T), int(MD_frame), Lang_coeff)
    vel_to_write += '%i\n%iK_%ifs_Lang%.4f_vels\n' % (N, int(T), int(MD_frame), Lang_coeff)
    xyz_to_write += ''.join(['%s\t%.6f\t%.6f\t%.6f\n' % (X, pos[k,0], pos[k,1], pos[k,2]) for k in range(N)])
    vels = atoms.get_momenta()/mass
    vel_to_write += ''.join(['%s\t%.6f\t%.6f\t%.6f\n' % (X, vels[k,0], vels[k,1], vels[k,2]) for k in range(N)])
    return xyz_to_write, vel_to_write

def xyzwrite_logs(xyz_to_write, vel_to_write):
    with open(startpath + xyzout + name + '_pos.xyz', 'a') as xyzfile:
        xyzfile.write(xyz_to_write)
    with open(startpath + velpath + name + '_vels.xyz', 'a') as velfile:
        velfile.write(vel_to_write)

def xyzwrite_lasts(xyzlast_to_write, vellast_to_write):
    with open(lastpath + xyzlast + name + '_pos_last.xyz', 'w') as xyzlastfile:
        xyzlastfile.write(xyzlast_to_write)
    with open(lastpath + vellast + name + '_vels_last.xyz', 'w') as vellastfile:
        vellastfile.write(vellast_to_write)

# Run MD
xyz_to_write = ''
vel_to_write = ''
if not bool(args['short_time']):
    for i in range(1,int(nsteps*1000000/record + 1)):
        # Run for a number of timesteps (while recording positions and velocities)
        dyn.run(steps=record)
        xyz_to_write, vel_to_write = xyz_acc(atoms.positions*1.0, xyz_to_write, vel_to_write)
        # Then, write positions and velocities
        if i % args['write_freq']*record == 0:
            xyzwrite_logs(xyz_to_write, vel_to_write)
            xyzlast_to_write, vellast_to_write = xyz_acc(atoms.positions*1.0, xyz_to_write='', vel_to_write='')
            xyzwrite_lasts(xyzlast_to_write, vellast_to_write)
            xyz_to_write = ''
            vel_to_write = ''
else:
    # Short time option for debugging purposes
    dyn.run(steps=1000)

#!/usr/bin/env python

from ase.calculators.calculator import Calculator
import Interp_CH4_Au_pot as interp
from numpy import round, tril, isnan, isinf

#Code currently set up for point particles
class InterpPot_CH4(Calculator):

	implemented_properties = ['energy', 'forces']
	default_parameters = {'epsilon': 1.0,
				'sigma': 1.0,
				'rc': None}

	nolabel = True

	def __init__(self, **kwargs):
		Calculator.__init__(self, **kwargs)

	def calculate(self, atoms=None, properties=['energy'],
			system_changes=['positions', 'numbers', 'cell', 'pbc', 'charges', 'magmoms']):

		Calculator.calculate(self, atoms, properties, system_changes)
        
		# Interpolated Potential Portion (CH4 to Surface)
		positions = self.atoms.positions*1.0 # Should work just as well as get_positions
		cell = self.atoms.cell*1.0 # Should work just as well as get_cell
		natoms = len(self.atoms)

		energy = sum(interp.energyf_v4(positions))
		forces = interp.force_CH4_v4(positions)

		# OPLS Potential Portion (CH4 - CH4 interactions)
		sigma = 3.73 # in Angstroms
		epsilon = 0.012749 # in eV
		rc = 3 * sigma
		e0 = 4 * epsilon * ((sigma / rc)**12 - (sigma / rc)**6)

		# Get distances (while account for pbcs)
		xdiffs_raw = positions[:,0].reshape(natoms,1) - positions[:,0].reshape(1, natoms) 
		ydiffs_raw = positions[:,1].reshape(natoms,1) - positions[:,1].reshape(1, natoms)
		xdiffs = xdiffs_raw - (cell[0,0]*round(xdiffs_raw/cell[0,0])) #np.round part is key
		ydiffs = ydiffs_raw - (cell[1,1]*round(ydiffs_raw/cell[1,1]))
		zdiffs = positions[:,2].reshape(natoms,1) - positions[:,2].reshape(1, natoms)
		d2s = xdiffs**2 + ydiffs**2 + zdiffs**2
		c6 = (sigma**2/d2s)**3
		c6[isinf(c6)] = 0.0
		c6[d2s > rc**2] = 0.0
		c12 = c6**2

		# Energy after taking into account LJ
		energy = energy*1.0
		energy -= e0 * (tril(c6) != 0.0).sum()
		energy += 4 * epsilon * (tril(c12 - c6)).sum()

		# Forces after taking into account LJ
		f_core = 24 * epsilon * ((2*c12-c6)) / d2s
		f_core[isnan(f_core)] = 0.0
		forces[:,0] = forces[:,0] - tril(f_core * xdiffs).sum(0) + tril(f_core * xdiffs).sum(1) # "First step", "Second step"
		forces[:,1] = forces[:,1] - tril(f_core * ydiffs).sum(0) + tril(f_core * ydiffs).sum(1)
		forces[:,2] = forces[:,2] - tril(f_core * zdiffs).sum(0) + tril(f_core * zdiffs).sum(1)
		self.results['energy'] = energy
		self.results['forces'] = forces

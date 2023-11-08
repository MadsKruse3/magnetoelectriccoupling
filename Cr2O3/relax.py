import numpy as np
import json
from ase.io import read
from ase.phonons import Phonons
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, MixerDif
from gpaw.occupations import FermiDirac
from gpaw.borncharges import borncharges
from ase.optimize import BFGS 
from ase.constraints import UnitCellFilter

structure = read('structure.json')

#structure.set_initial_magnetic_moments(magmoms)

calc = GPAW(mode=PW(600),
            xc='LDA',
            convergence={'forces': 1.0e-4},
            occupations=FermiDirac(0.001),
            #kpts=(8,8,8),
            symmetry='off',
            kpts={'size': (12, 12, 12), 'gamma': True},
            txt='relax.txt')

structure.set_calculator(calc)

uf = UnitCellFilter(structure, mask=[1, 1, 1, 1, 1, 1])
opt = BFGS(uf)
opt.run(fmax=0.001)

structure.get_potential_energy()
calc.write('Cr2O3_relaxed.gpw')

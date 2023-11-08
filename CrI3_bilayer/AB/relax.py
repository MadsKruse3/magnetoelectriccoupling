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
            symmetry='off',
            convergence={'forces': 1.0e-4},
            occupations=FermiDirac(0.001),
            #parallel={'domain': 1, 'band': 1},
            kpts=(12,12,1),
            txt='relax.txt')

structure.set_calculator(calc)

#uf = UnitCellFilter(structure, mask=[1, 1, 1, 1, 1, 1])
opt = BFGS(structure)
opt.run(fmax=0.001)

structure.get_potential_energy()
calc.write('CrI3_AB_relaxed.gpw')

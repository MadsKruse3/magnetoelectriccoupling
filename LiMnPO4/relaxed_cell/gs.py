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

#structure = read('LiMnPO4_springermat.cif')
calc = GPAW('LiMnPO4_relaxed.gpw')


#magmoms = [0.0, 0.0, 0.0, 0.0,
#           0.0, 0.0, 0.0, 0.0,
#           0.0, 0.0, 0.0, 0.0,
#           1.0, -1.0, -1.0, 1.0,
#           0.0, 0.0, 0.0, 0.0,
#           0.0, 0.0, 0.0, 0.0,
#           0.0, 0.0, 0.0, 0.0]

#structure.set_initial_magnetic_moments(magmoms)

calc = GPAW(mode = PW(800),
            xc = 'LDA',
            #kpts = (4, 4, 6),
            kpts={'size': (4, 4, 6), 'gamma': True},
            symmetry = 'off',
            parallel = {'domain': 1, 'band': 1},
            occupations = FermiDirac(0.001),
            maxiter = 5000,
            #experimental = {'magmoms': magmoms, 'soc': True},
            #nbands=120,
            eigensolver=Davidson(3),
            mixer = MixerDif(beta = 0.02, nmaxold=2, weight=100, beta_m=0.02, nmaxold_m=2, weight_m = 100),     
            convergence = {'density': 1.0e-9},
            txt = 'shifted_Ez%1.2f_madkru.txt' % i)            


structure.set_calculator(calc)

#uf = UnitCellFilter(structure, mask=[1, 1, 1, 1, 1, 1])
#opt = BFGS(uf)
#opt.run(fmax=0.001)

structure.get_potential_energy()
calc.write('LiMnPO4.gpw')

import numpy as np
import json
from ase.io import read
from ase.phonons import Phonons
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, MixerDif
from gpaw.occupations import FermiDirac
from gpaw.borncharges import borncharges 

structure = read('LiMnPO4_relaxed.cif')

calc = GPAW(mode=PW(600),
            xc='LDA',
            symmetry={'point_group': False},
            convergence={'forces': 1.0e-6},
            occupations=FermiDirac(0.001),
            parallel={'domain': 1, 'band': 1},
            kpts=(8,8,1),
            txt='phonons.txt')

structure.set_calculator(calc)
structure.get_potential_energy()

calc.write('LiMnPO4.gpw')

ph = Phonons(structure, calc, supercell=(1,1,1), delta=0.01)
ph.run()
ph.read(symmetrize=10, acoustic=True, born=True)

calc2 = GPAW('LiMnPO4.gpw', txt=None)
borncharges(calc2)


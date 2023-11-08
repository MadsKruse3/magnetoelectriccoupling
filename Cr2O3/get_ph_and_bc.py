import numpy as np
import json
from ase.io import read
from ase.phonons import Phonons
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, MixerDif
from gpaw.occupations import FermiDirac
from gpaw.borncharges import borncharges 

#structure = read('CrI3_AB_relaxed.gpw')
calc = GPAW('Cr2O3_relaxed.gpw')
structure = calc.atoms

#calc = GPAW(mode=PW(600),
#            xc='LDA',
            #symmetry={'point_group': False},
#            convergence={'forces': 1.0e-4},
#            occupations=FermiDirac(0.001),
#            parallel={'domain': 1, 'band': 1},
#            kpts=(4,4,1),
#            txt='phonons_CrI3_AB.txt')

#structure.set_calculator(calc)
#structure.get_potential_energy()
#calc.write('CrI3_AB.gpw')

borncharges(calc)

#ph = Phonons(structure, calc, supercell=(1,1,1), delta=0.01)
#ph.run()
#ph.read(symmetrize=10, acoustic=True)

#borncharges(calc)


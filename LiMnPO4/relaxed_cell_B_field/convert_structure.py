from ase.io import write
from gpaw import GPAW

calc = GPAW('LiMnPO4_relaxed.gpw')

atom = calc.atoms

atom.write('structure.json')

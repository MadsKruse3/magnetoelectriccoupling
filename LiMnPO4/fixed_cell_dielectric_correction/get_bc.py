import numpy as np
import json
from ase.io import read
from ase.phonons import Phonons
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, MixerDif
from gpaw.occupations import FermiDirac
from gpaw.borncharges import borncharges 


calc2 = GPAW('LiMnPO4.gpw', txt=None)

calc2.set(maxiter=1000)

borncharges(calc2)

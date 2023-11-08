import numpy as np
import json
from ase.io import read
from ase.phonons import Phonons
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, MixerDif, MixerSum
from gpaw.occupations import FermiDirac
from gpaw.borncharges import borncharges 

structure = read('LiMnPO4_relaxed.gpw')
calc = GPAW('LiMnPO4_relaxed.gpw', txt=None)

ph = Phonons(structure, calc, supercell=(1,1,1), delta=0.01)
#ph = Phonons(structure)
#ph.run()

ph.read(symmetrize=10, acoustic=True)

C_ij = ph.get_force_constant()[0]

eigen_v, eigenvec = np.linalg.eigh(C_ij)

Omega = np.diag(eigen_v[3:])

A = eigenvec[:, 3:]

C_inv_tilde = np.linalg.multi_dot([A, np.linalg.inv(Omega), A.T])

D = C_inv_tilde.shape[1]

with open('borncharges-0.01.json') as fd:
    Z_avv = json.load(fd)['Z_avv']
    
#allborncharges = np.asarray(Z_avv) # Born effective charge tensor 
                                   # - [Atom #Number, response in i, 
                                   #  given applied field in j] 
                                   # - Z^{a}_{i,j}

E_k = np.array([1,0,0])

E_max = 0.01

e_k = np.reshape(E_k / np.linalg.norm(E_k),(3,1)) ## ?? 

res = 8

pos_equilibrium = structure.get_positions()

magmoms = np.zeros((D//3, 3), float) # Do this smarter so its general
## Look up LiMnPO4 magnetic atoms.. 
magmoms[13] = [5.0, 0.0, 0.0]
magmoms[14] = [-5.0, 0.0, 0.0]
magmoms[15] = [-5.0, 0.0, 0.0]
magmoms[16] = [5.0, 0.0, 0.0]

#params = dict(mode=PW(600),
#              xc = 'LDA',
#              kpts=(4,4,6),
#              txt='relax.txt',
#              symmetry='off',
#              parallel={'domain': 1, 'band': 1},
#              occupations={'name': 'fermi-dirac','width': 0.001},
#              experimental={'magmoms': magmoms, 'soc': True})

F = C_inv_tilde.dot(np.reshape(Z_avv,(D,3)))

U = 3.0
Es = []
magmoms_x = []
magmoms_y = []
magmoms_z = []
for i in range(0, res + 1):
#for i in range(0, res + 1):
    E = (-E_max)*e_k + i*(2/res)*E_max*e_k
    ub_j = np.reshape((F.dot(E))[:,0] , (D//3, 3))
    pos_new = pos_equilibrium + ub_j
    
    calc = GPAW(mode = PW(600),
                xc = 'LDA',
                #kpts = (4, 4, 6),
                kpts={'size': (4, 4, 6), 'gamma': True},
                symmetry = 'off',
                parallel = {'domain': 1, 'band': 1},
                occupations = FermiDirac(0.001),
                maxiter = 5000,
                setups={'Mn': f':d,{U}'},
                experimental = {'magmoms': magmoms, 'soc': True},
                #mixer = MixerDif(0.05, 3, 100),   
                #mixer = MixerDif(0.05, 3, 20),   
                mixer = MixerDif(beta = 0.05, nmaxold=3, weight=100, beta_m=0.05, nmaxold_m=3, weight_m = 100),     
                convergence = {'density': 1.0e-9},
                txt = 'shifted_Ez%1.2f_madkru.txt' % i)            
    
    structure.set_positions(pos_new)
    structure.set_calculator(calc)
    structure.get_potential_energy()

    m_v, m_av = calc.density.estimate_magnetic_moments()
    #magmom  = calc.get_magnetic_moment()
    magmoms_x.append(m_v[0])
    magmoms_y.append(m_v[1])
    magmoms_z.append(m_v[2])
    Es.append(E)

np.savetxt('Ez.npy', Es)
np.savetxt('magmoms_x.npy', magmoms_x)
np.savetxt('magmoms_y.npy', magmoms_y)
np.savetxt('magmoms_z.npy', magmoms_z)

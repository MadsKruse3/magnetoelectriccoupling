import numpy as np
import json
from ase.io import read
from ase.phonons import Phonons
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, MixerDif
from gpaw.occupations import FermiDirac
from pprint import pprint

bulk = read('LiMnPO4_relaxed.gpw')

# Phonons
ph = Phonons(bulk)
ph.read(symmetrize=10, acoustic=True)
C_ij = ph.get_force_constant()[0]
c_n, c_in = np.linalg.eigh(C_ij)

# Trace out the acoustic modes and invert
ct_n = c_n[3:]
ct_in = c_in[:, 3:]
Ctinv_ij = ct_in.dot(np.diag(1 / ct_n)).dot(ct_in.T)

# Born Charges
with open('borncharges-0.01.json') as fd:
    Z_avv = json.load(fd)['Z_avv']
z_avv = np.zeros((len(Z_avv), 3, 3), float)
for a, Z_vv in enumerate(Z_avv):
    z_avv[a] = Z_vv
    #pprint(np.round(Z_vv, 2))

# Calculate displacement with E_x = 0.1 V / AA 
# C_abij * u_ai = E_i * Z_bij => u_ai = Cinv_abij * E_k * Z_bkj

D = Ctinv_ij.shape[1]

#E_x = 0.1
#ExZ_av = E_x * z_avv[:, 0] 
#r_av = np.reshape(np.dot(Ctinv_ij, np.reshape(ExZ_av, -1)), (10, 3))
#E_y = 0.1
#EyZ_av = E_y * z_avv[:, 1] 
#r_av = np.reshape(np.dot(Ctinv_ij, np.reshape(EyZ_av, -1)), (10, 3))
E_z = 0.1
EzZ_av = E_z * z_avv[:, 2] 
#r_av = np.reshape(np.dot(Ctinv_ij, np.reshape(EzZ_av, -1)), (10, 3))
r_av = np.reshape(np.dot(Ctinv_ij, np.reshape(EzZ_av, -1)), (D//3,3) )

# Calculate magnetization
pos = bulk.get_positions()
pos -= 1.2 * r_av

magmoms = np.zeros((D//3, 3), float)
magmoms[4] = [0.0, 0.0, -3.0]
magmoms[5] = [0.0, 0.0, 3.0]
magmoms[6] = [0.0, 0.0, 3.0]
magmoms[7] = [0.0, 0.0, -3.0]

magmoms[13] = [5.0, 0.0, 0.0]
magmoms[14] = [-5.0, 0.0, 0.0]
magmoms[15] = [-5.0, 0.0, 0.0]
magmoms[16] = [5.0, 0.0, 0.0]


for i in range(11):
    pos += 0.2 * r_av
    bulk.set_positions(pos)
    #E = -E_x + 0.2 * i * E_x
    #E = -E_y + 0.2 * i * E_y
    E = -E_z + 0.2 * i * E_z
    calc = GPAW(mode=PW(600),
                xc='LDA',
                symmetry='off',
                parallel={'band':1, 'domain':1},
                kpts=(4, 4, 6),
                maxiter=5000,
                experimental={'magmoms': magmoms, 'soc': True},
                mixer=MixerDif(0.05, 3, 100),
                convergence={'density': 1.0e-9},
                occupations=FermiDirac(width=0.001),
                #txt='shifted_Ex%1.2f.txt' % E,
                #txt='shifted_Ey%1.2f.txt' % E,
                txt='shifted_Ez%1.2f_tols.txt' % E,
                )
    bulk.set_calculator(calc)
    bulk.get_potential_energy()
    bulk.get_forces()

    m_v, m_av = calc.density.estimate_magnetic_moments()
    #f = paropen('magmoms_Ex.dat', 'a')
    #f = paropen('magmoms_Ey.dat', 'a')
    f = paropen('magmoms_Ez_tols.dat', 'a')
    print(E, m_v[0], m_v[1], m_v[2], file=f)
    f.close()

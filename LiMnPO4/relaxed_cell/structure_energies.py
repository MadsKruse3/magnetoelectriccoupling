import numpy as np 
from ase.io import read
import matplotlib.pyplot as plt
from ase.parallel import paropen
from ase.io import read

energies = np.loadtxt('energies.npy')
Ez = np.loadtxt('Ez.npy')
structure = read('LiMnPO4_relaxed.gpw')
V = structure.get_volume()

energies = energies - energies[-1]
energies = energies*1e3


plt.figure()
plt.plot(Ez, energies, '-ro')
plt.xlabel('Electric field along z [V/Å]')
plt.ylabel('Relative energies [meV]')
plt.tight_layout()
plt.savefig('structure_energies.png')
plt.close()

energies = energies/V

plt.figure()
plt.plot(Ez, energies, '-ro')
plt.xlabel('Electric field along z [V/Å]')
plt.ylabel('Relative energies pr. unit cell volume [meV/Å$^{3}$]')
plt.tight_layout()
plt.savefig('structure_energies_pr_unitcell.png')
plt.close()

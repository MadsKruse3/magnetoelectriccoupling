import numpy as np 
from ase.io import read
import matplotlib.pyplot as plt

mu_0 = 4*np.pi*10**(-7)

structure = read('LiMnPO4_relaxed.gpw')

cell = structure.get_cell()

pos = structure.get_positions()
c = pos[:, 2].max() - pos[:, 2].min()

#filename = ??
#magmoms_Ex = np.loadtxt(filename, usecols=(3,4,5))
#Ex = np.loadtxt(filename, usecols=(0))

#filename = ??
#magmoms_Ey = np.loadtxt(filename, usecols=(3,4,5))
#Ey = np.loadtxt(filename, usecols=(1))

filename = 'magmoms_Ez.dat'
Ez = np.loadtxt(filename, usecols=(0))
magmom_Ex = np.loadtxt(filename, usecols=(1))
magmom_Ey = np.loadtxt(filename, usecols=(2))
magmom_Ez = np.loadtxt(filename, usecols=(3))

V = np.linalg.norm(np.cross(cell[0,:], cell[1,:])) * c
magnetization = magmom_Ez / V
plt.figure()

plt.plot(Ez, magnetization, '-ro')
#plt.plot(Ez, magnetization[1], '-bo')
#plt.plot(Ez, magnetization[2], '-go')

plt.xlabel('Ez')
plt.ylabel('Magnetization')
plt.title('Magnetoelectric coupling')

plt.tight_layout()
plt.savefig('me_coupling_plot.png')

def estimate_coef(x,y):

    n = np.size(x)

    m_x, m_y = np.mean(x), np.mean(y)

    SS_xy = np.sum(x*y) - n*m_x*m_y

    SS_xx = np.sum(x*x) - n*m_x*m_x
    
    SS_yy = np.sum(y*y) - n*m_y*m_y

    alpha = SS_xy / SS_xx
    offset = m_y - alpha*m_x

    SS_R = np.sum(((alpha*x + offset) - m_y)**2)
    r_squared = SS_R / SS_yy

    return(alpha, r_squared)

Alpha = np.zeros([3,3])
R_squared = np.zeros([3,3])

#for i in range(0,3):
#    for j in range(0,3):
#        Alpha[i,j], R_squared[i,j] = estimate_coef(E[:,i], mu_0*magnetization[i,:,j])


#Alpha *= 9.274009994*10**8

#f = open('alpha_LiMnPO4.txt', 'a')
#f.write('alpha: ' + '\n')
#print(Alpha, file=f)
#f.write('R^2: ' + '\n')
#print(R_squared, file=f)
#f.close()




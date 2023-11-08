import numpy as np 
from ase.io import read
import matplotlib.pyplot as plt
from ase.parallel import paropen
from ase.io import read

mu_0 = 4*np.pi*10**(-7)

structure = read('Cr2O3_relaxed.gpw')

cell = structure.get_cell()
V = structure.get_volume()
V = V*1e-30

magmoms_Ex = np.loadtxt('magmoms_Ex.npy')
magmoms_Ey = np.loadtxt('magmoms_Ey.npy')
magmoms_Ez = np.loadtxt('magmoms_Ez.npy')
Ex = np.loadtxt('Ex.npy')
Ey = np.loadtxt('Ey.npy')
Ez = np.loadtxt('Ez.npy')

#filename = 'magmoms_Ez.txt'
#Ez = np.loadtxt(filename)

#magmom_Ex = np.loadtxt(filename, usecols=(1))
#magmom_Ey = np.loadtxt(filename, usecols=(2))
#magmom_Ez = np.loadtxt(filename, usecols=(3))
#magmom = np.loadtxt(filename, usecols=(1,2,3))

#magmom_Ex = magmoms[:,0]
#magmom_Ey = magmoms[:,1]
#magmom_Ez = magmoms[:,2]

magmoms_Ex = mu_0*magmoms_Ex/V
magmoms_Ey = mu_0*magmoms_Ey/V
magmoms_Ez = mu_0*magmoms_Ez/V

# Unit convergence
Ex = Ex*1e10
Ey = Ey*1e10
Ez = Ez*1e10

magmoms_Ex = (1e-12)*magmoms_Ex
magmoms_Ey = (1e-12)*magmoms_Ey
magmoms_Ez = (1e-12)*magmoms_Ez

#plt.figure()
#plt.plot(Ex, magmoms_Ex, '-ro')
#plt.plot(Ey, magmoms_Ey, '-bo')
#plt.plot(Ez, magmoms_Ez, '-go')
#plt.legend(['M$_{x}$','M$_{y}$','M$_{z}$'])
#plt.xlabel('Electric field along x [V/Å]')
#plt.ylabel('Magnetization pr. unit cell [M/Å$^{3}$]')
#plt.tight_layout()
#plt.savefig('me_coupling_plot.png')

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
    print(r_squared)
    return(alpha, r_squared)

#Alpha = np.zeros([3,3])
#R_squared = np.zeros([3,3])

#for i in range(0,3):
#    for j in range(0,3):
#        Alpha[i,j], R_squared[i,j] = estimate_coef(Ez[i], mu_0*magnetization[i,:,j])

Alpha = np.zeros([3,3])
R_squared = np.zeros([3,3])

Alpha[0][0], R_squared[0][0] = estimate_coef(Ex[:,0], magmoms_Ex[:,0])
Alpha[0][1], R_squared[0][1] = estimate_coef(Ex[:,0], magmoms_Ey[:,0])
Alpha[0][2], R_squared[0][2] = estimate_coef(Ex[:,0], magmoms_Ez[:,0])
Alpha[1][0], R_squared[1][0] = estimate_coef(Ey[:,1], magmoms_Ex[:,1])
Alpha[1][1], R_squared[1][1] = estimate_coef(Ey[:,1], magmoms_Ey[:,1])
Alpha[1][2], R_squared[1][2] = estimate_coef(Ey[:,1], magmoms_Ez[:,1])
Alpha[2][0], R_squared[2][0] = estimate_coef(Ez[:,2], magmoms_Ex[:,2])
Alpha[2][1], R_squared[2][1] = estimate_coef(Ez[:,2], magmoms_Ey[:,2])
Alpha[2][2], R_squared[2][2] = estimate_coef(Ez[:,2], magmoms_Ez[:,2])

print(Alpha)
print(R_squared)

#Alpha *= 9.274009994*10**8

#f = open('alpha_LiMnPO4.txt', 'a')
#f.write('alpha: ' + '\n')
#print(Alpha, file=f)
#f.write('R^2: ' + '\n')
#print(R_squared, file=f)
#f.close()




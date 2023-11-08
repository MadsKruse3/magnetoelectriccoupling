import numpy as np 
from ase.io import read
import matplotlib.pyplot as plt
from ase.parallel import paropen
from ase.io import read

mu_0 = 4*np.pi*10**(-7)

structure = read('CrI3_AB_prime_relaxed.gpw')

cell = structure.get_cell()
V = structure.get_volume()

#magmoms_x = np.loadtxt('magmoms_x.npy')
#magmoms_y = np.loadtxt('magmoms_y.npy')
#magmoms_z = np.loadtxt('magmoms_z.npy')
#Ez = np.loadtxt('Ez.npy')

#Ex = [[-0.03, 0. ,0.  ], 
#      [-0.0225, 0., 0.],
#      [-0.015, 0., 0.],
#      [-0.0075, 0., 0.],
#      [0., 0., 0.],
#      [0.0075, 0., 0.],
#      [0.015, 0., 0.],
#      [0.0225, 0., 0.],
#      [0.03, 0., 0.]]

Ex = [-0.03,   
      -0.0225,
      -0.015,
      -0.0075,
      0.,
      0.0075,
      0.015,
      0.0225,
      0.03]

Ex = np.array(Ex)
magmoms_z = [-2.2937035571660852e-14,
             4.1837952373391385e-14,
            1.3852234517861689e-14,
             -2.503852770193056e-14,
            4.959633747181691e-14,
            2.5007377218262336e-14,
            -4.204023061745929e-14,
            -3.467422882022855e-14,
            -2.161275038060559e-14]

#magmoms_z = [[0.0, 0.0, -2.2937035571660852e-14],
#           [0.0, 0.0, 4.1837952373391385e-14],
#            [0.0, 0.0, 1.3852234517861689e-14],
#            [0.0, 0.0, -2.503852770193056e-14],
#            [0.0, 0.0, 4.959633747181691e-14],
#            [0.0, 0.0, 2.5007377218262336e-14],
#            [0.0, 0.0, -4.204023061745929e-14],
#            [0.0, 0.0, -3.467422882022855e-14],
#            [0.0, 0.0, -2.161275038060559e-14]]



#filename = 'magmoms_Ez.txt'
#Ez = np.loadtxt(filename)

#magmom_Ex = np.loadtxt(filename, usecols=(1))
#magmom_Ey = np.loadtxt(filename, usecols=(2))
#magmom_Ez = np.loadtxt(filename, usecols=(3))
#magmom = np.loadtxt(filename, usecols=(1,2,3))

#magmom_Ex = magmoms[:,0]
#magmom_Ey = magmoms[:,1]
#magmom_Ez = magmoms[:,2]

#magmoms_x = mu_0*magmoms_x/V
#magmoms_y = mu_0*magmoms_y/V
magmoms_z = mu_0*np.array(magmoms_z)/V

plt.figure()
#plt.plot(Ez, magmoms_x, '-ro')
#plt.plot(Ez, magmoms_y, '-bo')
plt.plot(Ex, magmoms_z, '-go')
#plt.legend(['M$_{x}$','M$_{y}$','M$_{z}$'])
plt.xlabel('Electric field along x [V/Å]')
plt.ylabel('Magnetization pr. unit cell [M/Å$^{3}$]')

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
    print(r_squared)
    return(alpha, r_squared)

#Alpha = np.zeros([3,3])
#R_squared = np.zeros([3,3])

#for i in range(0,3):
#    for j in range(0,3):
#        Alpha[i,j], R_squared[i,j] = estimate_coef(Ez[i], mu_0*magnetization[i,:,j])

Alpha = np.zeros([3,1])
R_squared = np.zeros([3,1])

#Alpha[0][0], R_squared[0][0] = estimate_coef(Ez[:], magmoms_x[:])
#Alpha[1][0], R_squared[1][0] = estimate_coef(Ez[:], magmoms_y[:])
Alpha[2][0], R_squared[2][0] = estimate_coef(Ex[:], magmoms_z[:])

print(Alpha)
print(R_squared)

#Alpha *= 9.274009994*10**8

#f = open('alpha_LiMnPO4.txt', 'a')
#f.write('alpha: ' + '\n')
#print(Alpha, file=f)
#f.write('R^2: ' + '\n')
#print(R_squared, file=f)
#f.close()




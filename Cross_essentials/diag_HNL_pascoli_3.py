import numpy as np
import matplotlib.pyplot as plt
import math
import h5py as h5
import sys
from numpy import linalg as la
e_m = 0.001



y_nu = np.ones((3,1))
N_hnl=3 # nr of HNLS
N_d =2

m_n =1
v_H = 2.46 #10^2 GeV
v_phi = 2.40 

y_nu = np.ones((3,1)).T #times yukawa coupling
y_N =np.ones((N_hnl,1)).T #times yukawa coupling
m_D = v_H/np.sqrt(2)  *y_nu #(N,3) #one Md and filled with zeros
mu =np.ones((N_hnl,N_hnl))*m_n
L = v_phi/np.sqrt(2)*y_N #(1,N)
LR = v_phi/np.sqrt(2)*y_N #(1,N)

m_D[0,0] = 0.00950 #(N,3) #one Md and filled with zeros
m_D[0,1] = 0.278 #(N,3) #one Md and filled with zeros
m_D[0,2] = 0.190 #(N,3) #one Md and filled with zeros

L[0,0] = -2.39  *10  #MEV
L[0,1] = 19.0 *10
L[0,2] =0.00 *10

LR[0,0] = -2.39 *10
LR[0,1] = 19.0  *10 #MeV 
LR[0,2] =0.00 *10

mu_s = np.array([-0.0429,1.10,-1.10])*10**3
mu = np.diag(mu_s)

m_x=-1.21*10**2

A = np.zeros((3,3))
z = np.zeros(np.shape(m_D))


M = np.block([[0,0,0,z,0,0 ],
              [0,0,0,z,0,0 ],
              [0,0,0,m_D,0,0 ],
               [z.T,z.T,m_D.T,mu,L.T,LR.T ],
               [0,0,0,L,0,m_x],
               [0,0,0,LR,m_x,0]])

e,v =la.eig(M)

#input parameters theory

m_D = np.array([[0.00950, -0.0347, 0.0336, 0.120],
               [0.287, 1.98, -0.635, 6.72],
                [0.190, -3.89, -1.03, -11.4]]) #MeV

m_diag = np.array([[-0.0429, -0.0900, -0.0963, 0.206],
               [1.10, 6.00, 5.07, 6.00],
                   [-1.10, -18.0, -10.1, -18.0]]) *1000 #MeV

L_L = np.array([[-2.39, 3.75, 3.51, 15.0],
               [19.0, 24.0, 25.5, 33.0],
               [0.00, 0.00, 12.7, 0.00]]) *10 #MeV

L_R = np.array([[-2.39, -2.81, -4.04, -14.9],
               [19.0, 54.0, 44.1, -12.9],
                [0.00, 0.00, -38.1, 0.00]]) *10 #MeV

M_X = np.array([-1.21, 1.96, 1.56, 3.50]) *100 #MeV

masses = np.zeros((4,8))
VU = np.zeros(4) #magnitude VU
vectors = np.zeros((4,8,8))
M = np.zeros((8,8))

###--------------------------------------------------------------
### Constructing the axial matrix (could be improved? Block)#####
####                   and calc eig val and vectors        ######
#################################################################
for i in range(4):
    k =2 
    M[0,3:6] = m_D[:,i]
    M[1,3:6] = m_D[:,i]
    M[2,3:6] = m_D[:,i]

    M[3:6,0] = m_D[:,i]
    M[3:6,1] = m_D[:,i]
    M[3:6,2] = m_D[:,i]


    M[1+k:4+k,1+k:4+k] = np.diag(m_diag[:,i])
    M[4+k,1+k:4+k] = L_L[:,i]
    M[1+k:4+k,4+k] = L_L[:,i]
    M[5+k,1+k:4+k] = L_R[:,i]
    M[1+k:4+k,5+k] = L_R[:,i]
    M[4+k,5+k] = M_X[i]
    M[5+k,4+k] = M_X[i]

    #eigval calc, eig vector
    e, v = la.eig(M)
    e = np.sqrt(e*e)
    ind_sort = np.argsort(e)
    masses[i,:] = e[ind_sort] #store eigval, masses
    vectors[i,:,:] = v[:,ind_sort] #store vectors


    
hf = h5.File('HNL_mass_pascoli_3.hdf5', 'w')
hf.create_dataset('m_D', data=m_D)
hf.create_dataset('L_L', data=L_L)
hf.create_dataset('L_R', data=L_R)
hf.create_dataset('mu', data=m_diag)
hf.create_dataset('M_X', data=M_X)
hf.create_dataset('M_diag', data=masses)
hf.create_dataset('mix', data=vectors)
hf.close()






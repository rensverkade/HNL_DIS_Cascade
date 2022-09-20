#! /usr/bin/env python

## Python LHAPDF6 usage example

import lhapdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py as h5
from scipy import integrate as sint

lhapdf.setPaths(['/home/rverkade/.local/share/lhapdf'])
print(lhapdf.paths())
## Getting a PDF member object
p = lhapdf.mkPDF("CT18NLO/0")

s2_W =0.22290 #sin_W^2
## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g

hf = h5.File('/home/rverkade/Documents/BO_thesis/HNL_mass_pascoli_3.hdf5', 'r')
N_masses = np.array(hf.get('M_diag')) #neutrino masses from diag
mix = np.array(hf.get('mix'))  #mixing  matrix VU
hf.close()
#mix = mix[0,:,:]
V=np.zeros((8,8))
Vert = np.zeros((4,3))

for l in range(4):  #Benchmark
    for z in range(3):  #HNL 1, 2 ,3
        V=np.zeros((8,8))
        for i in range(8):
            for j in range(8):
                for d in range(6,8): #dark states
                    #print('D',d)
                    V[i,j] =+ mix[l,d,i]*mix[l,d,j]    
        Vert[l,z] = np.sum(mix[l,2,:3])*np.sum(V[:3,3+z])/np.sqrt(np.sum(mix[l,:3,:3]**2))
        

#----------------
############
# Time for #
# F1 F2 FU #
############

# LIST of PARAMETER values
# we need to make so you can input these or
#take them from the ones defined in DARKNEWS
e_q = np.array([1/3,2/3,1/3,2/3]) #el charge of d,u,s,c
G_F = 1.188*10**(-5) #GeV^-2
M_W =  80.38 #GeV
mp = 0.94 #GeV
m_mu = 0.105 #GeV muon
mN = 0.185 #GeV neutrino
R =0.3 #Ratio for F_L calc
R_1 = 1/R
G_F = 1.188*10**(-5) #GeV^-2 Fermi constant
M_W =  80.38 #GeV

M_z = 91 #GeV Z-boson mass
g_A = np.array([-0.5, 0.5,-0.5,0.5]) #coupling of quarks = isospin
g_V = g_A -2*e_q*s2_W
g_X_c = 0.32*4*np.pi #dark coupling
X_c = 0.0003 #10^-4 dark mixing
W = np.arcsin(np.sqrt(0.22290)) #Weinberg mixing angle
M_Z2= 1.25 #GeV Dark photon mass
M_rat_2 =( M_Z2/M_z)**2  #ratio of Z and Z' mass
v_f = 500 #MeV VEV of dark higgs
v_H = 265*10**3 #MeV Higgs VEV
g_g = M_z*2/v_H #=sqrt(g'^2+g^2)  

#calc of geometric vals necessary for L_int calcs
X = X_c
g_X = g_X_c
s_W = np.sin(W)
c_W = np.cos(W)
e = c_W/g_g *s_W #g*s_W
#print('e',e)
s_X = np.sin(X)
c_X = np.cos(X)
t_X = s_X/c_X
t_B =t_X*s_W/(1-M_rat_2) #meaning NC Z'=0    
c_B = 1/(np.sqrt(1+t_B**2))
s_B = t_B*c_B



# calculating the L_int vals for axial (A=True) and vectorial coupling
def CCF_2_3(Q_e,Q_L,Q_X,A): 

    if A== False:              
        T_1 = e/(2*s_W*c_W) * ( -Q_L*(c_B+s_B*s_W*t_X)+2*Q_e*s_W*(-s_B*t_X + c_B*s_W))
        T_2 = g_X*Q_X*s_B/c_X #Z

        T_3 = e/(2*s_W*c_W) * ( Q_L*(s_B-c_B*s_W*t_X)*c_B+2*Q_e*s_W*(c_B*t_X - s_B*s_W) )
        T_4 = g_X*Q_X*c_B/c_X #Z'

    elif A== True:
        T_1 =  e/(2*s_W*c_W) * (Q_L*(c_B+s_B*s_W*t_X))
        T_2 = 0
        T_3 = e/(2*s_W*c_W) * ( -Q_L*(s_B-c_B*s_W*t_X)*c_B )
        T_4 = 0 #Z'

    return(T_1 - T_2,T_3 -T_4) # Z and Z' int

flavs = p.flavors()[1:-2]#d,u,s,c excluding b and gluon
f_len = len(flavs)
f_half = int(f_len*0.5)

zeros = np.zeros(f_half)

################################################
##Calculating the couplings                   ##
## Vectorial and axial for Z (Z1) and Z' (Z2) ##
################################################
# Z-boson Z'-phot,charges: e,I3,Q_X, axiality
g_Z1_V,g_Z2_V =  CCF_2_3(e_q,g_A,zeros,A=False) 
g_Z1_A,g_Z2_A =  CCF_2_3(e_q,g_A,zeros,A=True)
    
g_Z1_V_D,g_Z2_V_D =  CCF_2_3(0,0,1,A=False)
g_Z1_A_D,g_Z2_A_D =  CCF_2_3(0,0,1,A=True)

g_Z1_V_H,g_Z2_V_H =  CCF_2_3(1,0,0,A=False)
g_Z1_A_H,g_Z2_A_H =  CCF_2_3(1,0,0,A=True)

# structure functions calc, likely need to add interference lvls
def struct_PDF(x,Q2):
    xf = np.array(p.xfxQ2(flavs, x, Q2)) #notice func of Q^2    quarks           anti-quarks
    #xf_g = np.array(p.xfxQ2(21, x, Q2))
    #F2 = np.sum(xf[f_half:]+xf[:f_half]) #sea qu
    F2_Z2 = np.sum((g_Z2_V**2+g_Z2_A**2) *(xf[f_half:]+xf[:f_half]))
    #F2_yZ2 = np.sum(2*e_q*g_Z2_V *(xf[f_half:]+xf[:f_half]))
    xF3_Z2 = np.sum(2*g_Z2_V*g_Z2_A *(xf[f_half:]-xf[:f_half]))
    #xF3_yZ2 = np.sum(2*e_q*g_Z2_A *(xf[f_half:]-xf[:f_half]))
    FL_Z2= (R_1+1)*F2_Z2 #WHERE am is supposeed to find disnshit
    
    x2F1_Z2 = F2_Z2*(1+4*mp**2 * x**2 /Q2)/(1+FL_Z2)
    return (x2F1_Z2,F2_Z2,xF3_Z2)   # (x2F1,



G_F = 1.188*10**(-5) #GeV^-2
mp = 0.94 #GeV
M_W =  80.38 #GeV

#differential cross section, likely what DN needs       
def ds_v(x,y,Vert,E,anti): #M is nucleon mass
    Q2 = x*y*2*mp*E
    if anti ==1:
        x2F1, F2,xF3 = struct_PDF(x,Q2)
    elif anti== -1:
        x2F1, F2,xF3 = struct_PDF(x,Q2)
        
    ds1 = 0.5*y**2 *x2F1 +(1-y-x*y*mp/(2*E))*F2 # +0.5*y**2 *x2F1 +
    ds2 = y*(1-y*0.5)*xF3
    couple = g_Z2_V_H*g_Z2_V_D
    #couple = e*g_X_c*c_B**2 *c_W *t_X/c_X
    #anti =1
    ds = (Vert*couple)**2  *mp*E/(np.pi*(M_Z2**2+Q2)**2) *(ds1+anti*ds2) 
    return ds
#---------------------

#integration and parsing of cross section
def Cross_v(Vert,mN,E,anti):
    x_max =1
    y_max = 1
    
    x_min = mN**2/(2*mp*E)
    print('min',x_min)
    #y_min = mN**4/(8*mp*x*E**3)
    
    sig = sint.dblquad(ds_v,x_min,x_max,lambda x: mN**4/(8*mp*x*E**3),y_max, args = (np.array([Vert,E,anti])))
    
    return sig #Gev^-2


n= 20 #nr of energy evaluations
E_nus = np.linspace(100,400,n)
sigs_E = np.zeros((4,3,n))
err_sigs_E = np.zeros((4,3,n))

#N_masses[j,k+3] are in MeV
for j in range(4):#benchmarks
    for k in range(3): #HNL 4,5,6
        
        for i in range(n):
            sigs_E[j,k,i], err_sigs_E[j,k,i] = Cross_v(Vert[j,k],N_masses[j,k+3]*10**(-3),E_nus[i],1)



    

hf = h5.File('nu_sig_DC_L.hdf5', 'w')
hf.create_dataset('E_nus',data=E_nus)
hf.create_dataset('sigs_E',data=sigs_E)
hf.create_dataset('err_sigs_E',data=err_sigs_E)
#hf.create_dataset('sigs_E_bar',data=sigs_E_bar)
#hf.create_dataset('err_sigs_E_bar',data=err_sigs_E_bar)
hf.close()


quit()        

## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

#! /usr/bin/env python

## Python LHAPDF6 usage example

import lhapdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py as h5

lhapdf.setPaths(['/home/rverkade/.local/share/lhapdf'])
print(lhapdf.paths())
## Getting a PDF member object
#p = lhapdf.mkPDF("CT10nlo", 0)
#p = lhapdf.mkPDF("CT10nlo/0")
p = lhapdf.mkPDF("CT18NLO/0")
pw = lhapdf.mkPDF("CT10wnlo/0")

s2_W =0.22290 #sin_W^2
## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g

W = np.arcsin(np.sqrt(0.22290))
M_z = 91*10**3 #MeV
M_z_prime = 1.25*10**3 
M_rat_2 =( M_z_prime/M_z)**2
#M_rat = 0
v_f = 500 #MeV
v_H = 265*10**3 #MeV
g_g = M_z*2/v_H #=sqrt(g'^2+g^2)  
g_X_c = 0.32*4*np.pi

xs_x = [x for x in np.logspace(-4, 0, 80)]
qs_Q = [q for q in np.logspace(0, 4, 80)]
xs_Q= np.array([10**-8,10**-6, 0.0001,0.01,0.1,0.2,0.5,0.8])
qs_x = np.array([10,50,100,200,500,1000,2000,5000])
qs_x2 = np.array([0.05,0.1,1,5,10,50,100,200,500])
xs_Q2= np.array([10**-10,10**-9,10**-8,10**-6, 0.0001,0.01,0.1])
gluon_xfs_x = np.zeros([len(qs_x), len(xs_x)])
gluon_xfs_Q = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
l=len(xs_Q)
print(l)
k = len(xs_x)
close = True
colours = np.array(['r','b','g','orange','purple','pink','black','cyan'])
#e_q = np.array([-1/3,2/3,-1/3,2/3,-1/3])
e_q = np.array([1/3,2/3,1/3,2/3])


#----------------
############
# Time for #
# F1 F2 FU #
############


G_F = 1.188*10**(-5) #GeV^-2
M_W =  80.38 #GeV
mp = 0.94 #GeV
g_A = np.array([-0.5, 0.5,-0.5,0.5])
g_V = g_A -2*e_q*s2_W
g_X_c = 0.32*4*np.pi
X_c = 0.0003 #10^-4
W = np.arcsin(np.sqrt(0.22290))


def CCF_2_3(Q_e,Q_L,Q_X,A):#ymu Z_mu for now only Left handed input
    #mu = 2*g_X*v_f/(g_g*v_H)
    X = X_c
    g_X = g_X_c
    s_W = np.sin(W)
    c_W = np.cos(W)
    e = c_W/g_g *s_W #g*s_W
    #print('e',e)
    s_X = np.sin(X)
    c_X = np.cos(X)
    t_X = s_X/c_X

    #t_B =t_X*s_W #meaning NC Z'=0
    #M_Z' = (1-swtx/tb) *M_Z_SM
    t_B =t_X*s_W/(1-M_rat_2) #meaning NC Z'=0
    
    c_B = 1/(np.sqrt(1+t_B**2))
    s_B = t_B*c_B
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

    return(T_1 - T_2,T_3 -T_4)




def struct_PDF(x,Q2):
    flavs = p.flavors()[1:-2]#d,u,s,c no excluding b and gluon
    f_len = len(flavs)
    f_half = int(f_len*0.5)
    xf = np.zeros(f_len)
    xf = np.array(p.xfxQ2(flavs, x, Q2)) #notice func of Q^2    quarks           anti-quarks
    #xf_g = np.array(p.xfxQ2(21, x, Q2))
    #F2 = np.sum(xf[f_half:]+xf[:f_half]) #sea qu
    g_Z1_V,g_Z2_V =  CCF_2_3(e_q,g_A,np.zeros(f_half),A=False)
    g_Z1_A,g_Z2_A =  CCF_2_3(e_q,g_A,np.zeros(f_half),A=True)
    F2_y = np.sum(e_q**2 *(xf[f_half:]+xf[:f_half])) #sea quarks+valence, ?c dominated??
    F2_Z2 = np.sum((g_Z2_V**2+g_Z2_A**2) *(xf[f_half:]+xf[:f_half]))
    F2_yZ2 = np.sum(2*e_q*g_Z2_V *(xf[f_half:]+xf[:f_half]))

    xF3_Z2 = np.sum(2*g_Z2_V*g_Z2_A *(xf[f_half:]-xf[:f_half]))
    xF3_yZ2 = np.sum(2*e_q*g_Z2_A *(xf[f_half:]-xf[:f_half]))
    xF3_y = 0  #valence quarks(sea cancels)#times 2 for neutrino
    FL=0 #WHERE am is supposeed to find disnshit
    x2F1_y = F2_y*(1+4*mp**2 * x**2 /Q2)/(1+FL)
    return (F2_y,F2_Z2,F2_yZ2,xF3_y,xF3_Z2,xF3_yZ2)   # (x2F1,

F1_x_y = np.zeros([len(qs_x), len(xs_x)])
F2_x_y = np.zeros([len(qs_x), len(xs_x)])
F3_x_y = np.zeros([len(qs_x), len(xs_x)])
F1_Q_y = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F2_Q_y = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F3_Q_y = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)

F1_x_yZ2 = np.zeros([len(qs_x), len(xs_x)])
F2_x_yZ2 = np.zeros([len(qs_x), len(xs_x)])
F3_x_yZ2 = np.zeros([len(qs_x), len(xs_x)])
F1_Q_yZ2 = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F2_Q_yZ2 = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F3_Q_yZ2 = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)

F1_x_Z2 = np.zeros([len(qs_x), len(xs_x)])
F2_x_Z2 = np.zeros([len(qs_x), len(xs_x)])
F3_x_Z2 = np.zeros([len(qs_x), len(xs_x)])
F1_Q_Z2 = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F2_Q_Z2 = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)
F3_Q_Z2 = np.zeros([len(xs_Q), len(qs_Q)]) #(8,80)

for i in range(l): #l 0-8
    for j in range(k): #0-80   p.xfxQ2 takes Q^2 as arg
        F2_Q_y[i,j], F2_Q_yZ2[i,j],F2_Q_Z2[i,j],F3_Q_y[i,j],F3_Q_yZ2[i,j],F3_Q_Z2[i,j] = struct_PDF(xs_Q[i], qs_Q[j])#taking Q^2 as argument
        F2_x_y[i,j] ,F2_x_yZ2[i,j] ,F2_x_Z2[i,j] ,F3_x_y[i,j],F3_x_yZ2[i,j],F3_x_Z2[i,j]  = struct_PDF(xs_x[j], qs_x[i])#taking Q^2 as argument
        


hf = h5.File('Struct_func_dark_ZZ.hdf5', 'w')
hf.create_dataset('xs_x',data=xs_x)
hf.create_dataset('qs_x',data=qs_x)
hf.create_dataset('xs_Q',data=xs_Q)
hf.create_dataset('qs_Q',data=qs_Q)

hf.create_dataset('F2_x_y', data=F2_x_y)
hf.create_dataset('F2_x_yZ2', data=F2_x_yZ2)
hf.create_dataset('F2_x_Z2', data=F2_x_Z2)

hf.create_dataset('F3_x_y', data=F3_x_y)
hf.create_dataset('F3_x_yZ2', data=F3_x_yZ2)
hf.create_dataset('F3_x_Z2', data=F3_x_Z2)

hf.create_dataset('F2_Q_y', data=F2_Q_y)
hf.create_dataset('F2_Q_yZ2', data=F2_Q_yZ2)
hf.create_dataset('F2_Q_Z2', data=F2_Q_Z2)

hf.create_dataset('F3_Q_y', data=F3_Q_y)
hf.create_dataset('F3_Q_yZ2', data=F3_Q_yZ2)
hf.create_dataset('F3_Q_Z2', data=F3_Q_Z2)

hf.close()
quit()



## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

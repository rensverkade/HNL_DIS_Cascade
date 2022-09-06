import numpy as np
import matplotlib.pyplot as plt
import math
import h5py as h5
import sys
from numpy import linalg as la
#hf = h5.File('HNL_mass_pascoli.hdf5',r)
#mu = np.array(hf.get('mu'))
#hf.close()
e_m = 0.001

e=1
W = np.arcsin(np.sqrt(0.22290))
M_z = 91*10**3 #MeV
M_z_prime = 1.25*10**3 
M_rat_2 =( M_z_prime/M_z)**2
#M_rat = 0
v_f = 500 #MeV
v_H = 265*10**3 #MeV
g_g = M_z*2/v_H #=sqrt(g'^2+g^2)  
#mu = 2*g_x*v_f/(g_g*v_H)

#X is free param
def CCF_1(e,Q_e): #ymu A_mu
    return -e*Q_e

#Q_x = 1 for nu_D 0 for rest?
#Q_e el charge, Q_L isospin
def CCF_2_3(W,X,g_X,Q_e,Q_L,Q_X):#ymu Z_mu for now only Left handed input
    #mu = 2*g_X*v_f/(g_g*v_H)
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
                  
    T_1 = e/(2*s_W*c_W) * ( -Q_L*(c_B+s_B*s_W*t_X)+2*Q_e*s_W*(-s_B*t_X + c_B*s_W))
    T_2 = g_X*Q_X*s_B/c_X #Z

    T_3 = e/(2*s_W*c_W) * ( Q_L*(s_B-c_B*s_W*t_X)*c_B+2*Q_e*s_W*(c_B*t_X - s_B*s_W) )
    T_4 = g_X*Q_X*c_B/c_X #Z'

    return(T_1 - T_2,T_3 -T_4)


def NCF_2_3(W,X,g_X,Q_e,Q_L,Q_X): #still need to change this shit
    #mu = 2*g_X*v_f/(g_g*v_H)
    s_W = np.sin(W)
    c_W = np.cos(W)
    e = c_W/g_g *s_W #g*s_W
    #print('e',e)
    s_X = np.sin(X)
    c_X = np.cos(X)
    t_X = s_X/c_X

    #M_Z' = 1-swtx/tb *M_Z_SM
    t_B =t_X*s_W/(1-M_rat_2) #meaning NC Z'=0
       
    c_B = 1/(np.sqrt(1+t_B**2))
    s_B = t_B*c_B

 
    if Q_X == 0:
        T1 = -g_g/(c_W*2) *(c_B +s_B*s_W*t_X)
        T2 = -g_g/(c_W*2) *(c_B -s_B*s_W*t_X)


    else:
        T1 = -s_B/c_X
        T2 = -c_B/c_X

    return T1,T2 #Z, Z'
X_max = 0.0003
#previously 0.03
print(X_max)
Xs = np.linspace(0,X_max,30)
g_Xs =np.linspace(0.03,1.5,40)
#g_X_c = 0.1

g_X_c = 0.32*4*np.pi
#CC                     W, X, g_x, Q_e,Q_L,Q_x
E_Z1_X,E_Z2_X = np.abs( CCF_2_3(W,Xs,g_X_c,1,-0.5,0))
E_Z1_g,E_Z2_g = np.abs( CCF_2_3(W,X_max,g_Xs,1,-0.5,0))
vE_Z1_X,vE_Z2_X =  np.abs(CCF_2_3(W,Xs,g_X_c,0,-0.5,0))
vE_Z1_g,vE_Z2_g =  np.abs(CCF_2_3(W,X_max,g_Xs,0,-0.5,0))
vD_Z1_X,vD_Z2_X =  np.abs(CCF_2_3(W,Xs,g_X_c,0,0,1))
vD_Z1_g,vD_Z2_g =  np.abs(CCF_2_3(W,X_max,g_Xs,0,0,1))
#NC

vE_Z1_X_NC,vE_Z2_X_NC = np.abs( NCF_2_3(W,Xs,g_X_c,0,0.5,0))
vE_Z1_g_NC,vE_Z2_g_NC =  np.abs(NCF_2_3(W,X_max,g_Xs,0,0.5,0))
vD_Z1_X_NC,vD_Z2_X_NC =  np.abs(NCF_2_3(W,Xs,g_X_c,0,0,1))
vD_Z1_g_NC,vD_Z2_g_NC =  np.abs(NCF_2_3(W,X_max,g_Xs,0,0,1))

#print('e_Z_g',E_Z1_g)
#print('e_Z2_g',E_Z2_g)
#print('ve_Z_g',vE_Z1_g)
#print('ve_Z2_g',vE_Z2_g)
#print('vD_Z_g',vD_Z1_g)
#print('vD_Z2_g',vD_Z2_g)

closs =True
#box1 = dict( facecolor='none', edgecolor = 'black')
half = int(len(Xs)/2)
#folder = 'figs'
#folder = 'SSB_figs'
folder = 'SSB_figs_gX_lim'
#folder = 'figs_gX'

plt.title(r'$\mathcal{L}_{int} \, Z$')
plt.plot(Xs, E_Z1_X, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$',color='k')
plt.text(0, E_Z1_X[half]-0.075, s=r'$\mathcal{L}^e_{int}\simeq$ '+str(round(E_Z1_X[half],2)), fontsize = 13, bbox=dict( facecolor='none',edgecolor = 'black') )
plt.plot(Xs, vE_Z1_X, label= r'$\nu_e$ $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
plt.text(0, vE_Z1_X[half]-0.075, r'$\mathcal{L}^{\nu_e}_{int}\simeq$ '+str(round(vE_Z1_X[half],2)), fontsize = 13, bbox = dict(facecolor='none',edgecolor = 'black',ls='dashed'))
plt.plot(Xs, vD_Z1_X, label= r'$\nu_D$ $Q_e, \, Q_l, \, Q_X =0, \, 0, \, 1$', ls='dotted',color='k')
plt.text(0, vD_Z1_X[half]+0.075, r'$\mathcal{L}^{\nu_D}_{int}\simeq$ '+str(round(vD_Z1_X[half],2)), fontsize = 13, bbox = dict(facecolor='none', edgecolor = 'black', ls='dotted'))
plt.legend()
plt.xlabel(r'$\chi$')
plt.ylabel(r'$\left|L_{Int}\right|$')
plt.savefig(folder+'/L_X_Z1.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()

plt.title(r'$\mathcal{L}_{int} \, Z^{\prime}$')
plt.plot(Xs, E_Z2_X, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$', ls='solid',color='k')
plt.text(0, E_Z2_X[half]+0.3, s=r'$\mathcal{L}^e_{int}\simeq$ '+str(round(E_Z2_X[half],2)), fontsize = 13, bbox=dict( facecolor='none',edgecolor = 'black') )
plt.plot(Xs, vE_Z2_X, label= r'$\nu_e$ $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
plt.text(0.025, vE_Z2_X[half]+0.3, r'$\mathcal{L}^{\nu_e}_{int}\simeq$ '+str(round(vE_Z2_X[half],2)), fontsize = 13, bbox = dict(facecolor='none',edgecolor = 'black',ls='dashed'))
plt.plot(Xs, vD_Z2_X, label= r'$\nu_D$ $Q_e, \, Q_l, \, Q_X =0, \, 0, \, 1$', ls='dotted',color='k')
plt.text(0, vD_Z2_X[half]-0.4, r'$\mathcal{L}^{\nu_D}_{int}\simeq$ '+str(round(vD_Z2_X[half],2)), fontsize = 13, bbox = dict(facecolor='none', edgecolor = 'black', ls='dotted'))
plt.legend()
plt.xlabel(r'$\chi$')
plt.ylabel(r'$\left|L_{Int}\right|$')
plt.savefig(folder+'/L_X_Z2.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()


plt.title(r'$\mathcal{L}_{int} \, Z$')
plt.plot(g_Xs, E_Z1_g,color='k', label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$')
plt.text(0, E_Z1_g[half]-0.1, s=r'$\mathcal{L}^e_{int}\simeq$ '+str(round(E_Z1_g[half],2)), fontsize = 13, bbox=dict( facecolor='none',edgecolor = 'black') )
plt.plot(g_Xs, vE_Z1_g,color='k', label= r' $\nu_e \, Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed')
plt.text(0, vE_Z1_g[half]-0.1, r'$\mathcal{L}^{\nu_e}_{int}\simeq$ '+str(round(vE_Z1_g[half],2)), fontsize = 13, bbox = dict(facecolor='none',edgecolor = 'black',ls='dashed'))
plt.plot(g_Xs, vD_Z1_g,color='k', label= r' $\nu_D \, Q_e, \, Q_l, \, Q_X =0, \, 0, \, 1$', ls='dotted')
plt.text(0, vD_Z1_g[half]+0.1, r'$\mathcal{L}^{\nu_D}_{int}\simeq$ '+str(round(vD_Z1_g[half],2)), fontsize = 13, bbox = dict(facecolor='none', edgecolor = 'black', ls='dotted'))
plt.xlabel(r'$g_X$')
plt.ylabel(r'$\left|L_{Int}\right|$')
#plt.yscale('log')
plt.legend()
plt.savefig(folder+'/L_g_Z1.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()

    

plt.title(r'$\mathcal{L}_{int} \, Z^{\prime}$')
plt.plot(g_Xs, E_Z2_g,color='k', label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$')
plt.text(g_Xs[half], E_Z2_g[half]+0.2, s=r'$\mathcal{L}^e_{int}\simeq$ '+str(round(E_Z2_g[half],2)), fontsize = 13, bbox=dict( facecolor='none',edgecolor = 'black') )
plt.plot(g_Xs, vE_Z2_g,color='k', label= r'$\nu_e \, Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed')
plt.text(g_Xs[half]+0.4, E_Z2_g[half]+0.2, r'$\mathcal{L}^{\nu_e}_{int}\simeq$ '+str(round(vE_Z2_g[half],2)), fontsize = 13, bbox = dict(facecolor='none',edgecolor = 'black',ls='dashed'))
plt.plot(g_Xs, vD_Z2_g,color='k', label= r'$ \nu_D \, Q_e, \, Q_l, \, Q_X =0, \,0, \, 1$', ls='dotted')
plt.text(0, vD_Z2_g[half]+0.1, r'$\mathcal{L}^{\nu_D}_{int}\simeq$ '+str(round(vD_Z2_g[half],2)), fontsize = 13, bbox = dict(facecolor='none', edgecolor = 'black', ls='dotted'))
plt.xlabel(r'$g_X$')
plt.ylabel(r'$\left|L_{Int}\right|$')
#plt.yscale('log')
plt.legend()
plt.savefig(folder+'/L_g_Z2.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()


    
quit()
#NC AS FUNC OF X







plt.plot(Xs, vE_Z1_X_NC, label= r'e $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
plt.plot(Xs, vD_Z1_X_NC, label= r'e $Q_e, \, Q_l, \, Q_X =0, \, 0, \, 1$', ls='dotted',color='k')
plt.legend()
plt.xlabel(r'$\chi$')
plt.ylabel(r'$\left|L_{Int}\right|$')
plt.savefig('figs/L_X_Z1_NC.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()

plt.plot(Xs, vE_Z2_X_NC, label= r'$\nu_e$ $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
plt.plot(Xs, vD_Z2_X_NC, label= r'$\nu_D$ $Q_e, \, Q_l, \, Q_X =0, \, 0, \, 1$', ls='dotted',color='k')
plt.legend()

plt.xlabel(r'$\chi$')
plt.ylabel(r'$L_{Int}$')

plt.savefig('figs/L_X_Z2_NC.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()

plt.plot(g_Xs, vE_Z1_g_NC, label= r'e $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
#plt.plot(g_Xs, vD_Z1_g_NC, label= r'e $Q_e, \, Q_l, \, Q_X =0, \, 0, \, 1$', ls='dotted',color='k')
plt.legend()
plt.xlabel(r'$g_X$')
plt.ylabel(r'$L_{Int}$')
plt.savefig('figs/L_g_Z1_NC.pdf',dpi=800)
plt.show()

plt.plot(g_Xs, vE_Z2_g_NC, label= r'$\nu_e$ $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
#plt.plot(g_Xs, vD_Z2_g_NC, label= r'$\nu_D$ $Q_e, \, Q_l, \, Q_X =0, \, 0, \, 1$', ls='dotted',color='k')
plt.legend()

plt.xlabel(r'g_X')
plt.ylabel(r'$L_{Int}$')

plt.savefig('figs/L_g_Z2_NC.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()





quit()









plt.plot(g_Xs, E_Z1_g, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$',color='k')
#plt.plot(g_Xs, vE_Z1_g, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$', ls='dashed')
#plt.plot(g_Xs, vD_Z1_g, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$', ls='dotted')
plt.xlabel(r'$g_X$')
plt.ylabel(r'$L_{Int}$')
#plt.yscale('log')
plt.legend()
plt.savefig('figs/L_ge_Z1.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()

#plt.plot(g_Xs, E_Z1_g, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$',color='k')
plt.plot(g_Xs, vE_Z1_g, label= r'$\nu_e$ $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$', ls='dashed',color='k')
#plt.plot(g_Xs, vD_Z1_g, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$', ls='dotted')
plt.xlabel(r'$g_X$')
plt.ylabel(r'$L_{Int}$')
#plt.yscale('log')
plt.legend()
plt.savefig('figs/L_gv_Z1.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()


#plt.plot(g_Xs, E_Z1_g, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$')
#plt.plot(g_Xs, vE_Z1_g, label= r'e $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
plt.plot(g_Xs, vD_Z1_g, label= r'$\nu_D$ $Q_e, \, Q_l, \, Q_X =0, \, 0, \,1$', ls='dotted',color='k')
plt.xlabel(r'$g_X$')
plt.ylabel(r'$L_{Int}$')
plt.legend()
plt.savefig('figs/L_gd_Z1.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()


plt.plot(g_Xs, E_Z2_g, label= r'e $Q_e, \, Q_l, \, Q_X =1, \, -1/2, \, 0$')
plt.plot(g_Xs, vE_Z2_g, label= r'$\nu_e$ $Q_e, \, Q_l, \, Q_X =0, \, -1/2, \, 0$', ls='dashed',color='k')
plt.plot(g_Xs, vD_Z2_g, label= r'$\nu_D$ $Q_e, \, Q_l, \, Q_X =0, \, 0, \,1$', ls='dotted',color='k')
plt.xlabel(r'$g_X$')
plt.ylabel(r'$L_{Int}$')
plt.legend()
plt.savefig('figs/L_g_Z2.pdf',dpi=800)
if closs == True:
    plt.close()
else:
    plt.show()










quit()
#Tau = e =1-, 



print(masses)
hf = h5.File('HNL_mass_pascoli.hdf5', 'w')
hf.create_dataset('m_D', data=m_D)
hf.create_dataset('L_L', data=L_L)
hf.create_dataset('L_R', data=L_R)
hf.create_dataset('mu', data=m_diag)
hf.create_dataset('M_X', data=M_X)
hf.create_dataset('M_diag', data=masses)
hf.close()






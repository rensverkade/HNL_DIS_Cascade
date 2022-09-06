#! /usr/bin/env python

## Python LHAPDF6 usage example

import lhapdf
import apfel
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py as h5
from scipy import integrate as sint

lhapdf.setPaths(['/home/rverkade/.local/share/lhapdf'])
print(lhapdf.paths())
## Getting a PDF member object
#p = lhapdf.mkPDF("CT10nlo", 0)
#p = lhapdf.mkPDF("CT10nlo/0")
p = lhapdf.mkPDF("CT18NLO/0")
#pw = lhapdf.mkPDF("CT10wnlo/0")

#apfel.SetPDFSet("CT18NLO")
# Initializes integrals on the grids
apfel.InitializeAPFEL();



quit()
s2_W =0.22290 #sin_W^2
## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g



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
mN = 0.105 #GeV
g_A = np.array([-0.5, 0.5,-0.5,0.5])
g_V = g_A -2*e_q*s2_W

def struct_PDF_Wm(x,Q2):
    flavs = p.flavors()[1:-2]#d,u,s,c no excluding b and gluon
    f_len = len(flavs)
    f_half = int(f_len*0.5)
    xf = np.zeros(f_len)
    xf = np.array(p.xfxQ2(flavs, x, Q2)) #notice func of Q^2    quarks           ant
    #                   u,        c,         d bar, s bar
    F2_W_m = 2*np.sum((xf[f_half+1]+xf[f_half+3]+xf[f_half-1]+xf[f_half-3]))
    xF3_W_m = 2*np.sum((xf[f_half]+xf[f_half+3]-xf[f_half-2]-xf[f_half-3]))
    FL=0 #WHERE am is supposeed to find disnshit
    #x2F1_W_m = F2_W*(1+4*mp**2 * x**2 /Q2)/(1+FL)
    x2F1_W_m = F2_W_m #in parton model
    return (x2F1_W_m,F2_W_m,xF3_W_m) #F2_y,xF3_y  # (x2F1,

def struct_PDF_Wp(x,Q2):
    flavs = p.flavors()[1:-2]#d,u,s,c no excluding b and gluon
    f_len = len(flavs)
    f_half = int(f_len*0.5)
    xf = np.zeros(f_len)
    xf = np.array(p.xfxQ2(flavs, x, Q2)) 
    F2_W_p = 2*np.sum((xf[f_half]+xf[f_half+2]+xf[f_half-2]+xf[f_half-4]))
    xF3_W_p = 2*np.sum((xf[f_half]+xf[f_half+2]-xf[f_half-2]-xf[f_half-4]))
    FL=0 #WHERE am is supposeed to find disnshit
    #x2F1_W = F2_W*(1+4*mp**2 * x**2 /Q2)/(1+FL)
    x2F1_W_p = F2_W_p #in parton model
    return (x2F1_W_p,F2_W_p,xF3_W_p) #F2_y,xF3_y  # (x2F1,

G_F = 1.188*10**(-5) #GeV^-2
mp = 0.94 #GeV
M_W =  80.38 #GeV
#M= mp
#def lower_x(x):
def ds_t(x,y,E): #M is nucleon mass
    Q2 = x*y*2*mp*E
    x2F1, F2,xF3 = struct_PDF(x,Q2)
    #ds1 = 0.5*y**2 *x2F1 +(1-y-Q2/4*E)*F2 # +0.5*y**2 *x2F1 +
    #ds2 = y*(1-y*0.5)*xF3
    anti =1
    ds = G_F**2 *mp*E*( 1/(np.pi*(1+Q2/M_W**2)**2)) *((1-y-x*y*mp/(2*E))*F2 +0.5*y**2 *x2F1+y*(1-y*0.5)*xF3)
    
    #ds = G_F**2 *mp*E*((1-y-Q2/4*E)*F2 +0.5*y**2 *x2F1+y*(1-y*0.5)*xF3)
    #ds =1/((1+Q2/M_W**2)**2)
    return ds
        
def ds_v(x,y,E,anti): #M is nucleon mass
    Q2 = x*y*2*mp*E
    if anti ==1:
        x2F1, F2,xF3 = struct_PDF_Wm(x,Q2)
    elif anti== -1:
        x2F1, F2,xF3 = struct_PDF_Wp(x,Q2)
        
    ds1 = 0.5*y**2 *x2F1 +(1-y-x*y*mp/(2*E))*F2 # +0.5*y**2 *x2F1 +
    ds2 = y*(1-y*0.5)*xF3

    #anti =1
    ds = G_F**2 *mp*E/(np.pi*(1+Q2/M_W**2)**2) *(ds1+anti*ds2)
    return ds
def Cross_v(E,anti):
    x_max =1
    y_max = 1
    
    x_min = mN**2/(2*mp*E)
    #y_min = mN**4/(8*mp*x*E**3)
    
    sig = sint.dblquad(ds_v,x_min,x_max,lambda x: mN**4/(8*mp*x*E**3),y_max, args = (np.array([E,anti])))
    #sig = sint.dblquad(ds_v,0.0001,1,0.00001,1, args = (np.array([E,anti])))
    return sig #Gev^-2
n= 20
E_nus = np.linspace(100,400,n)
sigs_E = np.zeros(n)
err_sigs_E = np.zeros(n)
sigs_E_bar = np.zeros(n)
err_sigs_E_bar = np.zeros(n)
sigs_a = np.zeros(n)
for i in range(n):
    sigs_E[i], err_sigs_E[i] = Cross_v(E_nus[i],1)
    sigs_E_bar[i], err_sigs_E_bar[i] = Cross_v(E_nus[i],-1)
    #sigs_E[i], err_sigs_E[i] = sint.dblquad(ds_v,0.0001,1,0.00001,1, args = (np.array([E_nus[i],1])))
    
    #a_sigs_E[i], a_err_sigs_E[i] = Cross_v(0.001,1,0.001,-1,E_nus[i])
    print(i)
plt.plot(E_nus,sigs_E*3.89379372/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, \mu^- +\,X$')
plt.plot(E_nus,sigs_E_bar*3.89379372/E_nus,label=r'$\overline{\nu}\, +\, p\,\rightarrow\, \mu^+ +\,X$')
#plt.plot(a_E_nus,a_sigs_a)
plt.xlabel('E (GeV)')
plt.ylabel(r'$\sigma / E 10^{-28} cm^2$')
plt.legend()
plt.savefig('cross.pdf',dpi=800)
plt.show()

hf = h5.File('nu_sig_CC.hdf5', 'w')
hf.create_dataset('E_nus',data=E_nus)
hf.create_dataset('sigs_E',data=sigs_E)
hf.create_dataset('err_sigs_E',data=err_sigs_E)
hf.create_dataset('sigs_E_bar',data=sigs_E_bar)
hf.create_dataset('err_sigs_E_bar',data=err_sigs_E_bar)
hf.close()


quit()        

## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

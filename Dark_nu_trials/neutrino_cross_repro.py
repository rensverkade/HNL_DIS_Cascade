#! /usr/bin/env python

## Python LHAPDF6 usage example

import lhapdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate as sint


G_F = 1.188*10**(-5) #GeV^-2
M_W =  80.38 #GeV
mp = 0.94 #GeV

flav =np.append( np.arange(-5,0,1,dtype=int), np.arange(1,6,1,dtype=int))
e_q = np.array([-1/2,2/3,-1/3,2/3,-1/3])
e_q2 = np.append(np.flip(e_q),e_q)
print(flav)
print(e_q)
print(e_q2)
lhapdf.setPaths(['/home/rverkade/.local/share/lhapdf'])
print(lhapdf.paths())
## Getting a PDF member object
#p = lhapdf.mkPDF("CT10nlo", 0)
#p = lhapdf.mkPDF("CT18NNLO/0")
p = lhapdf.mkPDF("CT10nlo/0")
print(p.flavors())
#quit()
xs_Q= np.array([10**-8,10**-6, 0.0001,0.01,0.1,0.2,0.5,0.8])
## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g'

def struct_PDF(x,Q2,flavs,M):
    f_len = len(flavs)
    f_half = int(f_len*0.5)
    xf = np.zeros(f_len)
    xf = np.array(p.xfxQ2(flavs, x, Q2)) #notice func of Q^2    quarks           anti-quarks
    xf_g = np.array(p.xfxQ2(21, x, Q2))
    #print('xf',xf)
    F2 = np.zeros(f_len)
    F2 = 2*np.sum(xf[f_half:]+xf[:f_half]) #sea quarks+valence, ?c dominated??
    xF3 = 2*np.sum(xf[f_half:]-xf[:f_half]) #valence quarks(sea cancels)
    FL=0 #WHERE am is supposeed to find disnshit
    x2F1 = F2*(1+4*M**2 * x**2 /Q2)/(1+FL)

    return (x2F1,F2,xF3,xf)

xs_x = [x for x in np.logspace(-4, 0, 100)]

x2F1 = np.zeros(len(xs_x))
F2 = np.zeros(len(xs_x))
xF3 = np.zeros(len(xs_x))
x2F1_2 = np.zeros(len(xs_x))
F2_2 = np.zeros(len(xs_x))
xF3_2 = np.zeros(len(xs_x))
xf_s = np.zeros((len(flav),len(xs_x)))
xf_s2 = np.zeros((len(flav),len(xs_x)))
for i in range(len(xs_x)):
    x2F1[i],F2[i],xF3[i],xf_s[:,i] = struct_PDF(xs_x[i],90,flav,mp)
    x2F1_2[i],F2_2[i],xF3_2[i],xf_s2[:,i] = struct_PDF(xs_x[i],10,flav,mp)

#plt.plot(xs_x,F2, label='90')
#plt.plot(xs_x,F2_2, label='3.5')
plt.plot(xs_x,xF3, label='90')
plt.plot(xs_x,xF3_2, label='3.5')
plt.xlabel(r'$x$')
plt.ylabel(r'$F2(x,Q)$')
plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.savefig('test_F3.pdf',dpi=800)
plt.close()
plt.plot(xs_x,xf_s[5,:], label=r'$90 GeV^2$, fl'+str(flav[5]))
#plt.plot(xs_x,xf_s2[5,:], label=r'$10 GeV^2$, fl'+str(flav[5]))
plt.plot(xs_x,xf_s[6,:], label=r'$90 GeV^2$, fl'+str(flav[6]))
#plt.plot(xs_x,xf_s2[6,:], label=r'$10 GeV^2$, fl'+str(flav[6]))
plt.xlabel(r'$x$')
plt.ylabel(r'$xf(x,Q)$')
plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.savefig('test_fx_90.pdf',dpi=800)
plt.close()
#quit()



#                     -1 for anti
def ds_v(x,y,M,E,flavs,anti): #M is nucleon mass
    Q2 = x*y*2*M*E
    x2F1,F2,xF3,fx = struct_PDF(x,Q2,flavs,M)
    ds1 = 0.5*y**2 *x2F1 +(1-y-Q2/4*E)*F2 # +0.5*y**2 *x2F1 +
    ds2 = y*(1-y*0.5)*xF3
    #print('ds',ds1,ds2)
    ds = G_F**2 *M/(np.pi*(1+Q2/M_W**2)**2) *(ds1+anti*ds2)
    #print('ds',ds)
    return ds
#             0      1    0,   1
def Cross_v(x_min,x_max,y_min,y_max,M,E,flavs,anti):
    sig = sint.dblquad(ds_v,x_min,x_max, y_min,y_max, args = (M,E,flavs,anti))
    return sig



#quit()
#print('cross', Cross_v(0.001,1,0.001,1,mp,200,flav,1))

#quit()
p = lhapdf.mkPDF("CT10wnlo", 0)

flav =np.append( np.arange(-5,0,1,dtype=int), np.arange(1,6,1,dtype=int))
E_nus = np.linspace(100,400,30)
sigs_E = np.zeros(30)
err_sigs_E = np.zeros(30)
sigs_a = np.zeros(30)
for i in range(30):
    #sigs_E[i], err_sigs_E[i] = Cross_v(0.001,1,0.001,1,mp,E_nus[i],flav,1)
    sigs_E[i] = ds_v(0.001,0.01,1,E_nus[i],flav,1)
    #sigs_a[i] = 5.53*E_nus[i]**(0.363-1) *10**2


print(sigs_E)
plt.plot(E_nus,sigs_E)
#plt.plot(E_nus,sigs_a)
plt.xlabel('E')
plt.ylabel('sig')
plt.savefig('cross.pdf',dpi=800)
plt.show()



    

    


#plt.close()


## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

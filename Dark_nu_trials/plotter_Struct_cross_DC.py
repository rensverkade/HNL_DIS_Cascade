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
#p = lhapdf.mkPDF("CT10nlo", 0)
#p = lhapdf.mkPDF("CT10nlo/0")
p = lhapdf.mkPDF("CT18NLO/0")
#pw = lhapdf.mkPDF("CT10wnlo/0")

s2_W =0.22290 #sin_W^2
## Basic all-flavour PDF querying at x=0.01, Q=M_Z
#pid particle flavours
# - is anti-particle I asssume
# 1, 2, 3, 4, 5, 6, 7, 8, 21
# d, u, s, c, b, t, b', t', g


hf = h5.File('nu_sig_DC_L.hdf5', 'r')
E_nus = np.array(hf.get('E_nus'))
sigs_E = np.array(hf.get('sigs_E'))
err_sigs_E = np.array(hf.get('err_sigs_E'))
#sigs_E_bar = np.array(hf.get('sigs_E_bar'))
#err_sigs_E_bar = np.array(hf.get('err_sigs_E_bar'))
hf.close()

hf = h5.File('nu_sig_CC_L.hdf5', 'r')
E_nus_CC = np.array(hf.get('E_nus'))
sigs_E_CC = np.array(hf.get('sigs_E'))
sigs_E_CC_bar = np.array(hf.get('sigs_E_bar'))
hf.close()

hf = h5.File('nu_sig_NC_L.hdf5', 'r')
E_nus_NC = np.array(hf.get('E_nus'))
sigs_E_NC = np.array(hf.get('sigs_E'))
#sigs_E_CC_bar = np.array(hf.get('sigs_E_bar'))
hf.close()
#print(sigs_E_NC)


#----------------
############
# Time for #
# F1 F2 FU #
############
'''

#G_F = 1.188*10**(-5) #GeV^-2
#M_W =  80.38 #GeV
#mp = 0.94 #GeV
#mN = 0.105 #GeV
#g_A = np.array([-0.5, 0.5,-0.5,0.5])
#g_V = g_A -2*e_q*s2_W
plt.plot(E_nus,sigs_E[0,0,:]*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, N_4 +\,X$',color='k',ls='dotted')
plt.plot(E_nus,sigs_E[0,1,:]*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, N_5 +\,X$',color='r',ls='dotted')
plt.plot(E_nus,sigs_E[0,2,:]*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, N_6 +\,X$',color='b',ls='dotted')
plt.plot(E_nus,sigs_E_NC[0,0,:]*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, N_4 +\,X$',color='k')
plt.plot(E_nus,sigs_E_NC[0,1,:]*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, N_5 +\,X$',color='r')
plt.plot(E_nus,sigs_E_NC[0,2,:]*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, N_6 +\,X$',color='b')

#plt.plot(E_nus,sigs_E_CC*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, \mu^- +\,X$',color='k')
#plt.plot(E_nus,sigs_E_CC_bar*3.89379372*10**(10)/E_nus,label=r'$\overline{\nu}\, +\, p\,\rightarrow\, \mu^+ +\,X$',color='k',ls='dashed')
plt.xlabel('E (GeV)')
plt.ylabel(r'$\sigma^{\nu}\, /\, E_{\nu} \,\, \left(10^{-38}\, cm^2 \, GeV^{-1}\right)$')
plt.legend()
plt.savefig('figs_cross_L/cross_DC_NC.pdf',dpi=800,bbox_inches='tight')
plt.close()



quit()        
'''

close = True


fig, axs = plt.subplots(2, 2, sharex=True)
l = 3
#                    r       k
colours = np.array(['k','b','r','orange','purple','black','pink','cyan'])
for i in range(l): #l 0-8        
    axs[0,0].plot(E_nus,(sigs_E[0,i,:]+sigs_E_NC[0,0,:])*3.89379372*10**(12)/E_nus,label=r'$N$'+str(i+4)+r'$+X$',color=colours[i],ls='solid')
    axs[1,0].plot(E_nus,(sigs_E[0,i,:])*3.89379372*10**(12)/E_nus,color=colours[i],ls='dashed')
    axs[0,1].plot(E_nus,(sigs_E_NC[0,i,:])*3.89379372*10**(12)/E_nus,color=colours[i],ls='dotted')

    if  i==0:     
        axs[1,1].plot(E_nus,(sigs_E[0,i,:]+sigs_E_NC[0,0,:])*3.89379372*10**(12)/E_nus,label=r'$\sigma_Z +\sigma_{Z\prime}$',color=colours[i],ls='solid')
        axs[1,1].plot(E_nus,(sigs_E[0,i,:])*3.89379372*10**(12)/E_nus,label=r'$\sigma_{Z \prime}$',color=colours[i],ls='dashed')
        axs[1,1].plot(E_nus,(sigs_E_NC[0,i,:])*3.89379372*10**(12)/E_nus,label=r'$\sigma_Z$',color=colours[i],ls='dotted')
        
axs[1,0].set_xlabel('E (GeV)')

axs[1,1].set_xlabel('E (GeV)')
axs[0,0].set_ylabel(r'$\sigma^{\nu}\, /\, E_{\nu} \,\, \left(10^{-40}\, cm^2 \, GeV^{-1}\right)$')
axs[1,0].set_ylabel(r'$\sigma^{\nu}\, /\, E_{\nu} \,\, \left(10^{-40}\, cm^2 \, GeV^{-1}\right)$')
#plt.xscale('log')
#axs[0,0].set_yscale('log')
#axs[0,1].set_yscale('log')
#axs[1,0].set_yscale('log')
#axs[1,1].set_yscale('log')
#plt.ylim(0,1.5))
fig.tight_layout()
fig.legend(loc='center right')
plt.subplots_adjust(hspace=0.1,right=0.75)
plt.savefig('Cross_figs_Dark_SSM/cross_NC_DC_A_2.pdf', dpi=800,bbox_inches = 'tight')
if close == True:
    plt.close()
else:
    plt.show()






## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...

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


hf = h5.File('nu_sig_CC_L.hdf5', 'r')
E_nus = np.array(hf.get('E_nus'))
sigs_E = np.array(hf.get('sigs_E'))
err_sigs_E = np.array(hf.get('err_sigs_E'))
sigs_E_bar = np.array(hf.get('sigs_E_bar'))
err_sigs_E_bar = np.array(hf.get('err_sigs_E_bar'))
hf.close()
print(err_sigs_E_bar)


#----------------
############
# Time for #
# F1 F2 FU #
############


#G_F = 1.188*10**(-5) #GeV^-2
#M_W =  80.38 #GeV
#mp = 0.94 #GeV
#mN = 0.105 #GeV
#g_A = np.array([-0.5, 0.5,-0.5,0.5])
#g_V = g_A -2*e_q*s2_W

plt.plot(E_nus,sigs_E*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, \mu^- +\,X$',color='k')
plt.plot(E_nus,sigs_E_bar*3.89379372*10**(10)/E_nus,label=r'$\overline{\nu}\, +\, p\,\rightarrow\, \mu^+ +\,X$',color='k',ls='dashed')
#plt.plot(a_E_nus,a_sigs_a)
plt.xlabel('E (GeV)')
plt.ylabel(r'$\sigma^{\nu}_{CC}\, /\, E_{\nu} \,\, \left(10^{-38}\, cm^2 \, GeV^{-1}\right)$')
plt.legend()
plt.savefig('figs_cross_L/cross_SM_CC.pdf',dpi=800,bbox_inches='tight')
plt.close()



quit()        

plt.plot(E_nus,sigs_E*3.89379372*10**(10)/E_nus,label=r'$\nu\, +\, p\,\rightarrow\, \mu^- +\,X$',color='k')
plt.plot(E_nus,sigs_E_bar*3.89379372*10**(10)/E_nus,label=r'$\overline{\nu}\, +\, p\,\rightarrow\, \mu^+ +\,X$',color='k',ls='dashed')
plt.fill_between(E_nus,(sigs_E-err_sigs_E)*3.89379372*10**(10)/E_nus, (sigs_E+err_sigs_E)*3.89379372*10**(10)/E_nus,color='silver')
#plt.fill_between(E_nus,(sigs_E_bar-err_sigs_E_bar)*3.89379372*10**(10)/E_nus, (sigs_E_bar+err_sigs_E_bar)*3.89379372*10**(10)/E_nus,color='silver')
#plt.plot(a_E_nus,a_sigs_a)
plt.xlabel('E_{\nu} (GeV)')
plt.ylabel(r'$\sigma^{\nu}_{CC}\, /\, E_{\nu} \,\, \left(10^{-38}\, cm^2 \, GeV^{-1}\right)$')
plt.legend()
plt.savefig('figs_cross_L/cross_SM_CC_err.pdf',dpi=800,bbox_inches='tight')
plt.show()

## Version info, search paths, and metadata
print(lhapdf.version())
print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
print(lhapdf.paths())
# ...
